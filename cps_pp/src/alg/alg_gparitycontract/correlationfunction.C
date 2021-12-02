#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <util/qcdio.h>
#ifdef PARALLEL
#include <comms/sysfunc_cps.h>
#endif
#include <comms/scu.h>
#include <comms/glb.h>

#include <util/lattice.h>
#include <util/time_cps.h>
#include <util/smalloc.h>

#include <util/command_line.h>

#include<unistd.h>
#include<config.h>

#include <util/gjp.h>
#include <alg/correlationfunction.h>
#include <util/omp_wrapper.h>
CPS_START_NAMESPACE

CorrelationFunction::CorrelationFunction(const char *_label, const ThreadType &thread): wick(NULL),wick_threaded(NULL),ncontract(0),threadtype(thread),global_sum_on_write(true), label(_label){
  time_size=GJP.Tnodes()*GJP.TnodeSites();
}
CorrelationFunction::CorrelationFunction(const char *_label, const int &n_contractions, const ThreadType &thread): wick(NULL),wick_threaded(NULL),ncontract(0),threadtype(thread),global_sum_on_write(true), label(_label){
  time_size=GJP.Tnodes()*GJP.TnodeSites();
  setNcontractions(n_contractions);
}
CorrelationFunction::CorrelationFunction(): label("NoLabel"), wick(NULL),wick_threaded(NULL),ncontract(0),threadtype(UNTHREADED),global_sum_on_write(true){
  time_size=GJP.Tnodes()*GJP.TnodeSites();
}

//Set the thread type. Cannot be changed once memory has been allocated
void CorrelationFunction::setThreadType(const ThreadType &thread){
  if(wick!=NULL || wick_threaded!=NULL){
    if(!UniqueID()) printf("CorrelationFunction::setThreadType(const ThreadType &thread) : Thread type cannot be changed once internal memory has been allocated\n");
    exit(-1);
  }
  threadtype = thread;
}


void CorrelationFunction::setNcontractions(const int &n){
  clear();

  ncontract = n;
  wick = (Rcomplex**)smalloc(n*sizeof(Rcomplex*)); //allocate array of pointers (I wish we could just use std::vector here!)
    
  //contiguous memory slot
  Rcomplex *stack = (Rcomplex*)smalloc(n*time_size*sizeof(Rcomplex));
    
  for(int i=0;i<n;i++){
    wick[i] = stack+i*time_size;
    for(int t=0;t<time_size;t++){ wick[i][t].real(0.0); wick[i][t].imag(0.0); }
  }
  
  if(threadtype == THREADED){
    max_threads = omp_get_max_threads();
    wick_threaded = (Rcomplex***)smalloc(max_threads*sizeof(Rcomplex**));
    for(int s=0;s<max_threads;s++){
      wick_threaded[s] = (Rcomplex**)smalloc(n*sizeof(Rcomplex*));

      //contiguous memory slot
      Rcomplex *stack = (Rcomplex*)smalloc(n*time_size*sizeof(Rcomplex));
    
      for(int i=0;i<n;i++){
	wick_threaded[s][i] = stack+i*time_size;
	for(int t=0;t<time_size;t++){ wick_threaded[s][i][t].real(0.0); wick_threaded[s][i][t].imag(0.0); }
      } 
    }
  }

}
CorrelationFunction::~CorrelationFunction(){
  clear();
}

void CorrelationFunction::clear(){
  if(ncontract!=0){
    sfree(wick[0]); //wick array was created as contiguous mem block
    sfree(wick);    
    if(threadtype == THREADED){
      for(int s=0;s<max_threads;s++){
	sfree(wick_threaded[s][0]);
	sfree(wick_threaded[s]);
      }
      sfree(wick_threaded);
    }
  }
  wick = NULL;
  wick_threaded = NULL;
  global_sum_on_write = true;
}


Rcomplex & CorrelationFunction::operator()(const int &contraction_idx, const int &t){
  //if(threadtype==THREADED){ ERR.General("CorrelationFunction","operator(int,int)","Cannot call this version of operator() in a threaded environment"); }
  return wick[contraction_idx][t];
}
Rcomplex & CorrelationFunction::operator()(const int &thread_idx, const int &contraction_idx, const int &t){
  //if(threadtype==UNTHREADED){ ERR.General("CorrelationFunction","operator(int,int,int)","Cannot call this version of operator() in a non-threaded environment"); }
  if(thread_idx >= max_threads){ ERR.General("CorrelationFunction","operator(int,int,int)","thread_idx %d is >= max_threads %d\n",thread_idx,max_threads); }

  return wick_threaded[thread_idx][contraction_idx][t];
}

void CorrelationFunction::sumThreads(){
  if(threadtype == THREADED){
    //do the thread sum to wick
    for(int i=0;i<ncontract;i++){
      for(int t=0;t<time_size;t++){
	wick[i][t] = 0.0;
	for(int s=0;s<max_threads;s++){
	  wick[i][t] += wick_threaded[s][i][t];
	}
      }
    }
  }
}


void CorrelationFunction::sumLattice(){
  sumThreads();

  //now do the lattice sum
  for(int i=0;i<ncontract;i++){
    for(int t=0;t<time_size;t++){
      slice_sum( (Float*)&wick[i][t], 2, 99); //2 for re/im, 99 is a *magic* number (we are abusing slice_sum here)
    }
  }
  setGlobalSumOnWrite(false); //disable further sums
}
void CorrelationFunction::write(const char *file){
  FILE *fp;
  if ((fp = Fopen(file, "w")) == NULL) {
    ERR.FileW("CorrelationFunction","write(const char *file)",file);
  }
  write(fp);
  Fclose(fp);
}
void CorrelationFunction::write(FILE *fp){
  if(global_sum_on_write) sumLattice(); //sum the correlation function over all nodes

  Fprintf(fp,"%s\n",label.c_str());
  Fprintf(fp,"%d contractions\n",ncontract);

  for(int c=0;c<ncontract;c++){
    Fprintf(fp,"Contraction %d\n",c);

    for(int t=0; t<time_size; t++)
      Fprintf(fp, "%d %.16e %.16e\n",t, wick[c][t].real(), wick[c][t].imag());
  }
}

CPS_END_NAMESPACE
