#ifndef __READNERSC_H__
#define __READNERSC_H__

// a slight modification on ReadLattice class (if not inheritance)
// to enable parallel reading/writing of "Gauge Connection Format" lattice data

#include <cstdlib>	// exit()
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <sys/time.h>
#include <util/gjp.h>
#include <util/vector.h>
#include <util/time_cps.h>
#include <util/lattice.h>
#include <util/data_types.h>
#include <util/verbose.h>
#include <util/error.h>
#include <util/qioarg.h>
#include <util/fpconv.h>
#include <util/iostyle.h>
#include <util/latheader.h>
#include <alg/do_arg.h>

CPS_START_NAMESPACE

template <class headertype, int _Ndim, typename dtype>
class  ReadNERSC : public QioControl
{
  // which determines parallel reading or serial reading

 private:
  const char *cname;
  headertype hd;
  int UseParIO;
  const int default_concur=16;
  int Ndim;
  int data_per_site;

 public:
  // ctor for 2-step loading
 ReadNERSC(int ndata)
   : QioControl(), cname("ReadNERSC"),Ndim(_Ndim),data_per_site(ndata)
    { 
      //CK set architecture-dependent default IO style rather than hardcoding it.
      setDefault();
    }

  // ctor invoking loading behavior
 ReadNERSC(int ndata,dtype *data , const char * filename): QioControl(), cname("ReadNERSC"),Ndim(_Ndim),data_per_site(ndata)
    {
      setDefault();
    //CK set architecture-dependent default IO style rather than hardcoding it.
    QioArg rd_arg(filename);
	rd_arg.ConcurIONumber=default_concur;
    read(data,rd_arg);
  }

  // ctor invoking loading behavior
 ReadNERSC(int ndata,dtype * data, const QioArg & rd_arg) 
   : QioControl(), cname("ReadNERSC"),Ndim(_Ndim),data_per_site(ndata)
    {
      //CK set architecture-dependent default IO style rather than hardcoding it.
      setDefault();

      read(data,rd_arg);
    }
  
  virtual ~ReadNERSC() {}

  void read(dtype * data, const char *filename){
    QioArg rd_arg(filename);
	rd_arg.ConcurIONumber=default_concur;
    read(data,rd_arg);
  }


 public:
  inline void setParallel() { UseParIO = 1; }
  inline void setSerial() { UseParIO = 0; }
  inline void setDefault() { UseParIO = 1; }
  inline int parIO() const { return UseParIO; }

//  void read(dtype * ptr, const QioArg & rd_arg);

#undef PROFILE
#define PROFILE

//template <class headertype, int _Ndim, typename dtype>
//void ReadNERSC<headertype,_Ndim,dtype>::read(dtype *data, const QioArg & rd_arg)
void read(dtype *data, const QioArg & rd_arg)
{
  const char * fname = "read()";
  VRB.Func(cname,fname);
  VRB.Result(cname,fname,"data=%p\n",data);
  
  char loginfo[strlen(rd_arg.FileName) + 10];
  sprintf(loginfo,"Load %s",rd_arg.FileName);
  startLogging(loginfo);

#ifdef PROFILE
  struct timeval start,end;
  gettimeofday(&start,NULL);
#endif

  io_good = false;
  int error = 0;

  if (isRoot()) { // commander, analyze file header
    
    std::ifstream input(rd_arg.FileName);
    if ( !input.good() )
      {
	if(!UniqueID()){ printf("%s:%s Could not open file ptr %p for input.\n",cname,fname,rd_arg.FileName); fflush(stdout); }
	if(!UniqueID()){ printf("%s:%s Filename %s\n",cname,fname,rd_arg.FileName); fflush(stdout); }
	error = 1;
      }

    if(!error) {
      hd.setHeader(data_per_site);
      hd.read(input);
      VRB.Result(cname,fname,"data_per_site=%d hd.data_per_site=%d\n",data_per_site,hd.data_per_site);
      assert(data_per_site == hd.data_per_site);
      input.close();
    }
  }
  if(synchronize(error) != 0)
    ERR.FileR(cname, fname, rd_arg.FileName);
  log();

  broadcast(&hd.data_start, sizeof(std::streamoff));



  // check all conditions between FILE and GJP
  int nx = rd_arg.Xnodes() * rd_arg.XnodeSites();
  int ny = rd_arg.Ynodes() * rd_arg.YnodeSites();
  int nz = rd_arg.Znodes() * rd_arg.ZnodeSites();
  int nt = rd_arg.Tnodes() * rd_arg.TnodeSites();

  if(isRoot()) {
    if(hd.dimension[0] != nx || hd.dimension[1] != ny || hd.dimension[2] != nz || hd.dimension[3] != nt) {
      VRB.Flow(cname,fname,"Dimensions in file DISAGREE with GlobalJobParameter!\n");
      VRB.Flow(cname,fname,"In File: %d x %d x %d x %d\n",
	       hd.dimension[0],hd.dimension[1], hd.dimension[2], hd.dimension[3]);
      VRB.Flow(cname,fname,"In GJP:  %d x %d x %d x %d\n",nx, ny, nz, nt);
      error = 1;
    }

  }

  if(synchronize(error) != 0)  
    ERR.General(cname, fname, "Wrong Parameters Specified\n");

  // see if file Floating Points is acceptable
  if(isRoot()) {
    fpconv.setFileFormat(hd.floating_point);
  }
  VRB.Flow(cname,fname,"FileFormat=%d",hd.floating_point);
  broadcastInt((int*)&fpconv.fileFormat);
  if(fpconv.fileFormat == FP_UNKNOWN) {
    ERR.General(cname,fname, "Data file Floating Point Format UNKNOWN\n");
  }
  
  VRB.Flow(cname,fname,"A copy of header info from file:\n");
  if(isRoot())  hd.show();

//  int data_per_site = hd.recon_row_3 ? 4*12 : 4*18;

  // read lattice data, using parallel style or serial (all on node 0) style
  unsigned int csum;

  //CK: removed hardcoded IO style. Replaced with default IO style in constructor
// #if TARGET != QCDOC
//   setSerial();
// #endif

  log();
  const size_t chars_per_site = data_per_site * sizeof(dtype);

  VRB.Flow(cname,fname,"Reading configuation to address: %p\n", rd_arg.StartConfLoadAddr);
  if(parIO()) {
    ParallelIO pario(rd_arg);
    if(!doGparityReconstructUstarField() ){
      if(!UniqueID()) printf("ReadNERSC is disabling reconstruction bit on the ParallelIO object\n");
      pario.disableGparityReconstructUstarField();
    }

    if(! pario.load((char*)data, data_per_site, chars_per_site,
		    hd, fpconv, 4, &csum))  
      ERR.General(cname,fname,"Load Failed\n");  // failed to load
  }
#if 1
  else {
    SerialIO serio(rd_arg);

    if(! serio.load((char*)data, data_per_site, chars_per_site,
		    hd, fpconv, 4, &csum))  
      ERR.General(cname,fname,"Load Failed\n");  // failed to load
  }
#endif

  log();

//  printf("Node %d: lattice read csum=%x\n",UniqueID(),csum);
  //  cout << "loader finish, csum = " << hex << csum << dec << endl << endl;
  //  cout << "loader done" << endl << endl;


  // After reading...
  // STEP 1: checksum
  if(rd_arg.Scoor() == 0)
    csum = globalSumUint(csum);
  else
    globalSumUint(0);

  if(isRoot()) {
    if( hd.checksum != csum ) {
      VRB.Flow(cname,fname, "CheckSUM error !! Header: %x  Host calc: %x\n",hd.checksum,csum);
      
      printf("Node %d: CheckSUM error !! Header: %x  Host calc: %x\n",UniqueID(),hd.checksum,csum);
      error = 1;
    }
    else
      VRB.Flow(cname,fname,"CheckSUM is ok\n");
  }

  if(synchronize(error) != 0) 
    ERR.General(cname, fname, "Checksum error\n");


  // STEP 2: reconstruct row 3
  int size_matrices( rd_arg.XnodeSites() * rd_arg.YnodeSites() 
		     * rd_arg.ZnodeSites() * rd_arg.TnodeSites() * 4); 


#ifdef PROFILE
  gettimeofday(&end,NULL);
  print_flops(cname,"read",0,&start,&end);
#endif

  io_good = true;

  log();
//  finishLogging();

  VRB.FuncEnd(cname,fname);
}

};

CPS_END_NAMESPACE
#endif
