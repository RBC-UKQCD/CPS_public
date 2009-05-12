#ifdef USE_QMP
/*! \file
  \brief  Definition of parallel transport definitions for QCDOC.
  
  $Id: pt_ordering.C,v 1.5 2009-05-12 21:50:01 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2009-05-12 21:50:01 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_base/qmp/pt_ordering.C,v 1.5 2009-05-12 21:50:01 chulwoo Exp $
//  $Id: pt_ordering.C,v 1.5 2009-05-12 21:50:01 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: pt_ordering.C,v $
//  $Revision: 1.5 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_base/qmp/pt_ordering.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
#include "asq_data_types.h"
#include "pt_int.h"

static unsigned long PEC = 0xb0000000;
static unsigned long PLB = 0xb0000000;

#define TESTING
#undef TESTING
#undef CPP
//CPS_START_NAMESPACE
//void dirac_cmv_jcw_agg_cpp( int sites, long chi, long u,long in, long out);

//External function definitions
extern "C"{

  //matrix multiply for checkerboarded fields
  void pt_cmm_cpp(int sites, long u, long in, long out, long gauge_field);
  void pt_cmm_dag_cpp(int sites, long u, long in, long out, long gauge_field);

  //------------------------------------------------------------------------
  //C++ routines
#ifdef CPP
  //Matrix multiplication for full matrix fields
  void cmm_agg_cpp(gauge_agg *chi, matrix *phi, matrix *result, int counter);
  void cmv_agg_cpp( int sites, long u,long in, long out);
  #define partrans_cmm_agg(A,B,C,D) cmm_agg_cpp(A,B,C,D)
  #define partrans_cmv_agg(A,B,C,D) cmv_agg_cpp(A,B,C,D)

  //matrix vector multiply for checkerboarded fields
  void pt_cmv_cpp(int sites, ind_agg *agg, double *gauge_field, double *src, double *dest);
  void pt_cmv_dag_cpp(int sites, ind_agg *agg, double *gauge_field, double *src, double *dest);
  void pt_cmv_pad_cpp(int sites, ind_agg *agg, double *gauge_field, double *src, double *dest);
  void pt_cmv_dag_pad_cpp(int sites, ind_agg *agg, double *gauge_field, double *src, double *dest);
  #define partrans_cmv(A,B,C,D,E) pt_cmv_cpp(A,B,C,D,E)
  #define partrans_cmv_dag(A,B,C,D,E) pt_cmv_dag_cpp(A,B,C,D,E)
  #define partrans_cmv_pad(A,B,C,D,E) pt_cmv_pad_cpp(A,B,C,D,E)
  #define partrans_cmv_dag_pad(A,B,C,D,E) pt_cmv_dag_pad_cpp(A,B,C,D,E)
  //--------------------------------------------------------------------------
  //Assembly Routines
#else
  //Matrix multiplication for full matrix fields
  void pt_cmm_agg(gauge_agg *chi, matrix *phi,matrix *result, int counter);
  //void cmm_agg(gauge_agg *chi, matrix *phi,matrix *result, int counter);
  void pt_asqtad_agg( int sites, long chi, long u,long in, long out);
  void pt_asqtad_agg_s( int sites, long chi, long u,long in, long out);
  #define partrans_cmm_agg(A,B,C,D) pt_cmm_agg(A,B,C,D)
  #define partrans_cmv_agg(A,B,C,D) pt_asqtad_agg(A,0,B,C,D)

  void pt_cmv(int count, ind_agg *ind, double *gauge, double *src, double *dest);
  void pt_cmv_pad(int count, ind_agg *ind, double *gauge, double *src, double *dest);
  void pt_cmv_dag(int count, ind_agg *ind, double *gauge, double *src, double *dest);
  void pt_cmv_dag_pad(int count, ind_agg *ind, double *gauge, double *src, double *dest);
  void pt_cmv_s(int count, ind_agg *ind, float *gauge, float *src, float *dest);
  void pt_cmv_pad_s(int count, ind_agg *ind, float *gauge, float *src, float *dest);
  void pt_cmv_dag_s(int count, ind_agg *ind, float *gauge, float *src, float *dest);
  void pt_cmv_dag_pad_s(int count, ind_agg *ind, float *gauge, float *src, float *dest);
  #define partrans_cmv(A,B,C,D,E) pt_cmv(A,B,C,D,E)
  #define partrans_cmv_dag(A,B,C,D,E) pt_cmv_dag(A,B,C,D,E)
  #define partrans_cmv_pad(A,B,C,D,E) pt_cmv_pad(A,B,C,D,E)
  #define partrans_cmv_dag_pad(A,B,C,D,E) pt_cmv_dag_pad(A,B,C,D,E)
#endif
  //--------------------------------------------------------------------------

  void pt_copy(int count, ind_agg *ind, double *src, double *dest);
  void pt_copy_pad(int count, ind_agg *ind, double *src, double *dest);
  void pt_copy_s(int count, ind_agg *ind, float *src, float *dest);
  void pt_copy_pad_s(int count, ind_agg *ind, float *src, float *dest);

  void pt_copy_buffer(int n, long src, long dest, long ptable);
  // Assembler copying routines
  void copy_matrix(IFloat *res, IFloat *src, int *length, 
		   unsigned long *res_ptr, unsigned long *src_ptr);
  void copy_gauge(IFloat *res, struct gauge_agg *src, int *length,
		  unsigned long *res_ptr);
  // This is perhaps overkill but gives a couple of extra flops
  // cross_look - all input fields are lookup and sum to result
  void cross_look(IFloat *result, Float *fac, const IFloat *chi, const IFloat *phi,  
	     int counter, unsigned long *src, unsigned long *dest);
  // cross_lin - one input field is linear and sum to result
  void cross_lin(IFloat *result, Float *fac, const IFloat *chi, const IFloat *phi,  
	     int counter, unsigned long *dest, unsigned long *dest2);
  // cross_over_look - all input fields are lookup and overwrite result
  void cross_over_look(IFloat *result, Float *fac, const IFloat *chi, const IFloat *phi,  
	     int counter, unsigned long *src, unsigned long *dest);
  // cross_over_lin - one input field is linear and overwrite result
  void cross_over_lin(IFloat *result, Float *fac, const IFloat *chi, const IFloat *phi,  
	     int counter, unsigned long *dest, unsigned long *dest2);
  //Copies a vectors from v to u
  void copy_vector(IFloat *u, IFloat *v, int *length, unsigned long *dest, unsigned long *src);

  //---------------------------------------------------------------------------
  
  //---------------------------------------------------------------------------

  void m1m2_lookup(matrix *result, matrix *m1, matrix *m2, int length,
		   unsigned long *dest, unsigned long *dest2, unsigned long *src);
  void m1m2_lookup_copy(matrix *result2, matrix *result, matrix *m1, matrix *m2, 
			int length, unsigned long *dest2,  
			unsigned long *dest, unsigned long *dest3, 
			unsigned long *src);
  void m1m2_lin_copy(matrix *result2, matrix *result, matrix *m1, matrix *m2, 
		     int length, unsigned long *dest2,
		     unsigned long *dest, unsigned long *dest3);
  
}
inline  void pt_cmm_agg_print(gauge_agg *chi, matrix *phi,matrix *result, int counter){
   printf("pt_cmm_agg(%p %p %p %d)\n",chi,phi,result,counter);
//    for(int i =0;i<2*counter;i++){
//      printf("%d: %d %d\n",i,chi[i].src,chi[i].dest);
//    }
   printf("pt_cmm_agg(%p %p %p %d) done \n",chi,phi,result,counter);
}

inline  void cross_over_lin_cpp(IFloat *result, Float *fac, const IFloat *chi,
 const IFloat *phi,  int counter, unsigned long *src, unsigned long *dest){
    printf("cross_over_lin(%p %0.4f %p %p %d %p %p)\n",
    result,*fac,chi,phi,counter,src,dest);
    for(int i =0;i<counter;i++){
      printf("%d: %d %d\n",i,src[i],dest[i]);
    }
    cross_over_lin(result,fac,chi,phi,counter,src,dest);
    printf("cross_over_lin(%p %0.4f %p %p %d %p %p) done\n");
}

inline  void cross_over_look_cpp(IFloat *result, Float *fac, const IFloat *chi,
 const IFloat *phi,  int counter, unsigned long *src, unsigned long *dest){
    printf("cross_over_look(%p %0.4f %p %p %d %p %p)\n",
    result,*fac,chi,phi,counter,src,dest);
    for(int i =0;i<counter;i++){
      printf("%d: %d %d\n",i,src[i],dest[i]);
    }
    cross_over_look(result,fac,chi,phi,counter,src,dest);
    printf("cross_over_look(%p %0.4f %p %p %d %p %p) done\n");
}

#ifdef ASQD_SINGLE
#define pt_asqtad_agg(A,B,C,D,E) pt_asqtad_agg_s(A,B,C,D,E)
#define pt_cmv(A,B,C,D,E) pt_cmv_s(A,B,C,D,E)
#define pt_cmv_dag(A,B,C,D,E) pt_cmv_dag_s(A,B,C,D,E)
#define pt_cmv_pad(A,B,C,D,E) pt_cmv_pad_s(A,B,C,D,E)
#define pt_cmv_dag_pad(A,B,C,D,E) pt_cmv_dag_pad_s(A,B,C,D,E)
#define pt_copy_pad(A,B,C,D) pt_copy_pad_s(A,B,C,D)
#define pt_copy(A,B,C,D) pt_copy_s(A,B,C,D)
#endif


//dest=src
void PT::cpy (IFloat *dest, IFloat *src){
  for(int i=0;i<18;i++)
    dest[i]=src[i];
}

//dest=src.Dagger()
void PT::dag_cpy (IFloat *dest, IFloat *src){
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++){
      dest[2*(3*i+j)]=src[2*(3*j+i)];
      dest[2*(3*i+j)+1]=-src[2*(3*j+i)+1];
    }
}

//Returns lexical index associated with coordinate x[4]
//where the 0th coordinate runs fastest, 3rd coordinate runs slowest
int PT::lex_xyzt(int *x){
//  printf("lex_xyzt(%d %d %d %d)\n",x[0],x[1],x[2],x[3]);
  int result = x[0] + size[0]*(x[1]+size[1]*(x[2]+size[2]*x[3]));
  return result;
}

//Returns checkerboard index associated with coordinate x[4]
int PT::lex_xyzt_cb_o(int *x){
//  printf("lex_xyzt_cb_o(%d %d %d %d)\n",x[0],x[1],x[2],x[3]);
  int result = x[0] + size[0]*(x[1]+size[1]*(x[2]+size[2]*x[3]));
  if ( (x[0]+x[1]+x[2]+x[3]+evenodd)%2 == 0) result = result/2+vol/2;
  else result = result/2;
  return result;
}

//Returns checkerboard index associated with coordinate x[4]
int PT::lex_xyzt_cb_e(int *x){
//  printf("lex_xyzt_cb_o(%d %d %d %d)\n",x[0],x[1],x[2],x[3]);
  int result = x[0] + size[0]*(x[1]+size[1]*(x[2]+size[2]*x[3]));
  if ( (x[0]+x[1]+x[2]+x[3]+evenodd)%2 == 1) result = result/2+vol/2;
  else result = result/2;
  return result;
}


//---------------------------------------------------------------------------
//Returns index for fields in the STAG storage order on lattice
//sites of a given parity
int PT:: lex_txyz_cb(int *x)
{
//  printf("lex_txyz_cb(%d %d %d %d)\n",x[0],x[1],x[2],x[3]);
  int result = x[3]+size[3]*(x[0]+size[0]*(x[1]+size[1]*x[2]));
  return result/2;
}
//---------------------------------------------------------------------------

//Returns index associated with x[4] for txyz ordering
int PT:: lex_txyz(int *x){
  return  (x[3] + size[3]*(x[0]+size[0]*(x[1]+size[1]*x[2])))/2 ;
}

//Returns first index associated with the surface x[3] = 0, on
//a checkerboarded lattice
int PT:: LexSurface(int *x){
  return  (x[0]+size[0]*(x[1]+size[1]*x[2]))/2 ;
}

//Returns index associated with gauge link in the mu direction 
//and coordinate x
int PT::lex_g_xyzt(int *x, int mu){
  int temp =  lex_xyzt(x);
  return (temp*NDIM + mu);
}

//---------------------------------------------------------------------------
//Returns block ordering for the gauge fields, where all directions
//are stored in one block and sites ordered txyz
int PT::lex_g_txyz(int *x, int mu){
  int temp = mu*vol+x[3] + size[3]*(x[0]+size[0]*(x[1]+size[1]*x[2]));
  return temp;
}

//Returns block ordering for the gauge fields, where all directions
//are stored in one block with sites checkerboarded txyz
int PT::lex_g_txyz_cb(int *x, int mu){
  int result = (x[3]+size[3]*(x[0]+size[0]*(x[1]+size[1]*x[2])))/2;
  return mu*vol + result + ((x[0]+x[1]+x[2]+x[3])%2)*vol/2;
}
//---------------------------------------------------------------------------

//Returns index associated with gauge link in the mu direction and 
//coordinate x for checkerboarded storage
int PT:: lex_g_xyzt_cb_o(int *x, int mu){
  int temp =  lex_xyzt_cb_o(x);
  return (mu*vol+temp);
}
int PT:: lex_g_xyzt_cb_e(int *x, int mu){
  int temp =  lex_xyzt_cb_e(x);
  return (mu*vol+temp);
}

// Calculate the required offset given the direction and hop
int PT::set_offset(int dir, int hop) {

  // if positive direction then start at 0
  if (dir%2 == 0) return 0;

  int temp=1;
  int offset=0;
  for(int i=0;i<dir/2+1;i++){
    offset = temp*(size[i]-hop);
    temp *= size[i];
  }
  return offset;

}

void PT::set_hop_pointer() {

  char *fname = "set_hop_pointer()";

//  VRB.Func("PT",fname);
  //Actual memory usage of vectors
  int vlen = VECT_LEN*sizeof(IFloat);
  int vlen2 =VECT_LEN_OUT*sizeof(IFloat);

  int x[NDIM], nei[NDIM];
  
  //Counts how many parallel transports of given length and direction are local
  //and non-local, respectively
  int hp_local_count[MAX_HOP][2*NDIM];
  int hp_non_local_count[MAX_HOP][2*NDIM];
  int hop, i;


  //Initialize local and non-local hop counters.
  for (hop=0; hop<MAX_HOP; hop++) {
    for (i=0; i<2*NDIM; i++) {
      hp_non_local_count[hop][i] = 0;
      hp_local_count[hop][i] = 0;
    }
  }
  
  //For a given length of the parallel transport
  for (hop = 1; hop <= MAX_HOP; hop++) {
    hop_pointer **h_l = hp_l[hop-1];
    hop_pointer **h_nl = hp_nl[hop-1];

    //Local and non-local counts for given length of the hop
    int *local_count = hp_local_count[hop-1];
    int *non_local_count = hp_non_local_count[hop-1];

    //Loop over all directions
    for (i=0; i<NDIM; i++) {

      //Total number of sites that require non-local communication
      int non_local_check = hop*non_local_chi[i*2];
      //Total number of sites where parallel transport can be done locally
      int local_check = vol - non_local_check;

      //Loop through all the sites on the lattice
      //nei represents the coordinates of the neighboring site.
      for(x[3]=0,nei[3]=0;x[3]<size[3];x[3]++,nei[3]++)
	for(x[2]=0,nei[2]=0;x[2]<size[2];x[2]++,nei[2]++)
	  for(x[1]=0,nei[1]=0;x[1]<size[1];x[1]++,nei[1]++)
	    for(x[0]=0,nei[0]=0;x[0]<size[0];x[0]++,nei[0]++){

	      //This is the parallel transport of the field in the 
	      //negative direction to another node
	      //"Positive hop" because the link variable points in the 
	      //positive direction, even though the resulting field is 
	      //"transported" in the negative direction
	      // positive direction

	      if((x[i] < hop) && (!local[i])){
		//This calculates the neighbor coordinate
		nei[i] = size[i]-hop+x[i];  

		//Sets the index for source and destination
		(h_nl[2*i]+non_local_count[2*i])->src = non_local_count[2*i]*vlen;
		(h_nl[2*i]+non_local_count[2*i])->dest = LexVector(nei)*vlen2;

		//Increments the non-local count
		non_local_count[i*2]++;

		//Make sure we haven't gone over the non non-local check
		if (non_local_count[i*2]>non_local_check)
		  fprintf(stderr,
			"%s:non_local_count[%d](%d)>non_local_check[%d](%d)\n",
			 fname,2*i,non_local_count[2*i],2*i,non_local_check);
		//The rest of the parallel transports in the local volume can 
		//be handled locally
	      } else {
		//Calculate the new coordinate
		nei[i] = (size[i]+x[i]-hop)%size[i];

		//if ( size[i] >2){
		//Calculate the index for the source and the destination
		(h_l[2*i]+local_count[2*i])->src = LexVector(x)*vlen;
		(h_l[2*i]+local_count[2*i])->dest = LexVector(nei)*vlen2;
                //}
		
		//Increment the local count
		local_count[i*2]++;
		//Make sure we haven't exceeded the number of local sites
		if (local_count[i*2]>local_check)
		  fprintf(stderr,"%s:local_count[%d](%d)>local_check[%d](%d)\n",
			      fname,2*i,local_count[2*i],2*i,local_check);
	      }
	      
	      //Consider hopping in the negative direction, which is parallel 
	      //transport in the positive direction
	      // negative direction
	      if( (x[i] >= (size[i]-hop)) && (!local[i])){
		//Calculate the non-local coordinate for this hop
		nei[i] = (x[i]+hop)%size[i];
		//Calculate source and destination indices
		(h_nl[2*i+1]+non_local_count[2*i+1])->src = non_local_count[2*i+1]*vlen;
		(h_nl[2*i+1]+non_local_count[2*i+1])->dest = LexVector(nei)*vlen2;

		//Increment the non-local count, check that bounds have not 
		//been exceeded
		non_local_count[i*2+1]++;
		if (non_local_count[i*2]>non_local_check)
		  fprintf(stderr,"%s:non_local_count[%d](%d)>non_local_check[%d](%d)\n",
			      fname,2*i,non_local_count[2*i],2*i,non_local_check);
	      } else {
		//Calculate the local coordinate for this hop
		nei[i] = (x[i]+hop)%size[i];
		//Calculate source and destination indices
		//if ( size[i] >2){
		(h_l[2*i+1]+local_count[2*i+1])->src = LexVector(x)*vlen;
		(h_l[2*i+1]+local_count[2*i+1])->dest = LexVector(nei)*vlen2;
		//}
		//Increment local count, check that bounds not exceeded
		local_count[i*2+1]++;
		if (local_count[i*2]>local_check)
		  fprintf(stderr,"%s:local_count[%d](%d)>local_check[%d](%d)\n",
			      fname,2*i,local_count[2*i],2*i,local_check);
	      }
	      // Need to reset the neighbour pointer
	      nei[i] = x[i];
	    }
    }
  }
//  VRB.Func("PT",fname);
//  exit(44);
}
//CPS_END_NAMESPACE
#endif
