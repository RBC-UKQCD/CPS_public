#ifdef USE_SSE
#include <config.h>
CPS_START_NAMESPACE
/*! \file
  \brief  Routine used internally in the DiracOpWilson class.
  
  $Id: wilson_init.C,v 1.4 2013-01-08 21:09:25 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2013-01-08 21:09:25 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/sse/wilson_init.C,v 1.4 2013-01-08 21:09:25 chulwoo Exp $
//  $Id: wilson_init.C,v 1.4 2013-01-08 21:09:25 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.4 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/sse/wilson_init.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
/****************************************************************************/
/* 10/16/97                                                                 */
/*                                                                          */
/* wilson_int:                                                              */
/*                                                                          */
/* This routine performs all initializations needed before wilson funcs     */
/* are called. It sets the addressing related arrays and reserves memory    */
/* for the needed temporary buffers. It only needs to be called             */
/* once at the begining of the program (or after a wilson_end call)         */
/* before any number of calls to wilson funcs are made.                     */
/*                                                                          */
/* WARNING:                                                                 */
/*                                                                          */
/* This set of routines will work only if the node sublattices have         */
/* even number of sites in each direction.                                  */
/*                                                                          */
/****************************************************************************/

/*--------------------------------------------------------------------------*/
/* Include header files                                                     */
/*--------------------------------------------------------------------------*/
CPS_END_NAMESPACE
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>

#include <util/data_types.h>
#include <util/wilson.h>
#include <util/gjp.h>
#include <util/smalloc.h>
#include <util/verbose.h>
#include <util/error.h>
//#include <omp.h>
#include <util/omp_wrapper.h>


CPS_START_NAMESPACE

//! try to read the CPU's frequency from /proc/cpuinfo
double read_cpu_MHz()
{
  return 2270;

  
  using namespace std;

  // default is proton, but with minus to show it uses the default value.
  double freq=-2270; 

  ifstream ifs("/proc/cpuinfo");
  string str;

  double max_freq=freq;
  if( ifs.is_open() )
    while(!ifs.eof()){
      getline(ifs, str);
      if( str.find("cpu MHz") != string::npos ){
	::printf( "wilson_init %s\n", str.c_str() );
	size_t pos1 = str.find(":");
	string f = str.substr( pos1+1, string::npos);
	freq = atof( f.c_str() );
	if( freq > max_freq )  max_freq = freq;
      }
    }
  return max_freq;
}

#undef USE_TAG
//! allocate communication buffers and prepare communication
void wilson_init_comm(int dir, int block, Wilson *wilson_p)
{

    int idx;
    int sflag=+1; //?check
    size_t len = sizeof(IFloat)*block;

    idx=dir; 
    
    wilson_p->send_buf[idx]=(Float *) smalloc(block*sizeof(Float));
    wilson_p->recv_buf[idx]=(Float *) smalloc(block*sizeof(Float));

    wilson_p->msgmem[idx][1] =
      QMP_declare_msgmem((void *)(wilson_p->send_buf[idx]), len);
    wilson_p->msgmem[idx][0] =
      QMP_declare_msgmem((void *)(wilson_p->recv_buf[idx]), len);
    wilson_p->msghandle[idx][1] =
#ifdef USE_TAG
      QMP_declare_send_relative_tag(wilson_p->msgmem[idx][1], dir, -sflag, 0,idx);
#else
      QMP_declare_send_relative(wilson_p->msgmem[idx][1], dir, -sflag, 0);
#endif
    wilson_p->msghandle[idx][0] =
#ifdef USE_TAG
      QMP_declare_receive_relative_tag(wilson_p->msgmem[idx][0], dir, +sflag, 0,idx);
#else
      QMP_declare_receive_relative(wilson_p->msgmem[idx][0], dir, +sflag, 0);
#endif
    wilson_p->multiple[idx] = QMP_declare_multiple(wilson_p->msghandle[idx], 2);
    
    idx=dir+4;

    wilson_p->send_buf[idx]=(Float *) smalloc(block*sizeof(Float));
    wilson_p->recv_buf[idx]=(Float *) smalloc(block*sizeof(Float));

    wilson_p->msgmem[idx][1] =
      QMP_declare_msgmem((void *)(wilson_p->send_buf[idx]), len);
    wilson_p->msgmem[idx][0] =
      QMP_declare_msgmem((void *)(wilson_p->recv_buf[idx]), len);
    wilson_p->msghandle[idx][1] =
#ifdef USE_TAG
      QMP_declare_send_relative_tag(wilson_p->msgmem[idx][1], dir, +sflag, 0,idx);
#else
      QMP_declare_send_relative(wilson_p->msgmem[idx][1], dir, +sflag, 0);
#endif
    wilson_p->msghandle[idx][0] =
#ifdef USE_TAG
      QMP_declare_receive_relative_tag(wilson_p->msgmem[idx][0], dir, -sflag, 0,idx);
#else
      QMP_declare_receive_relative(wilson_p->msgmem[idx][0], dir, -sflag, 0);
#endif
    wilson_p->multiple[idx] = QMP_declare_multiple(wilson_p->msghandle[idx], 2);
}
	       

void wilson_init(Wilson *wilson_p)  /* pointer to Wilson type structure    */
{
  char *cname = " ";
  char *fname = "wilson_init(Wilson*)";
  VRB.Func(cname,fname);

  int spinor_words;             /* size of the spinor field on the         */
				/* sublattice checkerboard                 */
  int size;


  if(  GJP.Snodes() != 1 )
    ERR.NotImplemented(cname,fname,"5-th direction must be local\n");


/*--------------------------------------------------------------------------*/
/* Reserve memory for the node sublattice sizes                             */
/*--------------------------------------------------------------------------*/
  size = 4*sizeof(int);
  wilson_p->ptr = (int *) smalloc(size);
  if( wilson_p->ptr == 0)
    ERR.Pointer(cname,fname, "ptr");
  VRB.Smalloc(cname,fname,
	      "ptr", wilson_p->ptr, size);

/*--------------------------------------------------------------------------*/
/* Set the node sublattice sizes                                            */
/*--------------------------------------------------------------------------*/
  wilson_p->ptr[0] = GJP.XnodeSites();
  wilson_p->ptr[1] = GJP.YnodeSites();
  wilson_p->ptr[2] = GJP.ZnodeSites();
  wilson_p->ptr[3] = GJP.TnodeSites();
  wilson_p->vol[0] = wilson_p->ptr[0] * wilson_p->ptr[1] *
                     wilson_p->ptr[2] * wilson_p->ptr[3] / 2;
  wilson_p->vol[1] = wilson_p->vol[0];

/*--------------------------------------------------------------------------*/
/* Reserve memory for 2  temporary spinors                                  */
/* Use the af[] array to pass them (for this routines af are not            */
/* spin projected half spinors but instead they are full 4 component        */
/* spinors.)                                                                */
/*--------------------------------------------------------------------------*/
  spinor_words = SPINOR_SIZE * wilson_p->vol[0];

  wilson_p->af[0] = (IFloat *) smalloc(spinor_words*sizeof(Float));
  if(wilson_p->af[0] == 0)
    ERR.Pointer(cname,fname, "af[0]");
  VRB.Smalloc(cname,fname,
	      "af[0]", wilson_p->af[0], spinor_words*sizeof(Float));

  wilson_p->af[1] = (IFloat *) smalloc(spinor_words*sizeof(Float));
  if(wilson_p->af[1] == 0)
    ERR.Pointer(cname,fname, "af[1]");
  VRB.Smalloc(cname,fname,
	      "af[1]", wilson_p->af[1], spinor_words*sizeof(Float));

  VRB.Debug("x = %d\n", wilson_p->ptr[0]);
  VRB.Debug("y = %d\n", wilson_p->ptr[1]);
  VRB.Debug("z = %d\n", wilson_p->ptr[2]);
  VRB.Debug("t = %d\n", wilson_p->ptr[3]);
  VRB.Debug("vol0 = %d\n", wilson_p->vol[0]);
  VRB.Debug("vol1 = %d\n", wilson_p->vol[1]);
  VRB.Debug("af0 = %x\n", wilson_p->af[0]);
  VRB.Debug("af1 = %x\n", wilson_p->af[1]);




  // FOR SSEOMP

  const int   lx = wilson_p->ptr[0];
  const int   ly = wilson_p->ptr[1];
  const int   lz = wilson_p->ptr[2];
  const int   lt = wilson_p->ptr[3];
  const int   vol = wilson_p->vol[0];

#pragma omp parallel
  wilson_p->num_threads= omp_get_num_threads();

#if 0
  /*
   * prep buffers for each omp thread
   */
  int num_threads ;
#pragma omp parallel
  num_threads= omp_get_num_threads();



  IFloat **omp_chi = (IFloat**) smalloc( num_threads * sizeof(IFloat*) );
  
  // length of T length, which each thread will do calculations
  int omp_size_t = lt/num_threads;
  if( lt%num_threads ) omp_size_t++;
  int omp_size_chi_f = omp_size_t * lx*ly*lz * SPINOR_SIZE;
  
#pragma omp parallel for schedule(static,1)  
  for(int i=0; i< num_threads; ++i)
    omp_chi[i] = (IFloat*)smalloc( vol*SPINOR_SIZE*sizeof( IFloat) );
  
#endif

  int block[4];
  block[0]=HALF_SPINOR_SIZE*ly*lz*lt/2;
  block[1]=HALF_SPINOR_SIZE*lx*lz*lt/2;
  block[2]=HALF_SPINOR_SIZE*lx*ly*lt/2;
  block[3]=HALF_SPINOR_SIZE*lx*ly*lz/2;


  if(GJP.Xnodes()!=1){
    int dir=0;     wilson_init_comm(dir, block[dir], wilson_p);
  }
  if(GJP.Ynodes()!=1){
    int dir=1;     wilson_init_comm(dir, block[dir], wilson_p);
  }
  if(GJP.Znodes()!=1){
    int dir=2;     wilson_init_comm(dir, block[dir], wilson_p);
  }
  if(GJP.Tnodes()!=1){
    int dir=3;     wilson_init_comm(dir, block[dir], wilson_p);
  }
  
#if 0
  // clear chi
#pragma omp parallel for
  for(int i=0;i< vol*SPINOR_SIZE; i+=2)
    _mm_stream_pd(chi_p_f +i, _mm_setzero_pd() ) ;
#endif
  

//#ifdef PROFILE
  int NITR= wilson_p-> NITR = 1; // iteration for better timing
  wilson_p -> CPU_GHZ = read_cpu_MHz()/1000;
  
  IFloat MultFlops = wilson_p-> MultFlops=NITR*(double)vol * 1320.0;
  wilson_p->MultFlops_bnd=0;
    if( GJP.Xnodes()!=1 ) wilson_p->MultFlops_bnd += (double)NITR*	\
			    double (block[0]/HALF_SPINOR_SIZE)*156;
    if( GJP.Ynodes()!=1 ) wilson_p->MultFlops_bnd +=(double)NITR*	\
			    double (block[1]/HALF_SPINOR_SIZE)*156;
    if( GJP.Znodes()!=1 ) wilson_p->MultFlops_bnd +=(double)NITR*	\
			    double (block[2]/HALF_SPINOR_SIZE)*156;
    if( GJP.Tnodes()!=1 ) wilson_p->MultFlops_bnd +=(double)NITR*	\
			    double (block[3]/HALF_SPINOR_SIZE)*156;
    
    wilson_p->MultFlops_blk=MultFlops - wilson_p->MultFlops_bnd;
//#endif  

}


CPS_END_NAMESPACE
#endif
