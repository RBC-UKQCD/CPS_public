#include <config.h>

//#define USE_BLAS

CPS_START_NAMESPACE
/*! \file
  \brief  Definition of DiracOp class Low Mode Approximation 

*/

CPS_END_NAMESPACE
#include <util/dirac_op.h>
#include <util/lattice.h>
#include <util/smalloc.h>
#include <util/vector.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/error.h>
#include <util/checksum.h>
#include <comms/glb.h>
#include <math.h>
#include <stdio.h>
//#include <qcdocos/gint.h>

#include <comms/sysfunc_cps.h>
#include <util/eigen_container.h>

#include<util/time_cps.h>

// Is STL available ?
#include <vector>

#ifdef USE_BLAS
#if TARGET == BGL
#define cblas_zdotc zdotc
#define cblas_zaxpy zaxpy

#include "essl.h"
#else
#include <util/qblas_extend.h>
#endif
#endif


CPS_START_NAMESPACE

int DiracOp::InvLowModeApprox(
			Vector *out, 
			Vector *in, 
			char* fname_eig_root,
			int neig,
			Float *true_res){

  time_elapse();
  
  int f_size_cb;     // Node checkerboard size of the fermion field
  int itr;                       // Current number of CG iterations
  int max_itr;                       // Max number of CG iterations
  Float stp_cnd;                   // Stop if residual^2 <= stp_cnd
  Float res_norm_sq_prv;          // The previous step |residual|^2
  Float res_norm_sq_cur;           // The current step |residual|^2 
  Float a;
  Float b;
  Float d;
  int i, j;
  char *fname = "InvLowModeApprox(V*,V*,C,I,F*)";

  Float mdagm_time=0.;
  Float gsum_time=0.;
  Float linalg_time=0.;

//------------------------------------------------------------------
// Initializations
//------------------------------------------------------------------

// Set the source vector pointer
//------------------------------------------------------------------
  Vector *src = in;
  IFloat *src_tmp = (IFloat *)src;

// Set the solution vector pointer
//------------------------------------------------------------------
  Vector *sol = out;

// Set the node checkerboard size of the fermion field
//------------------------------------------------------------------

  if(lat.Fclass() == F_CLASS_CLOVER) {
    f_size_cb = GJP.VolNodeSites() * lat.FsiteSize() / 2;
  } else {
    f_size_cb = GJP.VolNodeSites() * lat.FsiteSize() / (lat.FchkbEvl()+1);
  }


  const int n_fields =  GJP.SnodeSites();  //   *nk ; 
  const int f_size_per_site = lat.FsiteSize() / GJP.SnodeSites()  / (lat.FchkbEvl()+1);


  char fname_eig_root_bc[1024];
  snprintf(fname_eig_root_bc,1024,"%s.bc%d%d%d%d",fname_eig_root,
	   GJP.Bc(0),GJP.Bc(1),GJP.Bc(2),GJP.Bc(3));
	   


  //----------------------------------------------------------------------------------------------//
  // Dear Chulwoo :
  //  Here we allocate a new element of EigenCacheList
  //  In principle, it would be better to deallocate in this class, not to forget to do so.
  //
  //  But the cache list needes to live after the dirac operator destory, as the same eigen vector
  //  would be needed again (at least 12 times for the propagator, and more for other source locations)
  //  So cleanup  is  on the hand of user in either algorithm or in main.C. 
  //  For example, user could Cleanup(EigenCacheList) when one is done with the current configuration.
  //  Of course, forgetting to cleanup will cause a potentially huge memory leak... may force with you.
  //------------------------------------------------------------------------------------------------//
  

#define EIGEN_CACHE_IN_INV_LOWMODE_APPROX_YES

#ifdef EIGEN_CACHE_IN_INV_LOWMODE_APPROX_YES
  // search for eigen cache
  EigenCache* ecache;
  if( (ecache=EigenCacheListSearch( fname_eig_root_bc, neig ))== 0 ){

    ERR.General(cname,fname,"Eigenvector cache does not exist: neig %d name %s \n",neig,fname_eig_root_bc);
    //ecache = new EigenCache();  
    //EigenCacheList. push_back( ecache );
  }
#else
  EigenCache* ecache=0;
#endif
  
  Float mass = dirac_arg->mass;
  EigenContainer eigcon( lat, fname_eig_root_bc, neig, f_size_per_site, n_fields, ecache );

  Float* eval = eigcon.load_eval();
  //
  // make eigen values squared if needed
  //
  switch (dirac_arg->RitzMatOper) {
  case  MAT_HERM :
  case  MATPC_HERM : 
    for(int i=0;i<neig;++i)  eval[i]= eval[i]*eval[i];
    break;
  case  MATPCDAG_MATPC  :
  case  NEG_MATPCDAG_MATPC :
  case  MATDAG_MAT :
  case  NEG_MATDAG_MAT :
  case  MATDAG_MAT_NORM :
  case  NEG_MATDAG_MAT_NORM :
  case  MATPCDAG_MATPC_SHIFT :
    break;
  default:
    ERR.General(cname,fname,"RiztMatType %d not implemented\n",
                dirac_arg->RitzMatOper);
  }

  sol->VecZero( f_size_cb );

  for(int iev=0;iev<neig;++iev) {
    
    //time_elapse();
    Vector* evec = eigcon.nev_load( iev );
    //print_time("inv_lowmode_approx","loading", time_elapse());

#ifndef USE_BLAS
    Complex z = evec->CompDotProductGlbSum( src, f_size_cb );
#else
    Complex z;

#if TARGET != BGL
    cblas_zdotc_sub( f_size_cb / 2,
		     (double*)evec, 1, (double*)src, 1, (double*)&z);
#else
    *(complex<double>*)&z=cblas_zdotc( f_size_cb / 2,
		     (complex<double>*)evec, 1, (complex<double>*)src, 1);
#endif
    
    glb_sum((double*)&z);
    glb_sum((double*)&z+1);
#endif
    //print_flops("inv_lowmode_approx","compdotprod",f_size_cb*3, time_elapse());
    
    z /= eval[ iev ];

#ifndef USE_BLAS
    sol -> CTimesV1PlusV2(z, evec, sol, f_size_cb );
#else
#if TARGET == BGL
    cblas_zaxpy( f_size_cb / 2, *(complex<double>*)&z, (complex<double>*)evec, 1, (complex<double>*)sol,1);
#else
    cblas_zaxpy( f_size_cb / 2, (double*)&z, (Float*)evec, 1, (Float*)sol,1);
#endif    
#endif
    //print_flops("inv_lowmode_approx","ctimesv1pv2", f_size_cb*4, time_elapse());

}

  print_flops(cname,fname, f_size_cb*7*neig, time_elapse());
  
  *true_res = 0.0; // could compute the residue norm here, but saving a dirac operator for now
  return neig;
}


void DiracOp::InvLowModeProj( Vector *in, 
			char* fname_eig_root,
			int neig){

  time_elapse();
  
  int f_size_cb;     // Node checkerboard size of the fermion field
  int itr;                       // Current number of CG iterations
  int max_itr;                       // Max number of CG iterations
  Float stp_cnd;                   // Stop if residual^2 <= stp_cnd
  Float a;
  Float b;
  Float d;
  int i, j;
  char *fname = "InvLowModeApprox(V*,V*,C,I,F*)";

//------------------------------------------------------------------
// Initializations
//------------------------------------------------------------------

// Set the source vector pointer
//------------------------------------------------------------------
  Vector *src = in;
  IFloat *src_tmp = (IFloat *)src;

// Set the node checkerboard size of the fermion field
//------------------------------------------------------------------

  if(lat.Fclass() == F_CLASS_CLOVER) {
    f_size_cb = GJP.VolNodeSites() * lat.FsiteSize() / 2;
  } else {
    f_size_cb = GJP.VolNodeSites() * lat.FsiteSize() / (lat.FchkbEvl()+1);
  }


  const int n_fields =  GJP.SnodeSites();  //   *nk ; 
  const int f_size_per_site = lat.FsiteSize() / GJP.SnodeSites()  / (lat.FchkbEvl()+1);


  char fname_eig_root_bc[1024];
  snprintf(fname_eig_root_bc,1024,"%s.bc%d%d%d%d",fname_eig_root,
	   GJP.Bc(0),GJP.Bc(1),GJP.Bc(2),GJP.Bc(3));
	   


  //----------------------------------------------------------------------------------------------//
  // Dear Chulwoo :
  //  Here we allocate a new element of EigenCacheList
  //  In principle, it would be better to deallocate in this class, not to forget to do so.
  //
  //  But the cache list needes to live after the dirac operator destory, as the same eigen vector
  //  would be needed again (at least 12 times for the propagator, and more for other source locations)
  //  So cleanup  is  on the hand of user in either algorithm or in main.C. 
  //  For example, user could Cleanup(EigenCacheList) when one is done with the current configuration.
  //  Of course, forgetting to cleanup will cause a potentially huge memory leak... may force with you.
  //------------------------------------------------------------------------------------------------//
  

#define EIGEN_CACHE_IN_INV_LOWMODE_APPROX_YES

#ifdef EIGEN_CACHE_IN_INV_LOWMODE_APPROX_YES
  // search for eigen cache
  EigenCache* ecache;
  if( (ecache = EigenCacheListSearch( fname_eig_root_bc, neig ) ) == 0 ){
    ecache = new EigenCache();  

    EigenCacheList. push_back( ecache );
  }
#else
  EigenCache* ecache=0;
#endif
  
  Float mass = dirac_arg->mass;
  EigenContainer eigcon( lat, fname_eig_root_bc, neig, f_size_per_site, n_fields, ecache );

  Float* eval = eigcon.load_eval();
  //
  // make eigen values squared if needed
  //
  switch (dirac_arg->RitzMatOper) {
  case  MAT_HERM :
  case  MATPC_HERM : 
    for(int i=0;i<neig;++i)  eval[i]= eval[i]*eval[i];
    break;
  case  MATPCDAG_MATPC  :
  case  NEG_MATPCDAG_MATPC :
  case  MATDAG_MAT :
  case  NEG_MATDAG_MAT :
  case  MATDAG_MAT_NORM :
  case  NEG_MATDAG_MAT_NORM :
  case  MATPCDAG_MATPC_SHIFT :
    break;
  default:
    ERR.General(cname,fname,"RiztMatType %d not implemented\n",
                dirac_arg->RitzMatOper);
  }

  FILE *fp=Fopen("LowmodeSrcProj.dat", "o");

  for(int iev=0;iev<neig;++iev) {
    
    //time_elapse();
    Vector* evec = eigcon.nev_load( iev );
    //print_time("inv_lowmode_approx","loading", time_elapse());

#ifndef USE_BLAS
    Complex z = evec->CompDotProductGlbSum( src, f_size_cb );
#else
    Complex z;

#if TARGET != BGL
    cblas_zdotc_sub( f_size_cb / 2,
		     (double*)evec, 1, (double*)src, 1, (double*)&z);
#else
    *(complex<double>*)&z=cblas_zdotc( f_size_cb / 2,
		     (complex<double>*)evec, 1, (complex<double>*)src, 1);
#endif
    
    glb_sum((double*)&z);
    glb_sum((double*)&z+1);
#endif
    
    z /= eval[ iev ];

    if(!UniqueID())fprintf(fp,"%d %.14e %.14e\n",iev,z.real(),z.imag());

  }
  fclose(fp);
  return;

}

CPS_END_NAMESPACE
