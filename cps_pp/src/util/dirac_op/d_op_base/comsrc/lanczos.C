#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <util/dirac_op.h>
#include <util/lattice.h>
#include <util/smalloc.h>
#include <util/vector.h>
#include <util/gjp.h>
#include <comms/glb.h>
#include <util/verbose.h>
#include <util/error.h>
#include <math.h>
#include <util/qcdio.h>
//#include <alg/matrixpolynomial_arg.h>
#include <util/time_cps.h>

//#define TEST_MAT_YES
#define test_SIZE 50



// FIXME
//       Needs clean up
//    There is  DWF 4D precondition specific piece here
//    move this into DWF specific DiracOp  or Lattice ... ?
//

CPS_START_NAMESPACE

int tqli(float *alpha, float *beta, int m, float **q);
int tqli(Float *alpha, Float *beta, int m, Float **Q, int p);
int tqli(Float *alpha, Float *beta, int m, Float **Q, int p, Float *shifts);
int tqri(Float *alpha, Float *beta, int m, int mm, Float **Q);
int tqri(Float *alpha, Float *beta, Float* s, int m, Float **Q, int p);
int tqri2(Float *alpha, Float *beta, Float* s, int m, Float **Q, int p);
int tqri3(Float *alpha, Float *beta, Float* s, int m, Float **Q, int p);
void eigsrt(Float *d, Float **v, int n);
void eigsrt(Float *d, Float **v1, int n, Vector **v2, int n2);
void eigsrt(Float *d, Float *e, Float **v1, int n, Vector **v2, int n2);
void eigsrt(Float *d, Float *e, Float *tte, Float **v1, int n, Vector **v2, int n2);
void eigsrt_low(Float *d, Float **v1, int n, Vector **v2, int n2);
void eigsrt_low(Float *d,int n);
void eigsrt(Float *d,int n);
void testMat(Vector *Apsi, Vector *psi, int f_size);
void QRtrf(Float *d, Float *e, int nk, int n, Float **z, Float dsh, int kmin, int kmax);

void lanczos_GramSchm(Float *psi, Float **vec, int Nvec, int f_size, Float* alpha);
void lanczos_GramSchm_real(Float *psi, Float **vec, int Nvec, int f_size, Float* alpha);
void lanczos_GramSchm_test(Float *psi, Float **vec, int Nvec, int f_size, Float* alpha);



//#define PROFILE_LANCZOS
#ifdef  PROFILE_LANCZOS
#undef PROFILE_LANCZOS
#define PROFILE_LANCZOS(msg,a ...) do		\
    { if(!UniqueID())				\
	print_asctime(msg, ##a);		\
    }  while(0);

#else
#define time_elapse() 0
#define PROFILE_LANCZOS(msg,a ...) {}
#endif


//#define USE_BLAS
#ifndef USE_BLAS
#define MOVE_FLOAT( pa, pb, n )  moveFloat(pa, pb, n)
#define VEC_TIMESEQU_FLOAT(py, fact, n ) vecTimesEquFloat( py, fact, n)
#define AXPY(n, fact, px, py)  fTimesV1PlusV2(py, fact, px, py, n)
#define glb_DDOT(n, px, py, p_dot) { *(p_dot)=dotProduct(px,py,n); glb_sum((p_dot)); }
#else
#include <util/qblas_extend.h>
#define MOVE_FLOAT( pa, pb, n )  cblas_dcopy(n, pb, 1, pa, 1)
#define VEC_TIMESEQU_FLOAT(py, fact, n ) cblas_dscal( n,  fact, py,1 )
#define AXPY(n, fact, px, py)  cblas_daxpy(n, fact, px,1,py,1)
#define glb_DDOT(n, px, py, p_dot) { *(p_dot) = cblas_ddot(n,px,py); glb_sum((p_dot)); }
#endif


// Implicitly Restarted Lanczos Method Lehoucq and Sorensen (some SIAM Book)
// iterate Lanczos factorization and QL implicit factorization to diagonalize 
// a Hermitian dirac operator
//
//In this routine two different RitzMat will be used :
// 1.  RitzMat_lanczos, which is used in the kernel of lanczos process,
//      e.g. MATPCDAG_MATPC
// 2.  RitzMat_convcheck, which is used for check of convergence,
//     e.g. MATPC_HERM
//
// two matrixces has to be commutable, so that they shares
// the same set of eigenvectors.
//
// RitzMat_convcheck will be used for the eigenvalue computation
//  i.e. if RitzMat_convcheck==MATPC_HERM, then the eigenvalue file will be
//  MATPC_HERM's eigenvalue not the squared one.

//

int DiracOp::ImpResLanczos(Vector **V, //Lanczos vectors, eigenvectors of RitzMat on return
			   Float *alpha, // eigenvalues
			   LanczosArg* eig_arg
			   ) { 
  
  char *fname = "ImResLanczos(...)";
  
  VRB.Func(cname,fname);
  
  RitzMatType RitzMat_lanczos = eig_arg->RitzMat_lanczos;
  RitzMatType RitzMat_convcheck= eig_arg->RitzMat_convcheck;
  
  RitzMatType save_RitzMatOper =  dirac_arg->RitzMatOper;
  
  // set the RitzMatOper to the one for the convergence check
  dirac_arg->RitzMatOper = RitzMat_convcheck;
  
  int nk = eig_arg-> nk_lanczos_vectors; // number of wanted eigenvectors
  int np  = eig_arg-> np_lanczos_vectors; // extra, unwanted eigenvectors
  int MaxIters  = eig_arg-> maxiters; // maximum number of restarting
  Float StopRes  = eig_arg-> stop_residual;// stop when residual is smaller than this
    
  //arguments for the filtering polynomial
//  MatrixPolynomialArg *cheby_arg  = (MatrixPolynomialArg*) (eig_arg-> matpoly_arg);
  MatrixPolynomialArg *cheby_arg  = &(eig_arg-> matpoly_arg);


  // Set the node checkerboard size of the fermion field
  //------------------------------------------------------------------
  int f_size = RitzLatSize();
  
  // implicitly restarted lanczos
  int m = nk+np;
  // current residual vector
  Vector *r = (Vector *) smalloc(cname,fname, "r", f_size * sizeof(Float));
  // scratch vector
  Vector *Apsi = (Vector *) smalloc(cname,fname, "Apsi", f_size * sizeof(Float));
  // off-diagonal for the tridiagonal Lanczos matrix Tm
  Float *beta = (Float *) smalloc(cname,fname, "beta", m * sizeof(Float));
  Float *ttbeta = (Float *) smalloc(cname,fname, "ttbeta", m * sizeof(Float));
  Float *shifts = (Float *) smalloc(cname,fname, "shifts", m * sizeof(Float));

  // orthoganal matrix for similarity transform of Tm in QL factorization
  Float **Q = (Float **) smalloc(cname,fname, "Q", (m)* sizeof(Float*));
  Float* tmp_Qi = (Float *) smalloc(cname,fname,"tmp_Qi", (m*m) * sizeof(Float));
  for(int i=0;i<m;i++)  Q[i] = tmp_Qi+i*m;

  Float *scratch_2m = (Float*) smalloc(cname,fname,"scratch_m", m*2*sizeof(Float));

  // setup temporaly vectors for the the matrix polynomial
  cheby_arg->tmp1 = (Pointer) smalloc(cname,fname,"matrix_polynomial.tmp1", f_size *sizeof(Float));
  cheby_arg->tmp2 = (Pointer) smalloc(cname,fname,"matrix_polynomial.tmp2", f_size *sizeof(Float));  

  // save the eigenvalues so we don't have to computed them after last iter.
  //Float *savelambda = (Float *) smalloc(cname,fname, "savelambda", m * sizeof(Float));


  // initialize starting residual

#if 0

  for(int n=0;n<f_size;n+=2){
    //LRG.AssignGenerator(site);
#if 0
    *((Float*)r+n) = LRG.Grand();
    *((Float*)r+1+n) = LRG.Grand();
#else
    *((Float*)r+n) = 1.0;
    *((Float*)r+1+n) = 0.0;
#endif
  }
  
#else

  printf("norm v[0] = %g\n", V[0]->NormSqGlbSum(f_size));
  
  // use V[0] as the initial vector
  moveFloat( (Float*)r, (Float*)(V[0]), f_size);

#endif
  
  
  
#if 0
  // Initial filtering
  MatrixPolynomialArg cheby0_arg;  //arguments for the filtering polynomial
  
  cheby0_arg. Npol = 300 ; 
  cheby0_arg. params[0] = 2.40;
  cheby0_arg. params[1] = 0.015; // 0.030 is the 100th ?
  cheby0_arg. tmp1 = cheby_arg-> tmp1;
  cheby0_arg. tmp2 = cheby_arg-> tmp2;
  
  
  VRB.Result(cname,fname,"Initial Filter with Polynomial degree=%d\n",cheby0_arg.Npol);
  
  RitzMat(Apsi, r, &cheby0_arg);
  
  r->CopyVec(Apsi,f_size);
  
#endif
  
  //----------------------------------------------------
  
  
  // set eig vals to zero
  for(int j=0;j<m;j++) alpha[j]=0.0;
  
  //PROFILE_LANCZOS("dummy",time_elapse() );
  double time_in=dclock();
  
  // initial m-step lanczos
  lanczos(0, m, f_size, Apsi, r, V, alpha, beta, cheby_arg, RitzMat_lanczos);
  
  //PROFILE_LANCZOS("initial lanczos %e\n",time_elapse());
  if(!UniqueID()) printf("Initial Lanczos time(sec) %e\n", dclock()-time_in);
  
  // do until converged or hit maximum iterations
  int it=1;
  while(it<=MaxIters){
    
    //PROFILE_LANCZOS("dummy",time_elapse());
    double time_in=dclock();
    // get shifts
    moveFloat(shifts,alpha, m);
    moveFloat(ttbeta,beta, m);
    vecZero((IFloat*)Q[0],m*m);
    for(int j=0;j<m;j++) Q[j][j] = 1.0;
    
    int tqri_iters=tqri(shifts, ttbeta, m, m, Q);
    VRB.Result(cname,fname, "tqri converges at iter=%d\n",tqri_iters);
    
    eigsrt(shifts,m);
    //eigsrt_low(shifts,m);
    
    if( VRB.IsActivated(VERBOSE_DEBUG_LEVEL))
      for(int i=0; i< m; ++i)
	VRB.Debug(cname,fname, "shifts %d %.16e\n", i, shifts[i]);
    
    
    //np steps of QR factorization. 
    vecZero((IFloat*)Q[0],m*m);
    for(int j=0;j<m;j++) Q[j][j] = 1.0;
    
    
    if( VRB.IsActivated(VERBOSE_DEBUG_LEVEL))
      for(int j=0;j<m;j++) 
	VRB.Debug(cname,fname, "TDMAT Before QR alpha beta = : %d %d %.16e  %.16e\n",it,j,alpha[j],beta[j]);
    
    PROFILE_LANCZOS("vec zero etc %e\n",time_elapse());
    for(int i=0;i<np;i++){
      QRtrf(alpha,beta,m,m,Q,shifts[nk+i],0,m-1);
    }
    PROFILE_LANCZOS("QRtrf %e\n",time_elapse());
    
    //if( VRB.IsActivated(VERBOSE_DEBUG_LEVEL))
    for(int j=0;j<m;j++) 
      VRB.Result(cname,fname, "TDMAT After QR alpha beta = : %d %d %.16e  %.16e\n",it,j,alpha[j],beta[j]);
    
    //if(beta[nk-1] < 1.0e-150)break;
    
    
    
    // new Vk, without large scratch array, 
    // and without computing zero matrix elements in Q :
    //   Q[j][k]==0  for  j>np , 0<= k < j-np
    // e.g. for  nk=10, np=12
    //      Q[13][0]=0, Q[14][0]=Q[14][1]=0
    // FIXME:  would be slow due to non-localized memory access
    for(int n=0;n<f_size;n+=2){ 
      vecZero( (IFloat*)scratch_2m, 2*m);
      
      for(int j=0; j<np;++j){
	for(int k=0;k<m;++k){ 
	  scratch_2m[2*k]   += *((Float*)V[j]+n)   * Q[j][k];
	  scratch_2m[2*k+1] += *((Float*)V[j]+n+1) * Q[j][k];
	}
      }
      for(int j=np; j<m;++j){
	for(int k=j-np;k<m;++k){ 
	  scratch_2m[2*k]   += *((Float*)V[j]+n)   * Q[j][k];
	  scratch_2m[2*k+1] += *((Float*)V[j]+n+1) * Q[j][k];
	}
      }
      for(int k=0;k<m;++k) {
	*((Float*)V[k]+n) = scratch_2m[2*k];
	*((Float*)V[k]+n+1) = scratch_2m[2*k+1];
      }
    }

#if 0 
    for(int k=0;k<m;++k)
      for(int j=0;j<m;++j)
	printf("Q %d %d = %g\n",k,j,Q[k][j]);
#endif
    PROFILE_LANCZOS("diag Q %e\n",time_elapse());
    
    //
    // restart Lanczos with new residual r_k
    // first multiply r_m by Q(m,k)
    //
    if( VRB.IsActivated(VERBOSE_DEBUG_LEVEL)){
      VRB.Debug(cname,fname, "RESTVEC %e %e\n", Q[m-1][nk-1], beta[nk-1]);
      {Float* fp = (Float*)r;
	VRB.Debug(cname,fname,"RESTVEC %e %e %e %e\n", fp[0],fp[1],fp[2],fp[3]);
      }
      {Float* fp = (Float*)V[nk];
	VRB.Debug(cname,fname,"RESTVEC %e %e %e %e\n", fp[0],fp[1],fp[2],fp[3]);
      }
    }
    
    r->VecTimesEquFloat(Q[m-1][nk-1],f_size);
    if(nk<m){
      Apsi->CopyVec(V[nk],f_size);
      Apsi->VecTimesEquFloat(beta[nk-1], f_size);
      vecAddEquVec((Float*)r, (Float*)Apsi, f_size);
    }
    beta[nk-1]=sqrt(r->NormSqGlbSum(f_size));
    
    
    VRB.Debug(cname,fname, "RESTVEC new beta_k %e\n", beta[nk-1]);
    //if(beta[nk-1] < 1.0e-150)break;
    
    if( VRB.IsActivated(VERBOSE_DEBUG_LEVEL)){
      VRB.Debug(cname,fname,"new residual vector\n");
      for(int j=0;j<f_size;j++) 
	VRB.Debug(cname,fname, "%d %.16f\n",j,*((Float*)r+j));
    }
    
    PROFILE_LANCZOS("prep new resvec Q %e\n",time_elapse());

    if(it%eig_arg->conv_check==0){

      PROFILE_LANCZOS("before conv check %e\n",time_elapse());
      //
      // Convergence check 
      //   by monitoring the eigenV eq. for original operator
      //
      int Ndone(0); // Number of converged vector
      double time_conv=dclock();
      
      // get shifts
      moveFloat(shifts,alpha, m);
      moveFloat(ttbeta,beta, m);
      vecZero((IFloat*)Q[0],m*m);
      for(int j=0;j<m;j++) Q[j][j] = 1.0;
      
      //int exam_dim = m;   // full dimension, sometimes better ? 
      int exam_dim = nk; //only the wanted dimension
      
      tqri_iters=tqri(shifts, ttbeta, exam_dim, m, Q);
      VRB.Debug(cname,fname,  "tqri (convergence test) converges at iter=%d\n",tqri_iters);
      
      Vector* vtmp = (Vector*)(cheby_arg->tmp1); // reuse temp vector for polynomial here
      
      
      PROFILE_LANCZOS("before residual comp. %e\n",time_elapse());
      
      for(int k=0;k<nk;k++){
	
	vtmp->VecZero(f_size);
	for(int j=0;j<exam_dim;j++)
	  vtmp->FTimesV1PlusV2( Q[j][k], V[j], vtmp, f_size);  
	//moveFloat((Float*)vtmp,(Float*)V[k],f_size);
	
	// Monitor convergence for the original RitzMat (not the polynomial of it)
	RitzMat(Apsi, vtmp);
	
	// eigenvalue obtained from Rayleigh quotient
	Float bt=sqrt(vtmp->NormSqGlbSum(f_size));
	Float alp =  Apsi->ReDotProductGlbSum( vtmp, f_size) / bt;
	
	// calculate residue for eigenequation
	Apsi->FTimesV1PlusV2(-alp, vtmp, Apsi, f_size);  
	Float rnorm = sqrt(Apsi->NormSqGlbSum(f_size)) / bt;
	//Float rnorm = Apsi->NormSqGlbSum(f_size);
	if( rnorm < StopRes ) Ndone++;
	VRB.Result(cname,fname, "Residual %d %d %e lambda %e %e %e\n",it, k, 
		   rnorm, alp, shifts[k], alpha[k]);
	if(isnan(rnorm)) 
          ERR.General(cname,fname, "Residual is nan\n"); 
	// save eigenvalues/vectors so don't have to recompute in final conv. check
	if(it==MaxIters)alpha[k] = alp;
	//if(it==MaxIters)moveFloat((Float*)V[k],(Float*)vtmp,f_size);
      }
      PROFILE_LANCZOS("after residual comp %e\n",time_elapse());
      
      if(!UniqueID()) printf("# %d Lanczos: time of conv check(sec) %e\n",it,dclock()-time_conv);
      if(Ndone == nk) break;      
    }

    // np more steps of lanczos
    lanczos(nk, m, f_size, Apsi, r, V, alpha, beta, cheby_arg, RitzMat_lanczos);
    
    //PROFILE_LANCZOS("lanczos %e\n",time_elapse());
    if(!UniqueID()) printf("# %d Lanczos: time(sec) %e\n",it,dclock()-time_in);
    
    it++;
  }
  if(it==MaxIters+1)it--;
  //-------------------------------------------------------------------------

  PROFILE_LANCZOS("before the final  conv check %e\n",time_elapse());
  double time_conv=dclock();


  // Final diagonalization for V[j]
  // FIXEME: this is already done in the convergence check 
  // in the iteration.
  
  // get shifts
  moveFloat(shifts,alpha, m);
  moveFloat(ttbeta,beta, m);
  vecZero((IFloat*)Q[0],m*m);
  for(int j=0;j<m;j++) Q[j][j] = 1.0;
  
  //int exam_dim = m;   // full dimension, sometimes better ? 
  int exam_dim = nk; //only the wanted dimension
  
  int tqri_iters=tqri(shifts, ttbeta, exam_dim, m, Q);
  VRB.Debug(cname,fname,  "tqri (convergence test) converges at iter=%d\n",tqri_iters);
  
  PROFILE_LANCZOS("before new V %e\n",time_elapse());
  // new Vk, without large scratch array
  // FIXME:  would be slow due to non-localized memory access
#if 1
    int* tab;
    tab = (int*) smalloc(m*m*sizeof(int));
    int cnt=0;
    for(int k=0;k<m;++k){
      for(int j=0;j<m;++j){
	if(fabs(Q[j][k]) > 1e-16){
	  if(!UniqueID())printf("Qfinal %d %d = %g\n",j,k,Q[j][k]);
	  tab[cnt] = j+m*k;
	  cnt++;
	}
      }
    }
#endif
  for(int n=0;n<f_size;n+=2){ 
    vecZero( (IFloat*)scratch_2m, 2*m);
#if 0
    for(int k=0;k<m;++k){ 
      for(int j=0;j<m;++j){
	scratch_2m[2*k]   += *((Float*)V[j]+n)   * Q[j][k];
	scratch_2m[2*k+1] += *((Float*)V[j]+n+1) * Q[j][k];
      }
    }
#endif
    for(int i=0;i<cnt;++i){
      int j = tab[i]%m;
      int k = (tab[i]-j)/m;
      scratch_2m[2*k]   += *((Float*)V[j]+n)   * Q[j][k];
      scratch_2m[2*k+1] += *((Float*)V[j]+n+1) * Q[j][k];      
    }
    for(int k=0;k<m;++k) {
      *((Float*)V[k]+n) = scratch_2m[2*k];
      *((Float*)V[k]+n+1) = scratch_2m[2*k+1];
    }
  }
  
  PROFILE_LANCZOS("after new V %e\n",time_elapse());
  
  FILE* filep=0;
  if (eig_arg->results!=0){
    // Final matrix 
    filep=Fopen(eig_arg->results,"a");
    
    Fprintf( filep, "mass = %e\n",eig_arg->mass );
    Fprintf( filep, "nk=%d, np=%d\n", nk, np);
    Fprintf( filep, "number of restarting = %d\n",it);
  }
  PROFILE_LANCZOS("after final output %e\n",time_elapse());
  
  Complex *del = (Complex*)smalloc(cname,fname,"del", nk*nk*sizeof(Complex));
  
  for(int i=0;i<nk;++i){
    RitzMat(Apsi, V[i]);
    
    for(int j=0;j<nk;++j){ // can cut into half in the exact arithmetic
      del[i*nk+j] = Apsi->CompDotProductGlbSum(V[j],f_size);  // Apsi . V[i]
    }
    //Float alp =  Apsi->ReDotProductGlbSum( V[i], f_size);
    Float alp =  del[i*nk+i].real();
    alpha[i]=alp;
    Apsi->FTimesV1PlusV2(-alp, V[i], Apsi, f_size);  
    Float rnorm = sqrt(Apsi->NormSqGlbSum(f_size));
    Float norm = sqrt(V[i]->NormSqGlbSum(f_size));
    if( filep )
      Fprintf(filep, "Final True Residual %d %d %e norm %e lambda %.16e\n",it,i,rnorm,norm,alp);
    VRB.Result(cname,fname, "Final True Residual %d %d %e norm %e lambda %.16e\n",it,i,rnorm,norm,alp);
  }
  PROFILE_LANCZOS("after final residual comp. %e\n",time_elapse());
  
  if( filep ){
    Fprintf(filep, "\n");
    Fprintf(filep, "final matrix\n");
    for(int i=0;i<nk;++i)
      for(int j=0;j<nk;++j){
	Fprintf(filep, "%d %d %e %e\n",
		i, j, del[i*nk+j].real(), del[i*nk+j].imag());
      }
    Fclose(filep); // close file
  }
  sfree(cname,fname, "del",del);
  
  if(!UniqueID()) printf("# %d Lanczos: time of final conv check(sec) %e\n",it,dclock()-time_conv);
  
  PROFILE_LANCZOS("after the final  conv check %e\n",time_elapse());
  
#if 0
  for(int j=0;j<2;j++){
    Float* fp = (Float*)(V[j]);
    // phase convention
    Float phase = atan2( fp[1], fp[0] );
    Rcomplex rot_phase( cos( phase ), -sin( phase ) );
    
    for(int i=0;i<f_size;i+=2){
      Rcomplex C(fp[i], fp[i+1]);
      C *= rot_phase;
      printf("DEB %d %d %.16e %.16e %.16e\n", j,i, C.real(),C.imag(),C.norm());
    }
  }
#endif
  sfree(cname,fname, "tab", tab);
  sfree(cname,fname, "scratch_2m", scratch_2m);
  sfree(cname,fname, "tmp_Qi", tmp_Qi);
  sfree(cname,fname, "Q", Q);
  sfree(cname,fname, "shifts",shifts);
  sfree(cname,fname, "ttbeta",ttbeta);
  sfree(cname,fname, "beta",beta);
  sfree(cname,fname, "Apsi",Apsi);
  sfree(cname,fname, "r",r);
  sfree(cname,fname,"matrix_polynomial.tmp1", (Float *)cheby_arg->tmp1);
  sfree(cname,fname,"matrix_polynomial.tmp2", (Float *)cheby_arg->tmp2);

  // make sure we won't change the original RitMatOper.
  dirac_arg->RitzMatOper = save_RitzMatOper;

  return(it);
}

/* get m Lanczos vectors */
/*
k is the last Vector in an original k-step Lanczos factorization
m is the number of Lanczos vectors on return
f_size is the size of the Lanczos(eigen)vectors
psi is the kth Lanczos vector
Apsi is a scratch vector
r is the residual vector (on entry it is the kth lanczos vector)
V is the array of m Lanczos vectors
alpha: diagonal elements of the Lanczos tridiagonal matrix
beta: off-diagonal matrix elements of said tridiagonal matrix  
 */

void DiracOp::lanczos(int k0, int m, int f_size,
		      Vector *Apsi, Vector *r,
		      Vector **V, Float *alpha, Float* beta,
		      MatrixPolynomialArg* cheby_arg,
		      RitzMatType RitzMat_lanczos )
{

  RitzMatType save_RitzMatOper =  dirac_arg->RitzMatOper;
  dirac_arg->RitzMatOper = RitzMat_lanczos;
  
  //  if(k0==0)  
  {
    Float ff;
    ff=sqrt(r->NormSqGlbSum(f_size));
    //V[k0]->CopyVec(r, f_size);
    //V[k0]->VecTimesEquFloat(1.0/ff, f_size);
    MOVE_FLOAT( (Float*)(V[k0]), (Float*)r, f_size );
    VEC_TIMESEQU_FLOAT((Float*)(V[k0]), 1.0/ff, f_size );
  }
  
  for(int k = k0; k < m; ++k){
    
    // --- the first step ---
    if( k == 0) {

      double time_in=dclock();
      PROFILE_LANCZOS("before first mult %e\n",time_elapse());    
      RitzMat(r, V[k], cheby_arg);
      PROFILE_LANCZOS("first mult %e\n",time_elapse());    

      //alpha[k] = V[k]->ReDotProductGlbSum(r, f_size);
      glb_DDOT( f_size, (IFloat*)(V[k]), (IFloat*)r, alpha+k);

      //r->FTimesV1PlusV2(-alpha[k], V[k], r, f_size);
      AXPY(f_size, -alpha[k], (Float*)(V[k]), (Float*)r );
      
      //beta[k]=sqrt(r->NormSqGlbSum(f_size));
      glb_DDOT(f_size, (Float*)r, (Float*)r, beta+k);
      beta[k]=sqrt(beta[k]);

      

      //V[k+1]->CopyVec(r, f_size);
      MOVE_FLOAT((Float*)(V[k+1]), (Float*)r, f_size);
      //V[k+1]->VecTimesEquFloat(1.0/beta[k], f_size);
      VEC_TIMESEQU_FLOAT( (Float*)(V[k+1]), 1.0/beta[k], f_size);
      PROFILE_LANCZOS("first other linalg %e\n",time_elapse());

      if(!UniqueID()) printf("First Lanczos iter: time(sec) %e\n", dclock()-time_in);

    } else {
   // --- iteration step ---

      double time_in=dclock();
      PROFILE_LANCZOS("before mult %e\n",time_elapse());    
      RitzMat(r, V[k], cheby_arg);
      PROFILE_LANCZOS("mult %e\n",time_elapse());    

#if 0
      r->FTimesV1PlusV2(-beta[k-1], V[k-1], r, f_size);
      alpha[k] = V[k]->ReDotProductGlbSum(r, f_size); 
      r->FTimesV1PlusV2(-alpha[k], V[k], r, f_size);      
      beta[k]=sqrt(r->NormSqGlbSum(f_size));
#else
      AXPY(f_size, -beta[k-1], (Float*)V[k-1], (Float*)r);
      //alpha[k] = V[k]->ReDotProductGlbSum(r, f_size);
      glb_DDOT( f_size, (Float*)(V[k]), (Float*)r,alpha+k );
      AXPY(f_size, -alpha[k], (Float*)V[k], (Float*)r);
      //beta[k]=sqrt(r->NormSqGlbSum(f_size));
      glb_DDOT( f_size, (Float*)r, (Float*)r, beta+k);
      beta[k]=sqrt(beta[k]);
#endif

      PROFILE_LANCZOS("other lin alg %e\n",time_elapse());    

      //lanczos_GramSchm_real((Float*)r, (Float**)V, k+1, f_size, 0);
      lanczos_GramSchm_test((Float*)r, (Float**)V, k+1, f_size, 0);
      PROFILE_LANCZOS("gramschm %e\n",time_elapse());
      
      if( k+1 < m ){
	//V[k+1]->CopyVec(r, f_size);
	MOVE_FLOAT( (Float*) (V[k+1]), (Float*)r, f_size);
	//V[k+1]->VecTimesEquFloat(1.0/beta[k], f_size);
	VEC_TIMESEQU_FLOAT( (Float*)(V[k+1]), 1.0/beta[k], f_size);
      }
      PROFILE_LANCZOS("last linalg %e\n",time_elapse());
      if(!UniqueID()) printf("# %d Lanczos iter: time(sec) %e\n",k, dclock()-time_in);
    }
  }

  dirac_arg->RitzMatOper = save_RitzMatOper;
}


CPS_END_NAMESPACE
