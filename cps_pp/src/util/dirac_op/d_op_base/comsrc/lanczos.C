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
#ifdef USE_LAPACK
#include <lapacke.h>
#endif
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

CPS_START_NAMESPACE int tqli (float *alpha, float *beta, int m, float **q);
int tqli (Float * alpha, Float * beta, int m, Float ** Q, int p);
int tqli (Float * alpha, Float * beta, int m, Float ** Q, int p,
	  Float * shifts);
int tqri (Float * alpha, Float * beta, int m, int mm, Float ** Q);
int tqri (Float * alpha, Float * beta, Float * s, int m, Float ** Q, int p);
int tqri2 (Float * alpha, Float * beta, Float * s, int m, Float ** Q, int p);
int tqri3 (Float * alpha, Float * beta, Float * s, int m, Float ** Q, int p);
void eigsrt (Float * d, Float ** v, int n);
void eigsrt (Float * d, Float ** v1, int n, Vector ** v2, int n2);
void eigsrt (Float * d, Float * e, Float ** v1, int n, Vector ** v2, int n2);
void eigsrt (Float * d, Float * e, Float * tte, Float ** v1, int n,
	     Vector ** v2, int n2);
void eigsrt_low (Float * d, Float ** v1, int n, Vector ** v2, int n2);
void eigsrt_low (Float * d, int n);
void eigsrt (Float * d, int n);
void testMat (Vector * Apsi, Vector * psi, size_t f_size);
void QRtrf (Float * d, Float * e, int nk, int n, Float ** z, Float dsh,
	    int kmin, int kmax);

void lanczos_GramSchm (Float * psi, Float ** vec, int Nvec, size_t f_size,
		       Float * alpha);
void lanczos_GramSchm_real (Float * psi, Float ** vec, int Nvec, size_t f_size,
			    Float * alpha);
void lanczos_GramSchm_test (Float * psi, Float ** vec, int Nvec, size_t f_size,
			    Float * alpha);
void lanczos_GramSchm_test (Float * psi, float **vec, int Nvec, size_t f_size,
			    Float * alpha);
void double_to_float (Float * vec, size_t f_size);
void movefloattoFloat (Float * out, float *in, size_t f_size);
void moveFloattofloat (float *out, Float * in, size_t f_size);

#ifdef USE_LAPACK
void tqri_LAPACK (double *shifts, double *ttbeta, int NN, int m, double *Q)
{
//  int NN = m;
  double evals_tmp[NN];
  memset (evals_tmp, 0, sizeof (double) * NN);
  double evec_tmp[NN * NN];
  memset (evec_tmp, 0, sizeof (double) * NN * NN);
//      double AA[NN][NN];
  double DD[NN];
  double EE[NN];
  for (int i = 0; i < NN; i++)
    for (int j = i - 1; j <= i + 1; j++)
      if (j < NN && j >= 0) {
//                AA[i][j] = AH(i,j);
	if (i == j)
	  DD[i] = shifts[i];
//          if (i == j)
//            evals_tmp[i] = alpha[i];
	if (j == (i - 1))
	  EE[j] = ttbeta[j];
//              if (i<20 && j<20) 
//                              QDPIO:: cout << "AA["<<i<<"]["<<j<<"]="<<AA[i][j]<<endl;
      }
  int evals_found;
  int lwork = (18*NN);
//    ((18 * NN) > (1 + 4 * NN + NN * NN) ? (18 * NN) : (1 + 4 * NN + NN * NN));
  int liwork = 3 + NN * 10;
  int iwork[liwork];
  double work[lwork];
  int isuppz[2 * NN];
  char jobz = 'V';		// calculate evals & evecs
  char range = 'I';		// calculate all evals
  char uplo = 'U';		// refer to upper half of original matrix
  char compz = 'I';		// Compute eigenvectors of tridiagonal matrix
  int ifail[NN];
  int info;
  int total = QMP_get_number_of_nodes ();
  int node = QMP_get_node_number ();
  int interval = (NN / total) + 1;
  double vl = 0.0, vu = 0.0;
  int il = interval * node + 1, iu = interval * (node + 1);
  if (iu > NN)
    iu = NN;
  double tol = 0.0;
  int if_print= VRB.IsActivated(VERBOSE_DEBUG_LEVEL);
  if (il <= NN) {
    if (if_print) printf ("total=%d node=%d il=%d iu=%d\n", total, node, il, iu);
    LAPACK_dstegr (&jobz, &range, &NN, (double *) DD, (double *) EE, &vl, &vu, &il, &iu,	// these four are ignored if second parameteris 'A'
		   &tol,	// tolerance
		   &evals_found, evals_tmp, evec_tmp, &NN,
		   isuppz, work, &lwork, iwork, &liwork, &info);
    for (int i = iu - 1; i >= il - 1; i--) {
     if (if_print)  printf ("node=%d evals_found=%d evals_tmp[%d]  =%g evec_tmp[0][%d]= %g %g\n", node,
	      evals_found, i - (il - 1), evals_tmp[i - (il - 1)],
	       i - (il - 1), evec_tmp[(i - (il - 1))],evec_tmp[i-(il-1)+1]);
      evals_tmp[i] = evals_tmp[i - (il - 1)];
      if (il > 1)
	evals_tmp[i - (il - 1)] = 0.;
      for (int j = 0; j < NN; j++) {
	evec_tmp[i*NN  + j] = evec_tmp[(i - (il - 1))*NN + j];
	if (il > 1)
	  evec_tmp[(i - (il - 1))*NN + j] = 0.;
      }
    }
  }
  {
//      TIMER ("QMP_sum_double_array");
    QMP_sum_double_array (evals_tmp, NN);
    QMP_sum_double_array ((double *) evec_tmp, NN * NN);
  }
  for (int i = 0; i < NN; i++) {
    shifts[i] = evals_tmp[i];
    for (int j = 0; j < NN; j++) {
//      VRB.Result ("", "tqri_LAPACK", "evec_tmp[%d][%d]=%g\n", i, j, evec_tmp[i * NN + j]);
      Q[i * m + j] = evec_tmp[j * NN + i];
    }
  }
}
#endif

#undef PROFILE_LANCZOS
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

int DiracOp::ImpResLanczos (Vector ** V,	//Lanczos vectors, eigenvectors of RitzMat on return
			    Float * alpha,	// eigenvalues
			    LanczosArg * eig_arg)
{

  char *fname = "ImResLanczos(...)";

  VRB.Func (cname, fname);

  RitzMatType RitzMat_lanczos = eig_arg->RitzMat_lanczos;
  RitzMatType RitzMat_convcheck = eig_arg->RitzMat_convcheck;

  RitzMatType save_RitzMatOper = dirac_arg->RitzMatOper;

  // set the RitzMatOper to the one for the convergence check
  dirac_arg->RitzMatOper = RitzMat_convcheck;

  int nk = eig_arg->nk_lanczos_vectors;	// number of wanted +extra eigenvectors
  int nt = eig_arg->nt_lanczos_vectors;	// number of wanted eigenvectors
  int np = eig_arg->np_lanczos_vectors;	// extra, unwanted eigenvectors
  int MaxIters = eig_arg->maxiters;	// maximum number of restarting
  Float StopRes = eig_arg->stop_residual;	// stop when residual is smaller than this

  //arguments for the filtering polynomial
  MatrixPolynomialArg *cheby_arg = &eig_arg->matpoly_arg;
//    (MatrixPolynomialArg *) (eig_arg->matpoly_arg);


  // Set the node checkerboard size of the fermion field
  //------------------------------------------------------------------
  size_t f_size = RitzLatSize ();

  // implicitly restarted lanczos
  int m = nk + np;
  // current residual vector
  Vector *r = (Vector *) smalloc (cname, fname, "r", f_size * sizeof (Float));
  // scratch vector
  Vector *Apsi =
    (Vector *) smalloc (cname, fname, "Apsi", f_size * sizeof (Float));
  // off-diagonal for the tridiagonal Lanczos matrix Tm
  Float *beta = (Float *) smalloc (cname, fname, "beta", m * sizeof (Float));
  memset (beta, 0, m * sizeof (Float));
  Float *ttbeta =
    (Float *) smalloc (cname, fname, "ttbeta", m * sizeof (Float));
  Float *shifts =
    (Float *) smalloc (cname, fname, "shifts", m * sizeof (Float));

  // orthoganal matrix for similarity transform of Tm in QL factorization
  Float **Q = (Float **) smalloc (cname, fname, "Q", (m) * sizeof (Float *));
  Float *tmp_Qi =
    (Float *) smalloc (cname, fname, "tmp_Qi", (m * m) * sizeof (Float));
  for (int i = 0; i < m; i++)
    Q[i] = tmp_Qi + i * m;

  //Float *scratch_2m = (Float*) smalloc(cname,fname,"scratch_m", m*2*sizeof(Float));
  Float *scratch_2m =
    (Float *) smalloc (cname, fname, "scratch_m", m * f_size * sizeof (Float));

  // setup temporaly vectors for the the matrix polynomial
  cheby_arg->tmp1 =
    (Float *) smalloc (cname, fname, "matrix_polynomial.tmp1",
		       f_size * sizeof (Float));
  cheby_arg->tmp2 =
    (Float *) smalloc (cname, fname, "matrix_polynomial.tmp2",
		       f_size * sizeof (Float));

  // save the eigenvalues so we don't have to computed them after last iter.
  //Float *savelambda = (Float *) smalloc(cname,fname, "savelambda", m * sizeof(Float));


  // initialize starting residual

#if 0
  int nodes =
    GJP.Nodes (0) * GJP.Nodes (1) * GJP.Nodes (2) * GJP.Nodes (3) *
    GJP.Nodes (4);
  for (int n = 0; n < f_size; n += 2) {
    //LRG.AssignGenerator(site);
#if 0
    *((Float *) r + n) = LRG.Grand ();
    *((Float *) r + 1 + n) = LRG.Grand ();
#else
    *((Float *) r + n) = sqrt (2.0 / (Float) (nodes * f_size));
    *((Float *) r + 1 + n) = 0.0;
#endif
  }

#else

  // use V[0] as the initial vector
  movefloattoFloat ((Float *) r, (float *) V[0], f_size);
#endif
#if 0
  if (!UniqueID ()) {
    printf ("f_zise: %d\n", f_size);
  };
  for (int n = 0; n < f_size; n += 2) {
    printf ("evec0 r %d %d %e %e %e %e\n",
	    UniqueID (), n,
	    *((float *) V[0] + n), *((float *) V[0] + n + 1),
	    *((Float *) r + n), *((Float *) r + n + 1));
  }
#endif
  printf ("norm V[0] = %g\n", r->NormSqGlbSum (f_size));




#if 0
  // Initial filtering
  MatrixPolynomialArg cheby0_arg;	//arguments for the filtering polynomial

  cheby0_arg.Npol = 300;
  cheby0_arg.params[0] = 2.40;
  cheby0_arg.params[1] = 0.015;	// 0.030 is the 100th ?
  cheby0_arg.tmp1 = cheby_arg->tmp1;
  cheby0_arg.tmp2 = cheby_arg->tmp2;


  VRB.Result (cname, fname, "Initial Filter with Polynomial degree=%d\n",
	      cheby0_arg.Npol);

  RitzMat (Apsi, r, &cheby0_arg);

  r->CopyVec (Apsi, f_size);

#endif

  //----------------------------------------------------


  // set eig vals to zero
  for (int j = 0; j < m; j++)
    alpha[j] = 0.0;

  //PROFILE_LANCZOS("dummy",time_elapse() );
  double time_in = dclock ();

  // initial m-step lanczos
  lanczos (0, m, f_size, Apsi, r, V, alpha, beta, cheby_arg, RitzMat_lanczos);

  //PROFILE_LANCZOS("initial lanczos %e\n",time_elapse());
  if (!UniqueID ())
    printf ("Initial Lanczos time(sec) %e\n", dclock () - time_in);

  // do until converged or hit maximum iterations
  int it = 1;
  while (it <= MaxIters) {

    //PROFILE_LANCZOS("dummy",time_elapse());
    double time_in = dclock ();
    // get shifts
    moveFloat (shifts, alpha, m);
    moveFloat (ttbeta, beta, m);
    vecZero ((IFloat *) Q[0], m * m);
    for (int j = 0; j < m; j++)
      Q[j][j] = 1.0;


    int tqri_iters;
#ifndef USE_LAPACK
//#if 1
    tqri_iters = tqri (shifts, ttbeta, m, m, Q);
    VRB.Result (cname, fname, "tqri converges at iter=%d time=%g\n", tqri_iters,
		dclock () - time_in);
#else
    tqri_LAPACK (shifts, ttbeta, m, m, (double *) Q[0]);
    VRB.Result (fname, "tqri_LAPACK", " time(sec) %e\n", dclock () - time_in);
#endif

    eigsrt (shifts, m);
    //eigsrt_low(shifts,m);

    if (VRB.IsActivated (VERBOSE_RESULT_LEVEL))
      for (int i = 0; i < m; ++i)
	VRB.Result (cname, fname, "shifts %d %.16e\n", i, shifts[i]);


    //np steps of QR factorization. 
    vecZero ((IFloat *) Q[0], m * m);
    for (int j = 0; j < m; j++)
      Q[j][j] = 1.0;


    if (VRB.IsActivated (VERBOSE_DEBUG_LEVEL))
      for (int j = 0; j < m; j++)
	VRB.Debug (cname, fname,
		   "TDMAT Before QR alpha beta = : %d %d %.16e  %.16e\n", it, j,
		   alpha[j], beta[j]);

    PROFILE_LANCZOS ("vec zero etc %e\n", time_elapse ());
    for (int i = 0; i < np; i++) {
      QRtrf (alpha, beta, m, m, Q, shifts[nk + i], 0, m - 1);
    }
    PROFILE_LANCZOS ("QRtrf %e\n", time_elapse ());

    //if( VRB.IsActivated(VERBOSE_DEBUG_LEVEL))
    for (int j = 0; j < m; j++)
      VRB.Result (cname, fname,
		  "TDMAT After QR alpha beta = : %d %d %.16e  %.16e\n", it, j,
		  alpha[j], beta[j]);

    //if(beta[nk-1] < 1.0e-150)break;



    // new Vk, without large scratch array, 
    // and without computing zero matrix elements in Q :
    //   Q[j][k]==0  for  j>np , 0<= k < j-np
    // e.g. for  nk=10, np=12
    //      Q[13][0]=0, Q[14][0]=Q[14][1]=0
#if 1
    // FIXME:  would be slow due to non-localized memory access
    for (int n = 0; n < f_size; n += 2) {
      vecZero ((IFloat *) scratch_2m, 2 * m);

      for (int j = 0; j < np; ++j) {
	for (int k = 0; k < m; ++k) {
	  scratch_2m[2 * k] += *((float *) V[j] + n) * Q[j][k];
	  scratch_2m[2 * k + 1] += *((float *) V[j] + n + 1) * Q[j][k];
	}
      }
      for (int j = np; j < m; ++j) {
	for (int k = j - np; k < m; ++k) {
	  scratch_2m[2 * k] += *((float *) V[j] + n) * Q[j][k];
	  scratch_2m[2 * k + 1] += *((float *) V[j] + n + 1) * Q[j][k];
	}
      }
      for (int k = 0; k < m; ++k) {
	*((float *) V[k] + n) = (float) scratch_2m[2 * k];
	*((float *) V[k] + n + 1) = (float) scratch_2m[2 * k + 1];
      }
    }
#endif

    PROFILE_LANCZOS ("diag Q %e\n", time_elapse ());

    //
    // restart Lanczos with new residual r_k
    // first multiply r_m by Q(m,k)
    //
    if (VRB.IsActivated (VERBOSE_DEBUG_LEVEL)) {
      VRB.Debug (cname, fname, "RESTVEC %e %e\n", Q[m - 1][nk - 1],
		 beta[nk - 1]);
      {
	Float *fp = (Float *) r;
	VRB.Debug (cname, fname, "RESTVEC %e %e %e %e\n", fp[0], fp[1], fp[2],
		   fp[3]);
      }
      {
	float *fp = (float *) V[nk];
	VRB.Debug (cname, fname, "RESTVEC %e %e %e %e\n", fp[0], fp[1], fp[2],
		   fp[3]);
      }
    }

    r->VecTimesEquFloat (Q[m - 1][nk - 1], f_size);
    if (nk < m) {
      movefloattoFloat ((Float *) Apsi, (float *) V[nk], f_size);
      Apsi->VecTimesEquFloat (beta[nk - 1], f_size);
      vecAddEquVec ((Float *) r, (Float *) Apsi, f_size);
    }
    beta[nk - 1] = sqrt (r->NormSqGlbSum (f_size));


    VRB.Debug (cname, fname, "RESTVEC new beta_k %e\n", beta[nk - 1]);
    //if(beta[nk-1] < 1.0e-150)break;

    if (VRB.IsActivated (VERBOSE_DEBUG_LEVEL)) {
      VRB.Debug (cname, fname, "new residual vector\n");
      for (int j = 0; j < f_size; j++)
	VRB.Debug (cname, fname, "%d %.16f\n", j, *((Float *) r + j));
    }

    PROFILE_LANCZOS ("prep new resvec Q %e\n", time_elapse ());

    if (it % eig_arg->conv_check == 0) {

      PROFILE_LANCZOS ("before conv check %e\n", time_elapse ());
      //
      // Convergence check 
      //   by monitoring the eigenV eq. for original operator
      //
      int Ndone (0);		// Number of converged vector
      double time_conv = dclock ();

      // get shifts
      moveFloat (shifts, alpha, m);
      moveFloat (ttbeta, beta, m);
      vecZero ((IFloat *) Q[0], m * m);
      for (int j = 0; j < m; j++)
	Q[j][j] = 1.0;

      //int exam_dim = m;   // full dimension, sometimes better ? 
      int exam_dim = nk;	//only the wanted dimension

//      tqri_iters = tqri (shifts, ttbeta, exam_dim, m, Q);
      Float time_in = dclock ();
#ifndef USE_LAPACK
//#if 1
      tqri_iters = tqri (shifts, ttbeta, exam_dim, m, Q);
      VRB.Result (cname, fname,
		  "tqri (convergence test) converges at iter=%d time=%g\n",
		  tqri_iters, dclock () - time_in);
#else
      tqri_LAPACK (shifts, ttbeta, exam_dim, m, (double *) Q[0]);
      VRB.Result (fname, "tqri_LAPACK", " time(sec) %e\n",
		  dclock () - time_in);
#endif
      eigsrt (shifts , Q,  m);


      Vector *vtmp = (Vector *) (cheby_arg->tmp1);	// reuse temp vector for polynomial here
      Vector *vtmp2 = (Vector *) (cheby_arg->tmp2);	// reuse temp vector for polynomial here


      PROFILE_LANCZOS ("before residual comp. %e\n", time_elapse ());

      for (int k = 0; k < nk; k++) {

	vtmp->VecZero (f_size);
	for (int j = 0; j < exam_dim; j++) {
	  movefloattoFloat ((Float *) vtmp2, (float *) V[j], f_size);
	  vtmp->FTimesV1PlusV2 (Q[j][k], vtmp2, vtmp, f_size);
	}

	// Monitor convergence for the original RitzMat (not the polynomial of it)
	RitzMat (Apsi, vtmp);

	// eigenvalue obtained from Rayleigh quotient
	Float bt = sqrt (vtmp->NormSqGlbSum (f_size));
	Float alp = Apsi->ReDotProductGlbSum (vtmp, f_size) / bt;

	// calculate residue for eigenequation
	Apsi->FTimesV1PlusV2 (-alp, vtmp, Apsi, f_size);
	Float rnorm = sqrt (Apsi->NormSqGlbSum (f_size)) / bt;
	//Float rnorm = Apsi->NormSqGlbSum(f_size);
	if ((rnorm < StopRes) && (k<nt) )
	  Ndone++;
	VRB.Result (cname, fname, "Residual %d %d %e lambda %e %e %e\n", it, k,
		    rnorm, alp, shifts[k], alpha[k]);
	if (std::isnan (rnorm))
	  ERR.General (cname, fname, "Residual is nan\n");
	// save eigenvalues/vectors so don't have to recompute in final conv. check
	//if(it==MaxIters)alpha[k] = alp;
	//if(it==MaxIters)moveFloat((Float*)V[k],(Float*)vtmp,f_size);
      }
      PROFILE_LANCZOS ("after residual comp %e\n", time_elapse ());

      if (!UniqueID ())
	printf ("# %d Lanczos: time of conv check(sec) %e\n", it,
		dclock () - time_conv);
      if (Ndone == nt)
	break;
    }
    // np more steps of lanczos
    lanczos (nk, m, f_size, Apsi, r, V, alpha, beta, cheby_arg,
	     RitzMat_lanczos);

    //PROFILE_LANCZOS("lanczos %e\n",time_elapse());
    if (!UniqueID ())
      printf ("# %d Lanczos: time(sec) %e\n", it, dclock () - time_in);

    it++;
  }
  if (it == MaxIters + 1)
    it--;
  //-------------------------------------------------------------------------

  PROFILE_LANCZOS ("before the final  conv check %e\n", time_elapse ());
  double time_conv = dclock ();


  // Final diagonalization for V[j]
  // FIXEME: this is already done in the convergence check 
  // in the iteration.

  // get shifts
  moveFloat (shifts, alpha, m);
  moveFloat (ttbeta, beta, m);
  vecZero ((IFloat *) Q[0], m * m);
  for (int j = 0; j < m; j++)
    Q[j][j] = 1.0;

  //int exam_dim = m;   // full dimension, sometimes better ? 
  int exam_dim = nk;		//only the wanted dimension

  int tqri_iters;
//     = tqri (shifts, ttbeta, exam_dim, m, Q);
  time_in = dclock ();
#ifndef USE_LAPACK
//#if 1
  tqri_iters = tqri (shifts, ttbeta, exam_dim, m, Q);
  VRB.Result (cname, fname,
	      "tqri (convergence test) converges at iter=%d time=%g\n",
	      tqri_iters, dclock () - time_in);
#else
  tqri_LAPACK (shifts, ttbeta, exam_dim, m, (double *) Q[0]);
  VRB.Result (fname, "tqri_LAPACK", " time(sec) %e\n", dclock () - time_in);
#endif
  eigsrt (shifts , Q,  m);
#if 0
      for (int i = 0; i < m; i++){
         VRB.Result (cname, fname, "shifts[%d] = %g\n", i, shifts[i]);
	for (int j = 0; j < m; j++)
	  VRB.Result (cname, fname, "Q[%d][%d] = %g\n", j, i, Q[j][i]);
      }
      Float temp=0.; glb_sum(&temp);
//      exit(-43);
#endif

  PROFILE_LANCZOS ("before new V %e\n", time_elapse ());
  // new Vk, without large scratch array
  // FIXME:  would be slow due to non-localized memory access
#if 1
  int *tab;
  tab = (int *) smalloc (m * m * sizeof (int));
  int cnt = 0;
  for (int k = 0; k < m; ++k) {
    for (int j = 0; j < m; ++j) {
      if (fabs (Q[j][k]) > 1e-16) {
	if (!UniqueID ())
	  printf ("Qfinal %d %d = %g\n", j, k, Q[j][k]);
	tab[cnt] = j + m * k;
	cnt++;
      }
    }
  }
#endif

  for (int n = 0; n < f_size; n += 2) {
    vecZero ((IFloat *) scratch_2m, 2 * m);
#if 0
    for (int k = 0; k < m; ++k) {
      for (int j = 0; j < m; ++j) {
	scratch_2m[2 * k] += *((Float *) V[j] + n) * Q[j][k];
	scratch_2m[2 * k + 1] += *((Float *) V[j] + n + 1) * Q[j][k];
      }
    }
#endif
#if 1
    for (int i = 0; i < cnt; ++i) {
      int j = tab[i] % m;
      int k = (tab[i] - j) / m;
      scratch_2m[2 * k] += *((float *) V[j] + n) * Q[j][k];
      scratch_2m[2 * k + 1] += *((float *) V[j] + n + 1) * Q[j][k];
    }
#endif
    for (int k = 0; k < m; ++k) {
      *((float *) V[k] + n) = (float) scratch_2m[2 * k];
      *((float *) V[k] + n + 1) = (float) scratch_2m[2 * k + 1];
    }
  }

  PROFILE_LANCZOS ("after new V %e\n", time_elapse ());

  FILE *filep = 0;
  if (eig_arg->results != 0) {
    // Final matrix 
    filep = Fopen (eig_arg->results, "a");

    Fprintf (filep, "mass = %e\n", eig_arg->mass);
    Fprintf (filep, "nk=%d, np=%d\n", nk, np);
    Fprintf (filep, "number of restarting = %d\n", it);
  }
  PROFILE_LANCZOS ("after final output %e\n", time_elapse ());

  Complex *del =
    (Complex *) smalloc (cname, fname, "del", nk * nk * sizeof (Complex));

  Vector *vtmp = (Vector *) (cheby_arg->tmp1);	// reuse temp vector for polynomial here
  Vector *vtmp2 = (Vector *) (cheby_arg->tmp2);	// reuse temp vector for polynomial here

  for (int i = 0; i < nt; ++i) {

    movefloattoFloat ((Float *) vtmp, (float *) V[i], f_size);

    RitzMat (Apsi, vtmp);

    for (int j = 0; j < nt; ++j) {	// can cut into half in the exact arithmetic
      movefloattoFloat ((Float *) vtmp2, (float *) V[j], f_size);
      del[i * nk + j] = Apsi->CompDotProductGlbSum (vtmp2, f_size);	// Apsi . V[i]
    }
    //Float alp =  Apsi->ReDotProductGlbSum( V[i], f_size);
    Float alp = del[i * nk + i].real ();
    alpha[i] = alp;
    Apsi->FTimesV1PlusV2 (-alp, vtmp, Apsi, f_size);
    Float rnorm = sqrt (Apsi->NormSqGlbSum (f_size));
    //Float norm = sqrt(V[i]->NormSqGlbSum(f_size));
    Float norm = sqrt (vtmp->NormSqGlbSum (f_size));
    if (filep)
      Fprintf (filep, "Final True Residual %d %d %e norm %e lambda %.16e\n", it,
	       i, rnorm, norm, alp);
    VRB.Result (cname, fname,
		"Final True Residual %d %d %e norm %e lambda %.16e\n", it, i,
		rnorm, norm, alp);
  }
  PROFILE_LANCZOS ("after final residual comp. %e\n", time_elapse ());

  if (filep) {
    Fprintf (filep, "\n");
    Fprintf (filep, "final matrix\n");
    for (int i = 0; i < nk; ++i)
      for (int j = 0; j < nk; ++j) {
	Fprintf (filep, "%d %d %e %e\n",
		 i, j, del[i * nk + j].real (), del[i * nk + j].imag ());
      }
    Fclose (filep);		// close file
  }
  sfree (cname, fname, "del", del);

  if (!UniqueID ())
    printf ("# %d Lanczos: time of final conv check(sec) %e\n", it,
	    dclock () - time_conv);

  PROFILE_LANCZOS ("after the final  conv check %e\n", time_elapse ());

#if 0
  for (int j = 0; j < 2; j++) {
    Float *fp = (Float *) (V[j]);
    // phase convention
    Float phase = atan2 (fp[1], fp[0]);
    Rcomplex rot_phase (cos (phase), -sin (phase));

    for (int i = 0; i < f_size; i += 2) {
      Rcomplex C (fp[i], fp[i + 1]);
      C *= rot_phase;
      printf ("DEB %d %d %.16e %.16e %.16e\n", j, i, C.real (), C.imag (),
	      C.norm ());
    }
  }
#endif
  sfree (cname, fname, "tab", tab);
  sfree (cname, fname, "scratch_2m", scratch_2m);
  sfree (cname, fname, "tmp_Qi", tmp_Qi);
  sfree (cname, fname, "Q", Q);
  sfree (cname, fname, "shifts", shifts);
  sfree (cname, fname, "ttbeta", ttbeta);
  sfree (cname, fname, "beta", beta);
  sfree (cname, fname, "Apsi", Apsi);
  sfree (cname, fname, "r", r);
  sfree (cname, fname, "matrix_polynomial.tmp1", cheby_arg->tmp1);
  sfree (cname, fname, "matrix_polynomial.tmp2", cheby_arg->tmp2);

  // make sure we won't change the original RitMatOper.
  dirac_arg->RitzMatOper = save_RitzMatOper;

  return (it);
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

void DiracOp::lanczos (int k0, int m, size_t f_size,
		       Vector * Apsi, Vector * r,
		       Vector ** V, Float * alpha, Float * beta,
		       MatrixPolynomialArg * cheby_arg,
		       RitzMatType RitzMat_lanczos)
{

  const char *fname = "lanczos(i,i,i.....)";
  RitzMatType save_RitzMatOper = dirac_arg->RitzMatOper;
  dirac_arg->RitzMatOper = RitzMat_lanczos;

  //Vector* vtmp = (Vector*)(cheby_arg->tmp1); // reuse temp vector for polynomial here
  //Vector* vtmp2 = (Vector*)(cheby_arg->tmp2); // reuse temp vector for polynomial here
  Vector *vtmp = (Vector *) smalloc ("", "lanczos", "matrix_polynomial.tmp2",
				     f_size * sizeof (Float));
  Vector *vtmp2 = (Vector *) smalloc ("", "lanczos", "matrix_polynomial.tmp2",
				      f_size * sizeof (Float));

  //  if(k0==0)  
  {
    Float ff;
    ff = sqrt (r->NormSqGlbSum (f_size));
    //V[k0]->CopyVec(r, f_size);
    //V[k0]->VecTimesEquFloat(1.0/ff, f_size);
    MOVE_FLOAT ((Float *) vtmp, (Float *) r, f_size);
    //VEC_TIMESEQU_FLOAT((Float*)(V[k0]), 1.0/ff, f_size );
    VEC_TIMESEQU_FLOAT ((Float *) vtmp, 1.0 / ff, f_size);
    moveFloattofloat ((float *) V[k0], (Float *) vtmp, f_size);
  }

  for (int k = k0; k < m; ++k) {

    // --- the first step ---
    if (k == 0) {

      double time_in = dclock ();
      PROFILE_LANCZOS ("before first mult %e\n", time_elapse ());
      movefloattoFloat ((Float *) vtmp, (float *) V[k], f_size);
      RitzMat (r, vtmp, cheby_arg);
      PROFILE_LANCZOS ("first mult %e\n", time_elapse ());

      //alpha[k] = V[k]->ReDotProductGlbSum(r, f_size);
      //glb_DDOT( f_size, (IFloat*)(V[k]), (IFloat*)r, alpha+k);
      glb_DDOT (f_size, (IFloat *) vtmp, (IFloat *) r, alpha + k);

      //r->FTimesV1PlusV2(-alpha[k], V[k], r, f_size);
      //AXPY(f_size, -alpha[k], (Float*)(V[k]), (Float*)r );
      AXPY (f_size, -alpha[k], (Float *) vtmp, (Float *) r);

      //beta[k]=sqrt(r->NormSqGlbSum(f_size));
      glb_DDOT (f_size, (Float *) r, (Float *) r, beta + k);
      beta[k] = sqrt (beta[k]);


      //V[k+1]->CopyVec(r, f_size);
      //MOVE_FLOAT((Float*)(V[k+1]), (Float*)r, f_size);
      MOVE_FLOAT ((Float *) vtmp, (Float *) r, f_size);
      //V[k+1]->VecTimesEquFloat(1.0/beta[k], f_size);
      VEC_TIMESEQU_FLOAT ((Float *) vtmp, 1.0 / beta[k], f_size);
      moveFloattofloat ((float *) (V[k + 1]), (Float *) vtmp, f_size);
      PROFILE_LANCZOS ("first other linalg %e\n", time_elapse ());

//      if (!UniqueID ())
	VRB.Result (cname,fname,"First Lanczos iter: time(sec) %e\n", dclock () - time_in);

    } else {
      // --- iteration step ---

      double time_in = dclock ();
      PROFILE_LANCZOS ("before mult %e\n", time_elapse ());
      movefloattoFloat ((Float *) vtmp, (float *) V[k], f_size);
      //RitzMat(r, V[k], cheby_arg);
      RitzMat (r, vtmp, cheby_arg);
      PROFILE_LANCZOS ("mult %e\n", time_elapse ());

#if 0
      r->FTimesV1PlusV2 (-beta[k - 1], V[k - 1], r, f_size);
      alpha[k] = V[k]->ReDotProductGlbSum (r, f_size);
      r->FTimesV1PlusV2 (-alpha[k], V[k], r, f_size);
      beta[k] = sqrt (r->NormSqGlbSum (f_size));
#else
      movefloattoFloat ((Float *) vtmp2, (float *) V[k - 1], f_size);
      //AXPY(f_size, -beta[k-1], (Float*)V[k-1], (Float*)r);
      AXPY (f_size, -beta[k - 1], (Float *) vtmp2, (Float *) r);
      //alpha[k] = V[k]->ReDotProductGlbSum(r, f_size);
      //glb_DDOT( f_size, (Float*)(V[k]), (Float*)r,alpha+k );
      glb_DDOT (f_size, (Float *) vtmp, (Float *) r, alpha + k);
      //AXPY(f_size, -alpha[k], (Float*)V[k], (Float*)r);
      AXPY (f_size, -alpha[k], (Float *) vtmp, (Float *) r);
      //beta[k]=sqrt(r->NormSqGlbSum(f_size));
      glb_DDOT (f_size, (Float *) r, (Float *) r, beta + k);
      beta[k] = sqrt (beta[k]);
#endif

      PROFILE_LANCZOS ("other lin alg %e\n", time_elapse ());

      //lanczos_GramSchm_real((Float*)r, (Float**)V, k+1, f_size, 0);
      //lanczos_GramSchm_test((Float*)r, (Float**)V, k+1, f_size, 0);
      lanczos_GramSchm_test ((Float *) r, (float **) V, k + 1, f_size, 0);
      PROFILE_LANCZOS ("gramschm %e\n", time_elapse ());

      if (k + 1 < m) {
	//V[k+1]->CopyVec(r, f_size);
	//MOVE_FLOAT( (Float*) (V[k+1]), (Float*)r, f_size);
	MOVE_FLOAT ((Float *) vtmp, (Float *) r, f_size);
	//V[k+1]->VecTimesEquFloat(1.0/beta[k], f_size);
	VEC_TIMESEQU_FLOAT ((Float *) vtmp, 1.0 / beta[k], f_size);
	moveFloattofloat ((float *) (V[k + 1]), (Float *) vtmp, f_size);
      }

//      PROFILE_LANCZOS ("last linalg %e\n", time_elapse ());
//      if (!UniqueID ())
	VRB.Result (cname,fname,"# %d Lanczos iter: time(sec) %e\n", k, dclock () - time_in);
    }
  }

  dirac_arg->RitzMatOper = save_RitzMatOper;
  sfree (vtmp);
  sfree (vtmp2);
}

#if 0
#define NOINLINE_MACRO __attribute((noinline))
#else
#define NOINLINE_MACRO 
#endif


void moveFloattofloat NOINLINE_MACRO (float *out, Float * in, size_t f_size)
{
#if 1
  float flt;
  for (int i = 0; i < f_size; i++) {
    flt = (float) in[i];
    out[i] = flt;
  }
#endif
};

#if 1
void movefloattoFloat NOINLINE_MACRO (Float * out, float *in, size_t f_size)
{
  float flt;
  for (int i = 0; i < f_size; i++) {
    flt = in[i];
    out[i] = (Float) flt;
  }
};
#endif

CPS_END_NAMESPACE
