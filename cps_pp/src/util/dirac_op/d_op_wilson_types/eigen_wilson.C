#include <config.h>

CPS_START_NAMESPACE
/*! \file
 * \brief  Definition of DiracOpWilsonTypes class eigensolver methods.
 */
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson_types/eigen_wilson.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
/*  Compute the spectral flow using the Ritz functional minimization */
/*  routine, if desired in the Kalkreuter-Simma algorithm */

/*  Compute low lying eigenvalues of the hermitian (but not pos. def.) */
/*  matrix RitzEigMat using the Ritz functional minimization routine, */
/*  if desired in the Kalkreuter-Simma algorithm */

/*  psi		Eigenvectors			(Modify) */
/*  lambda_H	Eigenvalues			(Write) */
/*  valid_eig	flag for valid eigenvalues	(Write) */
/*  eig_arg	Argument structure		(Read) */
/*  NCG_tot	Total number of CG iter		(Return value) */

/*  Parameters in eig_arg */
/*  N_eig	Eigenvalue number		(Read) */
/*  N_min	Minimal CG iterations		(Read) */
/*  N_max	Maximal CG iterations		(Read) */
/*  Kalk_Sim	Use Kalkreuter-Simma criterion	(Read) */
/*  n_renorm	Renormalize every n_renorm iter.	(Read) */
/*  N_KS_max	Max number of Kalkreuter-Simma iterations	(Read) */
/*  RsdR_a	(absolute) residue		(Read) */
/*  RsdR_r	(relative) residue		(Read) */
/*  Rsdlam	relative accuracy of lambda	(Read) */
/*  Cv_fact	"Convergence factor" required	(Read) */
/*  n_KS	total number of Kalkreuter-Simma iterations	(Write) */
/*  ProjApsiP	flag for projecting A.psi	(Read) */

/*  n_valid	number of valid eigenvalues	(Write) */

CPS_END_NAMESPACE
#include <util/dirac_op.h>
#include <util/lattice.h>
#include <util/smalloc.h>
#include <util/vector.h>
#include <util/gjp.h>
#include <comms/glb.h>
#include <util/verbose.h>
#include <util/error.h>
#include <util/qcdio.h>
#include <math.h>
CPS_START_NAMESPACE
void MatHermElements( DiracOpWilsonTypes * dirac_op,
    Vector ** psi, 
    Vector *tmp,
    int n_vec, 
    int f_size, 
    Float * diag, 
    Complex * off_diag );



/*
 * A normalised Gram-Schmidt orthogonalisation 
 */

inline Float Norm( Vector* psi, int f_size )
{
  Float norm( sqrt(psi->NormSqGlbSum(f_size)) );
  psi->VecTimesEquFloat(1.0/norm,f_size);
  return norm;
}

/* CLAUDIO: now void function does not normalize vector! */
static void checkON( FILE *fp, Vector** psi, int N_eig, int f_size  )
{
  for (int i=0; i<N_eig; ++i ) {
    for (int j=i; j<N_eig; ++j) {
      const Complex xp(psi[i]->CompDotProductGlbSum(psi[j], f_size));
      //Fprintf(fp, "elemento [%d,%d]=(%e,%e)",i,j,xp.real(),xp.imag());
      if (xp.abs()>1.e-15 && i!=j)
        Fprintf(fp, "ORTHONORMALITY Error: elem (%d,%d)=%e\n",i,j,(IFloat)(xp.abs()));
      //Fprintf(fp, "\n");
        
    }
  }
}

/*
 * Stores the information needed to asses the error
 * on an eigenvector
 */
class EigErrInfo
{
  public:
    Float M2norm;
    Float Mnorm;

    EigErrInfo():
      M2norm(0),
    Mnorm(0)
    {;}
};

/*
 * passed eigenvector candidate, fill the 
 * EigErrInfo using Ritz
 */
inline void  EigErrRitz( DiracOpWilsonTypes *dirac_op,
    Vector* psi, 
    Vector* tmp,
    int f_size,
    EigErrInfo& inf )
{
  dirac_op->RitzMat(tmp, psi);
  inf.M2norm = tmp->ReDotProductGlbSum(tmp, f_size);
  inf.Mnorm  = psi->ReDotProductGlbSum(tmp, f_size);
}


/*
 * passed eigenvector candidate, fill the 
 * EigErrInfo with RitzEigMat
 */
inline void  EigErrRitzEig( DiracOpWilsonTypes *dirac_op,
    Vector* psi, 
    Vector* tmp,
    int f_size,
    EigErrInfo& inf )
{
  dirac_op->RitzEigMat(tmp, psi);
  inf.M2norm = tmp->ReDotProductGlbSum(tmp, f_size);
  inf.Mnorm  = psi->ReDotProductGlbSum(tmp, f_size);
}

static int CheckValidEigV(Float k, int *valid_eig, int N_eig, 
    Float *lambda_HH, Float *lambda_H, Complex *off_diag) 
{

  int i,N_dp;
  
  for(i=0; i<N_eig; ++i) {
    if ((lambda_H[i]*lambda_H[i])>(k*lambda_HH[i])) {
      valid_eig[i]=1;
    } else {
      break;
    }
  }
  N_dp=i;
  for(;i<N_eig; ++i){
    valid_eig[i]=0;
  }
  
  return N_dp;

}

static int CheckSignEigV(int *valid_eig, int N_eig, 
    Float *lambda_H, Float *lb, Float *up) 
{

  int i,N_dp;
  
  for(i=0; i<N_eig; ++i) {
    Float e2=lambda_H[i]*lambda_H[i];
    if (e2>=(up[i]-lb[i]) && e2<=up[i]) {
      valid_eig[i]=1;
    } else {
      break;
    }
  }
  N_dp=i;
  for(;i<N_eig; ++i){
    valid_eig[i]=0;
  }
  
  return N_dp;

}


static int OldCheckValidEigV(int *valid_eig, int N_eig, 
    Float *lambda_HH, Float *lambda_H, Complex *off_diag) 
{

  int i;
  
  /* 
   * start off considering all the eigenvectors
   * as possibly valid
   */
  for (i=0;i<N_eig;i++) { valid_eig[i] = 1; }

  /*
   * then iteratively remove the eigenvectors that
   * are invalid ( only count the elements of the
   * jacobi matrix between valid eigenvectors ) 
   */

  bool looping(true);
  while ( looping ) {
    looping = false;
    for (i=0;i<N_eig;i++) {
      if ( valid_eig[i] ) {
        int xy(0),x,y;
        const Float rEig2( 1/lambda_HH[i] );
        Float sum(lambda_H[i]*lambda_H[i] * rEig2);
        // loop over all off-diagonal element and pick up all those
        // of the correct row
        for (x=0;x<N_eig;x++) {
          for (y=0;y<x;y++) {
            if ( ( x==i && valid_eig[y] ) || ( y==i && valid_eig[x] ) ) {
              Float value( off_diag[xy].abs() );
              sum+=value*value*rEig2;
            }
            xy++;
          }
        }
        // sum should now have the sum of the squares of the elements
        // of row i of the matrix. If this deviates from 1 then we have 
        // some missing eigenvectors
        if ( sum < 0.98 ) {
          valid_eig[i] = 0;
          // tagged a new invalid eigenvalue -> need
          // a new iteration to take this into account
          looping = true;
        }
      } // if valid_eig
    } // loop over eignevalue
  } // while ( looping )
  
  int N_dp(0);
  for(i=0;i<N_eig;i++) {
    if ( valid_eig[i] == 1 ) N_dp++;
    else { break; }
  }
  // invalidate the rest of the eigenvectors
  //Fprintf(fp,"Number of valid eigenvectors: N_dp=%d\n",N_dp);
  for (i++;i<N_eig;i++) { valid_eig[i] = 0; }
  
  return N_dp;

}

//this assumes that psi are orthonormal
static void ProjEigV(int N_dp, int N_eig, Vector **psi, Vector *tmp, 
    Float *lambda_HH, DiracOpWilsonTypes *dirac_op, int f_size) 
{
  
  int corr=N_dp;
  
  while ( N_dp < N_eig ) {

    //ortogonalize respect previous eigenvector
    //this ortogonalization ensures that we do not take parallel 
    //vectors when we project later
    Complex c;
    for (int i=corr; i<N_dp; ++i){
      c=psi[i]->CompDotProductGlbSum(psi[N_dp], f_size);
      psi[N_dp]->CTimesV1PlusV2(-c, psi[i], psi[N_dp], f_size);
    }
    Norm(psi[N_dp], f_size);

    Float invlam(1./sqrt(fabs(lambda_HH[N_dp])));

    // tmp=  (HemiteMatrix) .  psi[N_dp] 
    dirac_op->RitzEigMat(tmp,psi[N_dp]);

    const Float proj(tmp->ReDotProductGlbSum(psi[N_dp], f_size));

    if ( proj<0. ) { invlam=-invlam; }

    /* CLAUDIO: 1/lam*Hpsi-psi (this is zero for eigenvectors of H) */
    tmp->FTimesV1PlusV2( invlam, tmp, psi[N_dp], f_size );

    //ortogonalize respect previous eigenvector
    for (int i=0; i<N_dp; ++i){
      c=psi[i]->CompDotProductGlbSum(psi[N_dp], f_size);
      psi[N_dp]->CTimesV1PlusV2(-c, psi[i], psi[N_dp], f_size);
    }
    Norm(psi[N_dp], f_size);
    
    N_dp++;
  }

}

/* non mi serve piu'...
void MatElements( DiracOpWilsonTypes * dirac_op,
    Vector ** psi,
    Vector * tmp,
    int n_vec, 
    int f_size, 
    Float * diag, 
    Complex * off_diag )
{
  int ij = 0;
  for(int i = 0; i < n_vec; i++) 
  {
    dirac_op->RitzMat(tmp, psi[i]);
    diag[i] = psi[i]->ReDotProductGlbSum(tmp, f_size);
    for(int j = 0; j < i; j++)
      off_diag[ij++] = psi[j]->CompDotProductGlbSum(tmp, f_size);
  }
}
*/

int DiracOpWilsonTypes::RitzEig(Vector **psi, Float lambda_H[], int valid_eig[], EigArg *eig_arg)
{
  char *fname = "RitzEig(V**,V**,I)";
  VRB.Func(cname,fname);

  // Initialise constants
  const Float RsdR_a  ( eig_arg->RsdR_a    );
  const Float RsdR_r  ( eig_arg->RsdR_r    );
  const Float Rsdlam  ( eig_arg->Rsdlam    );
  /* CLAUDIO: this is no longer a constant: if KS is used it returns the actual value of Cv_fact */
  Float Cv_fact ( eig_arg->Cv_fact   );
  const int N_eig     ( eig_arg->N_eig     );
  /* CLAUDIO: addedd eigacc */
  int N_eigacc  ( eig_arg->N_eigacc  );
  const int N_min     ( eig_arg->N_min     );
  const int N_max     ( eig_arg->N_max     );
  const int n_renorm  ( eig_arg->n_renorm  );
  const int MaxCG     ( eig_arg->MaxCG     );
  const int ProjApsiP ( eig_arg->ProjApsiP );

  int Kalk_Sim ( eig_arg->Kalk_Sim  );
  int N_KS_max ( eig_arg->N_KS_max  );

  // can't do a kal-sim if I only have one eigenvector
  if (N_eig == 1)    Kalk_Sim = 0;
  if (Kalk_Sim == 0) N_KS_max = 0;

  /* CLAUDIO: il using KS check if N_eigacc is >0 && < N_eig */
  /* CLAUDIO: if not using KS set it to N_eig */
  if (Kalk_Sim==0 || N_eigacc<1 || N_eigacc>N_eig) {
    N_eigacc=N_eig; /* if a non-correct value is specified, cosider all eigenvalues accurate */
  }

  // Print out input parameters
  VRB.Result(cname,fname,"mass = %g\n" ,IFloat(eig_arg->mass));
  /* CLAUDIO: added N_eigacc */
  VRB.Result(cname,fname,"N_eig = %d N_eigacc=%d fname=%s\n",N_eig,N_eigacc,eig_arg->fname);
  FILE * fp(Fopen(eig_arg->fname,"a"));
  Fprintf(fp,"neig = %i neigacc = %i\n",N_eig,N_eigacc);

  Float rsdl_sq = Rsdlam * Rsdlam;

  // Set the node checkerboard size of the fermion field
  //------------------------------------------------------------------
  const int f_size(RitzLatSize());

  // allocated a couple of eigenvector-sized vectors 
  Vector *tmp  = (Vector *) smalloc(f_size * sizeof(Float));
  if(tmp == 0)
    ERR.Pointer(cname,fname, "tmp");

  Vector *tmp2 = (Vector *) smalloc(f_size * sizeof(Float));
  if(tmp2 == 0)
    ERR.Pointer(cname,fname, "tmp2");

  Complex *off_diag (0x0);
  if (N_eig > 1)
  {
    off_diag = (Complex *) smalloc(cname,fname, "off_diag", N_eig*(N_eig-1)/2 * sizeof(Complex));
    if ( off_diag == 0x0 )
      ERR.Pointer(cname,fname,"off_diag");
    for(int i=0; i < N_eig*(N_eig-1)/2; ++i)
      off_diag[i] = 0.0;
  }
  //CLAUDIO: allocate memory for lower and upper bounds
  Float *lb=0, *ub=0, *eps=0;
  lb=(Float*)smalloc(cname,fname,"lb",3*N_eig*sizeof(Float));
  if(lb==0)
    ERR.Pointer(cname,fname,"lb");
  eps=lb+N_eig;
  ub=eps+N_eig;
  


  // Local vars 
  Float lambda_t, del_lamb;
  int j, n, ij, NCG_tot, n_KS;
  int i;
  int n_jacob;
  //Float * lambda = lambda_H;

  /*
   * we are going to pass lambda^2 back to the calling function
   * the calling function has lambda_H[N_eig...2N_eig-1] reserved so
   * we let lambda_old point to lambda_H[N_eig], and lambda^2 will 
   * be put there automatically
   */

  Float *lambda_HH = & lambda_H[N_eig]; 


  NCG_tot = 0;

  /* Note: psi will be normalized in Ritz! */

  for(n=0; n < N_eig; ++n)
  {
    lambda_HH[n] = 1.0;
    valid_eig [n] = 0;
  }

  del_lamb = 1.0;
  int n_valid=0;
  if (Kalk_Sim!=0) {
    Float oldlb2=0.; //for the stopping condition of psi[N_eigacc-1]
    for(n_KS = 1; n_KS<=N_KS_max && (del_lamb<0. || del_lamb>Rsdlam || n_valid<N_eigacc); n_KS++) {
      ij = 0;
      for(n = 1, i = 0; n <= N_eig; n++, i++) {
        /* CLAUDIO: reset Cv_fact to the requested value */
        Cv_fact = eig_arg->Cv_fact;
        //if (n<=N_eigacc+1)
        //  Cv_fact*=Cv_fact;
        
        //project all eigenvectors
        if (n_KS>3) {
          n_valid=CheckValidEigV(0.999, valid_eig, N_eig, lambda_HH, lambda_H, off_diag);
          ProjEigV(n_valid, N_eig, psi, tmp, lambda_HH, this, f_size); 
        }
        
        /* 
         * call the ritz. This call will return a negative number
         * in the case of an error. So far, the only error is
         * hitting the max_cg, in which case it will return -1;
         */
        int ncnt =  Ritz(psi, n, lambda_t, RsdR_a, RsdR_r, Rsdlam, 0.0 /*rsdl_sq*/,
            n_renorm, Kalk_Sim, N_min, N_max, Cv_fact, MaxCG,
            ProjApsiP);

        // free allocated memory and pass the error "up"
        if ( ncnt<0 ) { 
          if ( off_diag!=0x0 ) { sfree(off_diag); } 
          sfree(tmp);
          sfree(tmp2);
          //CLAUDIO: free memory for lower bounds
          sfree(lb);
          Fclose(fp);
          return ncnt; 
        }

        NCG_tot += ncnt;
        /* CLAUDIO: added Cv_fact to output */
        VRB.Result(cname,fname ,"PCNTINFO KS=%d n_eig=%d cgiters=%d Cv_fact=%e\n", n_KS, n, ncnt,(IFloat)Cv_fact);

        lambda_HH[i] = lambda_t;
        VRB.Debug(cname,fname,"yes KS: lambda[%d] = %g\n",i,(IFloat)lambda_HH[i]);

        //CLAUDIO: compute off-diagonal elements of <psi[i],HH psi[j]>
        RitzMat(tmp, psi[i]);

        for(j = 0; j < i; j++)
          off_diag[ij++] = psi[j]->CompDotProductGlbSum(tmp, f_size);
        
        //CLAUDIO: compute bounds before diagonalization for comparison
        tmp->FTimesV1PlusV2(-lambda_HH[i],psi[i],tmp,f_size);
        eps[i]=tmp->NormSqGlbSum(f_size);
        if(lambda_HH[i]<0.) Fprintf(fp,"ERROR:: lambda[%d] negative!!!\n",i); //CLAUDIO: rimuovere il check se non serve
        lb[i]=sqrt(eps[i]);
        if(lb[i]>lambda_HH[i]) lb[i]=lambda_HH[i]; //eigenvalues of HH are positive
        //ub[i]=lambda_HH[i]; //non needed before diag

      }
      
      VRB.Result(cname, fname, "NCG_tot=%d\n", NCG_tot);

      //CLAUDIO: compute bounds before diagonalization for comparison
      for(i=N_eig-1; i>0; --i) {
        Float delta=lambda_HH[i]-lb[i]-lambda_HH[i-1];
        if(delta>0.) {
          //compute new bound and compare it with the previous one
          Float newb=eps[i-1]/delta;
          if(newb<lb[i-1]) {
            lb[i-1]=newb;
          }
        }
      }
      Fprintf(fp,"BOUNDS before HH diagonalizations:\n");
      for(i=0; i<N_eig; ++i)
        Fprintf(fp,"KS=%d: lambda[%d]=%e err=%e relerr=%e\n",
            n_KS,i,(IFloat)(lambda_HH[i]),(IFloat)lb[i],(IFloat)(lb[i]/lambda_HH[i]));
      
      //CLAUDIO: diagonalize HH
      //checkON(fp, psi, N_eig, f_size);
      n_jacob = Jacobi(psi, N_eig, lambda_HH, off_diag, 1.e-16, 50);
      checkON(fp, psi, N_eig, f_size);

      //CLAUDIO: compute norm eps=|(HH-lambda) psi| for lower bounds
      //the order 1 lower bound is just lambda-eps
      for(i=0; i<N_eig; ++i) {
        RitzMat(tmp, psi[i]);
        tmp->FTimesV1PlusV2(-lambda_HH[i],psi[i],tmp,f_size);
        eps[i]=tmp->NormSqGlbSum(f_size);
        if(lambda_HH[i]<0.) Fprintf(fp,"ERROR:: lambda[%d] negative (%e)!!!\n",i,(IFloat)lambda_HH[i]); //CLAUDIO: rimuovere il check se non serve
        ub[i]=lambda_HH[i];
        lb[i]=sqrt(eps[i]);
        if(lb[i]>lambda_HH[i]) lb[i]=lambda_HH[i]; //eigenvalues of HH are positive
        if(i>0) if(lb[i]>(ub[i]-ub[i-1])) lb[i]=ub[i]-ub[i-1]; //force intervals to be separated.
      }
      //Check for interval separation
      //for(i=1;i<N_eig;++i) {
      //  if((ub[i]-lb[i])-ub[i-i]<-1.e-15)
      //    Fprintf(fp,"TEST INTERVALS1: %d lb=%e ub=%e\n",i,(IFloat)lb[i],(IFloat)ub[i-1]);
      //}
      
      //CLAUDIO: compute lower bounds
      //at second order we use the Kato-Temple bounds
      //we cannot use it on the last eigenvalue
      for(i=N_eig-1; i>0; --i) {
        Float delta=lambda_HH[i]-lambda_HH[i-1]-lb[i];
        if(delta>0.) {
          //compute new bound and compare it with the previous one
          Float newb=eps[i-1]/delta;
          if(newb<lb[i-1]) {
            lb[i-1]=newb;
          }
        }
      }
      Fprintf(fp,"BOUNDS after HH diagonalizations:\n");
      for(i=0; i<N_eig; ++i)
        Fprintf(fp,"KS=%d: lambda[%d]=%e err=%e relerr=%e\n",
            n_KS,i,(IFloat)(lambda_HH[i]),(IFloat)lb[i],(IFloat)(lb[i]/lambda_HH[i]));
      
      //Check for interval separation
      //for(i=1;i<N_eig;++i) {
      //  if((ub[i]-lb[i])-ub[i-i]<-1.e-15)
      //    Fprintf(fp,"TEST INTERVALS2: %d lb=%e ub=%e\n",i,(IFloat)lb[i],(IFloat)ub[i-1]);
      //}

      //CLAUDIO: check if eigenvectors of HH are also ev of H
      //compute matrix elements of H
      MatHermElements(this, psi, tmp, N_eig, f_size, lambda_H, off_diag); 
      //check for valid eigenvectors
      //n_valid=CheckValidEigV(0.80, valid_eig, N_eig, lambda_HH, lambda_H, off_diag);
      n_valid=CheckSignEigV(valid_eig, N_eig, lambda_H, lb, ub);
      Fprintf(fp,"VALID eigenvectors: %d\n",n_valid);
      //project non valid eigenvectors on eigenspaces of H
      //ProjEigV(n_valid, N_eig, psi, tmp, lambda_HH, this, f_size); 
      
      for(i=0; i<N_eig; ++i) {
        Float b1=(ub[i]==lb[i])?0.:sqrt(ub[i]-lb[i]);
        Float b2=sqrt(ub[i]);
        lb[i]=b1;
        ub[i]=b2;
        if (i<n_valid && lambda_H[i]<0.) {
          lb[i]=-b2;
          ub[i]=-b1;
        }
      }
      //Check for interval separation
      //for(i=1;i<N_eig;++i) {
      // if((lb[i]*lb[i])-(ub[i-i]*ub[i-1])<-1.e-15)
      //    Fprintf(fp,"TEST INTERVALS3: %d lb=%e ub=%e\n",i,(IFloat)lb[i],(IFloat)ub[i-1]);
      //}

      if (n_valid>0) {

        //now use only valid eigenvectors for bounds
        //do not recompute matrix elements
        //first convert lower and upper bound for H2 in corresponding
        //bounds for H
        
        Fprintf(fp,"BOUNDS ON H before H bounds:\n");
        for(i=0; i<N_eig; ++i){
          Float t=ub[i];
          if(i<n_valid && lambda_H[i]<0.) t=-lb[i];
          Fprintf(fp,"KS=%d: n=%d lb=%e ub=%e relerr=%e valid=%d\n",
              n_KS,i,(IFloat)(lb[i]),(IFloat)ub[i],(IFloat)((ub[i]-lb[i])/t),valid_eig[i]);
        }
        
        //sort eigenvalues of H from smallest to greatest
        int o[n_valid]; //use this array of indexes to remember the order
        int h=n_valid-1, l=0;
        for (i=n_valid;i>0; ) {
          --i;
          if (lambda_H[i]>0.) {
            o[h--]=i;
          } else {
            o[l++]=i;
          }
        }
        Fprintf(fp,"EigV intervals before H bounds:\n");
        for (i=0; i<n_valid; ++i){
          Fprintf(fp,"[+]%d (%e,%e,%e,%e)\n",i, (IFloat)lb[o[i]],
              (IFloat)lambda_H[o[i]],(IFloat)(sqrt(lambda_HH[o[i]])),(IFloat)ub[o[i]]);
        }
        for (i=n_valid; i<N_eig; ++i) {
          Fprintf(fp,"[-]%d (%e,%e,%e,%e)\n",i, (IFloat)lb[i],
              (IFloat)lambda_H[i],(IFloat)(sqrt(lambda_HH[i])),(IFloat)ub[i]);
        }
        //adjust upper bounds
        Float bound;
        //o[0] case
        if (N_eig==n_valid) {
          if(o[0]!=(N_eig-1)) {
            if(lambda_H[N_eig-1]<0.) Fprintf(fp,"ERROR: lambda_H deve essere positivo qui! (%e)\n",(IFloat)lambda_H[N_eig-1]);
            bound=-lb[N_eig-1];
          } else {
            //bound cannot be refined
            bound=1.e6; //put a big positive number here
          }

        } else { //n_valid<N_eig;
          bound=-lb[n_valid];
        }
        for(i=0; i<n_valid; ++i){
          int n=o[i];
          if (bound<lambda_H[n]) {
            //compute new upper bound and compare with the old one
            Float e2H=lambda_HH[n]-lambda_H[n]*lambda_H[n];
            if (e2H<0.) { Fprintf(fp,"ERROR:e2H[%d]=%e negative !!!\n",n,(IFloat)e2H); }
            else {
              bound=lambda_H[n]+e2H/(lambda_H[n]-bound);
              if (bound<ub[n]) ub[n]=bound;
            }
          } 
          bound=ub[n];
        }
        //adjust lower bounds
        //o[n_valid-1] case
        if(N_eig==n_valid) {
          if (o[n_valid-1]!=(N_eig-1)) {
            if(lambda_H[N_eig-1]>0.) Fprintf(fp,"ERRORE: lambda_H deve essere negativo qui! (%e)\n",(IFloat)lambda_H[N_eig-1]);
            bound=-ub[N_eig-1];
          } else {
            //bound cannot be refined
            bound=-1.e6; //put a big negative number
          }
        } else {
          bound=lb[n_valid];
        }
        for(i=n_valid; i>0; ) {
          --i;
          int n=o[i];
          if (bound>lambda_H[n]) {
            //compute new upper bound and compare with the old one
            Float e2H=lambda_HH[n]-lambda_H[n]*lambda_H[n];
            if (e2H<0.)  { Fprintf(fp,"ERRORE:e2H[%d]=%e negativo !!!\n",n,(IFloat)e2H); }
            else {
              bound=lambda_H[n]-e2H/(bound-lambda_H[n]);
              if (bound>lb[n]) lb[n]=bound;
            }
          } 
          bound=lb[n];
        }

        Fprintf(fp,"EigV intervals after H bounds:\n");
        for (i=0; i<n_valid; ++i){
          Fprintf(fp,"[+]%d (%e,%e,%e,%e)\n",i, (IFloat)lb[o[i]],
              (IFloat)lambda_H[o[i]],(IFloat)(sqrt(lambda_HH[o[i]])),(IFloat)ub[o[i]]);
        }
        for (i=n_valid; i<N_eig; ++i) {
          Fprintf(fp,"[-]%d (%e,%e,%e,%e)\n",i, (IFloat)lb[i],
              (IFloat)lambda_H[i],(IFloat)(sqrt(lambda_HH[i])),(IFloat)ub[i]);
        }
        
        Fprintf(fp,"BOUNDS ON H after H bounds:\n");
        for(i=0; i<N_eig; ++i){
          Float t=ub[i];
          if(i<n_valid && lambda_H[i]<0.) t=-lb[i];
          Fprintf(fp,"KS=%d: n=%d lb=%e ub=%e relerr=%e valid=%d\n",
              n_KS,i,(IFloat)(lb[i]),(IFloat)ub[i],(IFloat)((ub[i]-lb[i])/t),valid_eig[i]);
        }

        //Check for interval separation
        //for(i=1;i<N_eig;++i) {
        //  if((lb[i]*lb[i])-(ub[i-i]*ub[i-1])<-1.e-15)
        //    Fprintf(fp,"TEST INTERVALS4: %d lb=%e ub=%e\n",i,(IFloat)lb[i],(IFloat)ub[i-1]);
        //}
        
        //compute error
        i=0;
        if (N_eig==n_valid) {
          if(o[0]!=(N_eig-1)) {
            bound=-lb[N_eig-1];
          } else {
            //in this case we do not compute the error
            //go to next one 
            eps[o[0]]=0.;
            bound=ub[o[0]];
            i=1;
          }
        } else { //n_valid<N_eig;
          bound=-lb[n_valid];
        }
        for(; i<n_valid; ++i){
          int n=o[i];
          if (ub[n]!=bound) eps[n]=(ub[n]-lb[n])/(ub[n]-bound);
          else eps[n]=-1.;
          bound=ub[n];
        }
        
        i=n_valid;
        if(N_eig==n_valid) {
          if (o[n_valid-1]!=(N_eig-1)) {
            bound=-ub[N_eig-1];
          } else {
            //in this case we do not compute the error
            //go to next one 
            bound=lb[o[n_valid-1]];
            i=n_valid-1;
          }
        } else {
          bound=lb[n_valid];
        }
        for(; i>0; ) {
          --i;
          int n=o[i];
          Float err;
          if(lb[n]!=bound) err=(ub[n]-lb[n])/(bound-lb[n]);
          else err=-1;
          if (eps[n]<0.) {
            if (err<eps[n]) eps[n]=err;
          } else {
            if (err<0. || eps[n]<err) { eps[n]=err; }
          }
          bound=lb[n];
        }
        
      }
      
      //CLAUDIO: 
      //compute error of non valid eigv
      if (n_valid!=N_eig) { //compute error of max eigv
        for (i=N_eig; i>n_valid; ) {
          --i;
          eps[i]=0.;
          if (i!=N_eig-1) {
            if (lb[i+1]!=lb[i]) eps[i]=(ub[i]-lb[i])/(lb[i+1]-lb[i]);
            else eps[i]=-1.;
          }
          Float lv,err;
          if(i==n_valid) {
            if (n_valid>0) 
              lv=(lambda_H[n_valid-1]>0.)?ub[n_valid-1]:-lb[n_valid-1];
            else 
              lv=ub[i]; //we are considering the lower and are all nonvalid 
                        //we put here a number so that err<0.
          } else {
            lv=ub[i-1];
          }
          if (ub[i]!=lv) err=(ub[i]-lb[i])/(ub[i]-lv);
          else err=-1.;
          if (eps[i]<0.) {
            if (err<eps[i]) eps[i]=err;
          } else {
            if (err<0. || eps[i]<err) { eps[i]=err; }
          }
        }
      }
      //for the smallest we also require to be separated from zero
      del_lamb=ub[0]-lb[0];
      del_lamb/=(ub[0]>0.)?ub[0]:-lb[0];
      if(eps[0]>=0. && eps[0]<del_lamb) eps[0]=del_lamb;
      
      Fprintf(fp,"stopping conditions:\n");
      for(i=0; i<N_eig; ++i){
        Fprintf(fp,"n=%d eps=%e\n",i,eps[i]);
      }
      
      //now find max error of N_eigacc
      del_lamb=0.;
      //for(i=0; i<N_eigacc; ++i) {
      for(i=0; i<N_eigacc-1; ++i) {//check only N_eigacc-1 eigenvector if u want to use the other method below
        if (eps[i]<0.) {
          //set del_lamb negative and exit
          del_lamb=-1.;
          break;
        }
        if (eps[i]>del_lamb) del_lamb=eps[i];
      }
      //now make check on the last
      // this is somewhat faster but a little bit unsafer
      if(del_lamb>0.) {
        if(eps[N_eigacc-1]>0.) {
          Float lb2=lb[N_eigacc-1];
          //lb2*=lb2; //compute square
          Float ub2=ub[N_eigacc-1];
          //ub2*=ub2;
          //Float dl=(ub2-lb2)/(ub2-ub[0]*ub[0]);
          Float dl=(ub2-lb2)/(fabs(ub2)-fabs(ub[0]));
          //oldlb2=lb2;
          if (dl>Rsdlam || dl<0.) {
            if (eps[N_eigacc-1]>del_lamb) del_lamb=eps[N_eigacc-1];
            //n_valid=0;//make the loop continue
          }
        } else {
          del_lamb=-1.; 
        }
      }
      
      
      
      //Check for interval separation
      //for(i=1;i<N_eig;++i) {
      //  if((lb[i]*lb[i])-(ub[i-i]*ub[i-1])<-1.e-15)
      //    Fprintf(fp,"TEST INTERVALS5: %d lb=%e ub=%e\n",i,(IFloat)lb[i],(IFloat)ub[i-1]);
      //}

      //recheck validity requiring high precision
      //n_valid=CheckValidEigV(0.99, valid_eig, N_eig, lambda_HH, lambda_H, off_diag);
      
      /*
      Float oldup, up, lo, sw;
      up=ub[0]*ub[0];
      lo=lb[0]*lb[0];
      if(up<lo) { sw=lo; lo=up; up=sw; }
      oldup=up;
      del_lamb=(up-lo)/up;
      for(i=1; i<N_eigacc; ++i) {
        up=ub[i]*ub[i];
        lo=lb[i]*lb[i];
        if(up<lo) { sw=lo; lo=up; up=sw; }
        lambda_t=(up-lo)/(up-oldup);
        oldup=up;
        if (lambda_t<0.) Fprintf(fp,"ERRORE:lambda_t negativo !!!\n");
        if(lambda_t>del_lamb) del_lamb=lambda_t;
      }
      */

      //del_lamb=lb[0]/lambda_HH[0];
      //for(i=1; i<N_eigacc; ++i) {
      //  //lambda_t=lb[i]/lambda[i];
      //  lambda_t=lb[i]/(lambda_HH[i]-lambda_HH[i-1]);
      //  if(lambda_t>del_lamb) del_lamb=lambda_t;
      //}
      /* CLAUDIO: copy eigenvalues of HH in lambda_HH */
      // CLAUDIO: no longer needed
      //for(i=0; i<N_eig; ++i) lambda_HH[i] = lambda[i]; 

      //for(i=0; i<N_eig; ++i)
      //  VRB.Result(cname,fname,"KS=%d: lambda[%d]=%e err=%e relerr=%e (eps=%e)\n",
      //      n_KS,i,(IFloat)(lambda_HH[i]),(IFloat)lb[i],(IFloat)(lb[i]/lambda_HH[i]),(IFloat)sqrt(eps[i]));

    }

    // check to see if maxed out on KS steps.

    if (n_KS >= N_KS_max) {
      /* 
       * want to be able to recover from this in alg_eig (by
       * moving onto the next mass, so just return a negative
       * value for the iteration number and cope with the problem
       * higer up. Choose -2 so as to not clash with the error code
       * from the ritz call
       */
      VRB.Warn(cname,fname,"NONConversion_warning: n_KS = %d  N_KS_max = %d\n",
          n_KS,N_KS_max);

      if ( off_diag !=0x0 ) { sfree(off_diag); } 
      sfree(tmp); 
      sfree(tmp2); 
      //CLAUDIO: free memory for lb
      sfree(lb);
      Fclose(fp);
      return -2;
    }
  } else { /* Kalk_Sim == NO */
    ij = 0;
    for(n = 1, i = 0; n <= N_eig; n++, i++) {
      int nret(Ritz(psi, n, lambda_t, RsdR_a, RsdR_r, Rsdlam, 0.0 /*rsdl_sq*/,
            n_renorm, 0, 0, 0, Cv_fact, MaxCG, ProjApsiP));
      if ( nret < 0 ) {
        if ( off_diag !=0x0 ) { sfree(off_diag); } 
        sfree(tmp); 
        sfree(tmp2); 
        Fclose(fp);
        return nret;
      }
      NCG_tot += nret;
      lambda_HH[i] = lambda_t;
      MatHermElements(this, psi, tmp, N_eig, f_size, lambda_H, off_diag); 
      VRB.Debug(cname,fname,"no KS: lambda[%d] = %g\n",i,(IFloat)lambda_HH[i]);
    }
  }

  Fprintf(fp,"Eigenvalue intervals:\n");
  for(i=0; i<N_eig; ++i){
    Float t=ub[i];
    if(i<n_valid && lambda_H[i]<0.) t=-lb[i];
    lat.Gamma5(tmp, psi[i], GJP.VolNodeSites());
    Float chirality = psi[i]->ReDotProductGlbSum4D(tmp, f_size);
    Fprintf(fp,"[I] %d (%e,%e) %e relerr=%e valid=%d\n",
        i,(IFloat)(lb[i]*(4.+eig_arg->mass)),(IFloat)(ub[i]*(4.+eig_arg->mass)),
        (IFloat)chirality,(IFloat)((ub[i]-lb[i])/t),valid_eig[i]);
  }

  //workaround for spectral flow
  lambda_H[1]=(fabs(lb[1])<fabs(ub[1]))?lb[1]:ub[1];
  
  /* work out how good the eigenvectors are */
  Fprintf(fp,"M^2 eigenvector test (valid=%d)\n",n_valid);
  for (i=0;i<N_eig;i++) {
    EigErrInfo Inf;
    EigErrRitz( this, psi[i], tmp, f_size, Inf );
    Fprintf(fp,"%i %g %g %g\n",i,
        Inf.M2norm,
        Inf.Mnorm,
        1-Inf.M2norm/(Inf.Mnorm*Inf.Mnorm));
  }

  Fprintf(fp,"H eigenvector test\n");
  for (i=0;i<n_valid;i++) {
    EigErrInfo Inf;
    EigErrRitzEig( this, psi[i], tmp, f_size, Inf );
    Fprintf(fp,"%i %e %e %e\n",i,
        Inf.M2norm,
        Inf.Mnorm,
        1-Inf.M2norm/(Inf.Mnorm*Inf.Mnorm));
  }
  Fclose(fp);

  // deallocate memory
  if( off_diag!=0x0 ) sfree(off_diag);
  sfree(tmp);
  sfree(tmp2);
  //CLAUDIO: free memory for lb
  sfree(lb);

  // end of calculation, dump some "valuable" information.
  VRB.Result(cname,fname,"Final Jacobi: n_jacob=%d  NCG=%d\n", 
      n_jacob,NCG_tot); 

  for(n = 0; n < N_eig; n++)
    VRB.Debug(cname,fname,"lambda_HH[%d] = %g\n",n,(float)lambda_HH[n]);
  for(n = 0; n < N_eig; n++)
    VRB.Debug(cname,fname,"lambda_H[%d] = %g\n",n,(float)lambda_H[n]);

  return NCG_tot;

}



//-------------------------------------------------------------------
// Calculate the Matrix elements for the hermitian dirac operator.
// only the real part of the diagonal(the imaginary part should be 0.0
// and half of the complex off-diagonal elements are calculated.
//-------------------------------------------------------------------
void MatHermElements( DiracOpWilsonTypes * dirac_op,
    Vector ** psi,
    Vector * tmp,
    int n_vec, 
    int f_size, 
    Float * diag, 
    Complex * off_diag )
{
  int ij = 0;
  for(int i = 0; i < n_vec; i++) 
  {
    dirac_op->RitzEigMat(tmp, psi[i]);
    diag[i] = psi[i]->ReDotProductGlbSum(tmp, f_size);
    VRB.Debug("","MatHermElements","diag[%d] = %g\n",i,(float)diag[i]);
    for(int j = 0; j < i; j++)
      off_diag[ij++] = psi[j]->CompDotProductGlbSum(tmp, f_size);
  }
}



CPS_END_NAMESPACE
