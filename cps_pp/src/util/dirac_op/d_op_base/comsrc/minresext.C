/*!\file
  \brief  Definition of the Dirac operator MinResExt method.

  $Id: minresext.C,v 1.3 2004-09-02 16:56:44 zs Exp $
*/

#include <config.h>
#include <util/dirac_op.h>
#include <util/lattice.h>
#include <util/smalloc.h>
#include <util/vector.h>
#include <util/gjp.h>
#include <comms/glb.h>
#include <util/verbose.h>
#include <util/error.h>
#include <math.h>
CPS_START_NAMESPACE

/*!
  This computes the starting guess for the solver for use in the
  HMD force calculations using a  minimal residual chronological method.
  This computes the guess solution as a linear combination of a given number
  of previous solutions. 

  \param psi The chronological guess for the solution vector.
  \param phi The source vector in the equation to be solved.
  \param psi_old The previous solutions.
  \param vm The previous solutions.multiplied by \f$  M^\dagger M \f$, computed
  as a necessary by-product.
  \param degree The number of previous solutions.
  \return The residue of this guess.
*/  

Float DiracOp::MinResExt(Vector *psi, Vector *phi, Vector **psi_old, 
			Vector **vm, int degree) 
{

/*
   We want to find the best initial guess of the solution of
   A x = b, and we have N previous solutions x_i.
   The method goes something like this:

   1. Orthonormalise the x_i
   2. Form the matrix G_ij = x_i^dagger A x_j
   3. Form the vector B_i = x_i^dagger b
   4. solve A_ij a_j  = B_i
   5. x = a_i x_i
*/
    
  char *fname = "MinResExt(V*, V*, V**, V**, int)";
  VRB.Func(cname,fname);

  int f_size;

  if(lat.Fclass() == F_CLASS_CLOVER) {
    f_size = GJP.VolNodeSites() * lat.FsiteSize() / 2;
  } else {
    f_size = GJP.VolNodeSites() * lat.FsiteSize() / (lat.FchkbEvl()+1);
  }

  Vector **v = (Vector**) smalloc(degree * sizeof(Vector*));
  for (int j=0; j<degree; j++)
    *(v+j) = (Vector*) smalloc(f_size * sizeof(Float));
  Vector *r = (Vector*) smalloc(f_size *sizeof(Float));

  Float dot;
  Complex xp;
  
  // Array to hold the matrix elements
  Complex **G = (Complex**) smalloc(degree * sizeof(Complex*));
  for (int i=0; i<degree; i++) G[i] = (Complex*) smalloc(degree *sizeof(Complex));

  // Solution and source vectors
  Complex *a = (Complex*) smalloc(degree * sizeof(Complex));
  Complex *b = (Complex*) smalloc(degree * sizeof(Complex));

  // Orthonormalise the vector basis
  for (int i=0; i<degree; i++) v[i] -> CopyVec(psi_old[i],f_size);
  for (int i=0; i<degree; i++) {
    v[i] -> VecTimesEquFloat(1/sqrt(v[i]->NormSqGlbSum(f_size)),f_size);
    for (int j=i+1; j<degree; j++) {
      xp = v[i]->CompDotProductGlbSum(v[j], f_size);
      v[j]->CTimesV1PlusV2(-xp, v[i], v[j], f_size);
    }
  }

  // Perform sparse matrix multiplication and construct rhs
  for (int i=0; i<degree; i++) {
    b[i] = v[i] -> CompDotProductGlbSum(phi,f_size);
    MatPcDagMatPc(vm[i],v[i],&dot);
    G[i][i].real(dot);
    G[i][i].imag(0.0);
  }

  // Construct the matrix
  for (int j=0; j<degree; j++) {
    for (int k=j+1; k<degree; k++) {
      G[j][k] = v[j] -> CompDotProductGlbSum(vm[k],f_size);
      G[k][j] = conj(G[j][k]);
    }
  }

  // Gauss-Jordan elimination with partial pivoting
  for (int i=0; i<degree; i++) {

    // Perform partial pivoting
    int k = i;
    for (int j=i+1; j<degree; j++) if (abs(G[j][j]) > abs(G[k][k])) k = j;
    if (k != i) {
      xp = b[k];
      b[k] = b[i];
      b[i] = xp;
      for (int j=0; j<degree; j++) {
	xp = G[k][j];
	G[k][j] = G[i][j];
	G[i][j] = xp;
      }
    }

    // Convert matrix to upper triangular form
    for (int j=i+1; j<degree; j++) {
      xp = G[j][i]/G[i][i];
      b[j] -= xp * b[i];
      for (int k=0; k<degree; k++) G[j][k] -= xp * G[i][k];
    }
  }

  // Use Gaussian Elimination to solve equations and calculate initial guess
  psi -> VecTimesEquFloat(0.0,f_size);
  r -> VecTimesEquFloat(0.0,f_size);
  for (int i=degree-1; i>=0; i--) {
    a[i] = 0.0;
    for (int j=i+1; j<degree; j++) a[i] += G[i][j] * a[j];
    a[i] = (b[i]-a[i])/G[i][i];
    psi -> CTimesV1PlusV2(a[i],v[i],psi,f_size);
    r -> CTimesV1PlusV2(a[i],vm[i],r,f_size);
  }
  r -> CTimesV1MinusV2(1.0, phi, r, f_size);
  Float error = sqrt(r->NormSqGlbSum(f_size));

  VRB.Result(cname,fname,"Degree = %d, Residual = %e\n", degree, error);
  
  for (int j=0; j<degree; j++) {
    sfree(*(v+j));
    sfree(*(G+j));
  }
  sfree(v);
  sfree(r);
  sfree(G);
  sfree(a);
  sfree(b);

  return error;

}

CPS_END_NAMESPACE
