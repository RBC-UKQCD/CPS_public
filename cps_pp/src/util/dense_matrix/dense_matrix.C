#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//	$Author: zs $
//	$Date: 2004-08-18 11:57:47 $
//	$Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dense_matrix/dense_matrix.C,v 1.5 2004-08-18 11:57:47 zs Exp $
//	$Id: dense_matrix.C,v 1.5 2004-08-18 11:57:47 zs Exp $
//	$Name: not supported by cvs2svn $
//	$Locker:  $
//	$Revision: 1.5 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dense_matrix/dense_matrix.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include<config.h>
#include<util/dense_matrix.h>
#include<util/data_types.h>
#include<util/smalloc.h>
CPS_START_NAMESPACE

CPS_END_NAMESPACE
#include<util/error.h>
CPS_START_NAMESPACE


const IFloat MAX_ERROR = 1e-7;

//----------------------------------------------------------------------
// Local auxiliary functions
//----------------------------------------------------------------------
inline IFloat abs(IFloat x) { return x > 0.0 ? x : -x;
}

//--------------------------------------------------------------------------
// void mat_hrm_cmpr(IFloat *mat_out, const IFloat *mat_in, int mat_n);
//--------------------------------------------------------------------------
// Purpose:
//   compress a hermitian matrix by storing only the lower triangular
//   matrix elements in the same order.
// Arguments:
//   in:   n x n hermitian matrix allocated as IFloat[n*n*2]
//   out:  compressed matrix allocated as IFloat[n*n]
// Warning:
//   in == out OK but never overlapped or embedded otherwise
//--------------------------------------------------------------------------
void mat_hrm_cmpr(IFloat *mat_out, const IFloat *mat_in, int n)
{
  // copy mat_in and mat_out. TARTAN SEEMS SCREWED UP WHEN THE ARGUMENT
  // POINTERS ARE MODIFIED !!!
  //------------------------------------------------------------------------
  IFloat *out = mat_out;         
  const IFloat *in = mat_in;
  
  // copy 2*row+1 IFloating numbers per row.
  //------------------------------------------------------------------------
  for (int row = 0; row < n; ++row) {
    for (int real_number = 0; real_number <= 2*row; ++real_number) {
      *out++ = *in++;      
    }
    in += 2*(n-row) - 1;
  }
}


//--------------------------------------------------------------------------
// void mat_hrm_decm(IFloat *mat_out, const IFloat *mat_in, int mat_n);
//--------------------------------------------------------------------------
// Purpose:
//   restore a hermitian matrix from its lower triangular half.
// Arguments:
//   out:   n x n hermitian matrix allocated as IFloat[n*n*2]
//   in:    compressed matrix allocated as IFloat[n*n]
// Warning:
//   in == out OK, but never overlapped or embeded otherwise
//--------------------------------------------------------------------------
void mat_hrm_decm(IFloat *mat_out, const IFloat *mat_in, int n)
{
  // copy the lower trianglar half in a memory-perserving order.
  //------------------------------------------------------------------------
  {
    IFloat *out = mat_out+n*n*2-1;                   // end of mat_out
    const IFloat *in = mat_in+n*n-1;                 // end of mat_in
    for (int row = 0; row < n; ++row) {             // row is from bottom up
      out -= 2*row;                                 // *out is mat[i][i][im]
      *out-- = 0.0;                                 // set mat[i][i][im] = 0
      for (int element = 0; element < 2*(n-row)-1; ++element) {
	*out-- = *in--;                             // set mat[i][i][re] 
      }
    }
  }
  
  // set the upper triangular part from the lower part.
  //------------------------------------------------------------------------
  {
    for (int row = 1; row < n; ++row) {            // row is from top down
      IFloat *out = mat_out + 2*n*row;              // mat[row][0][re]
      for (int col = 0; col < row; ++col) {
	int delta = 2*(n-1)*(row-col);             // [i][j][0]-[j][i][0]
	*(out-delta) = *out;                       // A[j][i][0]
	++out;
	*(out-delta) = - *out;                     // A[j][i][1]
	++out;
      }
    }
  }

}


//--------------------------------------------------------------------------
// void mat_hrm_ldl(IFloat *L, IFloat *D, const IFloat *A, int n);
//--------------------------------------------------------------------------
// Purpose:
//     to decompose a hermitian n by n matrix A into L D L^dag 
// Attention:
//    1. All matrices are stored in compact format for time efficiency.
//    2. A, D and L must be different from each other.
// Arguments:
//    A:  lower triangular half of a hermitian matrix. 
//        Allocated as IFloat[n*n]
//    D:  real diagonal matrix. allocated as Float[n]
//    L:  unit lower-triangular matrix. allocated as Float[n*(n-1)]
//--------------------------------------------------------------------------
// Theorem:
//
// Reference:
//  1. Matrix Computations, chapter 4 by Gene Golub and Charles Van Loan
//  2. Molecular dynamics for full QCD simulations with an improved 
//     action, by Xiang-Qian Luo, hep-lat/9603021
// Algorithm:
//   A = L D L^dag  
//     ==> A_ij = L_ik D_k L_jk*             where  k = 0 .. j,   j <= i
//     ==> A_ii = L_ik D_k L_ik* + D_i       where  k = 0 .. i-1
//         A_ij = L_ik D_k L_jk* + L_ij D_j  where  k = 0 .. j-1, j < i
//     ==> D_i      = A_ii - L_ik D_k L_ik*  where  k = 0 .. i-1
//         L_ij D_j = A_ij - L_ik D_k L_jk*  where  k = 0 .. j-1, j < i
//--------------------------------------------------------------------------
void mat_hrm_ldl(IFloat *L, IFloat *D, const IFloat *A, int n)
{
  Float *Lp = (Float *)L;  
  Float *Dp = (Float *)D;
  const Float *Ap = (const Float *)A;
  

  Float *Lij = Lp;                              // *Lij = L[i][j][0] always
  const Float *Aij = Ap;                        // *Aij = A[i][j][0] always
  for (int i = 0; i < n; ++i) {               
    for (int j = 0; j < i; ++j) {               // calculate L[i][j]
      Lij[0] = *Aij++;
      Lij[1] = *Aij++;
      Float *Lik = Lp + i * (i-1);              // *Lik = L[i][k][0]
      Float *Ljk = Lp + j * (j-1);              // *Ljk = L[j][k][0]
      for (int k = 0; k < j; ++k) {
	Lij[0] -= Dp[k] * (Lik[0] * Ljk[0] + Lik[1] * Ljk[1]);
	Lij[1] -= Dp[k] * (Lik[1] * Ljk[0] - Lik[0] * Ljk[1]);
	Lik += 2;
	Ljk += 2;	
      }
      *Lij++ /= Dp[j];                          // L[i][j][re] /= D[j]
      *Lij++ /= Dp[j];                          // L[i][j][im] /= D[j]
    }
    Dp[i] = *Aij++;                             // calculate D[i][i][0]
    Float *Lik = Lp + i * (i-1);                // *Lik = L[i][k][0]
    for (int k = 0; k < i; ++k) {
      Dp[i] -= Dp[k] * (Lik[0] * Lik[0] + Lik[1] * Lik[1]);
      Lik += 2;
    }
  }
}


//--------------------------------------------------------------------------
// void mat_inv(IFloat *out, const IFloat *in, int n, MAT_INV_ALG alg, 
//              IFloat *error_p);  
//--------------------------------------------------------------------------
// Purpose:
//         B = 1/A.    in=out OK
// Arguments:
//     A,B:   n x n complex matrix          (allocated as Float[n][n][2]
//     alg:   algorithm specified to solve the inverse
//     err:   sum { |A * B - I| } 
//            (the absolute value of all the matrix elements)
//--------------------------------------------------------------------------
// Algorithm LDL:
//  A = L D L^dag     ==> A_ij = L_ik D_k L_jk*   where  k = 0..j,   j <= i  
//  AB = I            ==> B = L^dag^inv D^inv L^inv
//  To solve Y = 1/L:
//        LY = I
//        I[i][j] = L[i][k]Y[k][j]                  k = 0..n-1
//                = Y[i][j] + L[i][k]Y[k][j]        k = 0..i-1
//        Y[i][j] = delta_ij - L[i][k] Y[k][j]      k = 0..i-1
//  To solve B = Y^dag D^inv Y:
//        B[i][j] = Y[k][i]*  Y[k][j] / D[k]        bad!
//        L^dag B = 1/D Y                           instead!
//        ==>  L[k][i]*  B[k][j] = 1/D[i] Y[i][j]   k = 0..n-1
//        ==>  B[i][j] = Y[i][j] / D[i] - L[k][i]* B[k][j],  k = i+1..n-1
//--------------------------------------------------------------------------
static void mat_inv_ldl_cmpr(Float *out, const Float *in, int n)
{
  // allocate temporary buffer
  //------------------------------------------------------------------------
  Float *buf = (Float *)smalloc(n*(n+2)*sizeof(*buf));
  Float *const Y = buf;           
  Float *const D = Y + n*2;
  Float *const L = D + n;

  // decompose the hermitian matrix into in = L D L^dag
  //------------------------------------------------------------------------
  mat_hrm_ldl((IFloat *)L, (IFloat *)D, (IFloat *)in, n);

  // get the inverse matrix.
  //------------------------------------------------------------------------
  for (int j=0; j<n; ++j) {
    int i;
    // calculate Y[i,j] = delta_ij - L[i,k] Y[k,j]  where k = 0..i-1
    //----------------------------------------------------------------------
    for (i=0; i<n; ++i) {
      Float *Lik = L + i*(i-1); 
      Y[2*i] = (i==j ? 1.0 : 0.0);
      Y[2*i+1] = 0.0;
      for (int k = 0; k<i; ++k, Lik+=2) {
	Y[2*i] -= Lik[0] * Y[2*k] - Lik[1] * Y[2*k+1];
	Y[2*i+1] -= Lik[0] * Y[2*k+1] + Lik[1] * Y[2*k];
      }
    }
    // calculate B[i][j] = Y[i][j] / D[i] - L[k][i]* B[k][j],  k = i+1..n-1
    //----------------------------------------------------------------------
    for (i = n-1; i>=j; --i) {
      Float re = Y[2*i] / D[i];
      Float im = Y[2*i+1] / D[i];      
      for (int k=i+1; k<n; ++k) {
	Float *Bkj = out + k*k + 2*j;                
	Float *Lki = L + k*(k-1) + 2*i;
	re -= Lki[0] * Bkj[0] + Lki[1] * Bkj[1];
	im -= Lki[0] * Bkj[1] - Lki[1] * Bkj[0];
      }
      Float *Bij = out + i*i + 2*j;
      Bij[0] = re;
      if (i!=j) {
	Bij[1] = im;
      } 
    }
  }
  sfree(buf); 
}

//--------------------------------------------------------------------------
//void mat_inv(IFloat *out, const IFloat *in, int n, MAT_INV_ALG alg, IFloat *)
//--------------------------------------------------------------------------
void mat_inv(IFloat *out, const IFloat *in, int n, 
	     MAT_INV_ALG alg, IFloat *err) 
{
  const Float *A = (const Float*)in;
  Float *B = (Float *)out;  

  switch (alg) {
  case MAT_INV_ALG_LDL_CMPR:     
    mat_inv_ldl_cmpr(B, A, n);
    break;
  case MAT_INV_ALG_LDL:
    mat_hrm_cmpr(out, in, n);                
    mat_inv_ldl_cmpr(B, B, n);
    mat_hrm_decm(out, out, n);
    break;
  default:
    break;
  }

  if (err) {   
    Float & error = *(Float *)err;    
    error = 0.0;
    for (int i = 0; i < n; ++i) {
      for (int k = 0; k < n; ++k) {
	Float re = (i == k ? -1.0 : 0.0);
	Float im = 0.0;
	for (int j = 0; j < n; ++j) {
	  const Float *Bjk = B;
	  const Float *Aij = A;	  
	  switch (alg) {
	  case MAT_INV_ALG_LDL_CMPR:     
	    Bjk += j*j + 2*k;
	    Aij += i*i + 2*j;	    
	    break;
	  case MAT_INV_ALG_LDL:
	    Bjk += 2*(n*j+k);
	    Aij += 2*(n*i+j);	    
	    break;
	  default:
	    break;
	  }
	  re += Aij[0] * Bjk[0] - Aij[1] * Bjk[1];	  
	  im += Aij[0] * Bjk[1] + Aij[1] * Bjk[0];
	}
	error += abs(re) + abs(im);	  
      }
    }
    if (error > MAX_ERROR)
      ERR.General("", "", "Lapack failed to invert.\n");    
  }
}

CPS_END_NAMESPACE
