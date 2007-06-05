#include <config.h>
CPS_START_NAMESPACE
/*--------------------------------------------------------------------*/
/*! \file
  \brief  Routine used internally in the DiracOpWilson class.
*/
/*--------------------------------------------------------------------*/

/***************************************************************************/
/*                                                                         */
/* wilson_mdagm: It calculates chi = M^dag * M * psi, where M is the       */
/* Wilson fermion matrix. M is a function of the gauge fields u.           */
/* The sum |M*psi|^2 is also calulated and stored in mp_sq_p.              */
/* Kappa is the usual Wilson parameter and lx,ly,lz,lt is the lattice      */
/* size.                                                                   */
/*                                                                         */
/* This routine is to be used with scalar machines.                        */
/*                                                                         */
/* WARNING:                                                                */
/*                                                                         */
/* This set of routines will work only if the node sublattices have        */
/* even number of sites in each direction.                                 */
/*                                                                         */
/***************************************************************************/

CPS_END_NAMESPACE
#include <util/data_types.h>
#include <util/wilson.h>
#include <util/vector.h>
#include <util/dirac_op.h>
CPS_START_NAMESPACE


//! Access to the elements of the \e SU(3) matrix
/*!
  Gets the element of the \e SU(3) matrix \e u with row \e row,
  column \e col and complex component \e d
*/
#define U(r,row,col,d,n,cb) *(u+(r+2*(row+3*(col+3*(d+4*(n+vol[0]*(cb)))))))
//! Access to the elements of a spinor vector.
/*!
  Gets the element of the spinor \e psi with spin \e s,
  colour \e c and complex component \e r
*/
#define PSI(r,c,s,n)     *(psi+(r+2*(c+3*(s+4*(n)))))
//! As above, but the vector is called chi
#define CHI(r,c,s,n)     *(chi+(r+2*(c+3*(s+4*(n)))))
#define TMP1(r,c,s,n)    *(tmp1+(r+2*(c+3*(s+4*(n)))))
#define TMP2(r,c,s,n)    *(tmp2+(r+2*(c+3*(s+4*(n)))))


/*--------------------------------------------------------------------------*/

/*!
  The matrix-vector product   \f$ M^\dagger M\psi \f$ is computed
  on all sites with odd parity.

  \pre The gauge field and spinor vector should be in odd-even order.
  
  \param chi_f The matrix-vector product
  \param u_f The gauge field
  \param psi_f The vector to be mutliplied
  \param mp_sq_p_f If this is initially non-NULL, then \f$ |M\psi|^2 \f$ is
  computed and stored at this address.
  \param kappa_f The Wilson matrix hopping parameter
  \param wilson_p The Wilson multiplication data structure.
 */

void wilson_mdagm(IFloat *chi_f, 
		  IFloat *u_f, 
		  IFloat *psi_f, 
		  IFloat *mp_sq_p_f,
		  IFloat kappa_f,
		  Wilson *wilson_p)
{
  IFloat *tmp1_f;
  IFloat *tmp2_f;
  Float sum;
  int vol;
  int r, c, s, n;
  int i, size;

/*--------------------------------------------------------------------------*/
/* Initializations                                                          */
/*--------------------------------------------------------------------------*/
  vol =  wilson_p->vol[0];
  size = 24*vol;
  tmp1_f = wilson_p->spinor_tmp;
  tmp2_f = tmp1_f + SPINOR_SIZE * vol;

  Float *chi = (Float *) chi_f;
  Float *psi = (Float *) psi_f;
  Float *mp_sq_p = (Float *) mp_sq_p_f;
  Float mkappa2 = -Float(kappa_f) * Float(kappa_f);
  Float *tmp2 = (Float *) tmp2_f;


/*--------------------------------------------------------------------------*/
/* Dslash_E0                                                                */
/*--------------------------------------------------------------------------*/
  wilson_dslash(tmp1_f, u_f, psi_f, 1, 0, wilson_p);

/*--------------------------------------------------------------------------*/
/* Dslash_0E                                                                */
/*--------------------------------------------------------------------------*/
  wilson_dslash(tmp2_f, u_f, tmp1_f, 0, 0, wilson_p);

/*--------------------------------------------------------------------------*/
/* 1_OO - kappa^2 * Dslash_0E * Dslash_E0                                   */
/*    tmp2[i] = psi[i] + mkappa2 * tmp2[i]                                  */
/*--------------------------------------------------------------------------*/
#if 0
  fTimesV1PlusV2(tmp2, mkappa2, tmp2, psi, size);
  DiracOp::CGflops += vol*24*2;
  if(mp_sq_p != 0){
    *mp_sq_p = dotProduct(tmp2, tmp2, size);
    DiracOp::CGflops += vol*24*2;
  }
#else
  if(mp_sq_p != 0) {
    xaxpy_norm(&mkappa2,tmp2,psi,size/6,mp_sq_p);
    DiracOp::CGflops += vol*24*4;
  } else {
    xaxpy(&mkappa2,tmp2,psi,size/6);
    DiracOp::CGflops += vol*24*2;
  }
#endif


/*--------------------------------------------------------------------------*/
/* DslashDag_E0                                                             */
/*--------------------------------------------------------------------------*/
  wilson_dslash(tmp1_f, u_f, tmp2_f, 1, 1, wilson_p);

/*--------------------------------------------------------------------------*/
/* DslashDag_0E                                                             */
/*--------------------------------------------------------------------------*/
  wilson_dslash(chi_f, u_f, tmp1_f, 0, 1, wilson_p);

/*--------------------------------------------------------------------------*/
/* [1_OO - kappa * DslashDag_0E * DslashDag_E0] *                           */
/*                                 [1_OO - kappa^2 * Dslash_0E * Dslash_E0] */
/*    chi[i] = tmp2[i] + mkappa2 * chi[i]                                   */
/*--------------------------------------------------------------------------*/
#if 0
  fTimesV1PlusV2(chi, mkappa2, chi, tmp2, size);
#else
  xaxpy(&mkappa2,chi,tmp2,size/6);
#endif

  DiracOp::CGflops += vol*24*2;

}







CPS_END_NAMESPACE
