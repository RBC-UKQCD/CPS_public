#include<config.h>
/*! \file
  \brief Declarations of routine used internally in the DiracOpWilson class.

  $Id: wilson.h,v 1.11 2013-04-05 17:46:30 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2013-04-05 17:46:30 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/wilson.h,v 1.11 2013-04-05 17:46:30 chulwoo Exp $
//  $Id: wilson.h,v 1.11 2013-04-05 17:46:30 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: wilson.h,v $
//  $Revision: 1.11 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/wilson.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
/****************************************************************************/
/* 10/16/97                                                                 */
/*                                                                          */
/* wilson.h                                                                 */
/*                                                                          */
/* C header file for the fermion wilson library wilson.lib                  */
/*                                                                          */
/****************************************************************************/

#ifndef INCLUDED_WILSON_H
#define INCLUDED_WILSON_H

#ifdef USE_QMP
#include <qmp.h>
#endif
#include <util/data_types.h>

CPS_START_NAMESPACE

/*--------------------------------------------------------------------------*/
/* Definitions                                                              */
/*--------------------------------------------------------------------------*/
#define ND                 4      //!< Number of space-time dimensions.
#define SPINOR_SIZE        24     //!< Number of floating-point numbers in a Dirac spinor.
#define HALF_SPINOR_SIZE   12     //!< Number of floating-point numbers in half a Dirac spinor.
#if TARGET != BGL && TARGET != BGP
#define PAD_HALF_SPINOR_SIZE 16
#endif
#define BLOCK   HALF_SPINOR_SIZE  //!< Number of floating-point numbers in a two-spinor.
#define COLUMN_SPINOR_SIZE  6     //!< Number of floating-point numbers in a colour vector.
#define GAUGE_SIZE         72     //!< Number of floating-point numbers in a colour matrix.

/*--------------------------------------------------------------------------*/
/* External                                                                 */
/*--------------------------------------------------------------------------*/
extern int wfm_wire_map[];     
//!< Numbers of the wires corresponding to logical directions.
/*!< The wire numbers are stored in the order
  +X, -X, +Y, -Y, +Z, -Z, +T, -T
*/

extern int wfm_max_scu_poll;
//!< Some crazy optimisation thing.
/*!
 The maximum number of times a wire is polled
   before wfm_scu_wait exits with an error. Implemented only in
   d_op_wilson_opt_nos_hdw_diag and d_op_dwf_nos_hdw_diag.
*/


extern int wfm_scu_diag[];
//!< Some crazy optimisation thing.
/*!< Contains diagnostic info about an scu failure to complete for
   wfm_max_scu_poll polls. Implemented only in 
   d_op_wilson_opt_nos_hdw_diag and d_op_dwf_nos_hdw_diag.
   The elements of the array are initialized to -1 (0xffffffff) in
   wilson_init (called by the lattice constructor of all wilson type 
   fermions) except for element 0 that is initialized to 0.

   The elements of the array have the following meaning:
-   0      Number of times an scu transfer has been initiated.
           Is set to 0 in wilson_init and is incremented every time 
           scu_comm_forward or scu_comm_backward are called.
-   1      The value of the clock at exit. Set in wilson_scu_error.
-   2      The value of the scu status register 813040. Set in wfm_scu_wait.
-   3-10   The number of polls left for each of the 8 wires. Set in 
           wfm_scu_wait. It corresponds to hardware wires 0-7.
-   11-18  The value of the 8 scu error status registers 813060-813067
	   Set in wfm_scu_wait. It corresponds to hardware wires 0-7.
*/

/*--------------------------------------------------------------------------*/
/* Structures                                                               */
/*--------------------------------------------------------------------------*/
/* This field makes one word = numblk[10bits]-blklen[10bits]-stride[12bits] */

struct comm_params
{
  unsigned int stride : 12;
  unsigned int blklen : 10;
  unsigned int numblk : 10;
};


/*--------------------------------------------------------------------------*/
/* Type definitions                                                         */
/*--------------------------------------------------------------------------*/
//! A container of data relevant to the Wilson matrix multiplication.
/*!
  This contains parameters describing the  local lattice geometry and
  pointers to workspace arrays needed by the Wilson matrix multiplication
  routines.
 */
typedef struct
{
  int   *ptr;               //!< The dimensions of the local lattice.
  int   yztmax;             /* # of points of the y-z-t lattice per node   */
  int   offset;             /* communication addressing related            */
  int   comm_offset[ND];    /* communication addressing related            */
  int   comm_stride[ND];    /* communication addressing related            */
  int   comm_blklen[ND];    /* communication addressing related            */
  int   comm_numblk[ND];    /* communication addressing related            */
  struct  comm_params comm[ND];  
  int   vol[2];             //!< The local lattice volume
  int   padded_subgrid_vol[ND];
  IFloat *spinor_tmp;        //!< Workspace array 
  IFloat *af[ND];     //!< Array of spinors
  IFloat *ab[ND];      //!< Array of spinors  

  // For Communication in SSEOMP
  IFloat *send_buf[8];
  IFloat *recv_buf[8];
#ifdef USE_QMP
  QMP_msgmem_t msgmem[8][2];
  QMP_msghandle_t msghandle[8][2];
  QMP_msghandle_t multiple[8];
#endif

  // For Profile in SSEOMP
  IFloat MultFlops,MultFlops_bnd, MultFlops_blk;
  IFloat CPU_GHZ;
  int NITR, num_threads;

} Wilson;



/*--------------------------------------------------------------------------*/
/* Function declarations                                                    */
/*--------------------------------------------------------------------------*/


/*==========================================================================*/
//! Initialisation of parameters and memory used in the Wilson matrix multiplication.
/*!
  This needs to be called once before the Wilson matrix is used.
  \param wilson_p The Wilson matrix data to be initialised
  \post Memory is allocated for spinor arrays and parameters pertaining to
  the lattice geometry are initiialised.
*/
/*==========================================================================*/
void wilson_init(Wilson *wilson_p);     /* pointer to Wilson struct. */


/*==========================================================================*/
//! Frees memory reserved by wilson_init                    
/*!
  This should be called when the Wilson matrix is no longer needed.
  \param wilson_p The Wilson matrix data.
*/
/*==========================================================================*/
void wilson_end(Wilson *wilson_p);     /* pointer to Wilson struct. */



/*==========================================================================*/
/* The Wilson fermion matrix * vector routines                              */
/*==========================================================================*/

extern "C"{

//! Multiplication by the odd-even preconditioned Wilson matrix    
/*!
  Performs the multiplication \f$ \chi = M^\dagger M \psi \f$.
  
  \param chi The result vector
  \param u The gauge field
  \param psi The vector to be multiplied
  \param m_sq_p The value of \f$ |M \psi|^2 \f$
  \param kappa The Wilson matrix mass parameter
  \param wilson_p pointer to a Wilson struct.
*/
void wilson_mdagm(IFloat  *chi,    
		  IFloat  *u,      
		  IFloat  *psi,    
		  IFloat  *mp_sq_p,
		  IFloat  Kappa,   
		  Wilson *wilson_p);

//--------------------------------------------------------------------

//! Multiplication by the Wilson matrix hopping term
/*!
  The multiplication is performed on all sites of a single parity.
  The resulting vector is then defined on sites of the opposite parity.
  
  \pre The gauge field and spinor vector should be in odd-even order.  

  \param chi The resulting vector
  \param u The gauge field
  \param psi The vector bing multiplied
  \param cb Multiplies a vector on odd parity sites if this is 0, even
            parity sites if this is 1.
  \param dag Multiply by the hermitian conjugate hopping term if this is 1,
             or not if this is 0.
  \param wilson_p The Wilson multiplication data structure.
 */
    
void wilson_dslash(IFloat *chi, 
		   IFloat *u, 
		   IFloat *psi, 
		   int cb,
		   int dag,
		   Wilson *wilson_p);

void wilson_dslash_two(Float *chi0, Float *chi1,
		   Float *u, 
		   Float *psi0, Float *psi1, 
		   int cb0, int cb1,
		   int dag,
		   Wilson *wp);

//--------------------------------------------------------------------

//! Multiplication by the odd-even preconditioned Wilson matrix
/*!
  The matrix-vector multiplication is performed on all sites with odd
  parity.

  \pre The gauge field and spinor vector should be in odd-even order.
  
  \param chi The matrix-vector product
  \param u The gauge field
  \param psi The vector to be mutliplied
  \param kappa The Wilson matrix hopping parameter
  \param wilson_p The Wilson multiplication data structure.
 */
    
void wilson_m(IFloat *chi, 
	      IFloat *u, 
	      IFloat *psi, 
	      IFloat kappa,
	      Wilson *wilson_p);

//----------------------------------------------------------------------

//! Multiplication by the odd-even preconditioned Wilson matrix
/*!
  The matrix-vector multiplication is performed on all sites of odd
   parity.

  \pre The gauge field and spinor vector should be in odd-even order.
  
  \param chi The matrix-vector product
  \param u The gauge field
  \param psi The vector to be mutliplied
  \param kappa The Wilson matrix hopping parameter
  \param wilson_p The Wilson multiplication data structure.
 */
    
void wilson_mdag(IFloat *chi, 
		 IFloat *u, 
		 IFloat *psi, 
		 IFloat kappa,
		 Wilson *wilson_p);

}

#if TARGET == BGL || TARGET == BGP
void wilson_set_sloppy(bool sloppy);
#else
inline void wilson_set_sloppy(bool sloppy) {}
#endif

CPS_END_NAMESPACE
#endif
