#include<config.h>
CPS_START_NAMESPACE
/*! \file
  \brief Declarations of routine used internally in the DiracOpWilson class.

  $Id: wilson.h,v 1.3 2004-04-27 03:51:17 cwj Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: cwj $
//  $Date: 2004-04-27 03:51:17 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/wilson.h,v 1.3 2004-04-27 03:51:17 cwj Exp $
//  $Id: wilson.h,v 1.3 2004-04-27 03:51:17 cwj Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: wilson.h,v $
//  $Revision: 1.3 $
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

CPS_END_NAMESPACE
#include <util/data_types.h>
CPS_START_NAMESPACE

/*--------------------------------------------------------------------------*/
/* Definitions                                                              */
/*--------------------------------------------------------------------------*/
#define ND                 4      //!< Number of space-time dimensions.
#define SPINOR_SIZE        24     //!< Number of floating-point numbers in a Dirac spinor.
#define HALF_SPINOR_SIZE   12     //!< Number of floating-point numbers in half a Dirac spinor.
#define BLOCK   HALF_SPINOR_SIZE  //!< Number of floating-point numbers in a two-spinor.
#define COLUMN_SPINOR_SIZE  6     //!< Number of floating-point numbers in a colour vector.
#define GAUGE_SIZE         72     //!< Number of floating-point numbers in a colour matrix.
#define PAD_HALF_SPINOR_SIZE 16   /* We pad to a divisor of the PEC reg size*/
                                  /* to avoid compulsory reg misses         */

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
/* The maximun number of times a wire is polled
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
typedef struct	{
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
  IFloat *ab[ND];     //!< Array of spinors  
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

void wilson_mdagm(IFloat  *chi,        /* chi = MdagM(u) psi          */
		  IFloat  *u,          /* Gauge field                 */
		  IFloat  *psi,        /* chi = MdagM(u) psi          */
		  IFloat  *mp_sq_p,    /* pointer to Sum |M psi|^2    */
		  IFloat  Kappa,       /* Wilson's kappa parameter    */
		  Wilson *wilson_p);  /* pointer to a Wilson struct. */

//--------------------------------------------------------------------

//! Multiplication by the Wilson matrix hopping term
/*!
  The multiplication is performed on all sites of a single parity.
  The resulting vector is then defined on sites of the opposite parity.
  
  \pre The gauge field and spinor vector should be in odd-even order.  

  \param chi_p_f The resulting vector
  \param u_p_f The gauge field
  \param psi_p_f The vector bing multiplied
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

//--------------------------------------------------------------------

//! Multiplication by the odd-even preconditioned Wilson matrix
/*!
  The matrix-vector multiplication is performed on all sites with odd
  parity.

  \pre The gauge field and spinor vector should be in odd-even order.
  
  \param chi_f The matrix-vector product
  \param u_f The gauge field
  \param psi_f The vector to be mutliplied
  \param kappa_f The Wilson matrix hopping parameter
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
  
  \param chi_f The matrix-vector product
  \param u_f The gauge field
  \param psi_f The vector to be mutliplied
  \param kappa_f The Wilson matrix hopping parameter
  \param wilson_p The Wilson multiplication data structure.
 */
    
void wilson_mdag(IFloat *chi, 
		 IFloat *u, 
		 IFloat *psi, 
		 IFloat kappa,
		 Wilson *wilson_p);

}

#ifdef PPC440QCDOC

/*Comms routines*/
void wfm_scatter_face(Wilson *wilson_p, int pm, int cb);
void wfm_comm_init(Wilson *wp);
void wfm_comm_forward_start(Wilson *wp);
void wfm_comm_backward_start(Wilson *wp);
void wfm_comm_forward_complete(Wilson *wp);
void wfm_comm_backward_complete(Wilson *wp);
#endif

#endif

CPS_END_NAMESPACE
