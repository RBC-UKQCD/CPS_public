#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:52 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/wilson.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
//  $Id: wilson.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.4  2001/08/16 10:50:32  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.2  2001/06/19 18:13:20  anj
//  Serious ANSIfication.  Plus, degenerate double64.h files removed.
//  Next version will contain the new nga/include/double64.h.  Also,
//  Makefile.gnutests has been modified to work properly, propagating the
//  choice of C++ compiler and flags all the way down the directory tree.
//  The mpi_scu code has been added under phys/nga, and partially
//  plumbed in.
//
//  Everything has newer dates, due to the way in which this first alteration was handled.
//
//  Anj.
//
//  Revision 1.2  2001/05/25 06:16:09  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: wilson.h,v $
//  $Revision: 1.1.1.1 $
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
#include<util/data_types.h>
CPS_START_NAMESPACE

/*--------------------------------------------------------------------------*/
/* Definitions                                                              */
/*--------------------------------------------------------------------------*/
#define ND                 4      /* Space time dimension                   */
#define SPINOR_SIZE        24     /* Dirac spinor size                      */
#define HALF_SPINOR_SIZE   12     /* Half of the Dirac spinor size          */
#define BLOCK   HALF_SPINOR_SIZE  /* Size of spin projected 2 comp. spinor  */
#define COLUMN_SPINOR_SIZE  6     /* Size of one spinor component           */
#define GAUGE_SIZE         72     /* Gauge field size per site              */

/*--------------------------------------------------------------------------*/
/* External                                                                 */
/*--------------------------------------------------------------------------*/
extern int wfm_wire_map[];     
/* Set in wilson_int. For a given index 0-7 corresponding to
   X+, X-, Y+, Y-, Z+, Z-, T+, T-  it gives the corresponding wire. */

extern int wfm_max_scu_poll;
/* Set in wilson_init. The maximun number of times a wire is polled
   before wfm_scu_wait exits with an error. Implemented only in
   d_op_wilson_opt_nos_hdw_diag and d_op_dwf_nos_hdw_diag. */

extern int wfm_scu_diag[];
/* Contains diagnostic info about an scu failure to complete for
   wfm_max_scu_poll polls. Implemented only in 
   d_op_wilson_opt_nos_hdw_diag and d_op_dwf_nos_hdw_diag.
   The elements of the array are initialized to -1 (0xffffffff) in
   wilsonon_init (called by the lattice constructor of all wilson type 
   fermions) except for element 0 that is initialized to 0.
   The ellements of the array have the following meaning:
   0    :  Number of times an scu transfer has been initiated.
           Is set to 0 in wilson_init and is incremented every time 
           scu_comm_forward or scu_comm_backward are called.
   1    :  The value of the clock at exit. Set in wilson_scu_error.
   2    :  The value of the scu status register 813040. Set in wfm_scu_wait.
   3-10 :  The number of polls left for each of the 8 wires. Set in 
           wfm_scu_wait. It correspond to hardware wires 0-7.
   11-18:  The value of the 8 scu error status registers 813060-813067
           Set in wfm_scu_wait. It correspond to hardware wires 0-7. */

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
/* The Wilson structure typedef */
typedef struct{
   int   *ptr;               /* pointer to an array with addressing offsets */
   int   yztmax;             /* # of points of the y-z-t lattice per node   */
   int   offset;             /* communication addressing related            */
   int   comm_offset[ND];    /* communication addressing related            */
   int   comm_stride[ND];    /* communication addressing related            */
   int   comm_blklen[ND];    /* communication addressing related            */
   int   comm_numblk[ND];    /* communication addressing related            */
   struct  comm_params comm[ND];  
   int   vol[2];
   int   padded_subgrid_vol[ND];
   IFloat *spinor_tmp;        /* temp spinor needed by mdagm                 */
   IFloat *af[ND];     /* point. array to 4 fwd spin projected half spinors  */
   IFloat *ab[ND];     /* point. array to 4 bwd spin projected half spinors  */
} Wilson;



/*--------------------------------------------------------------------------*/
/* Function declarations                                                    */
/*--------------------------------------------------------------------------*/


/*==========================================================================*/
/* This routine performs all initializations needed before wilson library   */
/* funcs are called. It sets the addressing related arrays and reserves     */
/* memory for the needed temporary buffers. It only needs to be called only */
/* once at the begining of the program (or after a wilson_end call)         */
/* before any number of calls to wilson functions are made.                 */
/*==========================================================================*/
void wilson_init(Wilson *wilson_p);     /* pointer to Wilson struct. */


/*==========================================================================*/
/* This routine frees any memory reserved by wilson_init                    */
/*==========================================================================*/
void wilson_end(Wilson *wilson_p);     /* pointer to Wilson struct. */



/*==========================================================================*/
/* The Wilson fermion matrix * vector routines                              */
/* They are external "C"                                                    */
/*==========================================================================*/

extern "C"{

void wilson_mdagm(IFloat  *chi,        /* chi = MdagM(u) psi          */
		  IFloat  *u,          /* Gauge field                 */
		  IFloat  *psi,        /* chi = MdagM(u) psi          */
		  IFloat  *mp_sq_p,    /* pointer to Sum |M psi|^2    */
		  IFloat  Kappa,       /* Wilson's kappa parameter    */
		  Wilson *wilson_p);  /* pointer to a Wilson struct. */

void wilson_dslash(IFloat *chi, 
		   IFloat *u, 
		   IFloat *psi, 
		   int cb,
		   int dag,
		   Wilson *wilson_p);

void wilson_m(IFloat *chi, 
	      IFloat *u, 
	      IFloat *psi, 
	      IFloat kappa,
	      Wilson *wilson_p);

void wilson_mdag(IFloat *chi, 
		 IFloat *u, 
		 IFloat *psi, 
		 IFloat kappa,
		 Wilson *wilson_p);

}

#endif
CPS_END_NAMESPACE
