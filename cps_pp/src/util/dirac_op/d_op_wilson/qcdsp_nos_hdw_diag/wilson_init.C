#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:46 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_nos_hdw_diag/wilson_init.C,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
//  $Id: wilson_init.C,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.4  2001/08/16 10:50:27  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.2  2001/06/19 18:13:12  anj
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
//  Revision 1.2  2001/05/25 06:16:08  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: wilson_init.C,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_nos_hdw_diag/wilson_init.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
/****************************************************************************/
/* 10/6/96                                                                  */
/*                                                                          */
/* wilson_int:                                                              */
/*                                                                          */
/* This routine performs all initializations needed before wilson func      */
/* are called. It sets the addressing related arrays and reserves memory    */
/* for the needed temporary buffers. It only needs to be called             */
/* once at the begining of the program (or after a wilson_end call)         */
/* before any number of calls to wilson funcs are made.                     */
/*                                                                          */
/****************************************************************************/

/*--------------------------------------------------------------------------*/
/* Include header files                                                     */
/*--------------------------------------------------------------------------*/
CPS_END_NAMESPACE
#include <sysfunc.h>
#include <scu_dir_arg.h>
#include<util/wilson.h>
#include<util/gjp.h>
#include<util/smalloc.h>
#include<util/verbose.h>
#include<util/error.h>
CPS_START_NAMESPACE

/*--------------------------------------------------------------------------*/
/* External variables                                                       */
/*--------------------------------------------------------------------------*/
int wfm_wire_map[8];
int wfm_scu_diag[19];
int wfm_max_scu_poll;

/*--------------------------------------------------------------------------*/
/* Function declarations                                                    */
/*--------------------------------------------------------------------------*/
void wfm_sublatt_pointers(int slx, 
			  int sly, 
			  int slz, 
			  int slt, 
			  int slatt_eo, 
			  Wilson *wilson_p);

/*=========================================================================*/
/* wilson_init:                                                            */
/*=========================================================================*/

void wilson_init(Wilson *wilson_p)  /* pointer to Wilson type structure    */
{
  char *cname = " ";
  char *fname = "wilson_init(Wilson*)";
  VRB.Func(cname,fname);

  int spinor_words;             /* size of the spinor field on the         */
				/* sublattice checkerboard                 */

  int half_spinor_words;        /* size of the spin-projected "half_spinors*/
                                /* on the sublattice checkerboard including*/
                                /* the communications padding              */

  int slx;                          /* x-direction size of node sublattice */
  int sly;                          /* y-direction size of node sublattice */
  int slz;                          /* z-direction size of node sublattice */
  int slt;                          /* t-direction size of node sublattice */
  int slatt_eo;                     /* =0/1 if the sublattice is even/odd. */
  int i;
  int size;



/*--------------------------------------------------------------------------*/
/* Set sublattice direction sizes                                           */
/*--------------------------------------------------------------------------*/
  slx = GJP.XnodeSites();
  sly = GJP.YnodeSites();
  slz = GJP.ZnodeSites();
  slt = GJP.TnodeSites();

/*--------------------------------------------------------------------------*/
/* Determine if the sublattice is even or odd from the node coordinates     */
/* If (px,py,pz,pt) are the coordinates of the node and the node            */
/* mesh is of size (nx,ny,nz,nt) then a node is even/odd if its             */
/* lexicographical number =  px + nx * ( py + ny * ( pz + nz * ( pt )))     */
/* is even/odd.                                                             */
/*--------------------------------------------------------------------------*/
/* A runtime system function is needed here to determine (px,py,pz,pt) and  */
/* (nx,ny,nz,nt). For now we set slat_eo = 0 which is a safe choice if      */
/* slx,sly,slz,slt are all even.                                            */
/* ??? */
  slatt_eo = 0;

/*--------------------------------------------------------------------------*/
/* Reserve memory for the node sublattice pointers                          */
/*--------------------------------------------------------------------------*/
  size = 40*sly*slz*slt*sizeof(int);
  wilson_p->ptr = (int *) smalloc(size);
  if( wilson_p->ptr == 0)
    ERR.Pointer(cname,fname, "ptr");
  VRB.Smalloc(cname,fname,
	      "ptr", wilson_p->ptr, size);

/*--------------------------------------------------------------------------*/
/* Set the node sublattice pointers                                         */
/*--------------------------------------------------------------------------*/
  wfm_sublatt_pointers(slx, sly, slz, slt, slatt_eo, wilson_p);


/*--------------------------------------------------------------------------*/
/* Reserve memory for 1  temporary spinor (nedded by mdagm)                 */
/*--------------------------------------------------------------------------*/
  spinor_words = SPINOR_SIZE * wilson_p->vol[1];

  wilson_p->spinor_tmp = (IFloat *) smalloc(spinor_words*sizeof(IFloat));
  if(wilson_p->spinor_tmp == 0)
    ERR.Pointer(cname,fname, "spinor_tmp");
  VRB.Smalloc(cname,fname,
	      "spinor_tmp", wilson_p->spinor_tmp, spinor_words*sizeof(IFloat));
    

/*--------------------------------------------------------------------------*/
/* Reserve memory for the 4 forward and 4 backward spin projected half      */ 
/* spinors.                                                                 */
/*--------------------------------------------------------------------------*/
  for(i=0; i<4; i++){
    half_spinor_words = HALF_SPINOR_SIZE * wilson_p->padded_subgrid_vol[i];

    wilson_p->af[i] = (IFloat *) smalloc(half_spinor_words*sizeof(IFloat));
    if(wilson_p->af[i] == 0)
      ERR.Pointer(cname,fname, "af[i]");
    VRB.Smalloc(cname,fname,
		"af[i]", wilson_p->af[i], half_spinor_words*sizeof(IFloat));

    wilson_p->ab[i] = (IFloat *) smalloc(half_spinor_words*sizeof(IFloat));
    if(wilson_p->ab[i] == 0)
      ERR.Pointer(cname,fname, "ab[i]");
    VRB.Smalloc(cname,fname,
		"ab[i]", wilson_p->ab[i], half_spinor_words*sizeof(IFloat));
  }

/*--------------------------------------------------------------------------*/
/* Set the wire map                                                         */
/*--------------------------------------------------------------------------*/
  wfm_wire_map[0] = SCURemap( SCU_XP );
  wfm_wire_map[1] = SCURemap( SCU_XM );
  wfm_wire_map[2] = SCURemap( SCU_YP );
  wfm_wire_map[3] = SCURemap( SCU_YM );
  wfm_wire_map[4] = SCURemap( SCU_ZP );
  wfm_wire_map[5] = SCURemap( SCU_ZM );
  wfm_wire_map[6] = SCURemap( SCU_TP );
  wfm_wire_map[7] = SCURemap( SCU_TM );

/*--------------------------------------------------------------------------*/
/* Set the maximum number an scu wire is polled                             */
/*--------------------------------------------------------------------------*/
  wfm_max_scu_poll = 100000000;


/*--------------------------------------------------------------------------*/
/* Initialize wfm_scu_diag to -1                                            */
/*--------------------------------------------------------------------------*/
  for(i=0; i<19; i++){
    wfm_scu_diag[i] = -1;
  }

/*--------------------------------------------------------------------------*/
/* Initialize the number of times an scu transfer has been                  */
/* initiated. Is incremented every time scu_comm_forward or                 */
/* scu_comm_backward are called.                                            */
/*--------------------------------------------------------------------------*/
  wfm_scu_diag[0] = 0;

}
























CPS_END_NAMESPACE
