**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: chulwoo $
**  $Date: 2004-01-13 20:39:42 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_lcl_nos/wfm_cb_reg0.asm,v 1.2 2004-01-13 20:39:42 chulwoo Exp $
**  $Id: wfm_cb_reg0.asm,v 1.2 2004-01-13 20:39:42 chulwoo Exp $
**  $Name: not supported by cvs2svn $
**  $Locker:  $
**  $Log: not supported by cvs2svn $
**  Revision 1.1.1.1.10.1  2003/11/06 20:23:02  cwj
**  *** empty log message ***
**
**  Revision 1.1.1.1  2003/11/04 05:05:07  chulwoo
**
**  starting again
**
**
**  Revision 1.1.1.1  2003/06/22 13:34:46  mcneile
**  This is the cleaned up version of the Columbia Physics System.
**  The directory structure has been changed.
**  The include paths have been updated.
**
**
**  Revision 1.2  2001/06/19 18:12:55  anj
**  Serious ANSIfication.  Plus, degenerate double64.h files removed.
**  Next version will contain the new nga/include/double64.h.  Also,
**  Makefile.gnutests has been modified to work properly, propagating the
**  choice of C++ compiler and flags all the way down the directory tree.
**  The mpi_scu code has been added under phys/nga, and partially
**  plumbed in.
**
**  Everything has newer dates, due to the way in which this first alteration was handled.
**
**  Anj.
**
**  Revision 1.2  2001/05/25 06:16:07  cvs
**  Added CVS keywords to phys_v4_0_0_preCVS
**
**  $RCSfile: wfm_cb_reg0.asm,v $
**  $Revision: 1.2 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_lcl_nos/wfm_cb_reg0.asm,v $
**  $State: Exp $
**
**--------------------------------------------------------------------

****************************************************************************************
* wfm_cb_reg0:
*
* wfm_s_cb0 saves the contents of the circular buffer control 
* register 0 in memory (at address save_cb0_p).
*
* wfm_r_cb0 restores the contents of the circular buffer control 
* register 0 from memory address save_cb0_p.
*
* These routines do not respect c calling conventions. Also they
* change the DP. 
*
****************************************************************************************
	.version	30

	.def	wfm_s_cb0
	.def	wfm_r_cb0

	.ref	cb_cntrl_b



	.text
save_cb0_p	.word	save_cb0
save_cb0	.word	0


****************************************************************************************
* wfm_s_cb0
****************************************************************************************
wfm_s_cb0:
	LDP	@cb_cntrl_b
	LDI	@cb_cntrl_b,	AR0
	LDP	@save_cb0_p
	LDI	@save_cb0_p,	AR1
	LDI	*AR0,		R0
	STI	R0,		*AR1
	RETS


****************************************************************************************
* wfm_r_cb0
****************************************************************************************
wfm_r_cb0:
	LDP	@cb_cntrl_b
	LDI	@cb_cntrl_b,	AR0
	LDP	@save_cb0_p
	LDI	@save_cb0_p,	AR1
	LDI	*AR1,		R0
	STI	R0,		*AR0
	RETS




