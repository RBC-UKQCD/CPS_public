**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: chulwoo $
**  $Date: 2004-01-13 20:39:42 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp/wfm_scu_wait.asm,v 1.2 2004-01-13 20:39:42 chulwoo Exp $
**  $Id: wfm_scu_wait.asm,v 1.2 2004-01-13 20:39:42 chulwoo Exp $
**  $Name: not supported by cvs2svn $
**  $Locker:  $
**  $Log: not supported by cvs2svn $
**  Revision 1.1.1.1.10.1  2003/11/06 20:22:58  cwj
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
**  Revision 1.2  2001/06/19 18:12:52  anj
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
**  Revision 1.2  2001/05/25 06:16:06  cvs
**  Added CVS keywords to phys_v4_0_0_preCVS
**
**  $RCSfile: wfm_scu_wait.asm,v $
**  $Revision: 1.2 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp/wfm_scu_wait.asm,v $
**  $State: Exp $
**
**--------------------------------------------------------------------
****************************************************************************************
*---------------------------------------------------------------------------------------
*
* _wfm_scu_wait
*
* This routine waits until the Serial Communications Unit has completed 
* all communications
*
* If STAND_ALONE = 1 it can be called from C as:
*
* wfm_scu_wait();
*
* If STAND_ALONE = 0 it can not be called from C. Instead the following is expected
* to be set up before the routine is called:
*
* 1) all the .ref in the STAND_ALONE = 0  case below must be defined
*
*---------------------------------------------------------------------------------------
****************************************************************************************

	.version	30

	.include	"wfm_nga_struct.hasm"

	.def _wfm_scu_wait

*---------------------------------------------------------------------------------------
* References
*---------------------------------------------------------------------------------------
	.ref	_WilsonSCUComplete

*---------------------------------------------------------------------------------------
* definitions
*--------------------------------------------------------------------------------------
FP	.set	AR3			; use FP for arguments only
SCU	.set	AR1
TMP	.set	R0

****************************************************************************************
* _wfm_scu_wait
****************************************************************************************
	.sect 	"T:wfm0"
_wfm_scu_wait:

*---------------------------------------------------------------------------------------
* If STAND_ALONE = 1: C-callable
*---------------------------------------------------------------------------------------
	.if	STAND_ALONE = 1
	PUSH    FP
	LDI     SP,FP
*  Save all registers that are important to C
	PUSH    R4
	PUSH    R5
	RND	R6
	PUSHF   R6
	RND	R7
	PUSHF   R7
	PUSH    AR4
	PUSH    AR5
	PUSH    AR6
	PUSH    AR7
	PUSH    FP              ; Local frame pointer
	PUSH    DP
	.endif

*---------------------------------------------------------------------------------------
* Routine starts here
*---------------------------------------------------------------------------------------
	PUSH    DP

	CALL	_WilsonSCUComplete

	POP	DP
*-------------------------------------------------------------------------
*
*  Restore the registers before returning to the C program                
*
	.if	STAND_ALONE = 1
	POP     DP
	POP     FP              ;  This is our frame pointer
	POP     AR7
	POP     AR6
	POP     AR5
	POP     AR4
	POPF    R7
	POPF    R6
	POP     R5
	POP     R4

	POP     FP              ;  This is our caller's frame pointer
	.endif

	RETS
	.end

