**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: chulwoo $
**  $Date: 2004-01-13 20:39:48 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_nos_hdw_diag/wfm_trick_segment.asm,v 1.2 2004-01-13 20:39:48 chulwoo Exp $
**  $Id: wfm_trick_segment.asm,v 1.2 2004-01-13 20:39:48 chulwoo Exp $
**  $Name: not supported by cvs2svn $
**  $Locker:  $
**  $Log: not supported by cvs2svn $
**  Revision 1.1.1.1.10.1  2003/11/06 20:24:31  cwj
**  *** empty log message ***
**
**  Revision 1.1.1.1  2003/11/04 05:05:10  chulwoo
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
**  Revision 1.2  2001/06/19 18:13:12  anj
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
**  Revision 1.2  2001/05/25 06:16:08  cvs
**  Added CVS keywords to phys_v4_0_0_preCVS
**
**  $RCSfile: wfm_trick_segment.asm,v $
**  $Revision: 1.2 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_nos_hdw_diag/wfm_trick_segment.asm,v $
**  $State: Exp $
**
**--------------------------------------------------------------------
*****************************************************************************************
*---------------------------------------------------------------------------------------
*
* wfm_trick_segment: This file contains 3 routines. Each routine is a segment of the 
*                    trick algorithm. Each one of these routines can not be used
*                    by itsself but only as part of the trick algorithm.
*                    Since there are more than one routines using the trick
*                    algorithm this saves substantial CRAM space at the cost
*                    of few CALL instructions.
*
*---------------------------------------------------------------------------------------
*****************************************************************************************

	.version	30

	.def	_wfm_trick_segment_1
	.def	_wfm_trick_segment_2
	.def	_wfm_trick_segment_3

*---------------------------------------------------------------------------------------
* definitions
*---------------------------------------------------------------------------------------
SIGN	.set	R6

AB0	.set	AR0
AB1	.set	AR1
AB2	.set	AR2
AB3	.set	AR3

CHI	.set	AR4
OCHI	.set	AR5

CHI_R0	.set	R0
CHI_I0	.set	R1
CHI_R1	.set	R2
CHI_I1	.set	R3
CHI_R2	.set	R4
CHI_I2	.set	R5

IK2	.set	R7

*---------------------------------------------------------------------------------------
* wfm_trick_segment_1:
*---------------------------------------------------------------------------------------
	.sect	"T:wfm1"

_wfm_trick_segment_1:
	MPYF	*OCHI++(1),	IK2,	CHI_R0		; OCHI / kappa_sq
	MPYF	*OCHI++(1),	IK2,	CHI_I0		; OCHI / kappa_sq
	MPYF	*OCHI++(1),	IK2,	CHI_R1		; OCHI / kappa_sq
	MPYF	*OCHI++(1),	IK2,	CHI_I1		; OCHI / kappa_sq
	MPYF	*OCHI++(1),	IK2,	CHI_R2		; OCHI / kappa_sq
	MPYF	*OCHI++(1),	IK2,	CHI_I2		; OCHI / kappa_sq
	
	SUBF	*CHI++(1),	CHI_R0				; - CHI
	SUBF	*CHI++(1),	CHI_I0				; - CHI
	SUBF	*CHI++(1),	CHI_R1				; - CHI
	SUBF	*CHI++(1),	CHI_I1				; - CHI
	SUBF	*CHI++(1),	CHI_R2				; - CHI
	SUBF	*CHI++(1),	CHI_I2				; - CHI

	SUBF	*AB0++(1),	CHI_R0				; - AB0
	SUBF	*AB0++(1),	CHI_I0				; - AB0
	SUBF	*AB0++(1),	CHI_R1				; - AB0
	SUBF	*AB0++(1),	CHI_I1				; - AB0
	SUBF	*AB0++(1),	CHI_R2				; - AB0
	SUBF	*AB0++(1),	CHI_I2				; - AB0

	SUBF	*AB1++(1),	CHI_R0				; - AB1
	SUBF	*AB1++(1),	CHI_I0				; - AB1
	SUBF	*AB1++(1),	CHI_R1				; - AB1
	SUBF	*AB1++(1),	CHI_I1				; - AB1
	SUBF	*AB1++(1),	CHI_R2				; - AB1
	SUBF	*AB1++(1),	CHI_I2				; - AB1

	SUBF	*AB2++(1),	CHI_R0				; - AB2
	SUBF	*AB2++(1),	CHI_I0				; - AB2
	SUBF	*AB2++(1),	CHI_R1				; - AB2
	SUBF	*AB2++(1),	CHI_I1				; - AB2
	SUBF	*AB2++(1),	CHI_R2				; - AB2
	SUBF	*AB2++(1),	CHI_I2				; - AB2

	SUBF	*AB3++(1),	CHI_R0				; - AB3
	SUBF	*AB3++(1),	CHI_I0				; - AB3
	SUBF	*AB3++(1),	CHI_R1				; - AB3
	SUBF	*AB3++(1),	CHI_I1				; - AB3
	SUBF	*AB3++(1),	CHI_R2				; - AB3
	SUBF	*AB3++(1),	CHI_I2				; - AB3

	RETS

*---------------------------------------------------------------------------------------
* wfm_trick_segment_2:
*---------------------------------------------------------------------------------------
	.sect	"T:wfm1"

_wfm_trick_segment_2:
	ADDF	*AB1++(1),	CHI_R0				; + AB1
	ADDF	*AB1++(1),	CHI_I0				; + AB1
	ADDF	*AB1++(1),	CHI_R1				; + AB1
	ADDF	*AB1++(1),	CHI_I1				; + AB1
	ADDF	*AB1++(1),	CHI_R2				; + AB1
	ADDF	*AB1--(11),	CHI_I2				; + AB1

	SUBF	*AB2++(1),	CHI_I0				; - AB2
	ADDF	*AB2++(1),	CHI_R0				; + AB2
	SUBF	*AB2++(1),	CHI_I1				; - AB2
	ADDF	*AB2++(1),	CHI_R1				; + AB2
	SUBF	*AB2++(1),	CHI_I2				; - AB2
	ADDF	*AB2++(1),	CHI_R2				; + AB2

	ADDF	*AB3++(1),	CHI_R0				; + AB3
	ADDF	*AB3++(1),	CHI_I0				; + AB3
	ADDF	*AB3++(1),	CHI_R1				; + AB3
	ADDF	*AB3++(1),	CHI_I1				; + AB3
	ADDF	*AB3++(1),	CHI_R2				; + AB3
	ADDF	*AB3++(1),	CHI_I2				; + AB3

	RETS

*---------------------------------------------------------------------------------------
* wfm_trick_segment_3:
*---------------------------------------------------------------------------------------
	.sect	"T:wfm1"

_wfm_trick_segment_3:
	SUBF	*AB1++(1),	CHI_R0				; - AB1
	SUBF	*AB1++(1),	CHI_I0				; - AB1
	SUBF	*AB1++(1),	CHI_R1				; - AB1
	SUBF	*AB1++(1),	CHI_I1				; - AB1
	SUBF	*AB1++(1),	CHI_R2				; - AB1
	SUBF	*AB1++(7),	CHI_I2				; - AB1

	ADDF	*AB2++(1),	CHI_I0				; + AB2
	SUBF	*AB2++(1),	CHI_R0				; - AB2
	ADDF	*AB2++(1),	CHI_I1				; + AB2
	SUBF	*AB2++(1),	CHI_R1				; - AB2
	ADDF	*AB2++(1),	CHI_I2				; + AB2
	SUBF	*AB2++(1),	CHI_R2				; - AB2

	ADDF	*AB3++(1),	CHI_R0				; + AB3
	ADDF	*AB3++(1),	CHI_I0				; + AB3
	ADDF	*AB3++(1),	CHI_R1				; + AB3
	ADDF	*AB3++(1),	CHI_I1				; + AB3
	ADDF	*AB3++(1),	CHI_R2				; + AB3
	ADDF	*AB3++(1),	CHI_I2				; + AB3

	RETS