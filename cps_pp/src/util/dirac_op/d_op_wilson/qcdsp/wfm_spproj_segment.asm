**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: chulwoo $
**  $Date: 2004-01-13 20:39:42 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp/wfm_spproj_segment.asm,v 1.2 2004-01-13 20:39:42 chulwoo Exp $
**  $Id: wfm_spproj_segment.asm,v 1.2 2004-01-13 20:39:42 chulwoo Exp $
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
**  $RCSfile: wfm_spproj_segment.asm,v $
**  $Revision: 1.2 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp/wfm_spproj_segment.asm,v $
**  $State: Exp $
**
**--------------------------------------------------------------------
*****************************************************************************************
*---------------------------------------------------------------------------------------
*
* wfm_spproj_segment: This file contains a segment of the spproj algorithm
*                     as it appears in the wfm_trick_kappa_spproj and
*                     wfm_trick_kappa routines. This routine can not be used
*                     by itsself but only as part of the spproj algorithm.
*                     By reusing this segment substantial CRAM space is saved 
*                     the cost of a CALL instruction.
*
*---------------------------------------------------------------------------------------
*****************************************************************************************

	.version	30

	.def	_wfm_spproj_segment

*---------------------------------------------------------------------------------------
* definitions
*---------------------------------------------------------------------------------------
AF0	.set	AR0
AF1	.set	AR1
AF2	.set	AR2
AF3	.set	AR3

PSI	.set	AR6

TMP	.set	R0
P2R	.set	R1
P2I	.set	R2
P3R	.set	R3
P3I	.set	R4

*---------------------------------------------------------------------------------------
* wfm_spproj_segment:
*---------------------------------------------------------------------------------------
	.sect	"T:ram1"

_wfm_spproj_segment:
	ADDF	*PSI++(IR0),	P3I,		TMP		; TMP  = P0R + P3I
	.if	ROUND=1
	RND	TMP
	.endif
	ADDF	*PSI++,		P2I,		TMP		; TMP  = P1R + P2I
||	STF	TMP,		*AF0++(IR0)			; AF0  = TMP
	.if	ROUND=1
	RND	TMP
	.endif
	SUBF	P2R,		*PSI--(IR0),	TMP		; TMP  = P1I - P2R
||	STF	TMP,		*AF0++				; AF0 += TMP
	.if	ROUND=1
	RND	TMP
	.endif
	SUBF	P3R,		*PSI--,		TMP		; TMP  = P0I - P3R
||	STF	TMP,		*AF0--(IR0)			; AF0 += TMP
	.if	ROUND=1
	RND	TMP
	.endif
	ADDF	*PSI++(IR0),	P3R,		TMP		; TMP  = P0R + P3R
||	STF	TMP,		*AF0++				; AF0 += TMP
	.if	ROUND=1
	RND	TMP
	.endif
	SUBF	P2R,		*PSI++,		TMP		; TMP  = P1R - P2R
||	STF	TMP,		*AF1++(IR0)			; AF1  = TMP
	.if	ROUND=1
	RND	TMP
	.endif
	SUBF	P2I,		*PSI--(IR0),	TMP		; TMP  = P1I - P2I
||	STF	TMP,		*AF1++				; AF1 += TMP
	.if	ROUND=1
	RND	TMP
	.endif
	ADDF	*PSI--,		P3I,		TMP		; TMP  = P0I + P3I
||	STF	TMP,		*AF1--(IR0)			; AF1 += TMP
	.if	ROUND=1
	RND	TMP
	.endif
	ADDF	*PSI++(IR0),	P2I,		TMP		; TMP  = P0R + P2I
||	STF	TMP,		*AF1++				; AF1 += TMP
	.if	ROUND=1
	RND	TMP
	.endif
	SUBF	P3I,		*PSI++,		TMP		; TMP  = P1R - P3I
||	STF	TMP,		*AF2++(IR0)			; AF2  = TMP
	.if	ROUND=1
	RND	TMP
	.endif
	ADDF	*PSI--(IR0),	P3R,		TMP		; TMP  = P1I + P3R
||	STF	TMP,		*AF2++				; AF2 += TMP
	.if	ROUND=1
	RND	TMP
	.endif
	SUBF	P2R,		*PSI--,		TMP		; TMP  = P0I - P2R
||	STF	TMP,		*AF2--(IR0)			; AF2 += TMP
	.if	ROUND=1
	RND	TMP
	.endif
	SUBF	P2R,		*PSI++(IR0),	TMP		; TMP  = P0R - P2R
||	STF	TMP,		*AF2++				; AF2 += TMP
	.if	ROUND=1
	RND	TMP
	.endif
	SUBF	P3R,		*PSI++,		TMP		; TMP  = P1R - P3R
||	STF	TMP,		*AF3++(IR0)			; AF3  = TMP

	RETS

