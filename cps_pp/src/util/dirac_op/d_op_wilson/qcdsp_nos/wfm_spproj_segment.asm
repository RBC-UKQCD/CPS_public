**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: zs $
**  $Date: 2004-08-18 11:57:59 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_nos/wfm_spproj_segment.asm,v 1.4 2004-08-18 11:57:59 zs Exp $
**  $Id: wfm_spproj_segment.asm,v 1.4 2004-08-18 11:57:59 zs Exp $
**  $Name: not supported by cvs2svn $
**  $Locker:  $
**  $Revision: 1.4 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_nos/wfm_spproj_segment.asm,v $
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

