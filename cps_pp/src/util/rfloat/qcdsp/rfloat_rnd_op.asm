**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: zs $
**  $Date: 2004-08-18 11:58:08 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/rfloat/qcdsp/rfloat_rnd_op.asm,v 1.4 2004-08-18 11:58:08 zs Exp $
**  $Id: rfloat_rnd_op.asm,v 1.4 2004-08-18 11:58:08 zs Exp $
**  $Name: not supported by cvs2svn $
**  $Locker:  $
**  $Revision: 1.4 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/rfloat/qcdsp/rfloat_rnd_op.asm,v $
**  $State: Exp $
**
**--------------------------------------------------------------------
*------------------------------------------------------------------
*
* rfloat_rnd_op.asm
*
* These overloaded operators perform the actual rounding.
*
* More overloaded operators have been directly defined
* in rfloat.h but they all use the operators below.
*
*------------------------------------------------------------------

	.version 30



;=====================================================================
; rfloat& rfloat::operator+=(float a)				     |
;								     |
;   Total words of code    : 8					     |
;   Volatile registers used: R0,AR0,DP,BK			     |
;   Parameters             : AR0	holds TEMP21		     |
;			     R0		holds a			     |
;   Stack frame            : quick (AR3 points to some old fram)     |
;=====================================================================
	.def	___apl__6rfloatFf
	.def	___apl__6rfloatFf$LAJ

	.sect	"T:___apl__6rfloatFf"

___apl__6rfloatFf:
	POP	BK
___apl__6rfloatFf$LAJ:
	POP	AR0
	POPF	R0
	ADDI	2,SP

	ADDF	*AR0,R0

	BUD	BK
	RND	R0
	STF	R0,*AR0
	LDIU	AR0,R0


;=====================================================================
; rfloat& rfloat::operator-=(float a)				     |
;								     |
;   Total words of code    : 8					     |
;   Volatile registers used: R0,AR0,DP,BK			     |
;   Parameters             : AR0	holds TEMP23		     |
;			     R0		holds a			     |
;   Stack frame            : quick (AR3 points to some old frame)    |
;=====================================================================
	.def	___ami__6rfloatFf
	.def	___ami__6rfloatFf$LAJ

	.sect	"T:___ami__6rfloatFf"

___ami__6rfloatFf:
	POP	BK
___ami__6rfloatFf$LAJ:
	POP	AR0
	POPF	R0
	ADDI	2,SP

	SUBRF	*AR0,R0

	BUD	BK
	RND	R0
	STF	R0,*AR0
	LDIU	AR0,R0



;=====================================================================
; rfloat& rfloat::operator*=(float a)				     |
;								     |
;   Total words of code    : 8					     |
;   Volatile registers used: R0,AR0,DP,BK			     |
;   Parameters             : AR0	holds TEMP25		     |
;			     R0		holds a			     |
;   Stack frame            : quick (AR3 points to some old frame)    |
;=====================================================================
	.def	___amu__6rfloatFf
	.def	___amu__6rfloatFf$LAJ

	.sect	"T:___amu__6rfloatFf"

___amu__6rfloatFf:
	POP	BK
___amu__6rfloatFf$LAJ:
	POP	AR0
	POPF	R0
	ADDI	2,SP

	MPYF	*AR0,R0

	BUD	BK
	RND	R0
	STF	R0,*AR0
	LDIU	AR0,R0




;=====================================================================
; rfloat& rfloat::operator/=(float a)				     |
;								     |
;   Total words of code    : 12					     |
;   Volatile registers used: R0,R1,R2,R3,AR0,AR2,DP,BK		     |
;   Parameters             : AR0	holds TEMP27		     |
;			     R0		holds a			     |
;   Stack frame            : quick (AR3 points to some old frame)    |
;=====================================================================
	.def	___adv__6rfloatFf
	.ref	DEFALT
	.ref	ARTDIVF32UZ
	.def	___adv__6rfloatFf$LAJ

	.sect	"T:___adv__6rfloatFf"

___adv__6rfloatFf$LAJ:
	PUSH	BK
___adv__6rfloatFf:
	LDIU	SP,AR0
	LDFU	*-AR0(2),R1
	LDIU	*-AR0(1),AR2

	PUSH	DP

	LDP	DEFALT,DP
	RND	*AR2,R0
	CALL	ARTDIVF32UZ
	RND	R0
	STF	R0,*AR2

	LDIU	AR2,R0

	POP	DP

	RETSU	



;=====================================================================
; rfloat operator-(const rfloat& a)				     |
;								     |
;   Total words of code    : 10					     |
;   Parameters             : AR0	holds TEMP14		     |
;			     AR1	holds a			     |
;   Stack frame            : quick (AR3 points to some old frame)    |
;=====================================================================
	.def	___mi__FRC6rfloat
	.ref	___ct__6rfloatFf
	.def	___mi__FRC6rfloat$LAJ

	.sect	"T:___mi__FRC6rfloat"

___mi__FRC6rfloat$LAJ:
	PUSH	BK
___mi__FRC6rfloat:
	LDIU	SP,AR0
	LDIU	*-AR0(2),AR1
	LDIU	*-AR0(1),AR0


	NEGF	*AR1,R0
	RND	R0
	PUSHF	R0
	PUSH	AR0
	CALL	___ct__6rfloatFf
	SUBI	2,SP


	RETSU	

;=====================================================================
	.end

;=====================================================================

