**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: zs $
**  $Date: 2004-08-18 11:57:38 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_ghb/qcdsp/cmhb_kern_sup.asm,v 1.4 2004-08-18 11:57:38 zs Exp $
**  $Id: cmhb_kern_sup.asm,v 1.4 2004-08-18 11:57:38 zs Exp $
**  $Name: not supported by cvs2svn $
**  $Locker:  $
**  $Revision: 1.4 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_ghb/qcdsp/cmhb_kern_sup.asm,v $
**  $State: Exp $
**
**--------------------------------------------------------------------
*****************************************************************

	.asg	1,	ROUND
	.def	_log_cram
	.ref	DIV_F30_cram_c

	.text
*****************************************************

_log_cram:
	LDI     SP,AR0
        LDF     *-AR0(1),R2           ;x -> R2
	LDP	LOG_ADR 	      ;load data page
	LDI	@LOG_ADR,AR0	      ;load data address in AR0
        PUSHF   R2                    ;save x
;
; Extract the mantissa of x = f
;
        LDF     1.0,R3                ;1 -> R3
        LDI     R2,R3                 ;mantissa x -> R3
        POP     AR2                   ;x -> AR2
        ASH     -24,AR2               ;exponent x -> AR2
        ADDI    1,AR2                 ;exponent x += 1
        LDF     0.5,R2                ;0.5 -> R2
        MPYF    R2,R3                 ;f = mantissa x * 0.5
        CMPF    *AR0++,R3             ;compare mantissa x to sqrt(0.5)
	BLT	LOGB1		      ;if x < sqrt(0.5), branch to LOGB1

        SUBF    R2,R3,R0              ;f - 0.5 -> R0
	BD	LOGB2		      ;branch to LOGB2
        SUBF    R2,R0                 ;znum = R0 - 0.5
        MPYF    R2,R3,R1              ;znum / 2 -> R1
        ADDF    R2,R1                 ;zden = R1 + 0.5
LOGB1:  SUBF    R2,R3,R0              ;znum = f - 0.5
        MPYF    R2,R0,R1              ;znum * 0.5 -> R1
        ADDF    R2,R1                 ;zden = R1 + 0.5
        SUBI    1,AR2                 ;exponent x -= 1
LOGB2:  PUSH    AR0                   ;SAVE AR0 BEFORE CALLING DIV_F30
				      ;z = znum / zden
	 CALL DIV_F30_cram_c
        POP     AR0                   ;RESTORE AR0
	PUSHF	R0		      ;save z
	MPYF	R0,R0,R1	      ;w = z * z
;
; Compute R = w*a0 / (w + b0)
;
        MPYF3   *AR0++,R1,R0          ;A = -0.552707855 * R0
        ADDF    *AR0++,R1             ;B = w + -6.632718214
        PUSH    AR0                   ;SAVE AR0 BEFORE CALLING DIV_F30
				      ;A / B -> R0
	 CALL DIV_F30_cram_c
        POP     AR0                   ;RESTORE AR0
	POPF	R1		      ;z -> R0
	MPYF	R1,R0		      ;z * R0 -> R0
	ADDF	R1,R0		      ;ln(mantissa x) = z + R0
        FLOAT   AR2,R2                ;exponent x -> R2
        MPYF    *AR0++,R2,R1          ;R2 * -2.121944005 -> R1
	ADDF	R1,R0		      ;R0 + R1 -> R0
        MPYF    *AR0++,R2,R1          ;R2 * 0.693359375 -> R1
	ADDF	R1,R0		      ;ln(x) = ln(mant x) + exp x * ln(2)
        RETS                          ;return from routine

***********************************************************************
*  DEFINE CONSTANTS FOR LOG
***********************************************************************
LOG_ADR:        .word  LOG

LOG             .float  0.7071067811865475244
		.float -0.5527074855
		.float -6.632718214
		.float -2.1219444005e-4
		.float	0.693359375
		.float	0.4342944925
                .float  3.4028235e+38                ; HUGH_VAL


	.end

