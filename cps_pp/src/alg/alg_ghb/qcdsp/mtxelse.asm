**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: zs $
**  $Date: 2004-08-18 11:57:39 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_ghb/qcdsp/mtxelse.asm,v 1.4 2004-08-18 11:57:39 zs Exp $
**  $Id: mtxelse.asm,v 1.4 2004-08-18 11:57:39 zs Exp $
**  $Name: not supported by cvs2svn $
**  $Locker:  $
**  $Revision: 1.4 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_ghb/qcdsp/mtxelse.asm,v $
**  $State: Exp $
**
**--------------------------------------------------------------------

	.asg	1,	ROUND

FP	.set	AR3

MX1r	.set	AR0
MX2r	.set	AR1
MX3r	.set	AR5
MX1i	.set	AR2
MX2i	.set	AR4
MX3i	.set	AR6

MX1	.set	AR0
MX2	.set	AR1
MX3	.set	AR2
	
*************************************************************

	.sect	".mtxelse"
*
	.global	_m_conjugate
	.global _m_conjugate__FP6rfloat

_m_conjugate
_m_conjugate__FP6rfloat
*  This routine is called as m_conjugate( mtx1 )
*  so the stack should look like:
*
*  SP-2		mtx1
*  SP-1		(return address)
*  SP		AR3  	(here called FP)

*  Push the old value of AR3 onto the stack
	PUSH	FP
*  Load the Stack Pointer into FP
	LDI	SP, FP
*  Save all registers that are important to C

	LDI	*-FP(2), MX1

	LDI	2, IR0
	LDI	3, IR1
*	Do Imaginary parts of diagonal elements.
	LDF	*+MX1(9), R0
	NEGF	R0
	STF	R0, *+MX1(9)
	LDF	*+MX1(13), R0
	NEGF	R0
	STF	R0, *+MX1(13)
	LDF	*+MX1(17), R0
	NEGF	R0
	STF	R0, *+MX1(17)

*	Do real parts of off diagonal elements
	LDI	MX1, MX2
	ADDI	5, MX2

	LDF	*+MX1(1), R0
||	LDF	*+MX1(IR1), R1
	STF	R0, *+MX1(IR1)
||	STF	R1, *+MX1(1)
	LDF	*+MX1(IR0), R0
||	LDF	*+MX2(1), R1
	STF	R0, *+MX2(1)
||	STF	R1, *+MX1(IR0)
	LDF	*+MX2(0), R0
||	LDF	*+MX2(IR0), R1
	STF	R0, *+MX2(IR0)
||	STF	R1, *+MX2(0)

*	Do imaginary parts of off diagonal elements
	ADDI	10, MX1
	ADDI	9, MX2

	LDF	*+MX1(0), R0
||	LDF	*+MX1(IR0), R1
	NEGF	R0
	NEGF	R1
	STF	R0, *+MX1(IR0)
||	STF	R1, *+MX1(0)
	LDF	*+MX1(1), R0
||	LDF	*+MX2(1), R1
	NEGF	R0
	NEGF	R1
	STF	R0, *+MX2(1)
||	STF	R1, *+MX1(1)
	LDF	*+MX2(0), R0
||	LDF	*+MX2(IR0), R1
	NEGF	R0
	NEGF	R1
	STF	R0, *+MX2(IR0)
||	STF	R1, *+MX2(0)

*
*  Clean up to leave
        POP     FP              ;  This is our caller's frame pointer

        RETS

*********************************************************

	.global	_m_equal

_m_equal
*  This routine is called as m_equal( mtx2, mtx1 )
*  so the stack should look like:
*
*  SP-3		mtx1
*  SP-2		mtx2
*  SP-1		(return equalress)
*  SP		AR3  	(here called FP)

*  Push the old value of AR3 onto the stack
	PUSH	FP
*  Load the Stack Pointer into FP
	LDI	SP, FP
*  Save all registers that are important to C
*	PUSH	AR4
*        PUSH    DP

	LDI	*-FP(3), MX1	;	mtx1
	LDI	*-FP(2), MX2	;	mtx2 = mtx1

*	PUSH	MX1
*	CALL	_iprint
*	SUBI	1, SP
*	LDI	*-FP(2), MX2	;	mtx2 = mtx1
*	PUSH	MX2
*	CALL	_iprint
*	SUBI	1, SP
*	LDI	*-FP(3), MX1	;	mtx1
*	LDI	*-FP(2), MX2	;	mtx2 = mtx1

	LDF	*MX1++, R0
	RPTS	16
	LDF	*MX1++, R0
||	STF	R0, *MX2++
	STF	R0, *MX2++

*
*  Clean up to leave
*        POP     DP
*	POP	AR4
        POP     FP              ;  This is our caller's frame pointer

        RETS

************************************************************

	.global	_m_identity

_m_identity
*  This routine is called as m_identity( mtx1 )
*  so the stack should look like:
*
*  SP-2		mtx1
*  SP-1		(return equalress)
*  SP		AR3  	(here called FP)

*  Push the old value of AR3 onto the stack
	PUSH	FP
*  Load the Stack Pointer into FP
	LDI	SP, FP
*  Save all registers that are important to C
*        PUSH    DP

	LDI	*-FP(2), MX1

	LDF	0., R0
	LDF	1., R1

	RPTS	17
	STF	R0, *MX1++

	STF	R1, *-MX1(18)
	STF	R1, *-MX1(14)
	STF	R1, *-MX1(10)

*
*  Clean up to leave
*        POP     DP
        POP     FP              ;  This is our caller's frame pointer

        RETS

************************************************************


*  Define some handy names for these things:
FP	.set	AR3
	.global	_m_fetch

	.asg	AR0,	dma_base
	.asg	*-FP(4),	length
	.asg	*-FP(3),	source
	.asg	*-FP(2),	destination

	.asg	0,	Control
	.asg	4,	Source_addr
	.asg	6,	Destination_addr
	.asg	8,	Counter

_m_fetch
*  Push the old value of AR3 onto the stack
	PUSH	FP
*  Load the Stack Pointer into FP
	LDI	SP, FP
*  Save all registers that are important to C

**************************************************************

 	LDI	source, AR0
 	LDI	destination, AR1
 	LDI	length, AR2

	SUBI	1, AR2

	LDI	*AR0++, R0
	RPTS	AR2
	LDI	*AR0++,R0
||	STI	R0, *AR1++

*  Restore the registers before returning to the C program
*
fin:	POP     FP              ;  This is our caller's frame pointer

        RETS

**************************************************************
        .end


