**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: chulwoo $
**  $Date: 2004-01-13 20:38:59 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_ghb/qcdsp/mtxfast.asm,v 1.2 2004-01-13 20:38:59 chulwoo Exp $
**  $Id: mtxfast.asm,v 1.2 2004-01-13 20:38:59 chulwoo Exp $
**  $Name: not supported by cvs2svn $
**  $Locker:  $
**  $Log: not supported by cvs2svn $
**  Revision 1.1.1.1.10.1  2003/11/06 00:10:38  cwj
**  *** empty log message ***
**
**  Revision 1.1.1.1  2003/11/04 05:04:57  chulwoo
**
**  starting again
**
**
**  Revision 1.1.1.1  2003/06/22 13:34:45  mcneile
**  This is the cleaned up version of the Columbia Physics System.
**  The directory structure has been changed.
**  The include paths have been updated.
**
**
**  Revision 1.2  2001/06/19 18:11:24  anj
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
**  Revision 1.2  2001/05/25 06:15:59  cvs
**  Added CVS keywords to phys_v4_0_0_preCVS
**
**  $RCSfile: mtxfast.asm,v $
**  $Revision: 1.2 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_ghb/qcdsp/mtxfast.asm,v $
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
	
************************************************************
************************************************************

*	.def	src_mfast
*	.def	end_mfast

*	.sect	".mtxfast"
	.text
*	.label	src_mfast
************************************************************

	.def	_m_add
_m_add
*  This routine is called as m_add( mtx3, mtx1, mtx2 )
*  so the stack should look like:
*
*  SP-4		mtx2
*  SP-3		mtx1
*  SP-2		mtx3
*  SP-1		(return address)
*  SP		AR3  	(here called FP)

*  Push the old value of AR3 onto the stack
	PUSH	FP
*  Load the Stack Pointer into FP
	LDI	SP, FP
*  Save all registers that are important to C

	LDI	*-FP(3), MX1
	LDI	*-FP(4), MX2
	LDI	*-FP(2), MX3

	LDI	17, RC
	RPTB	addloop
	LDF	*MX1++, R0
||	LDF	*MX2++, R1
	ADDF3	R0, R1, R2
	.if	ROUND
	RND	R2
	.endif
addloop STF	R2, *MX3++


*
*  Clean up to leave
        POP     FP              ;  This is our caller's frame pointer

        RETS


************************************************************

	.def	_m_multiply3
_m_multiply3
*  This routine is called as m_multiply3( mtx3, mtx1, mtx2 )
*  so the stack should look like:
*
*  SP-4		mtx2
*  SP-3		mtx1
*  SP-2		mtx3
*  SP-1		(return address)
*  SP		AR3  	(here called FP)

*stoptm:	PUSH	AR0		;  Stop the timer.
*	PUSHF	R1
*	PUSHF	R0
*	LDI	808h, AR0
*	LSH	12, AR0
*	ADDI	20h, AR0
*	LDI	*AR0, R0
*	LDI	80h, R1
*	NOT	R1
*	AND	R1, R0
*	STF	R0, *AR0
*	POPF	R0
*	POPF	R1
*tmstop:	POP	AR0		;  Timer is stopped

*  Push the old value of AR3 onto the stack
	PUSH	FP
*  Load the Stack Pointer into FP
	LDI	SP, FP
*  Save all registers that are important to C
        PUSH    AR4
        PUSH    AR5
        PUSH    AR6
	
*
	LDI	3, IR0
*	
	LDI	*-FP(3), MX1r
	LDI	*-FP(3), MX1i
	ADDI	9, MX1i

	LDI	*-FP(2), MX3r
	LDI	*-FP(2), MX3i
	ADDI	9, MX3i

	LDI	8, RC
	RPTB	END_MULTIPLY_LOOP
BEGIN_MULTIPLY_LOOP

	CMPI	8, RC
	BZ	_row_zero_only
;  3 cycles lost if branch 
	CMPI	5, RC
	BZ	_row_one_or_two_only
;  3 cycles lost if branch 
	CMPI	2, RC
	BZ	_row_one_or_two_only
;  3 cycles lost if branch 

_all_others
	MPYF3	*MX1r, *++MX2r, R0
||	ADDF3	R1, R2, R2
	LDF	*++MX2r(IR0), R3
	.if	ROUND
	RND	R2
	.endif
	BD	_common
	MPYF3	*++MX1r, R3, R1
||	STF	R2, *MX3i++
;  1 cycle lost to execute only delay
	MPYF3	*++MX1r, *++MX2r(IR0), R0
||	ADDF3	R0, R1, R2
	MPYF3	*MX1i, *++MX2i, R1
||	ADDF3	R0, R2, R2	
	
_row_zero_only
	LDI	*-FP(4), MX2r
	LDI	*-FP(4), MX2i
	ADDI	9, MX2i

;  2 cycle register conflict
	MPYF3	*MX1r, *MX2r, R0
	BD	_common
	MPYF3	*++MX1r, *++MX2r(IR0), R1
	MPYF3	*++MX1r, *++MX2r(IR0), R0
||	ADDF3	R0, R1, R2
	MPYF3	*MX1i, *MX2i, R1
||	ADDF3	R0, R2, R2	


_row_one_or_two_only
	LDI	*-FP(4), MX2r
	LDI	*-FP(4), MX2i
	ADDI	9, MX2i

;  2 cycle register conflict
	MPYF3	*++MX1r(IR0), *MX2r, R0
||	ADDF3	R1, R2, R2
	LDF	*++MX2r(IR0), R3
	.if	ROUND
	RND	R2
	.endif
	BD	_common
	MPYF3	*++MX1r, R3, R1
||	STF	R2, *MX3i++
;  1 cycle lost to execute only delay
	MPYF3	*++MX1r, *++MX2r(IR0), R0
||	ADDF3	R0, R1, R2
	MPYF3	*++MX1i(IR0), *MX2i, R1
||	ADDF3	R0, R2, R2	



_common
	MPYF3	*++MX1i, *++MX2i(IR0), R0
||	SUBF3	R1, R2, R2
	MPYF3	*++MX1i, *++MX2i(IR0), R1
||	SUBF3	R0, R2, R2
	MPYF3	*MX1i, *MX2r, R0
||	SUBF3	R1, R2, R2
	LDF	*MX2i, R3
	.if	ROUND
	RND	R2
	.endif
	MPYF3	*MX1r, R3, R1
||	STF	R2, *MX3r++
;  1 cycle lost to execute only delay
	MPYF3	*--MX1i, *--MX2r(IR0), R0
||	ADDF3	R0, R1, R2
	MPYF3	*--MX1r, *--MX2i(IR0), R1
||	ADDF3	R0, R2, R2
	MPYF3	*--MX1i, *--MX2r(IR0), R0
||	ADDF3	R1, R2, R2
END_MULTIPLY_LOOP
	MPYF3	*--MX1r, *--MX2i(IR0), R1
||	ADDF3	R0, R2, R2


*
	ADDF3	R1, R2, R2
	.if	ROUND
	RND	R2
	.endif
	STF	R2, *MX3i++
*
*  Clean up to leave
        POP     AR6
        POP     AR5
        POP     AR4
 
        POP     FP              ;  This is our caller's frame pointer

        RETS

************************************************************

	.end
