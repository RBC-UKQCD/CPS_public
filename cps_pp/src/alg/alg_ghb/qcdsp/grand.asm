**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: chulwoo $
**  $Date: 2004-01-13 20:38:59 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_ghb/qcdsp/grand.asm,v 1.2 2004-01-13 20:38:59 chulwoo Exp $
**  $Id: grand.asm,v 1.2 2004-01-13 20:38:59 chulwoo Exp $
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
**  $RCSfile: grand.asm,v $
**  $Revision: 1.2 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_ghb/qcdsp/grand.asm,v $
**  $State: Exp $
**
**--------------------------------------------------------------------

	.asg	AR3, FP			; FP is AR3

	.def	_grand
	.ref	_fprint2
	.ref	_el_seed_p

*	.sect	".grand"
	.text
_grand:
	PUSH	DP

	LDP	@_el_seed_p		; set DP to _el_seed_p data page
	LDI	@_el_seed_p, AR0	; set AR0 to _el_seed_p

	LDI	*+AR0(1), R1
	LDI	*+AR0(0), R0
	ADDI	*+AR0(4), R0

	SUBRI	R1, R0
	LDIU	0, AR2			; new carry if no borrow 
	BNC	store_carry

	LDIU	1, AR2			; new carry if borrow
	SUBI	18, R0			; mod 2^32 - 18

store_carry:
	LDI	*+AR0(2), RC
	STI	R1, *+AR0(0)		; _el_seed11 = _el_seed12
	STI	AR2, *+AR0(4)
	STI	RC, *+AR0(1)		; _el_seed12 = _el_seed13
	STI	R0, *+AR0(2)		; _el_seed13 = s

*****
* next section based on MPY_I30 using:
*  69069*n = (3533*(n>>16) + 1*(n&0xFFFF))<<16 + 3533*(n&0xFFFF)
*****

	LDI	*+AR0(3), R1
	LDI	R1, AR1			; n
	LSH	-16, R1			; 16 MSBs of n
	AND	0FFFFh, AR1		; 16 LSBs of n
	LDI	AR1, R2			; 1*(n&0xFFFF)
	MPYI	3533, R1		; 3533*(n>>16)
	MPYI	3533, AR1		; 3533*(n&0xFFFF)
	ADDI	R2, R1
	LSH	16, R1			;
	ADDI	AR1, R1			;final result n = 69069*n

	LDI	*+AR0(5), R2
	ADDI	R2, R1			; n = 69069*n + lconst
	STI	R1, *+AR0(3)		; _el_seed21 = n
	ADDI	R1, R0			; combine generators
	XOR	R1, R1
	PUSH	R1
	POPF	R1
*	ABSI	R0, R0			; rem for [-1,1]
	OR	R0, R1
	NORM	R1, R0
	RND	R0, R0

	POP	DP
	RETS

	.end
