#include<config.h>
CPS_START_NAMESPACE
**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: mcneile $
**  $Date: 2003-06-22 13:34:47 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/various/overlay/test.asm,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
**  $Id: test.asm,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
**  $Name: not supported by cvs2svn $
**  $Locker:  $
**  $Log: not supported by cvs2svn $
**  Revision 1.2  2001/06/19 18:12:31  anj
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
**  Revision 1.2  2001/05/25 06:16:04  cvs
**  Added CVS keywords to phys_v4_0_0_preCVS
**
**  $RCSfile: test.asm,v $
**  $Revision: 1.1.1.1 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/various/overlay/test.asm,v $
**  $State: Exp $
**
**--------------------------------------------------------------------
	.version	30

	.def	_test
	.def	tmpxxx

	.ref 	_out_buf

FP	.set	AR3

	.sect	"T:wfm1"
tmpxxx	.word	tmp1
tmp1	.word	1


	.sect	"T:wfm1"
_test:

*---------------------------------------------------------------------------------------
	PUSH	FP
	LDI	SP, FP
*  Save all registers that are important to C
	PUSH    R4
	PUSH    R5
	PUSHF   R6
	PUSHF   R7
	PUSH    AR4
	PUSH    AR5
	PUSH    AR6
	PUSH    AR7
	PUSH    FP              ; Local frame pointer
	PUSH    DP
	PUSH	ST		; Save the status register

	LDI	1010,	R0

	LDP	_out_buf,	AR7
	LSH	16,AR7
	OR	_out_buf,	AR7

	LDI	444h,		R0
	STI	R0,		*AR7++
	LDI	333h,		R0
	STI	R0,		*AR7++

	LDP	@tmpxxx
	LDI	@tmpxxx,	R0
	STI	R0,		*AR7++


*---------------------------------------------------------------------------------------
*  Restore the registers before returning to the C program                
*---------------------------------------------------------------------------------------
	POP	ST		; Restore the status register
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
	RETS
	.end



CPS_END_NAMESPACE
