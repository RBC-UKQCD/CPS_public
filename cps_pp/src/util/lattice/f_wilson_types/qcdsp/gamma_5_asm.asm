**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: chulwoo $
**  $Date: 2004-01-13 20:39:51 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/f_wilson_types/qcdsp/gamma_5_asm.asm,v 1.2 2004-01-13 20:39:51 chulwoo Exp $
**  $Id: gamma_5_asm.asm,v 1.2 2004-01-13 20:39:51 chulwoo Exp $
**  $Name: not supported by cvs2svn $
**  $Locker:  $
**  $Log: not supported by cvs2svn $
**  Revision 1.1.1.1.10.1  2003/11/06 20:24:31  cwj
**  *** empty log message ***
**
**  Revision 1.1.1.1  2003/11/04 05:05:14  chulwoo
**
**  starting again
**
**
**  Revision 1.1.1.1  2003/06/22 13:34:47  mcneile
**  This is the cleaned up version of the Columbia Physics System.
**  The directory structure has been changed.
**  The include paths have been updated.
**
**
**  Revision 1.2  2001/06/19 18:13:24  anj
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
**  Revision 1.2  2001/05/25 06:16:09  cvs
**  Added CVS keywords to phys_v4_0_0_preCVS
**
**  $RCSfile: gamma_5_asm.asm,v $
**  $Revision: 1.2 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/f_wilson_types/qcdsp/gamma_5_asm.asm,v $
**  $State: Exp $
**
**--------------------------------------------------------------------
	.version 30
;=============================================================================
; This file is   : gamma_5.asm						     |
; generated from : gamma_5						     |
; using          : Tartan C/C++ Compiler for the TMS320C3x/4x, v2.1.1	     |
; on (M/D/Y H:M) : 04/09/1998 11:31					     |
;----------------------------------------------------------------------------|
;									     |
; GLOBAL ATTRIBUTES AND ASSERTIONS					     |
;   Wait states            : Code = 0, Stack = 0, Data page = 0, Heap = 0    |
;   Max loop iterations    : 2**23 (no huge loops)			     |
;   Code Span              : 24-bits of address				     |
;									     |
; GLOBAL OPTIONS							     |
;   Optimization heuristics: space					     |
;   Debug information      : off					     |
;   Profiling              : off					     |
;   Max RPTS allowed       : 2**31					     |
;   Data Page              : No_Page					     |
;   Data Page is in ROM    : false					     |
;   *ARn same mem as *ARm? : true					     |
;									     |
;=============================================================================

;
	.def	_gamma_5__FPfT1i
	.def	_gamma_5__FPfT1i$LAJ
;=============================================================================
; LOCAL OPTIONS								     |
;   Calling convention     : TI C stack parameters			     |
;									     |
; GENERATED CODE PROPERTIES						     |
;   Total words of code    : 24						     |
;   Volatile registers used: R0,R1,R2,R3,AR0,AR1,AR2,DP,IR1,BK,RS,RE,RC      |
;   Registers for locals   : R2		holds site			     |
;			     R3		holds comp			     |
;			     AR1	holds p_out			     |
;			     AR2	holds p_in			     |
;   Parameters             : AR0	holds v_out			     |
;			     IR1	holds v_in			     |
;			     R0		holds num_sites			     |
;   Stack frame            : quick (AR3 points to some old frame)	     |
;=============================================================================

	.sect	"T:RTN1"

_gamma_5__FPfT1i:
	POP	BK
_gamma_5__FPfT1i$LAJ:
	POP	AR0
	POP	IR1
	POP	R0
	ADDI	3,SP

	BRD	L7
	LDIU	AR0,AR1
	LDIU	IR1,AR2
	LDIU	0,R2

; begin loop 1-1
L1:	LDIU	0,R3
	BR	L8

; begin loop 2-2
L3:	ADDI	1,R3
	LDIU	*AR2++(1),R1
	STI	R1,*AR1++(1)
L8:	CMPI	11,R3
	BLE	L3
; end loop 2-2

	LDIU	11,RC
	RPTB	L16

; begin loop 3-2
	NEGF	*AR2++(1),R1
	RND	R1,R1		; GRF
L16:	STF	R1,*AR1++(1)
; end loop 3-2

	ADDI	1,R2
L7:	CMPI	R2,R0
	BGT	L1
; end loop 1-1


	BU	BK

;
	.def	_CWInitgamma_5$src
	.def	_CWInitgamma_5$src$LAJ
;=============================================================================
; LOCAL OPTIONS								     |
;   Calling convention     : Tartan register parameters			     |
;									     |
; GENERATED CODE PROPERTIES						     |
;   Total words of code    : 3						     |
;   Volatile registers used: R0,DP,BK					     |
;   Stack frame            : quick (AR3 points to some old frame)	     |
;=============================================================================

	.sect	"T:RTN2"

_CWInitgamma_5$src:
	POP	BK
_CWInitgamma_5$src$LAJ:

	LDIU	0,R0

	BU	BK

;


;=============================================================================
; Total words of code = 27						     |
; Total words of data = 0						     |
;=============================================================================

	.end

