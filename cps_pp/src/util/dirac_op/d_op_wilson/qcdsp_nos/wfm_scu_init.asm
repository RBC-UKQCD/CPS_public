**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: chulwoo $
**  $Date: 2004-01-13 20:39:45 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_nos/wfm_scu_init.asm,v 1.2 2004-01-13 20:39:45 chulwoo Exp $
**  $Id: wfm_scu_init.asm,v 1.2 2004-01-13 20:39:45 chulwoo Exp $
**  $Name: not supported by cvs2svn $
**  $Locker:  $
**  $Log: not supported by cvs2svn $
**  Revision 1.1.1.1.10.1  2003/11/06 20:24:30  cwj
**  *** empty log message ***
**
**  Revision 1.1.1.1  2003/11/04 05:05:09  chulwoo
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
**  Revision 1.2  2001/06/19 18:13:04  anj
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
**  Revision 1.2  2001/05/25 06:16:07  cvs
**  Added CVS keywords to phys_v4_0_0_preCVS
**
**  $RCSfile: wfm_scu_init.asm,v $
**  $Revision: 1.2 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_nos/wfm_scu_init.asm,v $
**  $State: Exp $
**
**--------------------------------------------------------------------
****************************************************************************************
*---------------------------------------------------------------------------------------
*
* _wfm_scu_init
*
* This routine initializes the stride block-length and number-of-blocks in
* the Serial Communications Unit. Since these values depend on the lattice
* size only, and since the SCU does not change them during operation, they
* only need to be set once.
*
* If STAND_ALONE = 1 it can be called from C as:
*
* wfm_scu_init(Wilson *wilson_p);          ; Wilson struct.
*
* If STAND_ALONE = 0 it can not be called from C. Instead the following is expected
* to be set up before the routine is called:
*
* 1) all the .ref in the STAND_ALONE = 0  case below must be defined
* 2) The following registers must be loaded as follows:
*    AR0 <-- Wilson pointer
*
*---------------------------------------------------------------------------------------
****************************************************************************************

	.version	30

	.include	"../../include/wilson.hasm"
	.include	"wfm_nga_struct.hasm"

	.def _wfm_scu_init

*---------------------------------------------------------------------------------------
* References
*---------------------------------------------------------------------------------------
	.ref	scu_b
	.ref	_wfm_wire_map
	.ref _gjp_local_axis

*---------------------------------------------------------------------------------------
* definitions
*--------------------------------------------------------------------------------------
FP		.set	AR3			; use FP for arguments only
WILSON_AD	.set	AR0
SCUB		.set	AR1
SCU		.set	AR2
WIRE_MAP	.set	AR4
LOCAL_AXIS	.set	AR5	
WIRE		.set	R1
AXIS_FLAG	.set	R2
TMP		.set	R0

*---------------------------------------------------------------------------------------
* allocations
*--------------------------------------------------------------------------------------
	.text
wmap	.word	_wfm_wire_map
laxis	.word	_gjp_local_axis

	
****************************************************************************************
* _wfm_scu_init
****************************************************************************************
	.text
_wfm_scu_init:

*---------------------------------------------------------------------------------------
* If STAND_ALONE = 1: C-callable
*---------------------------------------------------------------------------------------
	.if	STAND_ALONE = 1
	PUSH    FP
	LDI     SP,FP
*  Save all registers that are important to C
	PUSH    R4
	PUSH    R5
	RND	R6
	PUSHF   R6	
	RND	R7
	PUSHF   R7
	PUSH    AR4
	PUSH    AR5
	PUSH    AR6
	PUSH    AR7
	PUSH    FP              ; Local frame pointer
	PUSH    DP
*  Load arguments from stack to registers
	LDI     *-FP(2), 	WILSON_AD		; Wilson structure
	.endif

	PUSH    DP
*---------------------------------------------------------------------------------------
* Routine starts here
*---------------------------------------------------------------------------------------
	LDP	@wmap
	LDI	@wmap,				WIRE_MAP	; wire map 

	LDP	@laxis
	LDI	@laxis,				LOCAL_AXIS	; local axis 

	LDP	@scu_b
	LDI	@scu_b,				SCUB		; SCU base address

* Set maximum number of words before access to DRAM is relinquished
	LDI	4,				TMP
	STI	TMP,			*+SCUB(scu.max_dram_words_wire) 	; GRPLEN = 4


* Set up communication parameters by initializing the SCU DMA registers. 
* wilson_p->comm[mu] contains a 32 bit word created by packaging
* the stride, number of blocks, and block length: 
* numblk[10bits]-blklen[10bits]-stride[12bits]

	LDI	*+LOCAL_AXIS(0),		AXIS_FLAG
	BNZ	SETX
	LDI	*+WILSON_AD(Wilson.comm+0),	TMP
	LDI	*+WIRE_MAP(0),			WIRE
	ADDI	WIRE,	SCUB,			SCU
	STI	TMP,				*+SCU(scu.dma_cntrl0) ; +X
	LDI	*+WIRE_MAP(1),			WIRE
	ADDI	WIRE,	SCUB,			SCU
	STI	TMP,				*+SCU(scu.dma_cntrl0) ; -X
SETX:	NOP
	
	LDI	*+LOCAL_AXIS(1),		AXIS_FLAG
	BNZ	SETY
	LDI	*+WILSON_AD(Wilson.comm+1),	TMP
	LDI	*+WIRE_MAP(2),			WIRE
	ADDI	WIRE,	SCUB,			SCU
	STI	TMP,				*+SCU(scu.dma_cntrl0) ; +Y
	LDI	*+WIRE_MAP(3),			WIRE
	ADDI	WIRE,	SCUB,			SCU
	STI	TMP,				*+SCU(scu.dma_cntrl0) ; -Y
SETY:	NOP

	LDI	*+LOCAL_AXIS(2),		AXIS_FLAG
	BNZ	SETZ
	LDI	*+WILSON_AD(Wilson.comm+2),	TMP
	LDI	*+WIRE_MAP(4),			WIRE
	ADDI	WIRE,	SCUB,			SCU
	STI	TMP,				*+SCU(scu.dma_cntrl0) ; +Z
	LDI	*+WIRE_MAP(5),			WIRE
	ADDI	WIRE,	SCUB,			SCU
	STI	TMP,				*+SCU(scu.dma_cntrl0) ; -Z
SETZ:	NOP

	LDI	*+LOCAL_AXIS(3),		AXIS_FLAG
	BNZ	SETT
	LDI	*+WILSON_AD(Wilson.comm+3),	TMP
	LDI	*+WIRE_MAP(6),			WIRE
	ADDI	WIRE,	SCUB,			SCU
	STI	TMP,				*+SCU(scu.dma_cntrl0) ; +T
	LDI	*+WIRE_MAP(7),			WIRE
	ADDI	WIRE,	SCUB,			SCU
	STI	TMP,				*+SCU(scu.dma_cntrl0) ; -T
SETT:	NOP

	POP	DP
*-------------------------------------------------------------------------
*
*  Restore the registers before returning to the C program                
*
	.if	STAND_ALONE = 1
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
	.endif

	RETS
	.end
