**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: zs $
**  $Date: 2004-08-18 11:57:58 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_lcl_nos/wilson_dslash.asm,v 1.4 2004-08-18 11:57:58 zs Exp $
**  $Id: wilson_dslash.asm,v 1.4 2004-08-18 11:57:58 zs Exp $
**  $Name: not supported by cvs2svn $
**  $Locker:  $
**  $Revision: 1.4 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_lcl_nos/wilson_dslash.asm,v $
**  $State: Exp $
**
**--------------------------------------------------------------------
****************************************************************************************
*---------------------------------------------------------------------------------------
*
* _wilson_dslash
*
* The Wilson fermion Dslash routine:
*
* chi = Dslash * psi, 
*
* void wilson_dslash(float  *chi,          chi = Dslash psi          
* 		     float  *u,            Gauge field                 
* 		     float  *psi,          chi = Dslash psi          
* 		     int    cb,		   checker board 0/1 -> even/odd
* 		     int    dag,           dagger 0/1 -> Dslash/Dslash^dagger
* 		     Wilson *wilson_p);    pointer to a Wilson struct. 
*
*---------------------------------------------------------------------------------------
****************************************************************************************

	.version	30

	.include	"../../include/wilson.hasm"

	.def	_wilson_dslash

	.ref	_wfm_dslash
        .ref    wfm_s_cb0
        .ref    wfm_r_cb0

*---------------------------------------------------------------------------------------
* References
*---------------------------------------------------------------------------------------
	.ref	u0
	.ref	u1
	.ref	wilson_p
	.ref	af
	.ref	ab

*---------------------------------------------------------------------------------------
* definitions
*---------------------------------------------------------------------------------------
AFBPT	.set	AR4
FP	.set	AR3
WILSONP .set	AR2
AF	.set	AR0
AB	.set	AR1

*---------------------------------------------------------------------------------------
* Reserve space for local data
*---------------------------------------------------------------------------------------
	.sect	"T:wfm1"
chi	.space	1
psi	.space	1
cb	.space	1
sign	.space  1

* Cache stuff
cf	.word	0010000000000b		; cache freeze mask 
ce	.word	0100000000000b		; cache enable mask
cc	.word	1000000000000b		; cache clear mask

******************************************************
* FUNCTION DEF : _wilson_mdagm
******************************************************
	.text
_wilson_dslash:

*---------------------------------------------------------------------------------------
* C-calling conventions and initial argument manipulations
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
*  Load arguments from stack to registers, and do necessary manipulations
	LDP	@psi
	LDI     *-FP(2), 	R0	
	STI	R0,		@chi		; chi
	LDI     *-FP(4), 	R0
	STI	R0,		@psi		; psi
	LDI     *-FP(5), 	R0
	STI	R0,		@cb		; cb
	LDI	*-FP(6),	R1		; dag
;>>>>	if dag=0 => sign=+1; if cb=1 => sign=-1
	CMPI	0,		R1
	BEQ	arg0
	CMPI	1,		R1
	BEQ	arg1
arg0:	LDF	1.0,		R0	
	BU	arg2
arg1:	LDF	-1.0,		R0	
	BU	arg2
arg2:	STF	R0,		@sign
	LDI     *-FP(7), 	WILSONP		; wilson pointer
	STI	WILSONP,	@wilson_p
;>>>> 	   u_eo[0] = u;
;>>>> 	   u_eo[1] = u + GAUGE_SIZE * wilson_p->vol[0];
	LDP	@u0
	LDI	*-FP(3),	R0		; gauge field u
	STI	R0,		@u0
	LDI	GAUGE_SIZE,	R1
	MPYI	*+WILSONP(Wilson.vol),		R1
	ADDI	R0,		R1
	STI	R1,		@u1

*---------------------------------------------------------------------------------------
* Enable instruction cache
*---------------------------------------------------------------------------------------
	LDP	@psi				; set DP to CRAM
	ANDN	@cf, ST				; clear CF -- cache not frozen
	OR	@cc, ST				; set CC -- cache initially cleared
	OR	@ce, ST				; set CE -- cache enabled
;;	ANDN	@ce, ST				; set CE -- cache disabled

*---------------------------------------------------------------------------------------
* Set the pointers of the temporary arrays
*---------------------------------------------------------------------------------------
	LDP	@af
	LDI	@af,		AF
	LDI	WILSONP,	AFBPT
	ADDI	Wilson.af,	AFBPT
	LDI	ND,		RC
	SUBI	1,		RC
	RPTB	afdir
	LDI	*AFBPT++,	R0
afdir:	STI	R0,		*AF++		; set the half spinor af[i] pointer

	LDI	@ab,		AB
	LDI	WILSONP,	AFBPT
	ADDI	Wilson.ab,	AFBPT
	LDI	ND,		RC
	SUBI	1,		RC
	RPTB	abdir
	LDI	*AFBPT++,	R0
abdir:	STI	R0,		*AB++		; set the half spinor ab[i] pointer

*---------------------------------------------------------------------------------------
* call wfm_dslash
*---------------------------------------------------------------------------------------
        CALL    wfm_s_cb0                       ; save contents of cbuf reg 0
	LDP	@psi				; set DP to CRAM
	
;>>>> 	   _wfm_dslash(chi, u_eo, psi, cb, sign, wilson_p);
	LDP	@psi
	LDI	@chi,		R0
	LDI	@psi,		R1
	LDI	@cb,		R2
	LDF	@sign,		R3
	CALL	_wfm_dslash

	CALL    wfm_r_cb0                       ; restore contents of cbuf reg 0

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
