**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: zs $
**  $Date: 2004-08-18 11:58:01 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_nos_hdw_diag/wilson_mdag.asm,v 1.4 2004-08-18 11:58:01 zs Exp $
**  $Id: wilson_mdag.asm,v 1.4 2004-08-18 11:58:01 zs Exp $
**  $Name: not supported by cvs2svn $
**  $Locker:  $
**  $Revision: 1.4 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_nos_hdw_diag/wilson_mdag.asm,v $
**  $State: Exp $
**
**--------------------------------------------------------------------
****************************************************************************************
*---------------------------------------------------------------------------------------
*
* _wilson_mdag
*
* The Wilson fermion Mdag routine:
*
* chi = [M^dag] * psi, 
*
* void wilson_mdag(float  *chi,          chi = Mdag(u) psi          
* 		    float  *u,            Gauge field                 
* 		    float  *psi,          chi = Mdag(u) psi          
* 		    float  Kappa,         Wilson's kappa parameter    
* 		    Wilson *wilson_p);    pointer to a Wilson struct. 
*
*---------------------------------------------------------------------------------------
****************************************************************************************

	.version	30

	.include	"../../include/wilson.hasm"

	.def	_wilson_mdag

	.ref	_wfm_mdag
	.ref	INV_F30
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
	.ref	kappa_sq
	.ref	i_kappa_sq
	.ref	m_kappa_sq

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

* Cache stuff
cf	.word	0010000000000b		; cache freeze mask 
ce	.word	0100000000000b		; cache enable mask
cc	.word	1000000000000b		; cache clear mask

******************************************************
* FUNCTION DEF : _wilson_mdag
******************************************************
	.text
_wilson_mdag:

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
	STI	R0,		@chi
	LDI     *-FP(4), 	R0
	STI	R0,		@psi
;>>>> 	   float  Kappa_sq = Kappa * Kappa;
	LDF	*-FP(5),	R1
	MPYF	R1,		R1,		R0
	STF	R0,		@kappa_sq
	NEGF	R0,		R1
	STF	R1,		@m_kappa_sq
	CALL	INV_F30
	LDP	@psi				; set DP to CRAM
	RND	R0
	STF	R0,		@i_kappa_sq
	LDI     *-FP(6), 	WILSONP
	STI	WILSONP,	@wilson_p
;>>>> 	   u_eo[0] = u;
;>>>> 	   u_eo[1] = u + GAUGE_SIZE * wilson_p->vol[0];
	LDP	@u0
	LDI	*-FP(3),	R0
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
* call mdag
*---------------------------------------------------------------------------------------
        CALL    wfm_s_cb0                       ; save contents of cbuf reg 0
	LDP	@psi				; set DP to CRAM
	
;>>>> 	   _wfm_mdag(chi, u_eo, psi, af, ab, Kappa, wilson_p);
	LDP	@psi
	LDI	@chi,		R0
	LDI	@psi,		R1
	CALL	_wfm_mdag

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
