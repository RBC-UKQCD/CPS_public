**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: mcneile $
**  $Date: 2003-06-22 13:34:46 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_lcl_nos/wfm_mdag.asm,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
**  $Id: wfm_mdag.asm,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
**  $Name: not supported by cvs2svn $
**  $Locker:  $
**  $Log: not supported by cvs2svn $
**  Revision 1.2  2001/06/19 18:12:57  anj
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
**  $RCSfile: wfm_mdag.asm,v $
**  $Revision: 1.1.1.1 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_lcl_nos/wfm_mdag.asm,v $
**  $State: Exp $
**
**--------------------------------------------------------------------
****************************************************************************************
*---------------------------------------------------------------------------------------
* _wfm_mdag
* It computes computes chi = M^dag. * psi  for Wilson fermions.
*
*>>>>      float wfm_mdag(float *chi, float *u_eo[2], 
*>>>> 		           float *psi,
*>>>> 			   float *af[4], float *ab[4],
*>>>> 	                   float Kappa,
*>>>> 			   Wilson *wilson_p)
*
* If called with STAND_ALONE=0 then the following registers are expected
* to be loaded as follows:
*	R0 = chi
*	R1 = psi
*
*---------------------------------------------------------------------------------------
****************************************************************************************
	.version	30

	.include "../../include/wilson.hasm"

	.def 	_wfm_mdag

	.ref	_wfm_spproj
	.ref	_wfm_cmat_spproj
	.ref	_wfm_mat_trick
	.ref	_wfm_trick_spproj
	.ref	INV_F30
	.ref	_wfm_trick_kappa

	.if 	SCU_ON = 0
	.ref	_wfm_comm_backward_l
	.ref	_wfm_comm_forward_l
	.endif

	.if 	SCU_ON = 1
	.ref	_wfm_scu_init
	.ref	_wfm_scu_wait
	.ref	_wfm_comm_backward
	.ref	_wfm_comm_forward
	.endif


*---------------------------------------------------------------------------------------
* References
*---------------------------------------------------------------------------------------
	.ref	u0
	.ref	u1
	.ref	wilson_p
	.ref	af
	.ref	af0
	.ref	af1
	.ref	af2
	.ref	af3
	.ref	tas0
	.ref	tas0_p
	.ref	tpsi0_p
	.ref	c_u0_p
	.ref	c_u1_p
	.ref	c_u0
	.ref	c_u1
	.ref	tas1
	.ref	tas1_p
	.ref	tpsi1_p
	.ref	ab
	.ref	ab0
	.ref	ab1
	.ref	ab2
	.ref	ab3
	.ref	ab0_t
	.ref	ab1_t
	.ref	ab2_t
	.ref	ab3_t
	.ref	parameters
	.ref	kappa_sq
	.ref	i_kappa_sq
	.ref	m_kappa_sq
	.ref	direct
	.ref	bank0
	.ref	bank1a
	.ref	bank1b
	.ref	bank2a
	.ref	bank2b
	.ref	bank3a
	.ref	bank3b
	.ref	bank4a
	.ref	bank4b
	.ref	scu_b	
	.ref	cb_cntrl_b

*---------------------------------------------------------------------------------------
* definitions
*---------------------------------------------------------------------------------------
FP	.set	AR3

*---------------------------------------------------------------------------------------
* Local data and pointers
*---------------------------------------------------------------------------------------
	.sect	"T:wfm1"
chi	.space	1
psi	.space	1

****************************************************************************************
* _wfm_mdag
****************************************************************************************
	.text
_wfm_mdag:

*---------------------------------------------------------------------------------------
* If STAND_ALONE = 1: C-callable
*---------------------------------------------------------------------------------------
	.if	STAND_ALONE = 1
	PUSH	FP
	LDI	SP, FP
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
	LDP	@af0
	LDI     *-FP(2), 	R0	
	STI	R0,		@chi
	LDI     *-FP(3), 	AR0
	LDI	*AR0++,		R0
	STI	R0,		@u0
	LDI	*AR0++,		R0
	STI	R0,		@u1
	LDI     *-FP(4), 	R0
	STI	R0,		@psi
	LDI     *-FP(5), 	AR0
	STI	AR0,		@af
	LDI	*AR0++,		R0
	STI	R0,		@af0
	LDI	*AR0++,		R0
	STI	R0,		@af1
	LDI	*AR0++,		R0
	STI	R0,		@af2
	LDI	*AR0++,		R0
	STI	R0,		@af3
	LDI     *-FP(6), 	AR0
	STI	AR0,		@ab
	LDI	*AR0++,		R0
	STI	R0,		@ab0
	LDI	*AR0++,		R0
	STI	R0,		@ab1
	LDI	*AR0++,		R0
	STI	R0,		@ab2
	LDI	*AR0++,		R0
	STI	R0,		@ab3
	LDI     *-FP(8), 	R0
	STI	R0,		@wilson_p
	LDF	*-FP(7),	R1
	MPYF	R1,		R1,		R0
	RND	R0
	STF	R0,		@kappa_sq
	NEGF	R0,		R1
	RND	R1
	STF	R1,		@m_kappa_sq
	CALL	INV_F30
	RND	R0
	STF	R0,		@i_kappa_sq
	.else
	LDP	@chi
	STI	R0,		@chi
	STI	R1,		@psi
	.endif

*---------------------------------------------------------------------------------------
* Routine starts here
*---------------------------------------------------------------------------------------
;>>>> 	   _wfm_scu_init(wilson_p)
	.if 	SCU_ON = 1
	LDI	@wilson_p,	AR0
	CALL	_wfm_scu_init
	.endif

;>>>> 	   _wfm_spproj(af[0], af[1], af[2], af[3], 
;>>>> 		       psi,
;>>>> 	               +1, wilson_p, 1);
	LDI	@psi,		AR6			; psi
	LDF	-1.0,		R5			; -1.0
	LDI	@wilson_p,	AR0			; wilson_p
	LDI	1,		AR1			; cb = 1
	CALL	_wfm_spproj

;>>>> 	   wfm_comm_backward(af0, af1, af2, af3,  wilson_p);
	.if 	SCU_ON = 0
	LDI	@wilson_p,	R0
	PUSH	R0
	LDI	@af3,		R0
	PUSH	R0
	LDI	@af2,		R0
	PUSH	R0
	LDI	@af1,		R0
	PUSH	R0
	LDI	@af0,		R0
	PUSH	R0
	CALL	_wfm_comm_backward_l
	SUBI	5,		SP
	.else
	LDI	@wilson_p,	AR0			; wilson_p
	CALL	_wfm_comm_backward
	.endif

;>>>> 	   wfm_cmat_spproj(ab[0], ab[1], ab[2], ab[3], 
;>>>> 			  u_eo[1], psi, 
;>>>> 			  -1, wilson_p, 1);
	LDI	@psi,		AR6			; psi
	LDI	@u1,		AR7			; u[1]
	LDF	-1.0,		R6			; sign = -1.0
	LDI	@wilson_p,	AR2			; wilson_p
	LDI	1,		AR1			; cb = 1
	CALL	_wfm_cmat_spproj

*---------------------- Synchronization ----------------------------------
;>>>>	scu_wait(); 
	.if 	SCU_ON = 1
	CALL	_wfm_scu_wait
	.endif
*-------------------------------------------------------------------------

;>>>> 	   wfm_comm_forward(ab0, ab1, ab2, ab3, wilson_p);
	.if 	SCU_ON = 0
	LDI	@wilson_p,	R0
	PUSH	R0
	LDI	@ab3,		R0
	PUSH	R0
	LDI	@ab2,		R0
	PUSH	R0
	LDI	@ab1,		R0
	PUSH	R0
	LDI	@ab0,		R0
	PUSH	R0
	CALL	_wfm_comm_forward_l
	SUBI	5,		SP
	.else
	LDI	@wilson_p,	AR0			; wilson_p
	CALL	_wfm_comm_forward
	.endif

;>>>> 	   wfm_mat_trick(chi, u_eo[0], 
;>>>> 			af[0], af[1], af[2], af[3], 
;>>>> 			-1, wilson_p, 1);
	LDI	@chi,		AR6			; chi
	LDI	@u0,		AR7			; u[0]
	LDF	-1.0,		R6			; sign = -1.0
	LDI	@wilson_p,	AR2			; wilson_p
	LDI	1,		AR1			; cb = 1
	CALL	_wfm_mat_trick

*---------------------- Synchronization ----------------------------------
;>>>>	scu_wait(); 
	.if 	SCU_ON = 1
	CALL	_wfm_scu_wait
	.endif
*-------------------------------------------------------------------------

;>>>> 	   wfm_trick_spproj(chi, 
;>>>> 			    af[0], af[1], af[2], af[3], 
;>>>> 			    ab[0], ab[1], ab[2], ab[3], 
;>>>> 			    -1, wilson_p, 1, 0);
	LDI	@chi,		AR4			; chi
	LDF	-1.0,		R6			; sign = -1.0
	LDI	@wilson_p,	AR0			; wilson_p
	LDI	1,		AR7			; cb = 1
	LDI	0,		R7			; skip_spproj=0 
	CALL	_wfm_trick_spproj

;>>>> 	   wfm_comm_backward(af0, af1, af2, af3,  wilson_p);
	.if 	SCU_ON = 0
	LDI	@wilson_p,	R0
	PUSH	R0
	LDI	@af3,		R0
	PUSH	R0
	LDI	@af2,		R0
	PUSH	R0
	LDI	@af1,		R0
	PUSH	R0
	LDI	@af0,		R0
	PUSH	R0
	CALL	_wfm_comm_backward_l
	SUBI	5,		SP
	.else
	LDI	@wilson_p,	AR0			; wilson_p
	CALL	_wfm_comm_backward
	.endif

;>>>> 	   wfm_cmat_spproj(ab[0], ab[1], ab[2], ab[3], 
;>>>> 			  u_eo[0], chi, 
;>>>> 			  -1, wilson_p, 0);
	LDI	@chi,		AR6			; chi
	LDI	@u0,		AR7			; u[0]
	LDF	-1.0,		R6			; sign = -1.0
	LDI	@wilson_p,	AR2			; wilson_p
	LDI	0,		AR1			; cb = 0
	CALL	_wfm_cmat_spproj

*---------------------- Synchronization ----------------------------------
;>>>>	scu_wait(); 
	.if 	SCU_ON = 1
	CALL	_wfm_scu_wait
	.endif
*-------------------------------------------------------------------------

;>>>> 	   wfm_comm_forward(ab0, ab1, ab2, ab3, wilson_p);
	.if 	SCU_ON = 0
	LDI	@wilson_p,	R0
	PUSH	R0
	LDI	@ab3,		R0
	PUSH	R0
	LDI	@ab2,		R0
	PUSH	R0
	LDI	@ab1,		R0
	PUSH	R0
	LDI	@ab0,		R0
	PUSH	R0
	CALL	_wfm_comm_forward_l
	SUBI	5,		SP
	.else
	LDI	@wilson_p,	AR0			; wilson_p
	CALL	_wfm_comm_forward
	.endif

;>>>> 	   wfm_mat_trick(chi, u_eo[1], 
;>>>> 			af[0], af[1], af[2], af[3], 
;>>>> 			-1, wilson_p, 0);
	LDI	@chi,		AR6			; chi
	LDI	@u1,		AR7			; u[1]
	LDF	-1.0,		R6			; sign = -1.0
	LDI	@wilson_p,	AR2			; wilson_p
	LDI	0,		AR1			; cb = 0
	CALL	_wfm_mat_trick

*---------------------- Synchronization ----------------------------------
;>>>>	scu_wait(); 
	.if 	SCU_ON = 1
	CALL	_wfm_scu_wait
	.endif
*-------------------------------------------------------------------------

;>>>> 	   wfm_trick_kappa(chi, psi,
;>>>> 		           ab[0], ab[1], ab[2], ab[3], 
;>>>> 		           -Kappa_sq, 1.0/Kappa_sq,
;>>>> 		           wilson_p, 0);           /* 2 Nc Ns flops */
	LDI	@chi,		AR4			; chi
	LDI	@psi,		AR5			; psi
	LDI	@wilson_p,	AR0			; wilson_p
	LDI	0,		AR7			; cb = 0
	CALL	_wfm_trick_kappa

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


