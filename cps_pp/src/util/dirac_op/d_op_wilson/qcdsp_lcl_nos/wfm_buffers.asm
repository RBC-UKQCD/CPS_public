**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: zs $
**  $Date: 2004-08-18 11:57:57 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_lcl_nos/wfm_buffers.asm,v 1.4 2004-08-18 11:57:57 zs Exp $
**  $Id: wfm_buffers.asm,v 1.4 2004-08-18 11:57:57 zs Exp $
**  $Name: not supported by cvs2svn $
**  $Locker:  $
**  $Revision: 1.4 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_lcl_nos/wfm_buffers.asm,v $
**  $State: Exp $
**
**--------------------------------------------------------------------
****************************************************************************************
*---------------------------------------------------------------------------------------
*
* buffers.asm
*
* Several pointers are allocated memory and made globably available 
* with .def. These pointers do not have a leading _ so they can
* not be confused with externals referenced by a C program.
* These pointers are needed by the routines called from mdagm, m, mdag, dslash.
*
*---------------------------------------------------------------------------------------
****************************************************************************************
	.version	30

	.include	"../../include/wilson.hasm"

*
* Space for CRAM temporaries
*
;;; mdagm
	.sect	"T:wfm1"
	.def	u0
	.def	u1
	.def	mp_sq_p
	.def	wilson_p
u0	.space	1
u1	.space	1
mp_sq_p	.space	1
wilson_p .space	1

;;; spproj
	.sect	"T:wfm1"
* Reserve space for af
	.def	af
	.def	af0
	.def	af1
	.def	af2
	.def	af3
af	.word	af0
af0	.space	1
af1	.space	1
af2	.space	1
af3	.space	1
* Reserve space for the temporary CRAM af

	.sect	"T:wfm1"
	.def	tas0
	.def	tas0_p
	.def	tpsi0_p
tpsi0_p	.word	tas0
tas0_p	.word	tas0
tas0	.space	AB_SIZE

;;; cmat_spproj
	.sect 	"T:wfm1" 
* Reserve space for DMA'd U in block 0
	.def	c_u0_p
	.def	c_u1_p
	.def	c_u0
	.def	c_u1
	
c_u0_p	.word 	c_u0 
c_u1_p	.word	c_u1 
c_u0	.space 	U_SIZE
c_u1	.space 	U_SIZE

	.sect 	"T:wfm0" 
* Reserve space for tas1.  
	.def	tas1
	.def	tas1_p
	.def	tpsi1_p
tpsi1_p	.word	tas1
tas1_p	.word	tas1
tas1	.space	AB_SIZE

	.sect 	"T:wfm1" 
* Reserve space for ab
	.def	ab
	.def	ab0
	.def	ab1
	.def	ab2
	.def	ab3
ab	.word	ab0
ab0	.space	1
ab1	.space	1
ab2	.space	1
ab3	.space	1

;;;	"trick" routines
	.sect 	"T:wfm1" 
* Reserve space for the temporary ab0_t, ab1_t, ab2_t, ab3_t
* used in the "trick" routines
	.def	ab0_t
	.def	ab1_t
	.def	ab2_t
	.def	ab3_t
ab0_t	.space	1
ab1_t	.space	1
ab2_t	.space	1
ab3_t	.space	1


;;; trick_kappa_spproj and trick_kappa
	.sect 	"T:wfm1" 
* Researve space for the parameters
	.def	parameters
	.def	kappa_sq
	.def	i_kappa_sq
	.def	m_kappa_sq
parameters	.word	kappa_sq
kappa_sq	.space	1
i_kappa_sq	.space	1
m_kappa_sq	.space	1
