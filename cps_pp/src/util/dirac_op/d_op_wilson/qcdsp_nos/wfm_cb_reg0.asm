**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: chulwoo $
**  $Date: 2004-06-04 21:14:09 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_nos/wfm_cb_reg0.asm,v 1.3 2004-06-04 21:14:09 chulwoo Exp $
**  $Id: wfm_cb_reg0.asm,v 1.3 2004-06-04 21:14:09 chulwoo Exp $
**  $Name: not supported by cvs2svn $
**  $Locker:  $
**  $Revision: 1.3 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_nos/wfm_cb_reg0.asm,v $
**  $State: Exp $
**
**--------------------------------------------------------------------

****************************************************************************************
* wfm_cb_reg0:
*
* wfm_s_cb0 saves the contents of the circular buffer control 
* register 0 in memory (at address save_cb0_p).
*
* wfm_r_cb0 restores the contents of the circular buffer control 
* register 0 from memory address save_cb0_p.
*
* These routines do not respect c calling conventions. Also they
* change the DP. 
*
****************************************************************************************
	.version	30

	.def	wfm_s_cb0
	.def	wfm_r_cb0

	.ref	cb_cntrl_b



	.text
save_cb0_p	.word	save_cb0
save_cb0	.word	0


****************************************************************************************
* wfm_s_cb0
****************************************************************************************
wfm_s_cb0:
	LDP	@cb_cntrl_b
	LDI	@cb_cntrl_b,	AR0
	LDP	@save_cb0_p
	LDI	@save_cb0_p,	AR1
	LDI	*AR0,		R0
	STI	R0,		*AR1
	RETS


****************************************************************************************
* wfm_r_cb0
****************************************************************************************
wfm_r_cb0:
	LDP	@cb_cntrl_b
	LDI	@cb_cntrl_b,	AR0
	LDP	@save_cb0_p
	LDI	@save_cb0_p,	AR1
	LDI	*AR1,		R0
	STI	R0,		*AR0
	RETS




