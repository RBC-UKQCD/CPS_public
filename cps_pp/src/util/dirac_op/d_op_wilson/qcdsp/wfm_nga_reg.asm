**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: zs $
**  $Date: 2004-08-18 11:57:56 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp/wfm_nga_reg.asm,v 1.4 2004-08-18 11:57:56 zs Exp $
**  $Id: wfm_nga_reg.asm,v 1.4 2004-08-18 11:57:56 zs Exp $
**  $Name: not supported by cvs2svn $
**  $Locker:  $
**  $Revision: 1.4 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp/wfm_nga_reg.asm,v $
**  $State: Exp $
**
**--------------------------------------------------------------------
*
* Some NGA register definitions
*
        .sect 	"T:wfm0" 
	
* memory offsets for circular buffer operation
	.def	direct
	.def	bank0
	.def	bank1a
	.def	bank1b
	.def	bank2a
	.def	bank2b
	.def	bank3a
	.def	bank3b
	.def	bank4a
	.def	bank4b
direct	.word	0880000h
bank0	.word	0000000h
bank1a	.word	0900000h
bank1b	.word	0980000h
bank2a	.word	0A00000h
bank2b	.word	0A80000h
bank3a	.word	0B00000h
bank3b	.word	0B80000h
bank4a	.word	0C00000h
bank4b	.word	0C80000h

* The base address of the SCU control registers in the NGA.
	.def	scu_b	
scu_b .word 813000h

* The base address of the CIRBUF control registers in the NGA.
	.def	cb_cntrl_b
cb_cntrl_b .word 815800h
