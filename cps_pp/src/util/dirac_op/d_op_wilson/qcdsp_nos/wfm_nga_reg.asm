**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: chulwoo $
**  $Date: 2004-01-13 20:39:45 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_nos/wfm_nga_reg.asm,v 1.2 2004-01-13 20:39:45 chulwoo Exp $
**  $Id: wfm_nga_reg.asm,v 1.2 2004-01-13 20:39:45 chulwoo Exp $
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
**  Revision 1.2  2001/06/19 18:13:03  anj
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
**  $RCSfile: wfm_nga_reg.asm,v $
**  $Revision: 1.2 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/qcdsp_nos/wfm_nga_reg.asm,v $
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
