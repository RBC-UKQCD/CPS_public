**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: zs $
**  $Date: 2004-08-18 11:57:38 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_ghb/qcdsp/assign_cram0.asm,v 1.4 2004-08-18 11:57:38 zs Exp $
**  $Id: assign_cram0.asm,v 1.4 2004-08-18 11:57:38 zs Exp $
**  $Name: not supported by cvs2svn $
**  $Locker:  $
**  $Revision: 1.4 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_ghb/qcdsp/assign_cram0.asm,v $
**  $State: Exp $
**
**--------------------------------------------------------------------
*  Publish all of the variable names as globals.
	.def	_core_iscratch

	.text

_core_iscratch		.int	buffer02
buffer02		.space	3

	.end
