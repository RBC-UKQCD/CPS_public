**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: zs $
**  $Date: 2004-08-18 11:57:38 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_ghb/qcdsp/assign_cram1.asm,v 1.4 2004-08-18 11:57:38 zs Exp $
**  $Id: assign_cram1.asm,v 1.4 2004-08-18 11:57:38 zs Exp $
**  $Name: not supported by cvs2svn $
**  $Locker:  $
**  $Revision: 1.4 $
**  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_ghb/qcdsp/assign_cram1.asm,v $
**  $State: Exp $
**
**--------------------------------------------------------------------
*  Publish all of the variable names as globals.
	.def	_el_seed_p
	.def	_fast_sigma
	.def	_fast_link
	.def	_core_fscratch

	.text
_fast_sigma:			.int	buffer09
buffer09:		.space	18		;	1*18
_fast_link:			.int	buffer0a
buffer0a:		.space	18		;	1*18
_core_fscratch:		.int	buffer0d
buffer0d:		.space	59
_el_seed_p:		.int	_el_seed
_el_seed:		.int	1F123BB5h
_el_seed_:		.int	159A55E5h
_el_seed__:		.int	00F6A3D9h
_el_seed___:		.int	436CBAE9h
carry:			.int	1			; carry
lconst:			.int	1013904243		; odd const for LCG
	.end

