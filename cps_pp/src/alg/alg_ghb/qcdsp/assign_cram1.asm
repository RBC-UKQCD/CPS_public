**--------------------------------------------------------------------
**  CVS keywords
**
**  $Author: mcneile $
**  $Date: 2003-06-22 13:34:45 $
**  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_ghb/qcdsp/assign_cram1.asm,v 1.1.1.1 2003-06-22 13:34:45 mcneile Exp $
**  $Id: assign_cram1.asm,v 1.1.1.1 2003-06-22 13:34:45 mcneile Exp $
**  $Name: not supported by cvs2svn $
**  $Locker:  $
**  $Log: not supported by cvs2svn $
**  Revision 1.2  2001/06/19 18:11:23  anj
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
**  Revision 1.2  2001/05/25 06:15:59  cvs
**  Added CVS keywords to phys_v4_0_0_preCVS
**
**  $RCSfile: assign_cram1.asm,v $
**  $Revision: 1.1.1.1 $
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

