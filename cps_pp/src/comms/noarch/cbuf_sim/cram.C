#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2006-09-18 05:07:40 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/noarch/cbuf_sim/cram.C,v 1.5 2006-09-18 05:07:40 chulwoo Exp $
//  $Id: cram.C,v 1.5 2006-09-18 05:07:40 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.5 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/noarch/cbuf_sim/cram.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
CPS_END_NAMESPACE
#include<comms/nga_reg.h>
CPS_START_NAMESPACE
//Allocate space for the SCRATCH CRAM 
//Workstation version only
#ifndef _TARTAN
int CRAM_SCRATCH_INTS[CRAM_SCRATCH_SIZE] ;
unsigned long CRAM_SCRATCH_ADDR = (unsigned long)CRAM_SCRATCH_INTS ;
#endif
CPS_END_NAMESPACE
