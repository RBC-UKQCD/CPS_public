#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2006-12-14 17:53:43 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/bgl/cbuf/cram.C,v 1.2 2006-12-14 17:53:43 chulwoo Exp $
//  $Id: cram.C,v 1.2 2006-12-14 17:53:43 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/bgl/cbuf/cram.C,v $
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
unsigned int CRAM_SCRATCH_ADDR = (unsigned int)CRAM_SCRATCH_INTS ;
#endif
CPS_END_NAMESPACE
