#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2013-04-05 17:51:13 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/mpi/cbuf/cram.C,v 1.4 2013-04-05 17:51:13 chulwoo Exp $
//  $Id: cram.C,v 1.4 2013-04-05 17:51:13 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: cram.C,v $
//  $Revision: 1.4 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/mpi/cbuf/cram.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
CPS_END_NAMESPACE
//#include<comms/nga_reg.h>
CPS_START_NAMESPACE
//Allocate space for the SCRATCH CRAM 
//Workstation version only
#ifndef _TARTAN
int CRAM_SCRATCH_INTS[CRAM_SCRATCH_SIZE] ;
unsigned int CRAM_SCRATCH_ADDR = (unsigned int)CRAM_SCRATCH_INTS ;
#endif
CPS_END_NAMESPACE
