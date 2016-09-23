#include<config.h>
#if 0
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//
//
//--------------------------------------------------------------------
CPS_END_NAMESPACE
//#include<comms/nga_reg.h>
CPS_START_NAMESPACE
//Allocate space for the SCRATCH CRAM 
//Workstation version only
#ifndef _TARTAN
int CRAM_SCRATCH_INTS[CRAM_SCRATCH_SIZE] ;
unsigned long CRAM_SCRATCH_ADDR = (unsigned long)CRAM_SCRATCH_INTS ;
#endif
CPS_END_NAMESPACE
#endif
