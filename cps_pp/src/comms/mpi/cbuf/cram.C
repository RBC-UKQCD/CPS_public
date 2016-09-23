#include<config.h>
//#include<comms/nga_reg.h>
CPS_START_NAMESPACE
//Allocate space for the SCRATCH CRAM 
//Workstation version only
#ifndef _TARTAN
int CRAM_SCRATCH_INTS[CRAM_SCRATCH_SIZE] ;
unsigned int CRAM_SCRATCH_ADDR = (unsigned int)CRAM_SCRATCH_INTS ;
#endif
CPS_END_NAMESPACE
