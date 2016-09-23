#include<config.h>
CPS_START_NAMESPACE
//-------------------------------------------------------------------
/*!\file
  \brief  Definition of local_checksum and global_checksum outine.

*/
//--------------------------------------------------------------------
//  CVS keywords
//  $Source: /space/cvs/cps/cps++/src/comms/bgl/glb/check_sum.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//--------------------------------------------------------------
//  local_checksum
//
//  Calculating checksum on a set of data on one node
//
//  global_checksum
//
//  Calculating checksum on each node and compare against each other
//--------------------------------------------------------------


CPS_END_NAMESPACE
//#include <util/data_types.h>
#include<comms/glb.h>
CPS_START_NAMESPACE

unsigned int local_checksum(Float * float_p, int len) {
  unsigned int csum = 0;
  unsigned int * uint_p = (unsigned int *) float_p;
  for(int i=0;i < len * sizeof(Float) / sizeof(unsigned int); i++) {
//    csum += *uint_p;
    csum ^= *uint_p;
    uint_p ++;
  }
  return csum;
}

unsigned int global_checksum(Float * float_p, int len) {
 return local_checksum(float_p,len);
 // ERR.NotImplemented("","global_checksum(Float *, int)");
}

unsigned int test_checksum(Float * float_p, int len) {
 return local_checksum(float_p,len);
 //ERR.NotImplemented("","test_checksum(Float *, int)");
}

CPS_END_NAMESPACE
