#include <iostream>
#include <unistd.h>
#include <util/gjp.h>
#include <sys/time.h>
#include <iomanip>
#include <cstring>
#include <comms/glb.h>
using namespace std;

#include <util/qioarg.h>

CPS_START_NAMESPACE


// The following are NEW functions added to QioControl class 
// to enable message passing between parallel processors, based on QMP calls
// (some are pretty useful...)

int QioControl::synchronize(const int errorStatus)  const {
  const char * fname = "synchronize()";
  int error = errorStatus;
 
//  QMP_barrier();
//  printf("Node %d:synchronize(%d)\n",UniqueID(),errorStatus);
  if(NumNodes()>1) {
    error = globalSumInt(error);
    if(error > 0) {
      VRB.Flow(cname,fname,"Totally %d nodes reported error!\n",error);
    }
  }
//  QMP_barrier();
//  printf("Node %d:synchronize(%d)\n",UniqueID(),errorStatus);
  return error;
}

void QioControl::broadcastInt(int * data, int fromID)  const {
  if(NumNodes() > 1) {
    if(unique_id != fromID) {
      *data = 0;
    }

    *data = globalSumInt(*data);
  }
}

void QioControl::broadcastFloat(Float * data, int fromID) const {
  if(NumNodes() > 1) {
    if(unique_id != fromID) {
      * data = 0;
    }
    *data = globalSumFloat(*data);
  }
}

int QioControl::round(const Float fdata) const{
  int ndata = (int)fdata;
  if(fdata - ndata >= 0.5) ndata++;
  if(fdata - ndata < -0.5) ndata--;
  return ndata;
}

int QioControl::globalSumInt(const int data) const{
  int temp=data;
#ifndef UNIFORM_SEED_TESTING
//  QMP_barrier();
//  printf("Node %d: before QMP_sum_int(%d)\n",UniqueID(),data);
  QMP_sum_int(&temp);
//  QMP_barrier();
//  printf("Node %d: after QMP_sum_int(%d)\n",UniqueID(),temp);
#endif
  return temp;
}

void sum_uint(void *inout, void *in){
  unsigned int *u_inout = (unsigned int *)inout;
  unsigned int *u_in= (unsigned int *)in;
  *u_inout += *u_in;
}

unsigned int QioControl::globalSumUint(const unsigned int data) const{
#if 1
  //  Gsum64Ext  gsum;
  //  return gsum.Sum(data);
  int hfbits = sizeof(unsigned int) * 8 / 2;
  unsigned int mask = (1 << hfbits) - 1;

  unsigned int sumd = data;
  unsigned int hi = sumd >> hfbits;
  unsigned int lo = sumd & mask;
  hi = round(globalSumFloat(hi));
  lo = round(globalSumFloat(lo));
  sumd = (hi<<hfbits)+lo;
  return sumd;
#else
  unsigned int temp=data;
//  QMP_sum_int(&temp);
  QMP_binary_reduction(&temp,sizeof(unsigned int),sum_uint);
  return temp;
#endif
}

Float QioControl::globalSumFloat(const Float data) const {
  double temp=data;
#ifndef UNIFORM_SEED_TESTING
  QMP_sum_double(&temp);
#endif
  return (Float)temp;
}

int QioControl::globalMinInt(const int data) const{
  double temp=data;
#ifndef UNIFORM_SEED_TESTING
  QMP_min_double(&temp);
#endif
  return (int)temp;
}

CPS_END_NAMESPACE
