#include <iostream>
#include <unistd.h>
#include <util/gjp.h>
#include <sys/time.h>
#include <iomanip>
#include <cstring>
#include <comms/glb.h>
using namespace std;

#include <util/qioarg.h>

#if TARGET == QCDOC
#include <util/gsum64ext.h>
#endif

CPS_START_NAMESPACE


// The following are NEW functions added to QioControl class 
// to enable message passing between parallel processors, based on QMP calls
// (some are pretty useful...)

int QioControl::synchronize(const int errorStatus)  const {
  const char * fname = "synchronize()";
  int error = errorStatus;
 
  if(NumNodes()>1) {
    error = globalSumInt(error);
    if(error > 0) {
      VRB.Flow(cname,fname,"Totally %d nodes reported error!\n",error);
    }
  }
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
#if TARGET == QCDOC
  //  Gsum64Ext  gsum;
  //  return gsum.Sum(data);
  int hfbits = sizeof(unsigned int) * 8 / 2;
  unsigned int mask = (1 << hfbits) - 1;

  int sumd = data;
  int hi = sumd >> hfbits;
  int lo = sumd & mask;
  hi = round(globalSumFloat(hi));
  lo = round(globalSumFloat(lo));
  sumd = (hi<<hfbits)+lo;
  return sumd;
#else
  return data;
#endif
}

unsigned int QioControl::globalSumUint(const unsigned int data) const{
#if TARGET == QCDOC
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
  return data;
#endif
}

Float QioControl::globalSumFloat(const Float data) const {
#if TARGET == QCDOC
  //  Gsum64Ext  gsum;
  //  return gsum.Sum(data);
  Float sumdata = data;
  glb_sum_five(&sumdata);
  return sumdata;
#else
  return data;
#endif
}

int QioControl::globalMinInt(const int data) const{
#if TARGET == QCDOC
  Float fdata = data;
  glb_min(&fdata);
  int res = round(fdata);
  return res;
#else
  return data;
#endif
}

CPS_END_NAMESPACE
