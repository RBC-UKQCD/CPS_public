#ifndef __IOSTYLE__
#define __IOSTYLE__


// a slight modification on ReadLattice class (if not inheritance)
// to enable parallel reading/writing of "Gauge Connection Format" lattice data

#include <util/qioarg.h>
#include <util/fpconv.h>
#include <util/latheader.h>

CPS_START_NAMESPACE

/************************************************************************************/

class IoStyle : protected QioControl {
 protected:
  QioArg qio_arg;

  IoStyle(const QioArg & qarg) : QioControl(), qio_arg(qarg) {  }

 public:
  virtual int load(char * data, const int data_per_site, const int site_mem,
		   const LatHeaderBase & hd, const DataConversion & dconv,
		   const int dimension /* 4 or 5 */,
		   unsigned int * ptrcsum, unsigned int * ptrpdcsum = 0,
		   Float * rand_sum = 0, Float * rand_2_sum = 0) = 0;
    virtual int store(std::iostream & output,
		    char * data, const int data_per_site, const int site_mem,
		    LatHeaderBase & hd, const DataConversion & dconv,
		    const int dimension /* 4 or 5 */,
		    unsigned int * ptrcsum, unsigned int * ptrpdcsum = 0,
		    Float * rand_sum = 0, Float * rand_2_sum = 0) = 0;
};


class ParallelIO : public IoStyle {
 private:
  char * cname;
 public:
  ParallelIO(const QioArg & qarg) : IoStyle(qarg), cname("ParallelIO") { }

  virtual int load(char * data, const int data_per_site, const int site_mem,
		   const LatHeaderBase & hd, const DataConversion & dconv,
		   const int dimension /* 4 or 5 */,
		   unsigned int * ptrcsum, unsigned int * ptrpdcsum = 0,
		   Float * rand_sum = 0, Float * rand_2_sum = 0);
    virtual int store(std::iostream & output,
		    char * data, const int data_per_site, const int site_mem,
		    LatHeaderBase & hd, const DataConversion & dconv,
		    const int dimension /* 4 or 5 */,
		    unsigned int * ptrcsum, unsigned int * ptrpdcsum = 0,
		    Float * rand_sum = 0, Float * rand_2_sum = 0);

};



#if 1  // Serial IO mode only applicable to QCDOC

class SerialIO : public IoStyle {
 private:
  char * cname;
 public:
  SerialIO(const QioArg & qarg) : IoStyle(qarg), cname("SerialIO") { }

  virtual int load(char * data, const int data_per_site, const int site_mem,
		   const LatHeaderBase & hd, const DataConversion & dconv,
		   const int dimension /* 4 or 5 */,
		   unsigned int * ptrcsum, unsigned int * ptrpdcsum = 0,
		   Float * rand_sum = 0, Float * rand_2_sum = 0);
    virtual int store(std::iostream & output,
		    char * data, const int data_per_site, const int site_mem,
		    LatHeaderBase & hd, const DataConversion & dconv,
		    const int dimension /* 4 or 5 */,
		    unsigned int * ptrcsum, unsigned int * ptrpdcsum = 0,
		    Float * rand_sum = 0, Float * rand_2_sum = 0);

 private:
  inline int isNode0() const { 
    return qio_arg.Xcoor()==0 && qio_arg.Ycoor()==0 && qio_arg.Zcoor()==0 && qio_arg.Tcoor()==0
      && qio_arg.Scoor()==0; 
  }
  inline int isRow0() const { 
    return qio_arg.Ycoor()==0 && qio_arg.Zcoor()==0 && qio_arg.Tcoor()==0 && qio_arg.Scoor()==0; 
  }
  inline int isFace0() const { 
    return qio_arg.Zcoor()==0 && qio_arg.Tcoor()==0 && qio_arg.Scoor()==0; 
  }
  inline int isCube0() const { 
    return qio_arg.Tcoor()==0 && qio_arg.Scoor()==0; 
  }
  inline int isSdim0() const { 
    return qio_arg.Scoor()==0; 
  }

  // default shifting direction 3->2->1->0->(N-1)->...
  void xShiftNode(char * data, const int xblk, const int dir = -1) const;
  void yShift    (char * data, const int xblk, const int dir = -1) const;
  void zShift    (char * data, const int xblk, const int dir = -1) const;
  void tShift    (char * data, const int xblk, const int dir = -1) const;
  void sShift    (char * data, const int xblk, const int dir = -1) const;
  void sSpread   (char * data, const int xblk) const;

 public:
  // these two do tests on *Shift() functions, when you are on large rack, try these
  int backForthTest(); // move data piece to next line (next node when on x dir) then back
  int rotateTest();    // move data piece until they return to starting place

};

#endif  // TARGET == QCDOC



CPS_END_NAMESPACE
#endif 


