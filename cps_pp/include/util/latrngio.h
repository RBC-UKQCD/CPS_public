#ifndef __LAT_RNG_IO__
#define __LAT_RNG_IO__

#include <config.h>
#include <util/qioarg.h>
#include <util/intconv.h>

CPS_START_NAMESPACE
using namespace std;


class LatRngIO : public QioControl {
 public:
  LatRngIO () : QioControl(), io_good(false) { }
  virtual ~LatRngIO() { }

  bool good() { return io_good; }

 protected:
  bool io_good;
  IntConv intconv;
  int local_n[5];
  int global_n[5];
  int coor[5];
  void calcDim(const QioArg & io_arg);

  // used for position-dependent checksum to verify layout of data
  int uniqueSiteID5d(int local_id_5d);
  int uniqueSiteID4d(int local_id_4d); 
};


class LatRngRead : public LatRngIO {
 public:
  LatRngRead();
  virtual ~LatRngRead();

  void read(UGrandomGenerator * ugran, UGrandomGenerator * ugran_4d,
	    const QioArg & rd_arg);

 private:
  GCFheaderPar hd;
};


class LatRngWrite : public LatRngIO {
 public:
  LatRngWrite();
  virtual ~LatRngWrite();

  void write(UGrandomGenerator * ugran, UGrandomGenerator * ugran_4d,
	     const QioArg & wt_arg);

};


CPS_END_NAMESPACE

#endif

