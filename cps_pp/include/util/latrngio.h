#ifndef __LAT_RNG_IO__
#define __LAT_RNG_IO__

#include <config.h>
#include <util/qioarg.h>
#include <util/iostyle.h>
#include <util/intconv.h>
#include <util/latheader.h>

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

  LatRngHeader hd;

 private:
    bool UseParIO;
 public:
    inline void setParallel() { UseParIO = 1; }

    inline void setSerial() { 
#if TARGET == QCDOC
      UseParIO = 0; 
#else
      cout << "On non-QCDOC platform, setSerial() has no effect!" << endl;
#endif
    }

    inline int parIO() const { return UseParIO; }

};


class LatRngRead : public LatRngIO {
 public:
  LatRngRead();
  virtual ~LatRngRead();

  void read(UGrandomGenerator * ugran, UGrandomGenerator * ugran_4d,
	    const QioArg & rd_arg);

};

class LatRngReadSerial : public LatRngRead {
 public:
  LatRngReadSerial() : LatRngRead() {
    setSerial();
  }
};

class LatRngWrite : public LatRngIO {
 public:
  LatRngWrite();
  virtual ~LatRngWrite();

  void write(UGrandomGenerator * ugran, UGrandomGenerator * ugran_4d,
	     const QioArg & wt_arg);

};

class LatRngWriteSerial : public LatRngWrite {
 public:
  LatRngWriteSerial() : LatRngWrite() {
    setSerial();
  }
};


CPS_END_NAMESPACE

#endif

