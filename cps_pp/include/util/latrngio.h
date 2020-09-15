#ifndef __LAT_RNG_IO__
#define __LAT_RNG_IO__

#include <config.h>
#include <util/qioarg.h>
#include <util/iostyle.h>
#include <util/intconv.h>
#include <util/latheader.h>
#include <util/random.h>

CPS_START_NAMESPACE

class LatRngIO : public QioControl {
 private:
  char * cname;
 public:
  LatRngIO () : QioControl(), cname("LatRngIO") { 
    //CK: add default IO style for different architectures rather than hardcoding it in
#if TARGET != BGQ
    setParallel();
#else
    setSerial();
#endif

  }
  virtual ~LatRngIO() { }

 protected:
  IntConv intconv;

  LatRngHeader hd;

 private:
    bool UseParIO;
 public:
    inline void setParallel() { UseParIO = 1; }

    inline void setSerial() { 
#if 1
      UseParIO = 0; 
#else
      const char * fname = "setSerial()";
      VRB.Flow(cname,fname, "On non-QCDOC platform, setSerial() has no effect!\n");
#endif
    }

    inline int parIO() const { return UseParIO; }

};


class LatRngRead : public LatRngIO {
 private:
  char * cname;
 public:
  LatRngRead();
  virtual ~LatRngRead();

#ifndef USE_C11_RNG
  void read(UGrandomGenerator * ugran, UGrandomGenerator * ugran_4d,
	    const QioArg & rd_arg);
#else
  void read(RNGSTATE *mtran_dump,
	    const QioArg & rd_arg);
#endif

};

class LatRngWrite : public LatRngIO {
 private:
  char * cname;
 public:
  LatRngWrite();
  virtual ~LatRngWrite();

#ifndef USE_C11_RNG
  void write(UGrandomGenerator * ugran, UGrandomGenerator * ugran_4d,
	     const QioArg & wt_arg);
#else
  void write(RNGSTATE *mtran_dump,
	    const QioArg & rd_arg);
#endif

};



CPS_END_NAMESPACE

#endif

