#ifndef __READU1LATTPAR__
#define __READU1LATTPAR__


// a slight modification on ReadU1Lattice class (if not inheritance)
// to enable parallel reading/writing of "Gauge Connection Format" lattice data

#include <stdlib.h>	// exit()
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <map>
#include <string>

#include <util/gjp.h>
#include <util/vector.h>
#include <util/lattice.h>
#include <util/data_types.h>
#include <util/verbose.h>
#include <util/error.h>
#include <alg/do_arg.h>
#include <util/qioarg.h>
#include <util/fpconv.h>
#include <util/iostyle.h>
#include <util/latheader.h>

CPS_START_NAMESPACE
using namespace std;


class  ReadU1LatticeParallel : public QioControl{
  // which determines parallel reading or serial reading

 private:
  char *cname;
  LatticeHeader hd;

public:
  // ctor for 2-step loading
  ReadU1LatticeParallel()
    : QioControl(), cname("ReadU1LatticeParallel"), UseParIO(1)
    {  }

  // ctor invoking loading behavior
  ReadU1LatticeParallel(Lattice & lat, const char * filename, const Float chkprec = 0.01)
    : 
    QioControl(),
    cname("ReadU1LatticeParallel") , 
    UseParIO(1)
    {        
    QioArg rd_arg(filename,chkprec);
    read(lat,rd_arg);
  }

  // ctor invoking loading behavior
  ReadU1LatticeParallel(Lattice & lat, const QioArg & rd_arg) 
    : QioControl(), cname("ReadU1LatticeParallel"), UseParIO(1)
  {
    read(lat,rd_arg);
  }
  
  virtual ~ReadU1LatticeParallel() {}

  void read(Lattice & lat, const char * filename, const Float chkprec = 0.01){
    QioArg rd_arg(filename,chkprec);
    read(lat,rd_arg);
  }

  void read(Lattice & lat, const QioArg & rd_arg);

 private:
  bool CheckPlaqU1Linktrace(Lattice & lat, const QioArg & rd_arg,
			    const Float plaq_inheader, const Float linktrace_inheader);

 private:
    bool UseParIO;
 public:
    inline void setParallel() { UseParIO = 1; }
    inline void setSerial() { 
#if 1
      UseParIO = 0; 
#else
      const char * fname = "setSerial()";
      VRB.Result(cname,fname,"On non-QCDOC platform, setSerial() has no effect!\n");
      exit(-42);
#endif
    }
    inline int parIO() const { return UseParIO; }
};


class ReadU1LatticeSerial : public ReadU1LatticeParallel {
 private:
  char * cname;

 public:
  // ctor for 2-step loading
  ReadU1LatticeSerial() : ReadU1LatticeParallel(), cname("ReadU1LatticeSerial") {
    setSerial();
  }

  // ctor invoking loading behavior
  ReadU1LatticeSerial(Lattice & lat, const char * filename, const Float chkprec = 0.01) 
    : ReadU1LatticeParallel(), cname("ReadU1LatticeSerial") {
    setSerial();
    read(lat, filename, chkprec);
  }

  // ctor invoking loading behavior
  ReadU1LatticeSerial(Lattice & lat, const QioArg & rd_arg) 
    : ReadU1LatticeParallel(), cname("ReadU1LatticeSerial") {
    setSerial();
    read(lat, rd_arg);
  }
  
  virtual ~ReadU1LatticeSerial() {}
};


CPS_END_NAMESPACE
#endif 


