#ifndef __READLATTPAR__
#define __READLATTPAR__


// a slight modification on ReadLattice class (if not inheritance)
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

CPS_START_NAMESPACE
using namespace std;


class  ReadLatticeParallel : private QioControl {
 private:
  char *cname;
  GCFheaderPar hd ;
  FPConv  fpconv;
  bool load_good;

 public:
  ReadLatticeParallel(Lattice & lat, const char * filename, const Float chkprec = 0.01)
    : QioControl(), load_good(false),cname("ReadLatticeParallel")    {        
    QioArg rd_arg(filename,chkprec);
    read(lat,rd_arg);
  }

  ReadLatticeParallel(Lattice & lat, const QioArg & rd_arg) 
    : QioControl(), load_good(false),cname("ReadLatticeParallel")   {
    read(lat,rd_arg);
  }
  
  virtual ~ReadLatticeParallel() {}

  void read( Lattice & lat, const QioArg & rd_arg);

  inline bool good() { return load_good; }

 private:
  //  bool CheckSum(char *fpoint, int size_Floats, const int Scoor);
  bool CheckPlaqLinktrace(Lattice & lat, const QioArg & rd_arg,
			  const Float plaq_inheader, const Float linktrace_inheader);
};

CPS_END_NAMESPACE
#endif 


