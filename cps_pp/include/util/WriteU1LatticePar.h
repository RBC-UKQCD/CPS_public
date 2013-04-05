
#ifndef _WRITEU1LATTICEPAR_H_
#define _WRITEU1LATTICEPAR_H_
// Write the format {load,unload}_lattice 

#include <stdlib.h>	// exit()
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

#include <util/gjp.h>
#include <util/vector.h>
#include <util/lattice.h>
#include <util/verbose.h>
#include <util/error.h>

#include <util/qioarg.h>
#include <util/fpconv.h>
#include <util/iostyle.h>
#include <util/latheader.h>


CPS_START_NAMESPACE
using namespace std;

//----------------------------------------
// WriteU1LatticeParallel class
// A modification to "WriteU1Lattice" class to enable parallel writing

class WriteU1LatticeParallel : public QioControl
{
  // IoStyle provides a function  IoStyle::store() 
  // which determines Parallel or Serial storing

 private:
  //    FPConv fpconv;
    int csum_pos;
    bool recon_row_3;
    char *cname;

 public:
    LatticeHeader hd;

 public:
    // ctor for 2-step unloading
    WriteU1LatticeParallel()  
      : QioControl(), cname("WriteU1LatticeParallel"), UseParIO(1) {
    }

    // ctor containing unloading behavior
    WriteU1LatticeParallel(Lattice & lat, const char * filename,
			 const FP_FORMAT dataFormat = FP_AUTOMATIC, const int recon_row_3 = 0)
      : QioControl(), cname("WriteU1LatticeParallel"), UseParIO(1)  {
      QioArg  wt_arg(filename, dataFormat, recon_row_3);
      write(lat, wt_arg);
    }

    // ctor containing unloading behavior
    WriteU1LatticeParallel(Lattice & lat, const QioArg & wt_arg)
      : QioControl(), cname("WriteU1LatticeParallel"), UseParIO(1)   {
      write(lat, wt_arg);
    }

    ~WriteU1LatticeParallel() {}

    void setHeader(const char * EnsembleId, const char * EnsembleLabel,
		   const int SequenceNumber) {
      hd.setHeader(EnsembleId, EnsembleLabel, SequenceNumber);
    }

    void write(Lattice & lat, const char * filename,
	       const FP_FORMAT dataFormat = FP_AUTOMATIC, const int recon_row_3 = 0) {
      QioArg  wt_arg(filename, dataFormat, recon_row_3);
      write(lat, wt_arg);
    }
    void write(Lattice & lat, const QioArg & wt_arg);

 private:
    bool UseParIO;
 public:
    inline void setParallel() { UseParIO = 1; }
    inline void setSerial() { 
#if 1
      UseParIO = 0; 
#else
      const char * fname = "setSerial()";
      VRB.Flow(cname,fname,"On non-QCDOC platform, setSerial() has no effect!\n");
#endif
    }
    inline int parIO() const { return UseParIO; }

};



class WriteU1LatticeSerial : public WriteU1LatticeParallel {
 private:
  char * cname;

 public:
    // ctor for 2-step unloading
    WriteU1LatticeSerial()
      : WriteU1LatticeParallel(), cname("WriteU1LatticeParallel") {
      setSerial();
    }

    // ctor containing unloading behavior
    WriteU1LatticeSerial(Lattice & lat, const char * filename,
		       const FP_FORMAT dataFormat = FP_AUTOMATIC, const int recon_row_3 = 0)
      : WriteU1LatticeParallel(), cname("WriteU1LatticeParallel") {
      setSerial();
      write(lat, filename, dataFormat, recon_row_3);
    }

    // ctor containing unloading behavior
    WriteU1LatticeSerial(Lattice & lat, const QioArg & wt_arg)
      : WriteU1LatticeParallel(), cname("WriteU1LatticeParallel"){
      setSerial();
      write(lat, wt_arg);
    }

    ~WriteU1LatticeSerial() {}
};


CPS_END_NAMESPACE
#endif
