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

#include <util/parctrl.h>
#include <util/fpconv.h>

CPS_START_NAMESPACE
using namespace std;

// GCFheaderPar class
// header parser for parallel IO
// removed "exit()"'s, others same as class GCFheader
typedef map<string,string> GCFHMapParT;

class GCFheaderPar
{
private:


  GCFHMapParT headerMap;
  bool  prevFound;

public:
  
  inline bool found() { return prevFound; }

  bool    add( string key_eq_value );

  int     asInt   ( string key );
  Float   asFloat ( string key );
  string  asString( string key );

  void Show();
};


class  ReadLatticeParallel {
 protected:
  Matrix *lpoint; // pointer to the lattice
  bool allocated; // check if holding memory
  GCFheaderPar hd ;
  Float _plaq_inheader;
  Float _linktrace_inheader;
  int recon_row_3;

  int error;
  int data_start;  // start pos of data in file
  ParallelControl  pc;
  
  FPConv  fpconv;

 public:
  
  DoArg   do_arg; // do_arg for this lattice
  inline Matrix * GaugeField()  { return lpoint; }
  inline void setConcurIONumber(int num) {  pc.setConcurIONumber(num); }

 public:
  ReadLatticeParallel():
    allocated(false), error(0)
    { 
      recon_row_3 = 0;
    }
  
  //  ReadLatticeParallel(const char* file);
  
  ~ReadLatticeParallel() { if ( allocated) delete[] lpoint; }
  
  void read( const char* file );
  unsigned calc_csum(char *fpoint);
  Float plaqInHeader() { return _plaq_inheader; }
  Float linktraceInHeader() { return _linktrace_inheader; }
  void CheckLinktrace(Matrix * mat, Float chkprec);
  void CheckPlaqLinktrace(Lattice &lattice, Float chkprec) ;
  //    void ReconRow3Single(float* dat);
  //    void ReconRow3(Matrix *lpoint,float *fpoint,int size_matrices);

};

CPS_END_NAMESPACE
#endif 


