#ifndef __READLATT__
#define __READLATT__

//=======================================
// Read the format {load,unload}_lattice 
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
CPS_START_NAMESPACE
using namespace std;

//----------------------------------------
// GCFHeader class

typedef map<string,string> GCFHMapT;

class GCFheader
{
private:


  GCFHMapT headerMap;

public:

  bool    add( string key_eq_value );

  int     asInt   ( string key );
  Float   asFloat ( string key );
  string  asString( string key );

  void Show();
};


//----------------------------------------
// ReadLattice class

class ReadLattice
{
private:

  Matrix *lpoint; // pointer to the lattice
  bool allocated; // check if holding memory
  GCFheader hd ;
  Float _plaq_inheader;
  Float _linktrace_inheader;
  int recon_row_3;

public:

  DoArg   do_arg; // do_arg for this lattice

public:
    ReadLattice():
        allocated(false)
        { recon_row_3 = 0;}
  
    ReadLattice(const char* file);
    
    ~ReadLattice() { if ( allocated) delete[] lpoint; }
    
    void read( const char* file );
    unsigned calc_csum(float *fpoint);
    Float plaqInHeader() { return _plaq_inheader; }
    Float linktraceInHeader() { return _linktrace_inheader; }
    void CheckPlaqLinktrace(Lattice &lattice, Float chkprec) ;
//    void ReconRow3Single(float* dat);
//    void ReconRow3(Matrix *lpoint,float *fpoint,int size_matrices);
    void Copy(Matrix *lpoint,float *fpoint,int size_matrices, int recon_row_3);

};
CPS_END_NAMESPACE
#endif 

