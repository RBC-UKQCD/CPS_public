#ifndef _WRITELATTICEPAR_H_
#define _WRITELATTICEPAR_H_
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

#include "qioarg.h"
#include "fpconv.h"


CPS_START_NAMESPACE
using namespace std;

//----------------------------------------
// WriteLatticeParallel class
// A modification to "WriteLattice" class to enable parallel writing

class WriteLatticeParallel : private QioControl
{

 public:
    // header strings
    string hd_ensemble_id ;
    string hd_ensemble_label ;
    int hd_sequence_number ;
    string hd_creator ;
    string hd_creator_hardware ;
    string hd_creation_date ;
    string hd_archive_date ;


 public:
    WriteLatticeParallel(Lattice & lat, const char * filename,
			 const FP_FORMAT dataFormat = FP_AUTOMATIC, const int recon_row_3 = 0)
      : QioControl(), unload_good(false)    {
      QioArg  wt_arg(filename, dataFormat, recon_row_3);
      write(lat, wt_arg);
    }

    WriteLatticeParallel(Lattice & lat, const QioArg & wt_arg)
      : QioControl(), unload_good(false)    {
      write(lat, wt_arg);
    }

    ~WriteLatticeParallel() {}

    void writeHeader(ostream & fout, Float link_trace, Float plaq, const QioArg & wt_arg); 

    void write(Lattice & lat, const QioArg & wt_arg);

    inline bool good() { return unload_good; }
 private:
    FPConv fpconv;
    bool unload_good;
    int data_start;
    int csum_pos;
    bool recon_row_3;
};


CPS_END_NAMESPACE
#endif
