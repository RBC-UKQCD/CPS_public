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

#include <util/parctrl.h>
#include <util/fpconv.h>


CPS_START_NAMESPACE
using namespace std;

//----------------------------------------
// WriteLatticeParallel class
// A modification to "WriteLattice" class to enable parallel writing

class WriteLatticeParallel
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
    WriteLatticeParallel();
    ~WriteLatticeParallel() {}
    void write(Lattice & lat, char * file, enum FP_FORMAT dataFormat = FP_AUTOMATIC, const int recon_row_3 = 0);
    inline void setConcurIONumber(int set_concur) { pc.setConcurIONumber(set_concur); }

 protected:
    char *filename; // output filename
    ParallelControl  pc;
    int data_start;
    int error;

    FPConv fpconv;

};


CPS_END_NAMESPACE
#endif
