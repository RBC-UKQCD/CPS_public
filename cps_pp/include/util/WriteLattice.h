#ifndef _WRITELATTICE_H_
#define _WRITELATTICE_H_
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
CPS_START_NAMESPACE


using namespace std;

//----------------------------------------
// WriteLattice class

class WriteLattice
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
    WriteLattice(char* file);
    ~WriteLattice() {;}
    void write();

private:
    char *filename; // output filename
    
};

CPS_END_NAMESPACE
#endif
