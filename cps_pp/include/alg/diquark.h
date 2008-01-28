#ifndef INCLUDED_DIQUARK_H
#define INCLUDED_DIQUARK_H

#include <comms/glb.h>
#include "wilson_matrix.h"

CPS_START_NAMESPACE

enum ProjectType {PPAR    = 0, 
		  NPAR    = 1, 
		  PPAR_5Z = 2, 
		  NPAR_5Z = 3,
		  PPAR_5Y = 4, 
		  NPAR_5Y = 5,  
		  PPAR_5X = 6, 
		  NPAR_5X = 7,
		  PPAR_XY = 8,
		  NPAR_XY = 9,
                  PPAR_5  = 10,
                  NPAR_5  = 11
} ;

// The Diquark class.
//------------------------------------------------------------------

class Diquark {
  WilsonVector q[4][4];

 public:
  Diquark(){};

  void Project(WilsonVector& res, ProjectType P) ;
  void D_diquark(WilsonMatrix& Q1, WilsonMatrix& Q2, 
		 WilsonMatrix& Q3, WilsonMatrix& Q4, int spin, int color) ;
  void U_diquark(WilsonMatrix& Q1, WilsonMatrix& Q2, int spin, int color) ;
    
  ~Diquark(){} ;
  
} ;

CPS_END_NAMESPACE
#endif

