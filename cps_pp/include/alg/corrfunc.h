//------------------------------------------------------------------
//
// corrfunc.h
//
// Header file for the CorrFunc class 
// 
// April 2001
//
// Kostas Orginos
//
//------------------------------------------------------------------


#ifndef INCLUDED_CORRFUNC_H
#define INCLUDED_CORRFUNC_H

#include <stdlib.h> 
#include <stdio.h>
#include <util/data_types.h>
#include <comms/glb.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/error.h>
#include <util/qcdio.h>

CPS_START_NAMESPACE

class CorrFunc {
  Complex *func ;
  int Nt ;
  
 public:
  CorrFunc() ;
  // copy constructor
  CorrFunc(const CorrFunc& rhs) ;
  ~CorrFunc(){ sfree(func);}
  // "equal" operator for CorrFunc 
  CorrFunc& operator=(const CorrFunc& rhs) ;
  CorrFunc& operator*=(const CorrFunc& rhs) ;
  CorrFunc& operator+=(const CorrFunc& rhs) ;
  CorrFunc& operator-=(const CorrFunc& rhs) ;
  CorrFunc& operator*=(const Float& r) ;
  CorrFunc& operator*=(const Complex& c) ;
  
  int TimeSize() const  {return Nt;} 

  // Output functions
  void print(FILE *fp) const 
    {
      for(int t=0; t<Nt; t++)
	Fprintf(fp,"%i  %e  %e\n", t,func[t].real(),func[t].imag());
    }
  
  void print(FILE *fp,int t) const
    {
      Fprintf(fp," %e  %e ",func[t].real(),func[t].imag());
    }

  void GlobalSum(void)
    {
      slice_sum((Float*)func, 2*Nt, 99);
    }

  Complex& operator[](int i) { return func[i] ; }

  void Zero(){ for(int t = 0 ; t<Nt; t++) func[t] = 0.0 ;}

};

CPS_END_NAMESPACE

#endif //!INCLUDED_CORRFUNC_H


