#include<config.h>
CPS_START_NAMESPACE

//---------------------------------------------------------------------------
//
// alg_w_spect.h
//
// Header file for all alg classes relevant to Wilson-type fermion
// spectroscopy. The type of glue or fermion is given as
// an argument of type Lattice& to the constructor. If the 
// fermion type is not F_CLASS_WILSON or F_CLASS_CLOVER or 
// F_CLASS_DWF the constructors exit with a general error.
//
//
// SEE COMMENTS AT THE END ABOUT THE NAME SCHEME OF MESON DATA FILES (4/4/99)
//---------------------------------------------------------------------------



#ifndef INCLUDED_ALG_W_SPECT_H
#define INCLUDED_ALG_W_SPECT_H


CPS_END_NAMESPACE
#include <alg/alg_base.h>  
#include <alg/w_spect_arg.h>

CPS_START_NAMESPACE

//---------------------------------------------------------------------------
// Forward Declarations
//---------------------------------------------------------------------------
class Lattice;                      
class CommonArg;                    
class WspectArg;                    


//---------------------------------------------------------------------------
// AlgWspect is derived from Alg and is relevant to  
// spectroscopy with Wilson type fermions (i.e. Wilson, Clover, Dwf). 
// The type of fermion is determined by the argument to the 
// constructor. If the fermion type is not F_CLASS_WILSON or 
// F_CLASS_CLOVER or F_CLASS_DWF the constructors exit with a 
// general error.
// Xiaodong: Only normal mesons are calculated
//           To calculate normal+extended mesons, use AlgWspectExtMeson
//---------------------------------------------------------------------------

class AlgWspect : public Alg
{
 private:
  static char *d_class_name;
 protected:                 
  static int   d_counter;
  static int   d_count_step;

  WspectArg *  d_arg_p;          // the spectrum arg for each quark mass
  CgArg *      cg_arg_p;         // added by mflin to really pass CgArg to the class
  int          d_num_args;       // num of non-degenerate quark masses
  
 public:
  AlgWspect(Lattice & latt, 
	    CommonArg *c_arg, 
	    WspectArg *arg,
	    CgArg *cg,
	    int n_quark_masses = 1);   // Ping -- needed for heavy-light

  static void SetCounter(int counter, int step);
  
  static int  GetCounter()              { return d_counter; }

  virtual ~AlgWspect();

  void run(void);
};

//derived class
class AlgWspectExtMeson : public AlgWspect
{
 private:
  static char *d_class_name;
  //simply override the run() function
 public:
  AlgWspectExtMeson(Lattice & latt, 
            CommonArg *c_arg, 
            WspectArg *arg,
		    CgArg *cg,
            int n_quark_masses = 1);
  virtual ~AlgWspectExtMeson();
  
  void run(void);
};


#endif


CPS_END_NAMESPACE
