#include<config.h>
CPS_START_NAMESPACE
#ifndef INCLUDED_W_GINFO
#define INCLUDED_W_GINFO


CPS_END_NAMESPACE
#include <util/data_types.h>        // Float, Complex
#include <alg/w_spect_arg.h>        // SourceKind
CPS_START_NAMESPACE

//---------------------------------------------------------------------------
// Forward Declarations  -- classes defined in other translation units
//---------------------------------------------------------------------------
class Lattice;                      // defined in util/include/lattice.h
class CgArg;                        // defined in  alg/include/cg_arg.h
class CommonArg;                    // defined in  alg/include/common_arg.h
class AlgWspect;                    // defined in  alg/include/alg_w_spect.h
class WspectGinfo;                      
class WspectHyperRectangle;                     
class WspectQuark;                     
class WspectMomenta;                     

//---------------------------------------------------------------------------
// class WspectGinfo -- a support class with ONLY static data members
//---------------------------------------------------------------------------
class WspectGinfo
{
  friend class AlgWspect;
  
public: 
  // constants
  //-------------------------------------------------------------------------
  enum {LORENTZs = 4,                          // range of Lorentz index
	DIRACs   = 4,                          // range of dirac index
	COLORs   = 3,                          // range of color index
	COMPLEXs = 2,                          // size  of a complex
	SPINORs  = DIRACs * COLORs * COMPLEXs, // size  of a Wilson spinor
	MESONs   = DIRACs * DIRACs,            // num   of mesons
	DELTAs   = LORENTZs,                    // num   of deltas
	SINK_GAMMAS = DIRACs*DIRACs,            // number of possible Gamma-matrices at sink
	SOUR_GAMMAS = DIRACs*DIRACs          // number of possible Gamma-matrices at source
	
  };
  

  // Why "protected" instead of "public":
  // so that its static data member or member functions can be accesses
  // only via an instance of this class, which means the class CTOR 
  // [where the initialization is done] has to called at least once.
  //-------------------------------------------------------------------------
protected:    

  // static data -- to reduce string redundance
  //-------------------------------------------------------------------------
  static char *ctor_str;
  static char *dtor_str;
  static char *empty_str;
  static char *out_range_str;
  static char *wrong_type_str;  
  static char *inconsistent_str;
  

  // static data -- to reduce code redundance
  //-------------------------------------------------------------------------
  static int initialized;  
  static int glb_sites[LORENTZs];     // num of lat sites on the whole machine
  static int lcl_sites[LORENTZs];     // num of lat sites on the local node
  static int lcl_node [LORENTZs];     // local node coord among the machine
  static int lcl2glb_offset[LORENTZs];// lcl2glb_offset = glb_x - lcl_x
  static int bnd_cnd[LORENTZs];       // (-)1 for (anti)periodic. 
                    
  
  // CTOR & DTOR
  //-------------------------------------------------------------------------
  WspectGinfo();
  ~WspectGinfo() {}


  // MEMBER FUNCTIONS  (non-static, so CTOR has to be called at least once)
  //-------------------------------------------------------------------------
  int isOutOfRange(const int site[], const int frame[]) const;
  // return 1 if the site[] is outside the hypercube frame[]
  //        0 otherwise

  int siteOffset(const int lcl_site[], int excluded_dir = -1) const; 
  // return the offset associated with lcl_site[].
  // excluded_dir:
  //      The offset is then evaluated on the hyperplane instead of
  //      the hypercube.
  //      no effect if it is not a valid Lorentz index.

  int glb2lcl(int lcl_output[], const int glb_input[]) const;  
  // translate a global site into a local site, 
  // return 0/1 if off/on-node.
  // two arguments could be the same.
  void lcl2glb(const int lcl_input[], int glb_output[]) const;

  // FOR THE PURPOSE OF DEBUGGING ONLY
  void printSite(FILE *fp, const int site[]) const;  
  void printSpinor(FILE *fp, const IFloat *) const;
};



#endif // ! _INCLUDED_W_GINFO

CPS_END_NAMESPACE
