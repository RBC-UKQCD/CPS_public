#include<config.h>
CPS_START_NAMESPACE
#ifndef INCLUDED_W_MESONS
#define INCLUDED_W_MESONS


CPS_END_NAMESPACE
#include <util/data_types.h>        // Float, Complex
#include <alg/w_spect_arg.h>        // SourceKind
#include <alg/w_ginfo.h>
#include <alg/w_momenta.h>
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
// class WspectMesons
// 
//---------------------------------------------------------------------------
class WspectMesons : public WspectGinfo {

public:
  // CTORs
  //-------------------------------------------------------------------------
  WspectMesons(const WspectQuark &quark1, 
	       const WspectQuark &quark2,
	       const WspectHyperRectangle & whr, 
	       const WspectMomenta& mom);
  WspectMesons(const IFloat *quark1, 
	       const IFloat *quark2, 
	       const WspectHyperRectangle & whr,
	       const WspectMomenta& mom);
      
  // DTOR
  //-------------------------------------------------------------------------
  ~WspectMesons();

  // ACCESSOR
  //-------------------------------------------------------------------------
  void print(WspectOutput *w_spect_output) const;  
  void print_mp(char *filename) const;


  // FOR THE PURPOSE OF DEBUGGING 
  //-------------------------------------------------------------------------
  void dumpData(char *filename) const;  


  // PRIVATEs
  //-------------------------------------------------------------------------
 private:
  static char *  d_class_name;
  Complex*       d_data_p;  // Complex[d_num_mom][global_slices][MESONs]
  int            d_size;    // size of d_data_p in IFloats

  //PROTECTED
  //-------------------------------------------------------------------------
 protected:
  int            d_num_mom; // zero-momentum included

  // constant references
  const IFloat *                d_quark1_p;
  const IFloat *                d_quark2_p;
  const WspectHyperRectangle & d_whr;  
  const WspectMomenta&         d_mom; 

  // The following data members are here for computation efficiency.
  int                          d_prop_dir;  
  int                          d_lclMin[LORENTZs]; 
  int                          d_lclMax[LORENTZs];
  int                          d_glb_walls;

  // computation buffer -- the four dirac indexes are: [D1x][D1y][D2x][D2y]
  Complex d_mom_proj[WspectMomenta::MAX_NUM+1][DIRACs][DIRACs][DIRACs][DIRACs];


  // step 1 -- color algebra
  //        -- same for all 16 mesons
  //        -- result += ...
  void ColorAlgebra(int D1x, int D2x,
		    int D1y, int D2y,
		    const int a_local_site[LORENTZs],
		    Complex &result) const;      

  // step 2 -- calculate the whole array d_mom_proj[mom][D1x][D1y][D2x][D2y] 
  //        -- same for all 16 mesons
  void MomProject(int D1x, int D2x,
		  int D1y, int D2y,
		  int lclW);            

  // step 3 -- dirac algebra for all mesons 
  //        -- different for different meson
  void DiracAlgebra(int lclW); 

  //do *result= sign*Trace(gam1*G1*gam2*transpose(G2)) 
  void traceDirac(int sign, IFloat* gam1, IFloat* gam2,int mom, Complex *result_p);
  // finished
  void Everything();  

};

#endif // ! _INCLUDED_W_MESONS




CPS_END_NAMESPACE
