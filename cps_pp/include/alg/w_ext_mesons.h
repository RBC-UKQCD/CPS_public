#include<config.h>
CPS_START_NAMESPACE

//
// w_ext_mesons.h:
//      
//   WspectExtendedMesons  -- implemented in w_ext_mesons.C
//

#ifndef INCLUDED_W_EXT_MESONS_H
#define INCLUDED_W_EXT_MESONS_H


CPS_END_NAMESPACE
#include <util/data_types.h>        // Float, Complex
#include <alg/w_spect_arg.h>        // SourceKind
#include <alg/w_ginfo.h>
#include <alg/w_hyper_rect.h>
#include <alg/w_gamma_mat.h>
#include <alg/w_fuzzing.h>
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


//-------------------------------
// struct WMesonOpInfo
//-------------------------------
#define WMESON_OP_MAX_NUM_TERMS 6
struct WMesonOpInfo{
  int num_terms;
  //Each term, contains: [0]-sign, [1]-gamma matrix, [2]-deriv operator
  int terms[WMESON_OP_MAX_NUM_TERMS][3];
};

//------------------------------
//struct WMesonStateInfo
//------------------------------
struct WMesonStateInfo{
  char *stateName;
  WMesonOutputName mesonId;//the output meson id
  WMesonCategory category;
  int polarization;
  int measure;
  WMesonOpKind sinkOp;
  WMesonOpKind srcOp;//int mesonOp[2];
};

//added by Thomas and Xiaodong
//---------------------------------------------------------------------------
// class WspectExtendedMesons
//---------------------------------------------------------------------------
class WspectExtendedMesons : public WspectGinfo {

  

 private:
  static struct WMesonOpInfo WMesonOpTable[NUM_WMESON_OP_KIND];
  static struct WMesonStateInfo WMesonStateTable[NUM_WMESON_STATE];
  static int tableInitialized; //avoid duplicate table for several objects
  int map[NUM_WMESON_STATE]; /*map each meson(real name) to its entry in coor_data_p
			   *which is ordered assuming prop_dir=3
			   *So if prop_dir==3, mapToCoor[i]=i
			   */

  Complex*       coor_data_p;  // storage [global_slices][EXTMESONS]
  int            d_size;    // size of coor_data_p in Floats
  Complex*       coor_output_p; //store polarization avaraged correlators
  int            d_output_size;

 public:
  // Constructors
  //-------------------------------------------------------------------------
  //allocate mem=0/1 used by subclass to avoid duplicate memory allocation
  WspectExtendedMesons(WspectArg *w_arg_p, const WspectHyperRectangle &whr, int fuzzing_index, int allocatemem=1); //whr for prop_dir info
 

  
      
  // DTOR
  //-------------------------------------------------------------------------
  ~WspectExtendedMesons();

  //collect meson correlators
  void collect(const WspectQuark &q_l, WspectQuark &q_nl, WspectFuzzing *sink_fuzz_p);
  void finish(); /*do global sum of all nodes, must be called after all collect's 
		  *and before print!*/
  // ACCESSOR write correlators to files
  //-------------------------------------------------------------------------
  void print() const;  
  //stateTable lookup, used also by AlgWspect
 Float table(int state, int sour_gamma, int sour_op, int sink_gamma, int sink_int) const;
 int isInOpGroup(int op, int groupId) const;
 int matchSUMOp(int op, int sum_op) const;

  // FOR THE PURPOSE OF DEBUGGING 
  //-------------------------------------------------------------------------
  void dumpData(char *filename) const;  

  
  // PRIVATEs and PROTECTED
  //-------------------------------------------------------------------------
 private:
  static char *  d_class_name;

 protected:
  //shared and overriden by sub classes
  WspectArg *arg_p;

  
  const int fuzzing_c_index;//fuzzing parameter index

  // The following data members are here for computation efficiency.
  int                          d_prop_dir;  
  int                          d_lclMin[LORENTZs]; 
  int                          d_lclMax[LORENTZs];
  int                          d_glb_walls; //total number of 'time' slices
  
  // constant references
  const WspectHyperRectangle & d_whr;   

  // computation buffer -- the four dirac indexes are: [D1x][D1y][D2x][D2y]
  Complex d_zero_mom_proj[DIRACs][DIRACs][DIRACs][DIRACs];

  //state table functions
  void setWMesonOpTerm(int *term_p, int weight, WGammaMatrix gammaMat, DEVOperatorKind opKind);
  void initWMesonOpTable();
  void initWMesonStateTable(WspectArg *arg);
 
  //called after applying sink operator on non-local prop(passed as qnl_p)
  //Note: qnl_p is not a full propagator, it's only propagator on a single
  //      time slice!
  void doAllAlgebra(int lclw, const Float* ql_p, const Float* qnl_p, DEVOperatorKind src_op, DEVOperatorKind);

  // Copied from normal mesons, should be modified step 1 -- color algebra
  //        -- same for all mesons
  //        return complex result[lcl_site][D1x][D1y][D2x][D2y];
  void ColorAlgebra(const Float* q1_p, const Float* q2_p,
		    int D1x, int D2x,
		    int D1y, int D2y,
		    const int a_local_site[LORENTZs],
		    Complex &result) const;      
  
  // step 2 -- zero momentum projection at sink wall(perp. to prop_dir)
  //calculate d_zero_mom_proj[D1x][D1y][D2x][D2y] 
  //        -- same for all 16 mesons
  void MomProject(const Float* q1_p, const Float* q2_p,int D1x, int D2x,
		  int D1y, int D2y,
		  int lclW);            
  
  // step 3 -- dirac 
  //        --use d_zero_mom_proj[D1x][D1y][D2x][D2y]
  //        -- different for different meson
  void DiracAlgebra(const Float* qp1, const Float* qp2, int lclW,
		    int sour_op, int sink_op);

  //do *result= sign*Trace(gam1*G1*gam2*transpose(G2)) 
  void traceDirac(Float* gam1, Float* gam2, Complex &result);
  void getBinary(const int i, int &i1, int &i2, int &i3, int &i4) const;
  void testCombination(const int sour_op, const int sink_op, Float &weight_test) const;
};
//-------------------------end of declaration of class WspectExtendedMesons


#endif // ! _INCLUDED_W_EXT_MESONS










CPS_END_NAMESPACE
