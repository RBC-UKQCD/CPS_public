#include<config.h>
CPS_START_NAMESPACE
#ifndef INCLUDED_W_EXT_MESONS_BE
#define INCLUDED_W_EXT_MESONS_BE

/*class WspectExtendedMesonsBE 
 *
 * derived from WspectExtendedMesons
 *
 *Use local magnetic and electric fields to construct Hybrid/Exotic meson
 *operators
 *
 */


CPS_END_NAMESPACE
#include <util/data_types.h>        // Float, Complex
#include <alg/w_spect_arg.h>        // SourceKind
#include <alg/w_ginfo.h>
#include <alg/w_momenta.h>
#include <alg/w_ext_mesons.h>
CPS_START_NAMESPACE

//------------------------------
//Enumeration Types
//------------------------------

//-------------------------------
// Structures
// struct WMesonOpInfo
// struct WMesonBEStateInfo
//-------------------------------
#define WMESONBE_OP_MAX_NUM_TERMS 6
struct WMesonBEOpInfo{
  int num_terms;
  //Each term, contains: [0]-sign, [1]-gamma matrix, [2]-FieldTensorId
  int terms[WMESONBE_OP_MAX_NUM_TERMS][3];
};


struct WMesonBEStateInfo{
  char *stateName; //used in filename(statename+filetail(from w_spect_arg))
  WExtMesonBEOutputName mesonId;//the output meson id
  WExtMesonBECategory category;
  int polarization;
  int measure;  
  WExtMesonBEOp srcOp;
  WExtMesonBEOp sinkOp;
};


/*Calculation procedure
 *1. Constructor (arg: qprop, fuzzedlink_ptr, w_spect_arg, fuzzing_id)
 *   Allocate Memory, initialize operator and state tables,
 *   (tables are static to eliminate redundancy when multiple
 *     objects are created) 
 *   
 *   Calculate magnetic and electric fields from fuzzed links(if fuzz_ptr!=0)
 *   (use point-like construction: P_uv is plaq 
 *     FieldTensor G_uv(x)=1/8*[(P_{u,v} -P_{v,u})+(P_{-u,-v}-P_{-v,-u})+
 *                              (P_{v,-u} -P_{-u,v})+(P_{-v,u}-P_{u,-v})]
 *     B_i=espi_ijk*G_jk
 *     E_i=G_0i
 *
 *step 2,3 are in run() function
 *2. Do Algebras
 *   1. Color Algebra
 *     Sum_C1xC2xC1yC2y(q1F1q2F2)
 *     output: colorSum(x)
 *   2. Mom Projection
 *      loop over D1xD1yD2xD2y,F1,F2(field tensors)
 *      colorSum(x)=ColorAlgebra(q1,q2,x,D1x,D1y,D2x,D2y,F1,F2)
 *      Sum_x(colorSum)
 *     output:mom[F1][F2][D1x][D1y][D2x][D2y]
 *
 *   3. Dirac Algebra
 *      
 *   Loop over each state
 *   Do Color Algebra, Momentum Projection and Dirac Trace.
 *
 *3. Do global sum
 *   
 *4  print out states
 *   print()
 *
 *
 */


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
class WspectFuzzing;
class WspectExtendedMesons;
class WspectField;

//---------------------------------------------------------------------------
// class WspectExtMesonsBE
//---------------------------------------------------------------------------
class WspectExtendedMesonsBE : public WspectExtendedMesons {

  //+++++++++++++++++++ public ++++++++++++++++++
 public:
  // CTORs
  //---------------------------------------------
  WspectExtendedMesonsBE(WspectArg *w_arg_p,
			 const WspectHyperRectangle & whr,
			 const int fuzzing_id,
			 const WspectField *field_p);
  
  // DTOR
  //--------------------------------------------
  ~WspectExtendedMesonsBE();
  
  //collect meson correlators
  void collect(const WspectQuark &q_l, WspectQuark &q_nl);
  void finish(); /*do global sum of all nodes, must be called after all collect's 
		  *and before print!*/

  // ACCESSOR
  //---------------------------------------------
  void print() const;  


   //stateTable lookup, used also by AlgWspect
  Float table(int state, int sour_gamma, int sour_fld, int sink_gamma, int sink_fld) const;
  int isInOpGroup(int op, int groupId) const;
  int matchSUMOp(int op, int sum_op) const;
  
  // FOR THE PURPOSE OF DEBUGGING 
  // not implemented yet
  //--------------------------------------------
  void dumpData(char *filename) const;  


  // PRIVATEs
  //++++++++++++++++++++++++++++++++++++++++++++++++
 private:
  //private members
  //--------------------------------
  //static tables 
  static struct WMesonBEOpInfo WMesonOpTable[NUM_WEXTMESON_BE_OPS];
  static struct WMesonBEStateInfo WMesonStateTable[NUM_WEXTMESON_BE_STATES];
  static int tableInitialized;
  int map[NUM_WEXTMESON_BE_STATES];

  Complex*       coor_data_p;  // storage [global_slices][EXTMESONS]
  int            d_size;    // size of coor_data_p in Floats
  Complex*       coor_output_p; //store polarization avaraged correlators
  int            d_output_size;
  
  static char *  d_class_name;

  
  // computation buffer -- the four dirac indexes are: [D1x][D1y][D2x][D2y]
  //Complex d_zero_mom_proj[SinkFieldTensor][DIRACs][DIRACs][DIRACs][DIRACs];
  //NUM_FLDS+1 include UNIT field
  
  Complex d_zero_mom_proj[NUM_FLDS+1][DIRACs][DIRACs][DIRACs][DIRACs];

  const WspectField *field_p;
  
  //private functions
  //---------------------------------
  //state table functions
  void setWMesonOpTerm(int *term_p, int weight, WGammaMatrix gammaMat, FieldTensorId fieldId);
  void initWMesonOpTable();
  void initWMesonStateTable(WspectArg *arg);
  
  
  // override algebras
  
  void doAllAlgebra(int lclw, const Float* ql_p, const Float* qnl_p, DEVOperatorKind src_op, DEVOperatorKind);
  void DiracAlgebra(const Float* qp1, const Float* qp2, int lclW,
		    int sour_op, int sink_op);
};

#endif // ! _INCLUDED_W_EXT_MESONSBE
 










CPS_END_NAMESPACE
