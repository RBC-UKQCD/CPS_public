#ifndef ALG_GPARITY_CONTRACT_H
#define ALG_GPARITY_CONTRACT_H

#include<config.h>
#include <alg/alg_base.h>
#include <util/lattice.h>
#include <alg/correlationfunction.h>
#include <alg/gparity_contract_arg.h>
#include <alg/common_arg.h>
#include <util/spincolorflavormatrix.h>
#include <vector>
#include <utility>
#include <map>

CPS_START_NAMESPACE

//Class for determining allowed sink momenta
class QuarkMomCombination{
  const char* cname;
  std::vector< std::vector<std::vector<int> > > quark_mom;
  std::vector<std::pair<PropagatorContainer const*,bool> > props;

  void mom_comb_recurse(std::vector<std::vector<int> > &into, std::vector<int> &cur_branch, const int &depth) const;

  void calc_allowed_combinations();
  
  std::vector<std::vector<int> > allowed_combinations;
  bool allowed_comb_calculated;
  int chosen_momcomb;
public:
  QuarkMomCombination(): cname("QuarkMomCombination"), chosen_momcomb(-1), allowed_comb_calculated(false){};
  void reset();
  void add_prop(const QPropWcontainer &prop, const bool &conj);

  void set_desired_momentum(const std::vector<int> &what);
  void set_desired_momentum(const int *what);

 //see if momentum combination 'what' is a valid combination for this contraction
  bool contains(const std::vector<int> &what);
  bool contains(const int *what);

  std::vector<Float> get_p() const;

  //Calculate the sink complex phase factor. Must have previously set a desired momentum combination
  Rcomplex phase_factor(const int *pos) const;
  Rcomplex phase_factor(const int &site) const;

  const std::vector<std::vector<int> > &get_allowed_momenta() const;
};

class ContractionQuarkMomCombination{
  std::vector<QuarkMomCombination> momcomb;
  std::map<int,int> contraction_map;
public:
  void add_contraction(const int &contraction, const QuarkMomCombination &cmomcomb);
  void same(const int &contraction, const int &same_as_contraction);
  void set_desired_momentum(const std::vector<int> &what);
  void set_desired_momentum(const int *what);

  std::vector<Float> get_p(const int &contraction) const;

  Complex phase_factor(const int &contraction, const int* pos) const;
  Complex phase_factor(const int &contraction, const int& site) const;
};

//Contraction class
class AlgGparityContract : public Alg{
public:
  //Left/right multiply by a gamma matrix structure in QDP-style conventions:
  //\Gamma(n) = \gamma_1^n1 \gamma_2^n2  \gamma_3^n3 \gamma_4^n4    where ni are bit fields. 
  static void qdp_gl(WilsonMatrix &wmat, const int &gidx);
  static void qdp_gr(WilsonMatrix &wmat, const int &gidx);
  static void qdp_gl(SpinColorFlavorMatrix &wmat, const int &gidx);
  static void qdp_gr(SpinColorFlavorMatrix &wmat, const int &gidx);

  //Coefficient when matrix is transposed or conjugated
  static Float qdp_gcoeff(const int &gidx, const bool &transpose, const bool &conj);
  static Float pauli_coeff(const int &pidx, const bool &transpose, const bool &conj);

  bool binary_write; //enable binary write where implemented
private:
  char *cname;
  GparityContractArg *args;
  
  void meson_LL_std(QPropWcontainer &prop, const int* sink_mom, const int &gamma_idx_1, const int &gamma_idx_2, FILE *fp);
  void meson_LL_gparity(QPropWcontainer &prop, const int* sink_mom, const int &gamma_idx_1, const int &gamma_idx_2, FILE *fp);

  void meson_HL_gparity(QPropWcontainer &prop_H, QPropWcontainer &prop_L, const int* sink_mom, const int &gamma_idx_1, const int &gamma_idx_2, FILE *fp);
  void meson_HL_std(QPropWcontainer &prop_H, QPropWcontainer &prop_L, const int* sink_mom, const int &gamma_idx_1, const int &gamma_idx_2, FILE *fp);

  void contract_OVVpAA_gparity(const ContractionTypeOVVpAA &args, const int &conf_idx);
  void contract_OVVpAA_std(const ContractionTypeOVVpAA &args, const int &conf_idx);
public:
  AlgGparityContract(Lattice & latt, CommonArg& c_arg, GparityContractArg& arg);
  AlgGparityContract(Lattice & latt, CommonArg& c_arg);

  static void global_coord(const int &site, int *into_vec);

  //Momentum phase. Momenta are in units of 2pi/L
  static Rcomplex phase_factor(const Float *p, const int* global_pos);
  static Rcomplex phase_factor(const Float *p, const int &site);

  void set_args(GparityContractArg *to){ args = to; }

  void enable_binary_write(const bool &val = true){ binary_write = val; } //enable binary write where implemented

  //Run the measurement programme defined by the 'job' input. 
  //By default it uses the args pointer provided when this object is constructed or set using set_args, but user can choose a different run script.
  void run(const int &conf_idx);
  void run(const int &conf_idx, const GparityContractArg& job);

  //Measure the quantity specified by the GparityMeasurement input (essentially a factory)
  void spectrum(const GparityMeasurement &measargs,const int &conf_idx);

  void contract_LL_mesons(const ContractionTypeLLMesons &args, const int &conf_idx);
  void contract_HL_mesons(const ContractionTypeHLMesons &args, const int &conf_idx);

  void contract_OVVpAA(const ContractionTypeOVVpAA &args, const int &conf_idx);
  void contract_OVVpAA_gparity(CorrelationFunction &corrfunc, const ContractionTypeOVVpAA &args);

  //Calculate \sum_x tr( M1 A^dag(x) M2 B(x) ) for all M1,M2, where M1 and M2 are spin(-flavor) matrices and A,B are propagors specified in the arguments
  void contract_all_bilinears(const ContractionTypeAllBilinears &args, const int &conf_idx);
  void contract_all_wallsink_bilinears_specific_momentum(const ContractionTypeAllWallSinkBilinearsSpecificMomentum &args, const int &conf_idx);

  void contract_fourier_prop(const ContractionTypeFourierProp &args, const int &conf_idx);
  void contract_bilinear_vertex(const ContractionTypeBilinearVertex &args, const int &conf_idx);
  void contract_quadrilinear_vertex(const ContractionTypeQuadrilinearVertex &args, const int &conf_idx);

  void measure_topological_charge(const ContractionTypeTopologicalCharge &args, const int &conf_idx);
  void measure_mres(const ContractionTypeMres &args, const int &conf_idx);
  void measure_mres(const ContractionTypeMres &args, CorrelationFunction &pion, CorrelationFunction &j5_q);
  void measure_mres_gparity(const ContractionTypeMres &args, CorrelationFunction &pion, CorrelationFunction &j5_q);

  void measure_wilson_flow(const ContractionTypeWilsonFlow &args, const int &conf_idx);

  void contract_a2a_bilinear(const ContractionTypeA2ABilinear &args, CorrelationFunction &result);
  void contract_a2a_bilinear(const ContractionTypeA2ABilinear &args, const int &conf);
};



// class MesonField2;
// class MFsource;

// #define TESTING

// class Gparity_KtoPiPi{
// #ifdef TESTING
// public:
// #else
// private:
// #endif

//   A2APropbfm *prop_H;
//   A2APropbfm *prop_L;
//   Lattice *lattice;

//   MFsource *pion_source;
//   MFsource *kaon_source;

//   MesonField2 *wdagL_S2_vL_pi1;
//   MesonField2 *wdagL_S2_vL_pi2;
//   MesonField2 *wdagL_g5_vH;
//   MesonField2 *wdagL_wH;

//   bool gparity_use_transconv_props;

//   //G*idx are indices with the following mapping
//   //0  M_{0,V-A} = F0 g^mu(1-g5)
//   //1  M_{0,V+A} = F0 g^mu(1+g5)
//   //2  M_{1,V-A} = -F1 g^mu(1-g5)
//   //3  M_{1,V+A} = -F1 g^mu(1+g5)

//   //First index is the Gidx, second is mu  
//   SpinColorFlavorMatrix Gamma[4][4];

//   static SpinColorFlavorMatrix _F0;
//   static SpinColorFlavorMatrix _F1;
//   static SpinColorFlavorMatrix g5;
//   static SpinColorFlavorMatrix S2;
//   static SpinColorFlavorMatrix unit;
//   static SpinColorFlavorMatrix gmu[4];
  
//   int t_size;

//   bool setup_called;

//   //Storage for results individual contractions. There are 4(Gamma1)*4(Gamma2) results
//   //This maps to an index within the results array
//   static inline int con_map(const int &gidx1, const int &gidx2){
//     return gidx1 + 4*gidx2;
//   }
//   static inline void con_inv_map(int idx, int &gidx1, int &gidx2){
//     gidx1 = idx%4; idx/=4;
//     gidx2 = idx;
//   }
//   //n_contract is the number of contractions
//   //cidx_start is the index of the first contraction that will be written into this results vector
//   void setup_resultvec(const int &n_contract, const int &cidx_start,  std::vector<CorrelationFunction> &rvec);

// public:
//   void setup(const ContractionTypeKtoPiPi &args, Lattice &lat); //Should be called before doing contractions

//   //We assume the pipi total momentum is zero. Define +p as the source momentum of pi1 and thus -p is the momentum of pi2.
//   //There are 2 choices of sink momentum for pi1/pi2: pi1(+p)pi2(-p) and pi1(-p)pi2(+p). Thus we have an extra input parameter
//   //pi1_snkmomsgn = '+' or '-'  for the sign of the pi1 sink momentum relative to the sign of its source momentum
//   void pipi(const int &t_sep_pion, char pi1_snkmomsgn, std::vector<CorrelationFunction> &into, std::vector<int> *tpi1_vals);

//   //t_sep_pi_k is the time separatation between the kaon and the closest of the two pions
//   //if tK_vals = NULL it will sum over all kaon source timeslices, otherwise it will restrict the sum to the set provided
//   void type1(const int &t_sep_pi_k, const int &t_sep_pion, std::vector<CorrelationFunction> &into, std::vector<int> *tK_vals = NULL);

//   //This version uses g5-hermiticity on the strange quark prop. Should be more efficient
//   void type1_propHg5conj(const int &t_sep_pi_k, const int &t_sep_pion, std::vector<CorrelationFunction> &into, std::vector<int> *tK_vals = NULL);

//   void type2(const int &t_sep_pi_k, const int &t_sep_pion, std::vector<CorrelationFunction> &into, std::vector<int> *tK_vals = NULL);
  
//   void type3(const int &t_sep_pi_k, const int &t_sep_pion, std::vector<CorrelationFunction> &into, std::vector<int> *tK_vals = NULL);

//   void type4(const int &t_sep_pi_k, const int &t_sep_pion, std::vector<CorrelationFunction> &into, std::vector<int> *tK_vals = NULL);

//   void psvertex_type3(const int &t_sep_pi_k, const int &t_sep_pion, std::vector<CorrelationFunction> &into, std::vector<int> *tK_vals = NULL);
//   void psvertex_type4(const int &t_sep_pi_k, const int &t_sep_pion, std::vector<CorrelationFunction> &into, std::vector<int> *tK_vals = NULL);

//   //Storage for vector of all contractions within a type. c indexes the correlation function within the type (0 for the first)
//   static inline int result_map(const int &c, const int &gidx1, const int &gidx2){
//     return con_map(gidx1,gidx2) + 4*4*c;
//   }
//   static inline void result_inv_map(int idx, int &c, int &gidx1, int &gidx2){
//     static const int sz = 4*4;
//     c = idx / sz;
//     con_inv_map(idx % sz, gidx1, gidx2);
//   } 


//   void run(const ContractionTypeKtoPiPi &args, Lattice &lat, const int &config);

//   Gparity_KtoPiPi(): wdagL_S2_vL_pi1(NULL),wdagL_S2_vL_pi2(NULL),wdagL_g5_vH(NULL),pion_source(NULL),kaon_source(NULL),setup_called(false){
//   }

//   ~Gparity_KtoPiPi(){
// #define DELETE_IT(M) if(M!=NULL) delete M
//     DELETE_IT(wdagL_S2_vL_pi1);
//     DELETE_IT(wdagL_S2_vL_pi2);
//     DELETE_IT(wdagL_g5_vH);
//     DELETE_IT(wdagL_wH);
//     DELETE_IT(pion_source);
//     DELETE_IT(kaon_source);
// #undef DELETE_IT
//   }
// };



template <typename T>
struct _multimom_helper{
  static void add_momenta(T &to, MomArg *momenta, const int &sz){
    const static Float pi_const = M_PI;
    std::vector<Float> mom(3);
    for(int i=0; i<sz; i++){
      mom[0] = momenta[i].p[0] * pi_const;
      mom[1] = momenta[i].p[1] * pi_const;
      mom[2] = momenta[i].p[2] * pi_const;
      to.add_momentum(mom);
    }
  }
  static void add_momenta(T &to, MomPairArg *momenta, const int &sz){
    const static Float pi_const = M_PI;
    std::pair< std::vector<Float>,std::vector<Float> > mom; mom.first.resize(3); mom.second.resize(3);
    std::vector<Float> mom2(3);
    for(int i=0; i<sz; i++){
      for(int j=0;j<3;j++){
	mom.first[j] = momenta[i].p1[j] * pi_const;
	mom.second[j] = momenta[i].p2[j] * pi_const;
      }
      to.add_momentum(mom);
    }
  }
};

CPS_END_NAMESPACE
#endif

