#ifndef CK_A2A
#define CK_A2A

#include<util/lattice.h>

#ifdef USE_GRID
#include <util/lattice/fgrid.h>
#endif

#ifdef USE_BFM
#include<util/lattice/bfm_mixed_solver.h>
#include <util/lattice/bfm_evo.h>
#include <alg/eigen/Krylov_5d.h>
#endif

#include<alg/a2a_arg.h>
#include<alg/a2a/CPSfield.h>
#include<alg/a2a/CPSfield_utils.h>
#include<alg/a2a/scfvectorptr.h>
#include<alg/a2a/utils.h>

#include<alg/a2a/a2a_params.h>
#include<alg/a2a/a2a_dilutions.h>
#include<alg/a2a/field_operation.h>
#include<alg/a2a/a2a_policies.h>
#include<alg/a2a/evec_interface.h>
#include<alg/a2a/grid_cgne_m_high.h>
CPS_START_NAMESPACE 

//If using SIMD, we don't want to vectorize across the time direction
template<typename FieldInputParamType>
struct checkSIMDparams{
  inline static void check(const FieldInputParamType &p){}
};
#ifdef USE_GRID
template<int Dimension>
struct checkSIMDparams<SIMDdims<Dimension> >{
  inline static void check(const SIMDdims<Dimension> &p){
    assert(p[3] == 1);
  }
};
#endif

template< typename mf_Policies>
class A2AvectorVfftw;

template< typename mf_Policies>
class A2AvectorV: public StandardIndexDilution, public mf_Policies::A2AvectorVpolicies{
public:
  typedef mf_Policies Policies;
  typedef typename Policies::FermionFieldType FermionFieldType;
  typedef typename FermionFieldType::FieldSiteType FieldSiteType;
  typedef typename FermionFieldType::InputParamType FieldInputParamType;
private:
  std::vector<PtrWrapper<FermionFieldType> > v;

public:
  typedef StandardIndexDilution DilutionType;

  A2AvectorV(const A2AArg &_args): StandardIndexDilution(_args){    
    v.resize(nv);
    //When computing V and W we can re-use previous V solutions as guesses. Set default to zero here so we have zero guess when no 
    //previously computed solutions
    this->allocInitializeFields(v,NullObject());
  }
  
  A2AvectorV(const A2AArg &_args, const FieldInputParamType &field_setup_params): StandardIndexDilution(_args){
    checkSIMDparams<FieldInputParamType>::check(field_setup_params);
    v.resize(nv);
    this->allocInitializeFields(v,field_setup_params);
  }
  
  static double Mbyte_size(const A2AArg &_args, const FieldInputParamType &field_setup_params);

  inline const FermionFieldType & getMode(const int i) const{ return *v[i]; }
  inline FermionFieldType & getMode(const int i){ return *v[i]; }  

  //Get a mode from the low mode part
  FermionFieldType & getVl(const int il){ return *v[il]; }
  const FermionFieldType & getVl(const int il) const{ return *v[il]; }

  //Get a mode from the high-mode part
  FermionFieldType & getVh(const int ih){ return *v[nl+ih]; }
  const FermionFieldType & getVh(const int ih) const{ return *v[nl+ih]; }

  //Get a particular site/spin/color element of a given mode 
  const FieldSiteType & elem(const int mode, const int x3d, const int t, const int spin_color, const int flavor) const{
    int x4d = v[mode]->threeToFour(x3d,t);
    return  *(v[mode]->site_ptr(x4d,flavor) + spin_color);
  }
  //Get a particular site/spin/color element of a given *native* (packed) mode. For V this does the same as the above
  inline const FieldSiteType & nativeElem(const int i, const int site, const int spin_color, const int flavor) const{
    return *(v[i]->site_ptr(site,flavor)+spin_color);
  }

  void importVl(const FermionFieldType &vv, const int il){
    *v[il] = vv;
  }
  void importVh(const FermionFieldType &vv, const int ih){
    *v[nl+ih] = vv;
  }

  //Set each float to a uniform random number in the specified range.
  //WARNING: Uses only the current RNG in LRG, and does not change this based on site. This is therefore only useful for testing*
  void testRandom(const Float &hi = 0.5, const Float &lo = -0.5){
    for(int i=0;i<nv;i++) v[i]->testRandom(hi,lo);
  }
  
};


template< typename mf_Policies>
class A2AvectorVfftw: public StandardIndexDilution, public mf_Policies::A2AvectorVfftwPolicies{  
public:
  typedef mf_Policies Policies;
  typedef typename Policies::FermionFieldType FermionFieldType;
  typedef typename FermionFieldType::FieldSiteType FieldSiteType;
  typedef typename FermionFieldType::InputParamType FieldInputParamType;

  #define VFFTW_ENABLE_IF_MANUAL_ALLOC(P) typename my_enable_if<  _equal<typename P::A2AvectorVfftwPolicies::FieldAllocStrategy,ManualAllocStrategy>::value , int>::type
private:
  std::vector<PtrWrapper<FermionFieldType> > v;

public:
  typedef StandardIndexDilution DilutionType;

  A2AvectorVfftw(const A2AArg &_args): StandardIndexDilution(_args){
    v.resize(nv);
    this->allocInitializeFields(v,NullObject());
  }
  A2AvectorVfftw(const A2AArg &_args, const FieldInputParamType &field_setup_params): StandardIndexDilution(_args){
    checkSIMDparams<FieldInputParamType>::check(field_setup_params);
    v.resize(nv);
    this->allocInitializeFields(v,field_setup_params);
  }

  static double Mbyte_size(const A2AArg &_args, const FieldInputParamType &field_setup_params);

  inline const FermionFieldType & getMode(const int i) const{ return *v[i]; }
  inline const FermionFieldType & getMode(const int i, const modeIndexSet &i_high_unmapped) const{ return getMode(i); }

  inline FermionFieldType & getMode(const int i){ return *v[i]; }
  
  //Set this object to be the threaded fast Fourier transform of the input field
  //Can optionally supply an object that performs a transformation on each mode prior to the FFT. 
  //We can use this to avoid intermediate storage for the gauge fixing and momentum phase application steps
  void fft(const A2AvectorV<Policies> &from, fieldOperation<FermionFieldType>* mode_preop = NULL);

  //Same as the above but allocates Vfft modes and deallocates V along the way to minimize memory usage. Only defined for manual alloc policies
  template<typename P = Policies>
  void destructivefft(A2AvectorV<P> &from, fieldOperation<typename P::FermionFieldType>* mode_preop = NULL, VFFTW_ENABLE_IF_MANUAL_ALLOC(P) = 0);
  
  void inversefft(A2AvectorV<Policies> &to, fieldOperation<FermionFieldType>* mode_postop = NULL) const;

  template<typename P = Policies>
  void destructiveInversefft(A2AvectorV<P> &to, fieldOperation<typename P::FermionFieldType>* mode_postop = NULL, VFFTW_ENABLE_IF_MANUAL_ALLOC(P) = 0);
  
  //For each mode, gauge fix, apply the momentum factor, then perform the FFT and store the result in this object
  void gaugeFixTwistFFT(const A2AvectorV<Policies> &from, const int _p[3], Lattice &_lat){
    gaugeFixAndTwist<FermionFieldType> op(_p,_lat); fft(from, &op);
  }

  template<typename P=mf_Policies>
  void destructiveGaugeFixTwistFFT(A2AvectorV<Policies> &from, const int _p[3], Lattice &_lat, VFFTW_ENABLE_IF_MANUAL_ALLOC(P) = 0){
    gaugeFixAndTwist<FermionFieldType> op(_p,_lat); destructivefft(from, &op);
  }
  
  //Unapply the phase and gauge fixing to give back a V vector
  void unapplyGaugeFixTwistFFT(A2AvectorV<Policies> &to, const int _p[3], Lattice &_lat) const{
    reverseGaugeFixAndTwist<FermionFieldType> op(_p,_lat); inversefft(to, &op);
  }
  
  template<typename P=mf_Policies>
  void destructiveUnapplyGaugeFixTwistFFT(A2AvectorV<Policies> &to, const int _p[3], Lattice &_lat, VFFTW_ENABLE_IF_MANUAL_ALLOC(P) = 0){
    reverseGaugeFixAndTwist<FermionFieldType> op(_p,_lat); destructiveInversefft(to, &op);
  }

  //Use the relations between FFTs to obtain the FFT for a chosen quark momentum
  //With G-parity BCs there are 2 disjoint sets of momenta hence there are 2 base FFTs
  void getTwistedFFT(const int p[3], A2AvectorVfftw<Policies> const *base_p, A2AvectorVfftw<Policies> const *base_m = NULL);

  void shiftFieldsInPlace(const std::vector<int> &shift);

  //A version of the above that directly shifts the base Wfftw rather than outputting into a separate storage
  //Returns the pointer to the Wfftw acted upon and the *shift required to restore the Wfftw to it's original form* (use shiftFieldsInPlace to restore)
  
  static std::pair< A2AvectorVfftw<mf_Policies>*, std::vector<int> > inPlaceTwistedFFT(const int p[3], A2AvectorVfftw<mf_Policies> *base_p, A2AvectorVfftw<mf_Policies> *base_m = NULL);
  
  const FieldSiteType & elem(const int mode, const int x3d, const int t, const int spin_color, const int flavor) const{
    int site = v[mode]->threeToFour(x3d,t);
    return *(v[mode]->site_ptr(site,flavor) + spin_color);
  }
  //Get a particular site/spin/color element of a given 'native' (packed) mode. For V this does the same thing as the above
  inline const FieldSiteType & nativeElem(const int i, const int site, const int spin_color, const int flavor) const{
    return *(v[i]->site_ptr(site,flavor)+spin_color);
  }

  //i_high_unmapped is the index i unmapped to its high mode sub-indices (if it is a high mode of course!)
  inline SCFvectorPtr<FieldSiteType> getFlavorDilutedVect(const int i, const modeIndexSet &i_high_unmapped, const int p3d, const int t) const{
    const FermionFieldType &field = getMode(i);
    const int x4d = field.threeToFour(p3d,t);
    FieldSiteType const *f0 = field.site_ptr(x4d,0);
    return SCFvectorPtr<FieldSiteType>(f0,f0+field.flav_offset());
  }
  //Return the pointer stride for between 3d coordinates for a given mode index and flavor. Relies on the dimension policy implementing dimpol_site_stride_3d
  inline int siteStride3D(const int i, const modeIndexSet &i_high_unmapped, const int f) const{
    const FermionFieldType &field = getMode(i);
    return field.dimpol_site_stride_3d()*field.siteSize();
  }
  
  //Replace this vector with the average of this another vector, 'with'
  void average(const A2AvectorVfftw<Policies> &with, const bool &parallel = true){
    if( !paramsEqual(with) ) ERR.General("A2AvectorVfftw","average","Second field must share the same underlying parameters\n");
    for(int i=0;i<nv;i++) v[i]->average(with.v[i]);
  }
  //Set each float to a uniform random number in the specified range.
  //WARNING: Uses only the current RNG in LRG, and does not change this based on site. This is therefore only useful for testing*
  void testRandom(const Float &hi = 0.5, const Float &lo = -0.5){
    for(int i=0;i<nv;i++) v[i]->testRandom(hi,lo);
  }

  template<typename extPolicies>
  void importFields(const A2AvectorVfftw<extPolicies> &r){
    if( !paramsEqual(r) ) ERR.General("A2AvectorVfftw","importFields","External field-vector must share the same underlying parameters\n");
    for(int i=0;i<nv;i++) v[i]->importField(r.getMode(i));
  }  

};


template< typename mf_Policies>
class A2AvectorW: public FullyPackedIndexDilution, public mf_Policies::A2AvectorWpolicies{
public:
  typedef mf_Policies Policies;
  typedef typename Policies::FermionFieldType FermionFieldType;
  typedef typename Policies::ComplexFieldType ComplexFieldType;
  typedef typename my_enable_if< _equal<typename FermionFieldType::FieldSiteType, typename ComplexFieldType::FieldSiteType>::value,  typename FermionFieldType::FieldSiteType>::type FieldSiteType;
  typedef typename my_enable_if< _equal<typename FermionFieldType::InputParamType, typename ComplexFieldType::InputParamType>::value,  typename FermionFieldType::InputParamType>::type FieldInputParamType;
private:
  std::vector<PtrWrapper<FermionFieldType> > wl; //The low mode part of the W field, comprised of nl fermion fields
  std::vector<PtrWrapper<ComplexFieldType> > wh; //The high mode random part of the W field, comprised of nhits complex scalar fields. Note: the dilution is performed later

  //Generate the wh field. We store in a compact notation that knows nothing about any dilution we apply when generating V from this
  //For reproducibility we want to generate the wh field in the same order that Daiqian did originally. Here nhit random numbers are generated for each site/flavor
  void setWhRandom(const RandomType &type);
  
public:
  typedef FullyPackedIndexDilution DilutionType;

  A2AvectorW(const A2AArg &_args): FullyPackedIndexDilution(_args){
    wl.resize(nl); this->allocInitializeLowModeFields(wl,NullObject());
    wh.resize(nhits); this->allocInitializeHighModeFields(wh,NullObject());
  }
  A2AvectorW(const A2AArg &_args, const FieldInputParamType &field_setup_params): FullyPackedIndexDilution(_args){
    checkSIMDparams<FieldInputParamType>::check(field_setup_params);
    wl.resize(nl); this->allocInitializeLowModeFields(wl,field_setup_params);
    wh.resize(nhits); this->allocInitializeHighModeFields(wh,field_setup_params);
  }

  static double Mbyte_size(const A2AArg &_args, const FieldInputParamType &field_setup_params);
  
  const FermionFieldType & getWl(const int i) const{ return *wl[i]; }
  const ComplexFieldType & getWh(const int hit) const{ return *wh[hit]; }

  FermionFieldType & getWl(const int i){ return *wl[i]; }
  ComplexFieldType & getWh(const int hit){ return *wh[hit]; }
  
  void importWl(const FermionFieldType &wlin, const int i){
    *wl[i] = wlin;
  }
  void importWh(const ComplexFieldType &whin, const int hit){
    *wh[hit] = whin;
  }

#ifdef USE_GRID
  //Generic Grid VW compute interface that can use either Grid or BFM-computed eigenvectors

  //Compute the low mode part of the W and V vectors.
  void computeVWlow(A2AvectorV<Policies> &V, Lattice &lat, EvecInterface<Policies> &evecs, const Float mass);

  //Compute the high mode parts of V and W. 
  void computeVWhigh(A2AvectorV<Policies> &V, Lattice &lat, EvecInterface<Policies> &evecs, const Float mass, const Float residual, const int max_iter);
#endif

#if defined(USE_BFM_LANCZOS)
  //In the Lanczos class you can choose to store the vectors in single precision (despite the overall precision, which is fixed to double here)
  //Set 'singleprec_evecs' if this has been done
  void computeVWlow(A2AvectorV<Policies> &V, Lattice &lat, BFM_Krylov::Lanczos_5d<double> &eig, bfm_evo<double> &dwf, bool singleprec_evecs);

  //singleprec_evecs specifies whether the input eigenvectors are stored in single precision
  //You can optionally pass a single precision bfm instance, which if given will cause the underlying CG to be performed in mixed precision.
  //WARNING: if using the mixed precision solve, the eigenvectors *MUST* be in single precision (there is a runtime check)
  void computeVWhigh(A2AvectorV<Policies> &V, BFM_Krylov::Lanczos_5d<double> &eig, bool singleprec_evecs, Lattice &lat, bfm_evo<double> &dwf_d, bfm_evo<float> *dwf_fp = NULL);

  void computeVW(A2AvectorV<Policies> &V, Lattice &lat, BFM_Krylov::Lanczos_5d<double> &eig, bool singleprec_evecs, bfm_evo<double> &dwf_d, bfm_evo<float> *dwf_fp = NULL){
    computeVWlow(V,lat,eig,dwf_d,singleprec_evecs);
    computeVWhigh(V,eig,singleprec_evecs,lat,dwf_d,dwf_fp);
  }
#endif


#if defined(USE_GRID_LANCZOS)
  void computeVWlow(A2AvectorV<Policies> &V, Lattice &lat, const std::vector<typename Policies::GridFermionField> &evec, const std::vector<Grid::RealD> &eval, const double mass);

  void computeVWhigh(A2AvectorV<Policies> &V, Lattice &lat, const std::vector<typename Policies::GridFermionField> &evec, const std::vector<Grid::RealD> &eval, const double mass, const Float residual, const int max_iter);

  void computeVW(A2AvectorV<Policies> &V, Lattice &lat, const std::vector<typename Policies::GridFermionField> &evec, const std::vector<Grid::RealD> &eval, const double mass, const Float high_mode_residual, const int high_mode_max_iter){
    computeVWlow(V,lat,evec,eval,mass);
    computeVWhigh(V,lat,evec,eval,mass,high_mode_residual,high_mode_max_iter);
  }

  //Single-precision variants (use mixed_CG internally)
  void computeVWlow(A2AvectorV<Policies> &V, Lattice &lat, const std::vector<typename Policies::GridFermionFieldF> &evec, const std::vector<Grid::RealD> &eval, const double mass);

  void computeVWhigh(A2AvectorV<Policies> &V, Lattice &lat, const std::vector<typename Policies::GridFermionFieldF> &evec, const std::vector<Grid::RealD> &eval, const double mass, const Float residual, const int max_iter);

  void computeVW(A2AvectorV<Policies> &V, Lattice &lat, const std::vector<typename Policies::GridFermionFieldF> &evec, const std::vector<Grid::RealD> &eval, const double mass, const Float high_mode_residual, const int high_mode_max_iter){
    computeVWlow(V,lat,evec,eval,mass);
    computeVWhigh(V,lat,evec,eval,mass,high_mode_residual,high_mode_max_iter);
  }
#endif

  //Get the diluted source with StandardIndex high-mode index dil_id.
  //We use the same set of random numbers for each spin and dilution as we do not need to rely on stochastic cancellation to separate them
  //For legacy reasons we use different random numbers for the two G-parity flavors, although this is not strictly necessary
  //Here dil_id is the combined spin-color/flavor/hit/tblock index
  template<typename TargetFermionFieldType>
  void getDilutedSource(TargetFermionFieldType &into, const int dil_id) const;

  //When gauge fixing prior to taking the FFT it is necessary to uncompact the wh field in the spin-color index, as these indices are acted upon by the gauge fixing
  //(I suppose technically only the color indices need uncompacting; this might be considered as a future improvement)
  void getSpinColorDilutedSource(FermionFieldType &into, const int hit, const int sc_id) const;

  //The spincolor, flavor and timeslice dilutions are packed so we must treat them differently
  //Mode is a full 'StandardIndex', (unpacked mode index)
  const FieldSiteType & elem(const int mode, const int x3d, const int t, const int spin_color, const int flavor) const{
    static FieldSiteType zero(0.0);
    if(mode < nl){
      int site = getWl(mode).threeToFour(x3d,t);
      return *(getWl(mode).site_ptr(site,flavor) + spin_color);
    }else{
      int mode_hit, mode_tblock, mode_spin_color,mode_flavor;
      const StandardIndexDilution &dilfull = static_cast<StandardIndexDilution const&>(*this);
      dilfull.indexUnmap(mode-nl,mode_hit,mode_tblock,mode_spin_color,mode_flavor);
      //flavor and time block indices match those of the mode, the result is zero
      int tblock = (t+GJP.TnodeSites()*GJP.TnodeCoor())/args.src_width;
      if(spin_color != mode_spin_color || flavor != mode_flavor || tblock != mode_tblock) return zero;
      int site = getWh(mode_hit).threeToFour(x3d,t);
      return *(getWh(mode_hit).site_ptr(site,flavor)); //we use different random fields for each time and flavor, although we didn't have to
    }
  }
  //Get a particular site/spin/color element of a given *native* (packed) mode 
  inline const FieldSiteType & nativeElem(const int i, const int site, const int spin_color, const int flavor) const{
    return i < nl ? 
      *(wl[i]->site_ptr(site,flavor)+spin_color) :
      *(wh[i-nl]->site_ptr(site,flavor)); //we use different random fields for each time and flavor, although we didn't have to
  }

  //Set each float to a uniform random number in the specified range.
  //WARNING: Uses only the current RNG in LRG, and does not change this based on site. This is therefore only useful for testing*
  void testRandom(const Float &hi = 0.5, const Float &lo = -0.5){
    for(int i=0;i<nl;i++) wl[i]->testRandom(hi,lo);
    for(int i=0;i<nhits;i++) wh[i]->testRandom(hi,lo);
  }

};


template< typename mf_Policies>
class A2AvectorWfftw: public TimeFlavorPackedIndexDilution, public mf_Policies::A2AvectorWfftwPolicies{
public:
  typedef mf_Policies Policies;
  typedef typename Policies::FermionFieldType FermionFieldType;
  typedef typename Policies::ComplexFieldType ComplexFieldType;
  typedef typename FermionFieldType::FieldSiteType FieldSiteType;
  typedef typename my_enable_if< _equal<typename FermionFieldType::InputParamType, typename ComplexFieldType::InputParamType>::value,  typename FermionFieldType::InputParamType>::type FieldInputParamType;

#define WFFTW_ENABLE_IF_MANUAL_ALLOC(P) typename my_enable_if<  _equal<typename P::A2AvectorWfftwPolicies::FieldAllocStrategy,ManualAllocStrategy>::value , int>::type
private:

  std::vector<PtrWrapper<FermionFieldType> > wl;
  std::vector<PtrWrapper<FermionFieldType> > wh; //these have been diluted in spin/color but not the other indices, hence there are nhit * 12 fields here (spin/color index changes fastest in mapping)

  FieldSiteType zerosc[12];
public:
  typedef TimeFlavorPackedIndexDilution DilutionType;

  A2AvectorWfftw(const A2AArg &_args): TimeFlavorPackedIndexDilution(_args){
    wl.resize(nl); this->allocInitializeLowModeFields(wl,NullObject());
    wh.resize(12*nhits); this->allocInitializeHighModeFields(wh,NullObject());
    for(int i=0;i<12;i++) CPSsetZero(zerosc[i]);
  }
  A2AvectorWfftw(const A2AArg &_args, const FieldInputParamType &field_setup_params): TimeFlavorPackedIndexDilution(_args){
    checkSIMDparams<FieldInputParamType>::check(field_setup_params);
    wl.resize(nl); this->allocInitializeLowModeFields(wl,field_setup_params);
    wh.resize(12*nhits); this->allocInitializeHighModeFields(wh,field_setup_params);
    for(int i=0;i<12;i++) CPSsetZero(zerosc[i]);
  }

  static double Mbyte_size(const A2AArg &_args, const FieldInputParamType &field_setup_params);
  
  inline const FermionFieldType & getWl(const int i) const{ return *wl[i]; }
  inline const FermionFieldType & getWh(const int hit, const int spin_color) const{ return *wh[spin_color + 12*hit]; }

  inline const FermionFieldType & getMode(const int i) const{ return i < nl ? *wl[i] : *wh[i-nl]; }

  inline FermionFieldType & getWl(const int i){ return *wl[i]; }
  inline FermionFieldType & getWh(const int hit, const int spin_color){ return *wh[spin_color + 12*hit]; }

  inline FermionFieldType & getMode(const int i){ return i < nl ? *wl[i] : *wh[i-nl]; }

  //This version allows for the possibility of a different high mode mapping for the index i by passing the unmapped indices: for i>=nl the modeIndexSet is used to obtain the appropriate mode 
  inline const FermionFieldType & getMode(const int i, const modeIndexSet &i_high_unmapped) const{ return i >= nl ? getWh(i_high_unmapped.hit, i_high_unmapped.spin_color): getWl(i); }
  
  //Set this object to be the threaded fast Fourier transform of the input field
  //Can optionally supply an object that performs a transformation on each mode prior to the FFT. 
  //We can use this to avoid intermediate storage for the gauge fixing and momentum phase application steps
  void fft(const A2AvectorW<Policies> &from, fieldOperation<FermionFieldType>* mode_preop = NULL);

  template<typename P = mf_Policies>
  void destructivefft(A2AvectorW<mf_Policies> &from, fieldOperation<FermionFieldType>* mode_preop = NULL, WFFTW_ENABLE_IF_MANUAL_ALLOC(P) = 0);
  
  void inversefft(A2AvectorW<Policies> &to, fieldOperation<FermionFieldType>* mode_postop = NULL) const;

  template<typename P=mf_Policies>
  void destructiveInversefft(A2AvectorW<mf_Policies> &to, fieldOperation<FermionFieldType>* mode_postop = NULL, WFFTW_ENABLE_IF_MANUAL_ALLOC(P) = 0);
  
  //For each mode, gauge fix, apply the momentum factor, then perform the FFT and store the result in this object
  void gaugeFixTwistFFT(const A2AvectorW<Policies> &from, const int _p[3], Lattice &_lat){
    gaugeFixAndTwist<FermionFieldType> op(_p,_lat); fft(from, &op);
  }

  template<typename P=mf_Policies>
  void destructiveGaugeFixTwistFFT(A2AvectorW<Policies> &from, const int _p[3], Lattice &_lat, WFFTW_ENABLE_IF_MANUAL_ALLOC(P) = 0){
    gaugeFixAndTwist<FermionFieldType> op(_p,_lat); destructivefft(from, &op);
  }
  
  //Unapply the phase and gauge fixing to give back a V vector
  void unapplyGaugeFixTwistFFT(A2AvectorW<Policies> &to, const int _p[3], Lattice &_lat) const{
    reverseGaugeFixAndTwist<FermionFieldType> op(_p,_lat); inversefft(to, &op);
  }

  template<typename P=mf_Policies>
  void destructiveUnapplyGaugeFixTwistFFT(A2AvectorW<Policies> &to, const int _p[3], Lattice &_lat, WFFTW_ENABLE_IF_MANUAL_ALLOC(P) = 0){
    reverseGaugeFixAndTwist<FermionFieldType> op(_p,_lat); destructiveInversefft(to, &op);
  }

  //Use the relations between FFTs to obtain the FFT for a chosen quark momentum
  //With G-parity BCs there are 2 disjoint sets of momenta hence there are 2 base FFTs
  void getTwistedFFT(const int p[3], A2AvectorWfftw<Policies> const *base_p, A2AvectorWfftw<Policies> const *base_m = NULL);

  void shiftFieldsInPlace(const std::vector<int> &shift);

  //A version of the above that directly shifts the base Wfftw rather than outputting into a separate storage
  //Returns the pointer to the Wfftw acted upon and the *shift required to restore the Wfftw to it's original form* (use shiftFieldsInPlace to restore)
  
  static std::pair< A2AvectorWfftw<mf_Policies>*, std::vector<int> > inPlaceTwistedFFT(const int p[3], A2AvectorWfftw<mf_Policies> *base_p, A2AvectorWfftw<mf_Policies> *base_m = NULL);
  
  //The flavor and timeslice dilutions are still packed so we must treat them differently
  //Mode is a full 'StandardIndex', (unpacked mode index)
  const FieldSiteType & elem(const int mode, const int x3d, const int t, const int spin_color, const int flavor) const{
    static FieldSiteType zero(0.0);
    if(mode < nl){
      int site = getWl(mode).threeToFour(x3d,t);
      return *(getWl(mode).site_ptr(site,flavor) + spin_color);
    }else{
      int mode_hit, mode_tblock, mode_spin_color,mode_flavor;
      const StandardIndexDilution &dilfull = static_cast<StandardIndexDilution const&>(*this);
      dilfull.indexUnmap(mode-nl,mode_hit,mode_tblock,mode_spin_color,mode_flavor);
      //flavor and time block indices match those of the mode, the result is zero
      int tblock = (t+GJP.TnodeSites()*GJP.TnodeCoor())/args.src_width;
      if(flavor != mode_flavor || tblock != mode_tblock) return zero;

      int site = getWh(mode_hit,mode_spin_color).threeToFour(x3d,t);
      return *(getWh(mode_hit,mode_spin_color).site_ptr(site,flavor) +spin_color); //because we multiplied by an SU(3) matrix, the field is not just a delta function in spin/color
    }
  }
  //Get a particular site/spin/color element of a given *native* (packed) mode 
  inline const FieldSiteType & nativeElem(const int i, const int site, const int spin_color, const int flavor) const{
    return i < nl ? 
      *(wl[i]->site_ptr(site,flavor)+spin_color) :
      *(wh[i-nl]->site_ptr(site,flavor)+spin_color); //spin_color index diluted out.
  }



  //Replace this vector with the average of this another vector, 'with'
  void average(const A2AvectorWfftw<Policies> &with, const bool &parallel = true){
    if( !paramsEqual(with) ) ERR.General("A2AvectorWfftw","average","Second field must share the same underlying parameters\n");
    for(int i=0;i<nl;i++) wl[i]->average(*with.wl[i]);
    for(int i=0;i<12*nhits;i++) wh[i]->average(*with.wh[i]);
  }

  //Set each float to a uniform random number in the specified range.
  //WARNING: Uses only the current RNG in LRG, and does not change this based on site. This is therefore only useful for testing*
  void testRandom(const Float &hi = 0.5, const Float &lo = -0.5){
    for(int i=0;i<nl;i++) wl[i]->testRandom(hi,lo);
    for(int i=0;i<12*nhits;i++) wh[i]->testRandom(hi,lo);
  }

  //BELOW are for use by the meson field
  //Meson field W-type indices are described in terms of the timePacked dilution index , where flavor has been diluted out in the process of computing the meson fields
  //This method performs the flavor dilution 'in-place' (i.e. without actually unpacking into another object). 
  //'site' is a local canonical-ordered, packed four-vector
  //i_high_unmapped is the index i unmapped to its high mode sub-indices (if it is a high mode of course!)

  inline SCFvectorPtr<FieldSiteType> getFlavorDilutedVect(const int i, const modeIndexSet &i_high_unmapped, const int p3d, const int t) const{
    const FermionFieldType &field = i >= nl ? getWh(i_high_unmapped.hit, i_high_unmapped.spin_color): getWl(i);
    bool zero_hint[2] = {false,false};
    if(i >= nl) zero_hint[ !i_high_unmapped.flavor ] = true;

    const int x4d = field.threeToFour(p3d,t);
    return SCFvectorPtr<FieldSiteType>(zero_hint[0] ? &zerosc[0] : field.site_ptr(x4d,0), zero_hint[1] ? &zerosc[0] : field.site_ptr(x4d,1), zero_hint[0], zero_hint[1]);
  }
  //Return the pointer stride for between 3d coordinates for a given mode index and flavor. Relies on the dimension policy implementing dimpol_site_stride_3d
  inline int siteStride3D(const int i, const modeIndexSet &i_high_unmapped, const int f) const{ 
    const FermionFieldType &field = i >= nl ? getWh(i_high_unmapped.hit, i_high_unmapped.spin_color): getWl(i);
    bool zero_hint[2] = {false,false};
    if(i >= nl) zero_hint[ !i_high_unmapped.flavor ] = true;
    return zero_hint[f] ? 0 : field.dimpol_site_stride_3d()*field.siteSize();
  }

  template<typename extPolicies>
  void importFields(const A2AvectorWfftw<extPolicies> &r){
    if( !paramsEqual(r) ) ERR.General("A2AvectorWfftw","importFields","External field-vector must share the same underlying parameters\n");
    for(int i=0;i<nl;i++) wl[i]->importField(r.getWl(i));
    for(int i=0;i<12*nhits;i++) wh[i]->importField(r.getWh(i/12,i%12));
  }  

};






#include <alg/a2a/a2a_impl.h>

#ifdef USE_GRID
#include<alg/a2a/evec_interface_impl.h>
#endif

//Can do Lanczos in BFM or Grid, and A2A in BFM or Grid. I have a BFM Lanczos -> Grid interface

#if defined(USE_BFM_A2A)
# warning "Using BFM A2A"

# ifndef USE_BFM
#  error "Require BFM for USE_BFM_A2A"
# endif

# ifdef USE_GRID_LANCZOS
#  error "No Grid Lanczos -> BFM A2A interface implemented"
# endif

# include <alg/a2a/a2a_impl_vwbfm.h>

#elif defined(USE_GRID_A2A)
# warning "Using Grid A2A"

# ifndef USE_GRID
#  error "Require Grid for USE_GRID_A2A"
# endif

# if defined(USE_BFM_LANCZOS) && !defined(USE_BFM)
#  error "BFM Lanczos -> Grid A2A interface requires BFM!"
# endif

# include <alg/a2a/a2a_impl_vwgrid.h>

#else

# error "Need either BFM or Grid to compute A2A vectors"

#endif


CPS_END_NAMESPACE

#endif
