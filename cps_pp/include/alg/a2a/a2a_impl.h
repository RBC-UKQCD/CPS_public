//Implementations of methods in a2a.h

template<typename VWtype>
inline double VW_Mbyte_size(const A2AArg &_args, const typename VWtype::FieldInputParamType &field_setup_params){
  typedef typename VWtype::DilutionType DilutionType;
  typedef typename VWtype::FermionFieldType FermionFieldType;
  DilutionType dil(_args); const int sz = dil.getNmodes();
  double field_size = double(FermionFieldType::byte_size(field_setup_params))/(1024.*1024.);
  return sz * field_size;
}


template< typename mf_Policies>
double A2AvectorV<mf_Policies>::Mbyte_size(const A2AArg &_args, const FieldInputParamType &field_setup_params){
  return VW_Mbyte_size<A2AvectorV<mf_Policies> >(_args,field_setup_params);
}
template< typename mf_Policies>
double A2AvectorVfftw<mf_Policies>::Mbyte_size(const A2AArg &_args, const FieldInputParamType &field_setup_params){
  return VW_Mbyte_size<A2AvectorVfftw<mf_Policies> >(_args,field_setup_params);
}
template< typename mf_Policies>
double A2AvectorW<mf_Policies>::Mbyte_size(const A2AArg &_args, const FieldInputParamType &field_setup_params){
  FullyPackedIndexDilution dil(_args);
  double ffield_size = double(FermionFieldType::byte_size(field_setup_params))/(1024.*1024.);
  double cfield_size = double(ComplexFieldType::byte_size(field_setup_params))/(1024.*1024.);
  return dil.getNl() * ffield_size + dil.getNhits() * cfield_size;
}
template< typename mf_Policies>
double A2AvectorWfftw<mf_Policies>::Mbyte_size(const A2AArg &_args, const FieldInputParamType &field_setup_params){
  return VW_Mbyte_size<A2AvectorWfftw<mf_Policies> >(_args,field_setup_params);
}

struct VFFTfieldPolicyBasic{
  template<typename T>
  static inline void actionOutputMode(T &v, const int i){}
  template<typename T>
  static inline void actionInputMode(T &v, const int i){}
};
struct VFFTfieldPolicyAllocFree{
  template<typename T>
  static inline void actionOutputMode(T &v, const int i){
    v.allocMode(i);
  }
  template<typename T>
  static inline void actionInputMode(T &v, const int i){
    v.freeMode(i);
  }
};


template<typename OutputType, typename InputType, typename FFTfieldPolicy>
struct _V_fft_impl{
  typedef typename InputType::FermionFieldType FermionFieldType;
  
  static inline void fft(OutputType &to, InputType &from, fieldOperation<FermionFieldType>* mode_preop){
    if(!UniqueID()){ printf("Doing V FFT\n"); fflush(stdout); }
    typedef typename FermionFieldType::InputParamType FieldParamType;
    FieldParamType field_setup = from.getMode(0).getDimPolParams();  
    FermionFieldType tmp(field_setup);
  
    Float preop_time = 0;
    Float fft_time = 0;

    const bool fft_dirs[4] = {true,true,true,false};
  
    for(int mode=0;mode<from.getNmodes();mode++){
      FermionFieldType const* init_gather_from = &from.getMode(mode);
      if(mode_preop != NULL){
	Float dtime = dclock();
	(*mode_preop)(from.getMode(mode),tmp);
	init_gather_from = &tmp;
	preop_time += dclock()-dtime;
      }
      Float dtime = dclock();
      
      FFTfieldPolicy::actionOutputMode(to, mode); //alloc
      
#ifndef MEMTEST_MODE
      cps::fft_opt(to.getMode(mode), *init_gather_from, fft_dirs);
#endif
      fft_time += dclock() - dtime;

      FFTfieldPolicy::actionInputMode(from, mode); //free
    }
    if(!UniqueID()){ printf("Finishing V FFT\n"); fflush(stdout); }
    print_time("A2AvectorVfftw::fft","Preop",preop_time);
    print_time("A2AvectorVfftw::fft","FFT",fft_time);
  }
};



//Set this object to be the fast Fourier transform of the input field
//Can optionally supply an object mode_preop that performs a transformation on each mode prior to the FFT
template< typename mf_Policies>
void A2AvectorVfftw<mf_Policies>::fft(const A2AvectorV<mf_Policies> &from, fieldOperation<FermionFieldType>* mode_preop){
  _V_fft_impl<A2AvectorVfftw<mf_Policies>, const A2AvectorV<mf_Policies>, VFFTfieldPolicyBasic>::fft(*this,from,mode_preop);
}

template< typename mf_Policies>
template<typename P>
void A2AvectorVfftw<mf_Policies>::destructivefft(A2AvectorV<P> &from, fieldOperation<typename P::FermionFieldType>* mode_preop, VFFTW_ENABLE_IF_MANUAL_ALLOC(P) ){
  _V_fft_impl<A2AvectorVfftw<P>, A2AvectorV<P>, VFFTfieldPolicyAllocFree>::fft(*this,from,mode_preop);
}



template<typename OutputType, typename InputType, typename FFTfieldPolicy>
struct _V_invfft_impl{
  typedef typename InputType::FermionFieldType FermionFieldType;

  static inline void inversefft(OutputType &to, InputType &from, fieldOperation<FermionFieldType>* mode_postop){
    if(!UniqueID()){ printf("Doing V inverse FFT\n"); fflush(stdout); }
    typedef typename FermionFieldType::InputParamType FieldParamType;
    FieldParamType field_setup = from.getMode(0).getDimPolParams();  
    FermionFieldType tmp(field_setup);
  
    Float postop_time = 0;
    Float fft_time = 0;

    const bool fft_dirs[4] = {true,true,true,false};
    for(int mode=0;mode<from.getNmodes();mode++){
      //if(!UniqueID()) printf("Mode %d, memory before output alloc\n",mode);
      //printMem();
      
      FFTfieldPolicy::actionOutputMode(to, mode); //alloc

      //if(!UniqueID()) printf("Mode %d, memory after output alloc\n",mode);
      //printMem();
      
      FermionFieldType* out = mode_postop == NULL ? &to.getMode(mode) : &tmp;
    
      Float dtime = dclock();
#ifndef MEMTEST_MODE
      cps::fft_opt(*out, from.getMode(mode), fft_dirs, true);
#endif

      //if(!UniqueID()) printf("Mode %d, memory before input free\n",mode);
      //printMem();
      
      FFTfieldPolicy::actionInputMode(from, mode); //alloc

      //if(!UniqueID()) printf("Mode %d, memory after input free\n",mode);
      //printMem();
      
      if(mode_postop != NULL){
	Float dtime = dclock();
	(*mode_postop)(tmp,to.getMode(mode));
	postop_time += dclock()-dtime;
      }
      fft_time += dclock() - dtime;
      //printMem();
    }
    if(!UniqueID()){ printf("Finishing V invert FFT\n"); fflush(stdout); }
    print_time("A2AvectorVfftw::inversefft","FFT",fft_time);
    print_time("A2AvectorVfftw::inversefft","Postop",postop_time);
  }
};


template< typename mf_Policies>
void A2AvectorVfftw<mf_Policies>::inversefft(A2AvectorV<Policies> &to, fieldOperation<FermionFieldType>* mode_postop) const{
  _V_invfft_impl<A2AvectorV<Policies>, const A2AvectorVfftw<mf_Policies>, VFFTfieldPolicyBasic>::inversefft(to,*this,mode_postop);
}

template< typename mf_Policies>
template<typename P>
void A2AvectorVfftw<mf_Policies>::destructiveInversefft(A2AvectorV<P> &to, fieldOperation<typename P::FermionFieldType>* mode_postop, VFFTW_ENABLE_IF_MANUAL_ALLOC(P) ){
  _V_invfft_impl<A2AvectorV<Policies>, A2AvectorVfftw<mf_Policies>, VFFTfieldPolicyAllocFree>::inversefft(to,*this,mode_postop);
}


struct WFFTfieldPolicyBasic{
  template<typename T>
  static inline void actionOutputLowMode(T &v, const int i){}
  template<typename T>
  static inline void actionOutputHighMode(T &v, const int i){}
  
  template<typename T>
  static inline void actionInputLowMode(T &v, const int i){}
  template<typename T>
  static inline void actionInputHighMode(T &v, const int i){}
};
struct WFFTfieldPolicyAllocFree{
  template<typename T>
  static inline void actionOutputLowMode(T &v, const int i){
    v.allocLowMode(i);
  }
  template<typename T>
  static inline void actionOutputHighMode(T &v, const int i){
    v.allocHighMode(i);
  }
  
  template<typename T>
  static inline void actionInputLowMode(T &v, const int i){
    v.freeLowMode(i);
  }
  template<typename T>
  static inline void actionInputHighMode(T &v, const int i){
    v.freeHighMode(i);
  }
};


template<typename OutputType, typename InputType, typename FFTfieldPolicy>
struct _W_fft_impl{
  typedef typename InputType::FermionFieldType FermionFieldType;

  inline static void fft(OutputType &to, InputType &from, fieldOperation<FermionFieldType>* mode_preop){
    if(!UniqueID()){ printf("Doing W FFT\n"); fflush(stdout); }
    typedef typename FermionFieldType::InputParamType FieldParamType;
    FieldParamType field_setup = from.getWh(0).getDimPolParams();  
    FermionFieldType tmp(field_setup), tmp2(field_setup);

    Float preop_time = 0;
    Float fft_time = 0;

    const bool fft_dirs[4] = {true,true,true,false};
  
    //Do wl
    for(int mode=0;mode<from.getNl();mode++){
      FermionFieldType const* init_gather_from = &from.getWl(mode);
      if(mode_preop != NULL){
	Float dtime = dclock();
	(*mode_preop)(from.getWl(mode),tmp);
	init_gather_from = &tmp;
	preop_time += dclock()-dtime;
      }
      FFTfieldPolicy::actionOutputLowMode(to, mode); //alloc
      Float dtime = dclock();
#ifndef MEMTEST_MODE
      cps::fft_opt(to.getWl(mode), *init_gather_from, fft_dirs);
#endif
      fft_time += dclock() - dtime;
      FFTfieldPolicy::actionInputLowMode(from, mode); //free
    }
    //Do wh. First we need to uncompact the spin/color index as this is acted upon by the operator
    for(int hit=0;hit<from.getNhits();hit++){
      for(int sc=0;sc<12;sc++){ //spin/color dilution index
	from.getSpinColorDilutedSource(tmp2,hit,sc);
	FermionFieldType* init_gather_from = &tmp2;
	if(mode_preop != NULL){
	  Float dtime = dclock();
	  (*mode_preop)(tmp2,tmp);
	  init_gather_from = &tmp;
	  preop_time += dclock()-dtime;
	}
	Float dtime = dclock();
	FFTfieldPolicy::actionOutputHighMode(to, sc+12*hit); //alloc
#ifndef MEMTEST_MODE
	cps::fft_opt(to.getWh(hit,sc), *init_gather_from, fft_dirs);
#endif
	fft_time += dclock()-dtime;
      }
      FFTfieldPolicy::actionInputHighMode(from, hit); //free
    }
    if(!UniqueID()){ printf("Finishing W FFT\n"); fflush(stdout); }
    print_time("A2AvectorWfftw::fft","Preop",preop_time);
    print_time("A2AvectorWfftw::fft","FFT",fft_time);
  }
};



//Set this object to be the fast Fourier transform of the input field
//Can optionally supply an object mode_preop that performs a transformation on each mode prior to the FFT
template< typename mf_Policies>
void A2AvectorWfftw<mf_Policies>::fft(const A2AvectorW<mf_Policies> &from, fieldOperation<FermionFieldType>* mode_preop){
  _W_fft_impl<A2AvectorWfftw<mf_Policies>, const A2AvectorW<mf_Policies>, WFFTfieldPolicyBasic>::fft(*this,from,mode_preop);
}

template< typename mf_Policies>
template<typename P>
void A2AvectorWfftw<mf_Policies>::destructivefft(A2AvectorW<mf_Policies> &from, fieldOperation<FermionFieldType>* mode_preop, WFFTW_ENABLE_IF_MANUAL_ALLOC(P) ){
  _W_fft_impl<A2AvectorWfftw<mf_Policies>, A2AvectorW<mf_Policies>, WFFTfieldPolicyAllocFree>::fft(*this,from,mode_preop);
}

template<typename OutputType, typename InputType, typename FFTfieldPolicy>
struct _W_invfft_impl{
  typedef typename InputType::FermionFieldType FermionFieldType;

  static inline void inversefft(OutputType &to, InputType &from, fieldOperation<FermionFieldType>* mode_postop){
    if(!UniqueID()){ printf("Doing W inverse FFT\n"); fflush(stdout); }
    typedef typename FermionFieldType::InputParamType FieldParamType;
    FieldParamType field_setup = from.getWh(0,0).getDimPolParams();  
    FermionFieldType tmp(field_setup), tmp2(field_setup);

    Float postop_time = 0;
    Float fft_time = 0;

    const bool fft_dirs[4] = {true,true,true,false};
  
    //Do wl
    for(int mode=0;mode<from.getNl();mode++){
      FFTfieldPolicy::actionOutputLowMode(to, mode); //alloc
      FermionFieldType * unfft_to = mode_postop == NULL ? &to.getWl(mode) : &tmp;

      Float dtime = dclock();
#ifndef MEMTEST_MODE
      cps::fft_opt(*unfft_to, from.getWl(mode), fft_dirs, true);
#endif
      fft_time += dclock() - dtime;

      if(mode_postop != NULL){
	Float dtime = dclock();
	(*mode_postop)(tmp,to.getWl(mode));
	postop_time += dclock()-dtime;
      }
      FFTfieldPolicy::actionInputLowMode(from, mode); //free
    }
    //Do wh. First we need to uncompact the spin/color index as this is acted upon by the operator
    for(int hit=0;hit<from.getNhits();hit++){
      FFTfieldPolicy::actionOutputHighMode(to, hit); //alloc
      typename InputType::ComplexFieldType & to_hit = to.getWh(hit);
    
      const int sc = 0;
      FermionFieldType * compress = mode_postop == NULL ? &tmp2 : &tmp;
      Float dtime = dclock();
#ifndef MEMTEST_MODE
      cps::fft_opt(tmp2, from.getWh(hit,sc), fft_dirs, true);
#endif
      fft_time += dclock()-dtime;

      if(mode_postop != NULL){
	Float dtime = dclock();
	(*mode_postop)(tmp2,tmp);
	postop_time += dclock()-dtime;
      }
      //Should give a multiple of the 12-component unit vector with 1 on index sc
#pragma omp parallel for
      for(int i=0;i<to_hit.nfsites();i++)
	*(to_hit.fsite_ptr(i)) = *(compress->fsite_ptr(i) + sc);

      for(int ssc=0;ssc<12;ssc++) FFTfieldPolicy::actionInputHighMode(from, ssc + 12*hit); //free for all sc
      
    }
    if(!UniqueID()){ printf("Finishing W inverse FFT\n"); fflush(stdout); }
    print_time("A2AvectorWfftw::fftinverse","FFT",fft_time);
    print_time("A2AvectorWfftw::fftinverse","Postop",postop_time);
  }
};
  

template< typename mf_Policies>
void A2AvectorWfftw<mf_Policies>::inversefft(A2AvectorW<mf_Policies> &to, fieldOperation<FermionFieldType>* mode_postop) const{
  _W_invfft_impl<A2AvectorW<mf_Policies>, const A2AvectorWfftw<mf_Policies>, WFFTfieldPolicyBasic>::inversefft(to,*this,mode_postop);
}

template< typename mf_Policies>
template<typename P>
void A2AvectorWfftw<mf_Policies>::destructiveInversefft(A2AvectorW<mf_Policies> &to, fieldOperation<FermionFieldType>* mode_postop, WFFTW_ENABLE_IF_MANUAL_ALLOC(P)){
  _W_invfft_impl<A2AvectorW<mf_Policies>, A2AvectorWfftw<mf_Policies>, WFFTfieldPolicyAllocFree>::inversefft(to,*this,mode_postop);
}

//Generate the wh field. We store in a compact notation that knows nothing about any dilution we apply when generating V from this
//For reproducibility we want to generate the wh field in the same order that Daiqian did originally. Here nhit random numbers are generated for each site/flavor
template<typename ComplexFieldType, typename complex_class>
struct _set_wh_random_impl{};

template<typename ComplexFieldType>
struct _set_wh_random_impl<ComplexFieldType, complex_double_or_float_mark>{
  static void doit(std::vector<PtrWrapper<ComplexFieldType> > &wh, const RandomType &type, const int nhits){
    typedef typename ComplexFieldType::FieldSiteType FieldSiteType;
    LRG.SetInterval(1, 0);
    int sites = wh[0]->nsites(), flavors = wh[0]->nflavors();
    
    for(int i = 0; i < sites*flavors; ++i) {
      int flav = i / sites;
      int st = i % sites;
      
      LRG.AssignGenerator(st,flav);
      for(int j = 0; j < nhits; ++j) {
	FieldSiteType* p = wh[j]->site_ptr(st,flav);
	RandomComplex<FieldSiteType>::rand(p,type,FOUR_D);
      }
    }
  }
};


template< typename mf_Policies>
void A2AvectorW<mf_Policies>::setWhRandom(const RandomType &type){
  _set_wh_random_impl<typename mf_Policies::ComplexFieldType, typename ComplexClassify<typename mf_Policies::ComplexFieldType::FieldSiteType>::type>::doit(wh,type,nhits);
}

//Get the diluted source with index id.
//We use the same set of random numbers for each spin and dilution as we do not need to rely on stochastic cancellation to separate them
//For legacy reasons we use different random numbers for the two G-parity flavors, although this is not strictly necessary
template< typename mf_Policies>
template<typename TargetFermionFieldType>
void A2AvectorW<mf_Policies>::getDilutedSource(TargetFermionFieldType &into, const int dil_id) const{
  typedef FieldSiteType mf_Complex;
  typedef typename TargetFermionFieldType::FieldSiteType TargetComplex;
  const char* fname = "getDilutedSource(...)";
  int hit, tblock, spin_color, flavor;
  StandardIndexDilution stdidx(getArgs());  
  stdidx.indexUnmap(dil_id,hit,tblock,spin_color,flavor);
    
  VRB.Result("A2AvectorW", fname, "Generating random wall source %d = (%d, %d, %d, %d).\n    ", dil_id, hit, tblock, flavor, spin_color);
  int tblock_origt = tblock * args.src_width;

  into.zero();

  if(tblock_origt / GJP.TnodeSites() != GJP.TnodeCoor()){
    VRB.Result("A2AvectorW", fname, "Not on node\n    ");
    return;
  }

  int tblock_origt_lcl = tblock_origt % GJP.TnodeSites();
    
  int src_size = GJP.VolNodeSites()/GJP.TnodeSites() * args.src_width; //size of source in units of complex numbers
#pragma omp parallel for
  for(int i=0;i<src_size;i++){
    int x[4];
    int rem = i;
    x[0] = rem % GJP.XnodeSites(); rem /= GJP.XnodeSites();
    x[1] = rem % GJP.YnodeSites(); rem /= GJP.YnodeSites();
    x[2] = rem % GJP.ZnodeSites(); rem /= GJP.ZnodeSites();
    x[3] = tblock_origt_lcl + rem;

    TargetComplex *into_site = (TargetComplex*)(into.site_ptr(x,flavor) + spin_color);
    mf_Complex const* from_site = (mf_Complex*)wh[hit]->site_ptr(x,flavor); //note same random numbers for each spin/color!
    *into_site = *from_site;
  }
}

//When gauge fixing prior to taking the FFT it is necessary to uncompact the wh field in the spin-color index, as these indices are acted upon by the gauge fixing
//(I suppose technically only the color indices need uncompacting; this might be considered as a future improvement)
template< typename mf_Policies>
void A2AvectorW<mf_Policies>::getSpinColorDilutedSource(FermionFieldType &into, const int hit, const int sc_id) const{
  const char* fname = "getSpinColorDilutedSource(...)";
  
  into.zero();

#pragma omp parallel for
  for(int i=0;i<wh[hit]->nfsites();i++){ //same mapping, different site_size
    FieldSiteType &into_site = *(into.fsite_ptr(i) + sc_id);
    const FieldSiteType &from_site = *(wh[hit]->fsite_ptr(i));
    into_site = from_site;
  }
}


template<typename mf_Policies, typename my_enable_if<_equal<typename ComplexClassify<typename mf_Policies::ComplexType>::type, complex_double_or_float_mark>::value, int>::type = 0>
void randomizeVW(A2AvectorV<mf_Policies> &V, A2AvectorW<mf_Policies> &W){
  typedef typename mf_Policies::FermionFieldType FermionFieldType;
  typedef typename mf_Policies::ComplexFieldType ComplexFieldType;
  
  int nl = V.getNl();
  int nh = V.getNh(); //number of fully diluted high-mode indices
  int nhit = V.getNhits();
  assert(nl == W.getNl());
  assert(nh == W.getNh());
  assert(nhit == W.getNhits());
  

  std::vector<FermionFieldType> wl(nl);
  for(int i=0;i<nl;i++) wl[i].setUniformRandom();
  
  std::vector<FermionFieldType> vl(nl);
  for(int i=0;i<nl;i++) vl[i].setUniformRandom();
  
  std::vector<ComplexFieldType> wh(nhit);
  for(int i=0;i<nhit;i++) wh[i].setUniformRandom();
  
  std::vector<FermionFieldType> vh(nh);
  for(int i=0;i<nh;i++) vh[i].setUniformRandom();
    
  for(int i=0;i<nl;i++){
    V.importVl(vl[i],i);
    W.importWl(wl[i],i);
  }

  for(int i=0;i<nh;i++)
    V.importVh(vh[i],i);
  
  for(int i=0;i<nhit;i++)
    W.importWh(wh[i],i);
}

//Ensure this generates randoms in the same order as the scalar version
template<typename mf_Policies, typename my_enable_if<_equal<typename ComplexClassify<typename mf_Policies::ComplexType>::type, grid_vector_complex_mark>::value, int>::type = 0>
void randomizeVW(A2AvectorV<mf_Policies> &V, A2AvectorW<mf_Policies> &W){
  typedef typename mf_Policies::FermionFieldType::FieldDimensionPolicy::EquivalentScalarPolicy ScalarDimensionPolicy;
  
  typedef CPSfermion4D<typename mf_Policies::ScalarComplexType, ScalarDimensionPolicy, DynamicFlavorPolicy, StandardAllocPolicy> ScalarFermionFieldType;
  typedef CPScomplex4D<typename mf_Policies::ScalarComplexType, ScalarDimensionPolicy, DynamicFlavorPolicy, StandardAllocPolicy> ScalarComplexFieldType;
  
  int nl = V.getNl();
  int nh = V.getNh(); //number of fully diluted high-mode indices
  int nhit = V.getNhits();
  assert(nl == W.getNl());
  assert(nh == W.getNh());
  assert(nhit == W.getNhits());

  ScalarFermionFieldType tmp;
  ScalarComplexFieldType tmp_cmplx;
  
  for(int i=0;i<nl;i++){
    tmp.setUniformRandom();
    W.getWl(i).importField(tmp);
  }
  for(int i=0;i<nl;i++){
    tmp.setUniformRandom();
    V.getVl(i).importField(tmp);
  }
  for(int i=0;i<nhit;i++){
    tmp_cmplx.setUniformRandom();
    W.getWh(i).importField(tmp_cmplx);
  }
  for(int i=0;i<nh;i++){
    tmp.setUniformRandom();
    V.getVh(i).importField(tmp);
  }
}

template< typename FieldType>
FieldType const * getBaseAndShift(int shift[3], const int p[3], FieldType const *base_p, FieldType const *base_m){
  //With G-parity base_p has momentum +1 in each G-parity direction, base_m has momentum -1 in each G-parity direction.
  //Non-Gparity directions are assumed to have momentum 0

  //Units of momentum are 2pi/L for periodic BCs, pi/L for antiperiodic and pi/2L for Gparity
  FieldType const * out = GJP.Gparity() ? NULL : base_p;
  for(int d=0;d<3;d++){
    if(GJP.Bc(d) == BND_CND_GPARITY){
      //Type 1 : f_{p=4b+1}(n) = f_+1(n+b)     // p \in {.. -7 , -3, 1, 5, 9 ..}
      //Type 2 : f_{p=4b-1}(n) = f_-1(n+b)     // p \n  {.. -5, -1, 3, 7 , 11 ..}
      if( (p[d]-1) % 4 == 0 ){
	//Type 1
	int b = (p[d]-1)/4;
	shift[d] = -b;  //shift f_+1 backwards by b
	if(out == NULL) out = base_p;
	else if(out != base_p) ERR.General("","getBaseAndShift","Momentum (%d,%d,%d) appears to be invalid because momenta in different G-parity directions do not reside in the same set\n",p[0],p[1],p[2]);
	
      }else if( (p[d]+1) % 4 == 0 ){
	//Type 2
	int b = (p[d]+1)/4;
	shift[d] = -b;  //shift f_-1 backwards by b
	if(out == NULL) out = base_m;
	else if(out != base_m) ERR.General("","getBaseAndShift","Momentum (%d,%d,%d) appears to be invalid because momenta in different G-parity directions do not reside in the same set\n",p[0],p[1],p[2]);
	
      }else ERR.General("","getBaseAndShift","Momentum (%d,%d,%d) appears to be invalid because one or more components in G-parity directions are not allowed\n",p[0],p[1],p[2]);
    }else{
      //f_b(n) = f_0(n+b)
      //Let the other directions decide on which base to use if some of them are G-parity dirs ; otherwise the pointer defaults to base_p above
      shift[d] = -p[d];
    }
  }
  if(!UniqueID()) printf("getBaseAndShift for p=(%d,%d,%d) determined shift=(%d,%d,%d) from ptr %c\n",p[0],p[1],p[2],shift[0],shift[1],shift[2],out == base_p ? 'p' : 'm');
  assert(out != NULL);
  
  return out;
}




//Use the relations between FFTs to obtain the FFT for a chosen quark momentum
//With G-parity BCs there are 2 disjoint sets of momenta hence there are 2 base FFTs
template< typename mf_Policies>
void A2AvectorWfftw<mf_Policies>::getTwistedFFT(const int p[3], A2AvectorWfftw<Policies> const *base_p, A2AvectorWfftw<Policies> const *base_m){
  Float time = -dclock();
  
  std::vector<int> shift(3);
  A2AvectorWfftw<mf_Policies> const* base = getBaseAndShift(&shift[0], p, base_p, base_m);
  if(base == NULL) ERR.General("A2AvectorWfftw","getTwistedFFT","Base pointer for twist momentum (%d,%d,%d) is NULL\n",p[0],p[1],p[2]);

  wl = base->wl;
  wh = base->wh;
  
  int nshift = 0;
  for(int i=0;i<3;i++) if(shift[i]) nshift++;

  if(nshift > 0){
    for(int i=0;i<this->getNmodes();i++)
      shiftPeriodicField( this->getMode(i), base->getMode(i), shift);
  }
  time += dclock();
  print_time("A2AvectorWfftw::getTwistedFFT","Twist",time);
}


template< typename mf_Policies>
void A2AvectorWfftw<mf_Policies>::shiftFieldsInPlace(const std::vector<int> &shift){
  Float time = -dclock();
  int nshift = 0;
  for(int i=0;i<3;i++) if(shift[i]) nshift++;
  if(nshift > 0){
    for(int i=0;i<this->getNmodes();i++)
      shiftPeriodicField( this->getMode(i), this->getMode(i), shift);
  }
  print_time("A2AvectorWfftw::shiftFieldsInPlace","Total",time + dclock());
}

//A version of the above that directly shifts the base Wfftw rather than outputting into a separate storage
//Returns the pointer to the Wfftw acted upon and the *shift required to restore the Wfftw to it's original form*
template< typename mf_Policies>
std::pair< A2AvectorWfftw<mf_Policies>*, std::vector<int> > A2AvectorWfftw<mf_Policies>::inPlaceTwistedFFT(const int p[3], A2AvectorWfftw<mf_Policies> *base_p, A2AvectorWfftw<mf_Policies> *base_m){
  Float time = -dclock();
  
  std::vector<int> shift(3);
  A2AvectorWfftw<mf_Policies>* base = getBaseAndShift(&shift[0], p, base_p, base_m);
  if(base == NULL) ERR.General("A2AvectorWfftw","getTwistedFFT","Base pointer for twist momentum (%d,%d,%d) is NULL\n",p[0],p[1],p[2]);

  base->shiftFieldsInPlace(shift);

  for(int i=0;i<3;i++) shift[i] = -shift[i];
  
  time += dclock();
  print_time("A2AvectorWfftw::inPlaceTwistedFFT","Twist",time);

  return std::pair< A2AvectorWfftw<mf_Policies>*, std::vector<int> >(base,shift);
}
  



template< typename mf_Policies>
void A2AvectorVfftw<mf_Policies>::getTwistedFFT(const int p[3], A2AvectorVfftw<Policies> const *base_p, A2AvectorVfftw<Policies> const *base_m){
  Float time = -dclock();
  
  std::vector<int> shift(3);
  A2AvectorVfftw<mf_Policies> const* base = getBaseAndShift(&shift[0], p, base_p, base_m);
  if(base == NULL) ERR.General("A2AvectorVfftw","getTwistedFFT","Base pointer for twist momentum (%d,%d,%d) is NULL\n",p[0],p[1],p[2]);
  
  v = base->v;
  
  int nshift = 0;
  for(int i=0;i<3;i++) if(shift[i]) nshift++;

  if(nshift > 0){
    for(int i=0;i<this->getNmodes();i++)
      shiftPeriodicField( this->getMode(i), base->getMode(i), shift);
  }
  time += dclock();
  print_time("A2AvectorVfftw::getTwistedFFT","Twist",time);
}

template< typename mf_Policies>
void A2AvectorVfftw<mf_Policies>::shiftFieldsInPlace(const std::vector<int> &shift){
  Float time = -dclock();
  int nshift = 0;
  for(int i=0;i<3;i++) if(shift[i]) nshift++;
  if(nshift > 0){
    for(int i=0;i<this->getNmodes();i++)
      shiftPeriodicField( this->getMode(i), this->getMode(i), shift);
  }
  print_time("A2AvectorVfftw::shiftFieldsInPlace","Total",time + dclock());
}

//A version of the above that directly shifts the base Wfftw rather than outputting into a separate storage
//Returns the pointer to the Wfftw acted upon and the *shift required to restore the Wfftw to it's original form*
template< typename mf_Policies>
std::pair< A2AvectorVfftw<mf_Policies>*, std::vector<int> > A2AvectorVfftw<mf_Policies>::inPlaceTwistedFFT(const int p[3], A2AvectorVfftw<mf_Policies> *base_p, A2AvectorVfftw<mf_Policies> *base_m){
  Float time = -dclock();
  
  std::vector<int> shift(3);
  A2AvectorVfftw<mf_Policies>* base = getBaseAndShift(&shift[0], p, base_p, base_m);
  if(base == NULL) ERR.General("A2AvectorWfftw","getTwistedFFT","Base pointer for twist momentum (%d,%d,%d) is NULL\n",p[0],p[1],p[2]);

  base->shiftFieldsInPlace(shift);

  for(int i=0;i<3;i++) shift[i] = -shift[i];
  
  time += dclock();
  print_time("A2AvectorVfftw::inPlaceTwistedFFT","Twist",time);

  return std::pair< A2AvectorVfftw<mf_Policies>*, std::vector<int> >(base,shift);
}
