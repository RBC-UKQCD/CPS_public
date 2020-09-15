#ifndef _A2A_SOURCES_H
#define _A2A_SOURCES_H

CPS_START_NAMESPACE

//Spatial source structure in *momentum-space*. Should assign the same value to both flavors if G-parity

//3D complex field. Defined for a *single flavor* if GPBC
template<typename mf_Complex,typename DimensionPolicy = SpatialPolicy, typename FieldAllocPolicy = StandardAllocPolicy, typename my_enable_if<DimensionPolicy::EuclideanDimension == 3, int>::type = 0>
class A2Asource{
public:
  typedef CPSfield<mf_Complex,1,DimensionPolicy,OneFlavorPolicy,FieldAllocPolicy> FieldType;  
protected:
  FieldType *src;
public:
  A2Asource(const typename FieldType::InputParamType &params){
    setup(params);
  }
  inline void setup(const typename FieldType::InputParamType &params){
    src = new FieldType(params);
  }
  
  A2Asource(): src(NULL){}
  A2Asource(const A2Asource &cp): src(cp.src == NULL ? NULL : new FieldType(*cp.src)){}
  
  ~A2Asource(){
    if(src != NULL) delete src;
  }

  
  inline const mf_Complex & siteComplex(const int site) const{ return *src->site_ptr(site); }
  inline const int nsites() const{ return src->nsites(); }

  template< typename extComplexType, typename extDimPol, typename extAllocPol>
  void importSource(const A2Asource<extComplexType,extDimPol,extAllocPol> &from){
    src->importField(*from.src);
  }
  FieldType & getSource(){ return *src; } //For testing

  //Periodic modulus operation
  inline static int pmod(const int x, const int Lx){
    return x <= Lx/2 ? x : Lx-x; //   0 ... L/2-1,  L/2, L/2-1, ... 1
  }
  //Periodic coordinate relative to boundary
  inline static int pcoord(const int x, const int Lx){
    return (x + Lx/2) % Lx - Lx/2; //   0 ... L/2-1, -L/2, 1-L/2, ... -1   includes sign. Note convention for sign at L/2
  }

  //Radial coordinate
  static Float pmodr(const int x[3], const int L[3]){
    Float ssq =0.;
    for(int i=0;i<3;i++){
      int sr = pmod(x[i],L[i]);
      ssq += sr*sr;
    }
    return sqrt(ssq);
  }
  
  //Spherical coordinates that know about the periodicity of the lattice  
  static void pmodspherical(Float &r, Float &theta, Float &phi, const int x[3], const int L[3]){
    int xp[3];
    Float ssq = 0.;    
    for(int i=0;i<3;i++){
      xp[i] = pcoord(x[i],L[i]);
      ssq += xp[i]*xp[i];
    }
    r = sqrt(ssq);
    theta = acos(xp[2]/r);
    phi = atan(xp[1]/xp[0]);
  }

  
};


//Use CRTP for 'setSite' method which should be specialized according to the source type
template<typename FieldPolicies, typename Child>
class A2AsourceBase: public A2Asource<typename FieldPolicies::ComplexType, typename FieldPolicies::DimensionPolicy, typename FieldPolicies::AllocPolicy>{
public:
  typedef FieldPolicies Policies;
  typedef typename A2Asource<typename FieldPolicies::ComplexType, typename FieldPolicies::DimensionPolicy, typename FieldPolicies::AllocPolicy>::FieldType::InputParamType FieldParamType;
  
  A2AsourceBase(const FieldParamType &p): A2Asource<typename FieldPolicies::ComplexType, typename FieldPolicies::DimensionPolicy, typename FieldPolicies::AllocPolicy>(p){};
  A2AsourceBase(): A2Asource<typename FieldPolicies::ComplexType, typename FieldPolicies::DimensionPolicy, typename FieldPolicies::AllocPolicy>(){}; //SOURCE IS NOT SETUP
  A2AsourceBase(const A2AsourceBase &r): A2Asource<typename FieldPolicies::ComplexType, typename FieldPolicies::DimensionPolicy, typename FieldPolicies::AllocPolicy>(r){}
  
  void fft_source(){
    assert(this->src != NULL);
    int glb_size[3]; for(int i=0;i<3;i++) glb_size[i] = GJP.Nodes(i)*GJP.NodeSites(i);

    //Generate a global 4d source
    CPSglobalComplexSpatial<cps::ComplexD,OneFlavorPolicy> glb; //always of this type
    glb.zero();
         
#pragma omp_parallel for
    for(int i=0;i<glb.nsites();i++){
      int x[3]; glb.siteUnmap(i,x); 
      *glb.site_ptr(i) = static_cast<Child const*>(this)->value(x,glb_size);
    }
    //Perform the FFT and pull out this nodes subvolume
    glb.fft();
    glb.scatter<typename FieldPolicies::ComplexType, typename FieldPolicies::DimensionPolicy, typename FieldPolicies::AllocPolicy>(*this->src);
  }
};


template<typename FieldPolicies, typename Derived>
class A2AhydrogenSourceBase: public A2AsourceBase<FieldPolicies, Derived >{
protected:
  Float radius;

public:
  typedef FieldPolicies Policies;
  typedef typename A2AsourceBase<FieldPolicies, Derived >::FieldParamType FieldParamType;
  typedef typename Policies::ComplexType ComplexType;
      
  A2AhydrogenSourceBase(const Float _radius, const FieldParamType &field_params): radius(_radius), A2AsourceBase<FieldPolicies, Derived >(field_params){
    this->fft_source();
  }

  A2AhydrogenSourceBase(): radius(0.), A2AsourceBase<FieldPolicies, Derived >(){} //src is not setup

  A2AhydrogenSourceBase(const A2AhydrogenSourceBase &r): radius(r.radius), A2AsourceBase<FieldPolicies, Derived >(r){}
  
  //Setup the source if the default constructor was used
  void setup(const Float _radius, const FieldParamType &field_params){
    this->A2AsourceBase<FieldPolicies, Derived >::setup(field_params);
    radius = _radius;
    this->fft_source();
  }
  void setup(const Float radius){
    return setup(radius, NullObject());
  }
    
  inline void siteFmat(FlavorMatrixGeneral<typename Policies::ComplexType> &out, const int site) const{
    out(0,0) = out(1,1) = this->siteComplex(site);
    out(0,1) = out(1,0) = typename Policies::ComplexType(0);    
  }
};

  
//Exponential (hydrogen wavefunction) source
//SrcParams is just a Float for the radius
template<typename FieldPolicies = StandardSourcePolicies>
class A2AexpSource: public A2AhydrogenSourceBase<FieldPolicies, A2AexpSource<FieldPolicies> >{
public:
  typedef FieldPolicies Policies;
  typedef typename A2AhydrogenSourceBase<FieldPolicies, A2AexpSource<FieldPolicies> >::FieldParamType FieldParamType;
  typedef typename Policies::ComplexType ComplexType;

  inline cps::ComplexD value(const int site[3], const int glb_size[3]) const{
    Float v = this->pmodr(site,glb_size)/this->radius;
    v = exp(-v)/(glb_size[0]*glb_size[1]*glb_size[2]);
    return ComplexD(v,0);
  }
    
  A2AexpSource(const Float _radius, const FieldParamType &field_params): A2AhydrogenSourceBase<FieldPolicies, A2AexpSource<FieldPolicies> >(_radius,field_params){ }
  A2AexpSource(const Float _radius): A2AhydrogenSourceBase<FieldPolicies, A2AexpSource<FieldPolicies> >(_radius, NullObject()){ }

  A2AexpSource(): A2AhydrogenSourceBase<FieldPolicies, A2AexpSource<FieldPolicies> >(){} //src is not setup
  A2AexpSource(const A2AexpSource &r): A2AhydrogenSourceBase<FieldPolicies, A2AexpSource<FieldPolicies> >(r){}
};

//General s-wave hydrogen wavefunction source 
template<typename FieldPolicies = StandardSourcePolicies>
class A2AhydrogenSource: public A2AhydrogenSourceBase<FieldPolicies, A2AhydrogenSource<FieldPolicies> >{
  int n, l, m;
public:
  typedef FieldPolicies Policies;
  typedef typename A2AhydrogenSourceBase<FieldPolicies, A2AhydrogenSource<FieldPolicies> >::FieldParamType FieldParamType;
  typedef typename Policies::ComplexType ComplexType;

  inline static cps::ComplexD zexp(const Float phi){
    return cps::ComplexD( cos(phi), sin(phi) );
  }
    
  
  inline cps::ComplexD value(const int site[3], const int glb_size[3]) const{
    assert(n>=0 && n <= 3 &&
	   l>=0 && l <= n-1 &&
	   abs(m) <= l);
    
    Float a0 = this->radius;
    Float na0 = a0 * n;
    Float r, theta, phi;
    if(l==0) r = this->pmodr(site,glb_size);   //don't need theta and phi for s-wave
    else this->pmodspherical(r,theta,phi,site,glb_size);
    
    cps::ComplexD v( exp(-r/na0)/(glb_size[0]*glb_size[1]*glb_size[2]), 0.);

    switch(n){
    case 1:
      break;
    case 2:
      switch(l){
      case 0:
	v *= (2. - r/a0); break;
      case 1:
	v *= r/a0;
	switch(m){
	case 0:
	  v *= cos(theta); break;
	case 1:
	case -1:
	  v *= sin(theta) * zexp(phi*m); break;
	}//m
      }//l
      break;
    case 3:
      switch(l){
      case 0:
	v *= (27. - 18.*r/a0 + 2.*r*r/a0/a0); break;
      case 1:
	v *= (6. -r/a0)*r/a0;
	switch(m){
	case 0:
	  v *= cos(theta); break;
	case 1:
	case -1:
	  v *= sin(theta) * zexp(phi*m); break;
	}//m
      case 2:
	v *= r*r/a0/a0;
	switch(m){
	case 0:
	  v *= 3.*cos(theta)*cos(theta) - 1.; break;
	case 1:
	case -1:
	  v *= sin(theta) * cos(theta) * zexp(phi*m); break;
	case 2:
	case -2:
	  v *= sin(theta) * sin(theta) * zexp(phi*m); break;	  
	}//m	
      }//l
      break;     
    }//n
    return v;      
  }

  A2AhydrogenSource(const int _n, const int _l, const int _m, const Float _radius, const FieldParamType &field_params): n(_n), l(_l), m(_m), A2AhydrogenSourceBase<FieldPolicies, A2AhydrogenSource<FieldPolicies> >(_radius,field_params){ }
  A2AhydrogenSource(const int _n, const int _l, const int _m, const Float _radius): n(_n), l(_l), m(_m), A2AhydrogenSourceBase<FieldPolicies, A2AhydrogenSource<FieldPolicies> >(_radius, NullObject()){ }

  A2AhydrogenSource(): A2AhydrogenSourceBase<FieldPolicies, A2AhydrogenSource<FieldPolicies> >(){} //src is not setup

  A2AhydrogenSource(const A2AhydrogenSource &r): n(r.n), l(r.l), m(r.m), A2AhydrogenSourceBase<FieldPolicies, A2AhydrogenSource<FieldPolicies> >(r){}
  
  void setup(const int _n, const int _l, const int _m, const Float _radius, const FieldParamType &field_params){
    n = _n; l=_l; m=_m;
    this->A2AhydrogenSourceBase<FieldPolicies, A2AhydrogenSource<FieldPolicies> >::setup(_radius,field_params);
  }
  void setup(const int _n, const int _l, const int _m, const Float radius){
    return setup(_n,_l,_m, radius, NullObject());
  }
};





//Box source. Unflavored so ignore second flav
//SrcParams is std::vector<Float> for the extents x,y,z . *These must be even numbers* (checked)
template<typename FieldPolicies = StandardSourcePolicies>
class A2AboxSource: public A2AsourceBase<FieldPolicies, A2AboxSource<FieldPolicies> >{
  int box_size[3];

  void box_setup_fft(const int _box_size[3]){    
    for(int i=0;i<3;i++){
      if(_box_size[i] % 2 == 1){
	ERR.General("A2AboxSource","A2AboxSource","box size must be multiple of 2");
      }
      box_size[i] = _box_size[i];
    }
    this->fft_source();
  }
  
public:
  typedef FieldPolicies Policies;
  typedef typename A2AsourceBase<FieldPolicies, A2AboxSource<FieldPolicies> >::FieldParamType FieldParamType;
  typedef typename Policies::ComplexType ComplexType;
  
  cps::ComplexD value(const int site[3], const int glb_size[3]) const{
    bool inbox = true;
    int V = glb_size[0]*glb_size[1]*glb_size[2];
    for(int i=0;i<3;i++){ 
      int bdist = this->pmod(site[i],glb_size[i]);
      
      if(bdist > box_size[i]){
	inbox = false; break;
      }
    }
    if(inbox)
      return cps::ComplexD(1./V);
  }
  
  A2AboxSource(const int _box_size[3],const FieldParamType &field_params): A2AsourceBase<FieldPolicies, A2AboxSource<FieldPolicies> >(field_params){
    this->box_setup_fft(_box_size);
  }
  A2AboxSource(const int _box_size[3]): A2AsourceBase<FieldPolicies, A2AboxSource<FieldPolicies> >(NullObject()){
    this->box_setup_fft(_box_size);
  }//syntatic sugar to avoid creating a NullObject
  A2AboxSource(const A2AboxSource &r):  A2AsourceBase<FieldPolicies, A2AboxSource<FieldPolicies> >(r){ memcpy(box_size,r.box_size,3*sizeof(int)); }
  
  
  void setup(const int _box_size[3], const FieldParamType &field_params = NullObject()){
    this->A2AsourceBase<FieldPolicies, A2AboxSource<FieldPolicies> >::setup(field_params);
    this->box_setup_fft(_box_size);
  }

  
  inline void siteFmat(FlavorMatrixGeneral<typename Policies::ComplexType> &out, const int site) const{
    out(0,0) = out(1,1) = this->siteComplex(site);
    out(0,1) = out(1,0) = typename Policies::ComplexType(0);    
  }
};

//Splat a cps::ComplexD onto a SIMD type. Just a plain copy for non-SIMD complex types
#ifdef USE_GRID
template<typename ComplexType>
inline void SIMDsplat(ComplexType &to, const cps::ComplexD &from, typename my_enable_if< _equal<  typename ComplexClassify<ComplexType>::type, grid_vector_complex_mark  >::value, int>::type = 0){
  vsplat(to,from);
}
#endif
template<typename ComplexType>
inline void SIMDsplat(ComplexType &to, const cps::ComplexD &from, typename my_enable_if< !_equal<  typename ComplexClassify<ComplexType>::type, grid_vector_complex_mark  >::value, int>::type = 0){
  to = ComplexType(from.real(),from.imag());
}
  

//Daiqian's original implementation sets the (1 +/- sigma_2) flavor projection on G-parity fields to unity when the two fermion fields coincide.
//I'm not sure this is actually necessary, but I need to be able to reproduce his numbers
//Derived class should setup the sources
template<typename SourceType>
class A2AflavorProjectedSource: public SourceType{
public:
  typedef typename SourceType::FieldParamType FieldParamType;
  typedef typename SourceType::Policies::ComplexType ComplexType;
protected:
  int sign;
  ComplexType *val000;
  virtual void dummy() = 0; //make sure this class can't be instantiated directly

  void setup_projected_src_info(const int p[3]){
    sign = getProjSign(p);
    int zero[3] = {0,0,0}; int L[3] = {GJP.NodeSites(0)*GJP.Nodes(0), GJP.NodeSites(1)*GJP.Nodes(1), GJP.NodeSites(2)*GJP.Nodes(2) };
    cps::ComplexD v = this->value(zero,L);
    SIMDsplat(*val000,v);    
  }
public:

  A2AflavorProjectedSource(): val000( (ComplexType*)memalign(128,sizeof(ComplexType)) ), SourceType(){}
  ~A2AflavorProjectedSource(){ free(val000); }
  
  A2AflavorProjectedSource(const A2AflavorProjectedSource &r): val000( (ComplexType*)memalign(128,sizeof(ComplexType)) ), sign(r.sign), SourceType(r){
    *val000 = *r.val000;
  }
  
  //Assumes momenta are in units of \pi/2L, and must be *odd integer* (checked)
  inline static int getProjSign(const int p[3]){
    if(!GJP.Gparity()){ ERR.General("A2AflavorProjectedSource","getProjSign","Requires GPBC in at least one direction\n"); }

    //Sign is exp(i\pi n_p)
    //where n_p is the solution to  p_j = \pi/2L( 1 + 2n_p )
    //Must be consistent for all j

    int np;
    for(int j=0;j<3;j++){
      if(GJP.Bc(j)!=BND_CND_GPARITY) continue;

      if(abs(p[j]) %2 != 1){ ERR.General("A2AflavorProjectedSource","getProjSign","Component %d of G-parity momentum (%d,%d,%d) is invalid as it is not an odd integer!\n",j,p[0],p[1],p[2]); }
      int npj = (p[j] - 1)/2;
      if(j == 0) np = npj;
      else if(abs(npj)%2 != abs(np)%2){ 
	ERR.General("A2AflavorProjectedSource","getProjSign","Momentum component %d of G-parity momentum (%d,%d,%d) is invalid because it doesn't differ from component 0 by multiple of 2pi (4 in these units). Got np(0)=%d, np(j)=%d\n",j,p[0],p[1],p[2],np,npj); 
      }
    }
    int sgn = (abs(np) % 2 == 0 ? 1 : -1); //exp(i\pi n_p) = exp(-i\pi n_p)  for all integer n_p
    if(!UniqueID()){ printf("A2AflavorProjectedSource::getProjSign got sign %d (np = %d) for p=(%d,%d,%d)pi/2L\n",sgn,np,p[0],p[1],p[2]); fflush(stdout); }

    return sgn;
  }

  inline void siteFmat(FlavorMatrixGeneral<ComplexType> &out, const int site) const{
    //Matrix is FFT of  (1 + [sign]*sigma_2) when |x-y| !=0 or 1 when |x-y| == 0
    //It is always 1 on the diagonals
    const ComplexType &val = this->siteComplex(site);
    
    out(0,0) = out(1,1) = val;
    //and has \pm i on the diagonals with a momentum structure that is computed by omitting site 0,0,0
    out(1,0) = multiplySignTimesI(sign,val - *val000);
    out(0,1) = -out(1,0); //-1 from sigma2
  }

  //Can change momentum sign without redoing FFT
  void setMomentum(const int p[3]){
    sign = getProjSign(p);
  }
};


template<typename FieldPolicies = StandardSourcePolicies>
class A2AflavorProjectedExpSource : public A2AflavorProjectedSource<A2AexpSource<FieldPolicies> >{
  void dummy(){}
public:
  typedef typename A2AflavorProjectedSource<A2AexpSource<FieldPolicies> >::FieldParamType FieldParamType;
  typedef typename A2AflavorProjectedSource<A2AexpSource<FieldPolicies> >::ComplexType ComplexType;
  
  A2AflavorProjectedExpSource(const Float radius, const int p[3], const FieldParamType &src_field_params = NullObject()){
    this->A2AexpSource<FieldPolicies>::setup(radius,src_field_params);
    this->A2AflavorProjectedSource<A2AexpSource<FieldPolicies> >::setup_projected_src_info(p);
  }
  A2AflavorProjectedExpSource(): A2AflavorProjectedSource<A2AexpSource<FieldPolicies> >(){}

  A2AflavorProjectedExpSource(const A2AflavorProjectedExpSource &r): A2AflavorProjectedSource<A2AexpSource<FieldPolicies> >(r){}
  
  void setup(const Float radius, const int p[3], const FieldParamType &src_field_params = NullObject()){
    this->A2AexpSource<FieldPolicies>::setup(radius,src_field_params);
    this->A2AflavorProjectedSource<A2AexpSource<FieldPolicies> >::setup_projected_src_info(p);
  }
  
};

template<typename FieldPolicies = StandardSourcePolicies>
class A2AflavorProjectedHydrogenSource : public A2AflavorProjectedSource<A2AhydrogenSource<FieldPolicies> >{
  void dummy(){}
public:
  typedef typename A2AflavorProjectedSource<A2AhydrogenSource<FieldPolicies> >::FieldParamType FieldParamType;
  typedef typename A2AflavorProjectedSource<A2AhydrogenSource<FieldPolicies> >::ComplexType ComplexType;
  
  A2AflavorProjectedHydrogenSource(const int _n, const int _l, const int _m, const Float _radius, const int p[3], const FieldParamType &src_field_params = NullObject()){
    this->A2AhydrogenSource<FieldPolicies>::setup(_n,_l,_m,_radius,src_field_params);
    this->A2AflavorProjectedSource<A2AhydrogenSource<FieldPolicies> >::setup_projected_src_info(p);
  }
  A2AflavorProjectedHydrogenSource(): A2AflavorProjectedSource<A2AhydrogenSource<FieldPolicies> >(){}

  A2AflavorProjectedHydrogenSource(const A2AflavorProjectedHydrogenSource &r): A2AflavorProjectedSource<A2AhydrogenSource<FieldPolicies> >(r){}
  
  void setup(const int _n, const int _l, const int _m, const Float _radius, const int p[3], const FieldParamType &src_field_params = NullObject()){
    this->A2AhydrogenSource<FieldPolicies>::setup(_n,_l,_m,_radius,src_field_params);
    this->A2AflavorProjectedSource<A2AhydrogenSource<FieldPolicies> >::setup_projected_src_info(p);
  }
  
};



define_test_has_enum(nSources); //a test for multisrc types (all should have enum nSources)

template<typename SourceList>
class A2AmultiSource{
public:
  typedef ListStruct<SourceList> SourceListStruct;
private:
  SourceListStruct sources;
public:
  enum { nSources = getSizeOfListStruct<SourceListStruct>::value };

  //Accessors for sources  (call like  src.template get<Idx>() )
  template<int i>
  typename getTypeFromList<SourceListStruct,i>::type & getSource(){ return getElemFromListStruct<SourceListStruct,i>::get(sources); }
  template<int i>
  const typename getTypeFromList<SourceListStruct,i>::type & getSource() const{ return getConstElemFromListStruct<SourceListStruct,i>::get(sources); }
};



CPS_END_NAMESPACE

#endif
