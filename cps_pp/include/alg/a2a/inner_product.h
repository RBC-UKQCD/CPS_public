#ifndef _CK_INNER_PRODUCT_H
#define _CK_INNER_PRODUCT_H

#include<alg/a2a/scfvectorptr.h>
#include<alg/a2a/a2a_sources.h>
#include<alg/a2a/gsl_wrapper.h>
#include<alg/a2a/conj_zmul.h>
#include<alg/a2a/inner_product_spincolorcontract.h>
#include<alg/a2a/inner_product_fmatspincolorcontract.h>

CPS_START_NAMESPACE

//Classes that perform the inner product of two spin-color-flavor vectors on a given (momentum-space) site
//Class must have an 
//complex<double> operator()(const SCFvectorPtr<mf_Complex> &l, const SCFvectorPtr<mf_Complex> &r, const int &p, const int &t) const
//p is the *local* 3-momentum coordinate in canonical ordering 
//t is the local time coordinate
//mf_Complex is the base complex type for the vectors
//Output should be *double precision complex* even if the vectors are stored in single precision. Do this to avoid finite prec errors on spatial sum

inline void doAccum(std::complex<double> &to, const std::complex<double> &from){
  to += from;
}
#ifdef USE_GRID
inline void doAccum(std::complex<double> &to, const Grid::vComplexD &from){
  to += Reduce(from);
}
inline void doAccum(Grid::vComplexD &to, const Grid::vComplexD &from){
  to += from;
}
#endif

//Simple inner product of a momentum-space scalar source function and a constant spin matrix
//Assumed diagonal matrix in flavor space if G-parity
template<typename mf_Complex, typename SourceType, bool conj_left = true, bool conj_right=false>
class SCmatrixInnerProduct{
  const WilsonMatrix &sc;
  const SourceType &src;
  bool conj[2];
public:
  typedef SourceType InnerProductSourceType;
  
  SCmatrixInnerProduct(const WilsonMatrix &_sc, const SourceType &_src): sc(_sc), src(_src){ }

  void operator()(std::complex<double> &into, const SCFvectorPtr<mf_Complex> &l, const SCFvectorPtr<mf_Complex> &r, const int p, const int t) const{
    std::complex<double> out(0.0,0.0);
    for(int f=0;f<1+GJP.Gparity();f++){
      // Mr
      std::complex<double> tvec[4][3];
      for(int s1=0;s1<4;s1++)
	for(int c1=0;c1<3;c1++){
	  tvec[s1][c1] = std::complex<double>(0.0,0.0);
	  for(int s2=0;s2<4;s2++)
	    for(int c2=0;c2<3;c2++)
	      tvec[s1][c1] += sc(s1,c1,s2,c2) * ( conj_right ? std::conj(r(s2,c2,f)) : r(s2,c2,f) );	    
	}      
      //l.(Mr)
      std::complex<double> outf(0.0,0.0);
      for(int s1=0;s1<4;s1++)
	for(int c1=0;c1<3;c1++)
	  outf += ( conj_left ? std::conj(l(s1,c1,f)) : l(s1,c1,f) ) * tvec[s1][c1];
	
      out += outf;
    }
    //Multiply by momentum-space source structure
    into += out * src.siteComplex(p);
  }
};


//Optimized gamma^5 inner product with unit flavor matrix
template<typename mf_Complex, typename SourceType, bool conj_left = true, bool conj_right=false>
class SCg5InnerProduct{
  const SourceType &src;
public:
  typedef SourceType InnerProductSourceType;
  
  SCg5InnerProduct(const SourceType &_src): src(_src){ }
    
  void operator()(std::complex<double> &into, const SCFvectorPtr<mf_Complex> &l, const SCFvectorPtr<mf_Complex> &r, const int p, const int t) const{
    std::complex<double> out(0.0,0.0);
    for(int f=0;f<1+GJP.Gparity();f++)
      out += OptimizedSpinColorContract<mf_Complex,conj_left,conj_right>::g5(l.getPtr(f),r.getPtr(f));    
    into += out * src.siteComplex(p);
  }
};

//Optimized inner product for general spin matrix
//Spin matrix indexed in range (0..15) following QDP convention: integer representation of binary(n4 n3 n2 n1) for spin structure  gamma1^n1 gamma2^n2 gamma3^n3 gamma4^n4
  //Thus   idx   matrix
  //       0      Unit
  //       1      gamma1
  //       2      gamma2
  //       3      gamma1 gamma2
  //       4      gamma3
  //       5      gamma1 gamma3
  //       6      gamma2 gamma3
  //       7      gamma1 gamma2 gamma3        =  gamma5 gamma4
  //       8      gamma4
  //       9      gamma1 gamma4
  //       10     gamma2 gamma4
  //       11     gamma1 gamma2 gamma4        = -gamma5 gamma3
  //       12     gamma3 gamma4
  //       13     gamma1 gamma3 gamma4        =  gamma5 gamma2
  //       14     gamma2 gamma3 gamma4        = -gamma5 gamma1
  //       15     gamma1 gamma2 gamma3 gamma4 =  gamma5
template<int smatidx,typename mf_Complex, typename SourceType, bool conj_left = true, bool conj_right=false>
class SCspinInnerProduct{
  const SourceType &src;
public:
  typedef SourceType InnerProductSourceType;
  
  SCspinInnerProduct(const SourceType &_src): src(_src){ }
    
  //When running with a multisrc type this returns the number of meson fields per timeslice = nSources
  template<typename S=SourceType>
  inline typename my_enable_if<has_enum_nSources<S>::value, int>::type mfPerTimeSlice() const{ return SourceType::nSources; }

  template<typename S=SourceType>
  inline typename my_enable_if<!has_enum_nSources<S>::value, int>::type mfPerTimeSlice() const{ return 1; }

  template<typename AccumType, typename ComplexType = mf_Complex>  
  void operator()(AccumType &into, const SCFvectorPtr<ComplexType> &l, const SCFvectorPtr<ComplexType> &r, const int p, const int t) const{
    assert(!GJP.Gparity());
    ComplexType out = GeneralSpinColorContractSelect<smatidx,ComplexType,conj_left,conj_right, typename ComplexClassify<ComplexType>::type>::doit(l.getPtr(0),r.getPtr(0));
    doAccum(into,out * src.siteComplex(p));
  }  
};


//Constant spin-color-flavor matrix source structure with position-dependent flavor matrix from source
// l M N r    where l,r are the vectors, M is the constant matrix and N the position-dependent
//For use with GPBC
template<typename mf_Complex, typename SourceType, bool conj_left = true, bool conj_right=false>
class SCFfmatSrcInnerProduct{
  const SourceType &src;
  const SpinColorFlavorMatrix &scf;
public:
  typedef SourceType InnerProductSourceType;
  
  SCFfmatSrcInnerProduct(const SpinColorFlavorMatrix &_scf, const SourceType &_src): 
    scf(_scf), src(_src){ 
    if(!GJP.Gparity()) ERR.General("SCFfmatSrcInnerProduct","SCFfmatSrcInnerProduct","Only for G-parity BCs");
  }

  void operator()(std::complex<double> &into, const SCFvectorPtr<mf_Complex> &l, const SCFvectorPtr<mf_Complex> &r, const int p, const int t) const{
    //Get source flavor matrix structure for this momentum site
    FlavorMatrix N = src.siteFmat(p);
    
    //Nr
    std::complex<double> rvec[4][3][2];
    for(int s1=0;s1<4;s1++)
      for(int c1=0;c1<3;c1++)
	for(int f1=0;f1<2;f1++){
	  rvec[s1][c1][f1] = std::complex<double>(0.0,0.0);
	  for(int f2=0;f2<2;f2++){
	    std::complex<double> rr = ( conj_right ? std::conj(r(s1,c1,f2)) : r(s1,c1,f2) );
	    rvec[s1][c1][f1] += N(f1,f2) * rr;	    
	  }
	}  
    //lM
    std::complex<double> lvec[4][3][2];
    for(int s1=0;s1<4;s1++)
      for(int c1=0;c1<3;c1++)
	for(int f1=0;f1<2;f1++){
	  lvec[s1][c1][f1] = std::complex<double>(0.0,0.0);
	  for(int s2=0;s2<4;s2++)
	    for(int c2=0;c2<3;c2++)
	      for(int f2=0;f2<2;f2++){
		std::complex<double> ll = ( conj_left ? std::conj(l(s2,c2,f2)) : l(s2,c2,f2) );
		lvec[s1][c1][f1] += ll * scf(s2,c2,f2,s1,c1,f1);
	      }
	}     
    std::complex<double> out(0.0,0.0);
    for(int s1=0;s1<4;s1++)
      for(int c1=0;c1<3;c1++)
	for(int f1=0;f1<2;f1++)
	  out += lvec[s1][c1][f1] * rvec[s1][c1][f1];
    into += out;
  }
};





template<typename SourceType, typename mf_Complex, int Remaining, int Idx=0>
struct _siteFmatRecurse{
  template<typename AccumVtype>
  static inline void doit(AccumVtype &into, const SourceType &src, const FlavorMatrixType sigma, const int p, const FlavorMatrixGeneral<mf_Complex> &lMr){
    FlavorMatrixGeneral<mf_Complex> phi;
    src.template getSource<Idx>().siteFmat(phi,p);
    phi.pl(sigma);
    
    doAccum(into[Idx], TransLeftTrace(lMr, phi));
    _siteFmatRecurse<SourceType,mf_Complex,Remaining-1,Idx+1>::doit(into,src,sigma,p,lMr);
  }
};
template<typename SourceType, typename mf_Complex, int Idx>
struct _siteFmatRecurse<SourceType,mf_Complex,0,Idx>{
  template<typename AccumVtype>
  static inline void doit(AccumVtype &into, const SourceType &src, const FlavorMatrixType sigma, const int p, const FlavorMatrixGeneral<mf_Complex> &lMr){}
};



//All of the inner products for G-parity can be separated into a part involving only the spin-color structure of the source and a part involving the flavor and smearing function.
//This case class implements the flavor/smearing function part and leaves the spin-color part to the derived class
template<typename mf_Complex, typename SourceType, typename SpinColorContractPolicy>
class GparityInnerProduct: public SpinColorContractPolicy{
  const SourceType &src;
  FlavorMatrixType sigma;
public:
  typedef SourceType InnerProductSourceType;

  GparityInnerProduct(const FlavorMatrixType &_sigma, const SourceType &_src): sigma(_sigma),src(_src){ }

  //When running with a multisrc type this returns the number of meson fields per timeslice = nSources
  template<typename S=SourceType>
  inline typename my_enable_if<has_enum_nSources<S>::value, int>::type mfPerTimeSlice() const{ return SourceType::nSources; }

  template<typename S=SourceType>
  inline typename my_enable_if<!has_enum_nSources<S>::value, int>::type mfPerTimeSlice() const{ return 1; }
  
  //Single source
  template<typename AccumType, typename ComplexType = mf_Complex>
  inline void operator()(AccumType &out, const SCFvectorPtr<ComplexType> &l, const SCFvectorPtr<ComplexType> &r, const int p, const int t) const{
#ifndef MEMTEST_MODE
    FlavorMatrixGeneral<ComplexType> lMr; //is vectorized
    this->spinColorContract(lMr,l,r);
    
    //Compute   lMr[f1,f3] s3[f1,f2] phi[f2,f3]  =   lMr^T[f3,f1] s3[f1,f2] phi[f2,f3] 
    FlavorMatrixGeneral<ComplexType> phi;
    src.siteFmat(phi,p);
    phi.pl(sigma);

    //Do the sum over the SIMD vectorized sites
    doAccum(out, TransLeftTrace(lMr, phi));
#endif
  }  
  
  //Multi source
  //Does out += op(l,r,p,t);
  template<typename AccumVtype, typename ComplexType = mf_Complex>
  inline void operator()(std::vector<AccumVtype> &out, const SCFvectorPtr<ComplexType> &l, const SCFvectorPtr<ComplexType> &r, const int p, const int t) const{
#ifndef MEMTEST_MODE
    FlavorMatrixGeneral<ComplexType> lMr; //is vectorized
    this->spinColorContract(lMr,l,r);

    _siteFmatRecurse<SourceType,ComplexType,SourceType::nSources>::doit(out,src,sigma,p,lMr);
#endif
  }
  
};

template<int smatidx, typename mf_Complex, typename SourceType, bool conj_left = true, bool conj_right=false>
class SCFspinflavorInnerProduct: public GparityInnerProduct<mf_Complex, SourceType, flavorMatrixSpinColorContract<smatidx,mf_Complex,conj_left,conj_right> >{
public:
  typedef SourceType InnerProductSourceType;
  
  SCFspinflavorInnerProduct(const FlavorMatrixType &_sigma, SourceType &_src):
    GparityInnerProduct<mf_Complex, SourceType, flavorMatrixSpinColorContract<smatidx,mf_Complex,conj_left,conj_right> >(_sigma,_src){}
};



template<typename SourceType, typename mf_Complex, int Remaining, int Idx=0>
struct _siteFmatRecurseShift{

  template<typename AccumVtype>
  static inline void doit(AccumVtype &into, const std::vector<SourceType*> &shifted_sources, const FlavorMatrixType sigma, const int p, const FlavorMatrixGeneral<mf_Complex> &lMr){
    FlavorMatrixGeneral<mf_Complex> phi;
    for(int i=0;i<shifted_sources.size();i++){
      shifted_sources[i]->template getSource<Idx>().siteFmat(phi,p);
      phi.pl(sigma);
      doAccum(into[Idx+SourceType::nSources*i], TransLeftTrace(lMr, phi));
    }
    _siteFmatRecurseShift<SourceType,mf_Complex,Remaining-1,Idx+1>::doit(into,shifted_sources,sigma,p,lMr);
  }
};
template<typename SourceType, typename mf_Complex, int Idx>
struct _siteFmatRecurseShift<SourceType,mf_Complex,0,Idx>{
  template<typename AccumVtype>
  static inline void doit(AccumVtype &into, const std::vector<SourceType*> &shifted_sources, const FlavorMatrixType sigma, const int p, const FlavorMatrixGeneral<mf_Complex> &lMr){}
};


template<typename SourceType,int Remaining, int Idx=0>
struct _shiftRecurse{
  static void inline doit(SourceType &what, const std::vector<int> &shift){
    shiftPeriodicField(  what.template getSource<Idx>().getSource(), what.template getSource<Idx>().getSource(), shift);
    _shiftRecurse<SourceType,Remaining-1,Idx+1>::doit(what,shift);
  }
};
template<typename SourceType, int Idx>
struct _shiftRecurse<SourceType,0,Idx>{
  static void inline doit(SourceType &what, const std::vector<int> &shift){}
};


template<typename mf_Complex, typename SourceType, typename SpinColorContractPolicy>
class GparitySourceShiftInnerProduct: public SpinColorContractPolicy{
  const SourceType &src;
  FlavorMatrixType sigma;
  std::vector< std::vector<int> > shifts; //should be a set of 3-vectors
  std::vector<SourceType*> shifted_sources; 
  std::vector<int> cur_shift;

  template<typename S=SourceType>
  inline typename my_enable_if<!has_enum_nSources<S>::value, void>::type
  shiftTheSource(SourceType &what, const std::vector<int> &shift){ shiftPeriodicField(what.getSource(),what.getSource(), shift); }

  template<typename S=SourceType>
  inline typename my_enable_if<has_enum_nSources<S>::value, void>::type
  shiftTheSource(SourceType &what, const std::vector<int> &shift){ _shiftRecurse<S,S::nSources>::doit(what, shift); }
  
  void shiftSource(SourceType &what, const std::vector<int> &shift){
    std::vector<int> actual_shift(shift);
    for(int i=0;i<3;i++) actual_shift[i] -= cur_shift[i]; //remove current shift in process
    shiftTheSource(what,actual_shift);
  }
  
public:
  typedef SourceType InnerProductSourceType;
  
  GparitySourceShiftInnerProduct(const FlavorMatrixType &_sigma, const SourceType &_src): sigma(_sigma),src(_src), shifts(0), cur_shift(3,0){ }
  GparitySourceShiftInnerProduct(const FlavorMatrixType &_sigma, const SourceType &_src, const std::vector< std::vector<int> > &_shifts): sigma(_sigma),src(_src), shifts(_shifts), cur_shift(3,0){ }

  ~GparitySourceShiftInnerProduct(){ for(int i=0;i<shifted_sources.size();i++) delete shifted_sources[i]; }
  
  //When running with a multisrc type this returns the number of meson fields per timeslice = nSources * nshift
  template<typename S=SourceType>
  inline typename my_enable_if<has_enum_nSources<S>::value, int>::type mfPerTimeSlice() const{ return shifts.size() * SourceType::nSources; } //indexed by  source_idx + nSources*shift_idx

  //When running with a single src type this returns the number of meson fields per timeslice = nshift
  template<typename S=SourceType>
  inline typename my_enable_if<!has_enum_nSources<S>::value, int>::type mfPerTimeSlice() const{ return shifts.size(); }

  void setShifts(const std::vector< std::vector<int> > &to_shifts){
    shifts = to_shifts;
    for(int i=0;i<shifted_sources.size();i++) delete shifted_sources[i];
    shifted_sources.resize(shifts.size());

    for(int i=0;i<shifts.size();i++){
      shifted_sources[i] = new SourceType(src);
      shiftSource(*shifted_sources[i], shifts[i]);
    }    
  }

  //Single src, output vector indexed by src shift index
  template<typename AccumVtype, typename ComplexType = mf_Complex, typename S = SourceType>
  inline typename my_enable_if<!has_enum_nSources<S>::value, void>::type
  operator()(AccumVtype &out, const SCFvectorPtr<ComplexType> &l, const SCFvectorPtr<ComplexType> &r, const int p, const int t) const{
#ifndef MEMTEST_MODE
    FlavorMatrixGeneral<ComplexType> lMr;
    this->spinColorContract(lMr,l,r);

    FlavorMatrixGeneral<ComplexType> phi;
    for(int i=0;i<shifts.size();i++){
      shifted_sources[i]->siteFmat(phi,p);
      phi.pl(sigma);
      doAccum(out[i],TransLeftTrace(lMr, phi));
    }
#endif
  }  

  //Multi src. output indexed by source_idx + nSources*shift_idx
  template<typename AccumVtype, typename ComplexType = mf_Complex, typename S = SourceType>
  inline typename my_enable_if<has_enum_nSources<S>::value, void>::type
  operator()(AccumVtype &out, const SCFvectorPtr<ComplexType> &l, const SCFvectorPtr<ComplexType> &r, const int p, const int t) const{
#ifndef MEMTEST_MODE
    FlavorMatrixGeneral<ComplexType> lMr;
    this->spinColorContract(lMr,l,r);

    _siteFmatRecurseShift<S,ComplexType,SourceType::nSources>::doit(out,shifted_sources,sigma,p,lMr);
#endif
  }
  
};




CPS_END_NAMESPACE

#endif
