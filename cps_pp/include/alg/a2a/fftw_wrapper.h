#ifndef FFTW_WRAPPER
#define FFTW_WRAPPER

//#ifdef USE_FFTW
#include <fftw3.h>
//#endif

CPS_START_NAMESPACE 

//Basic, slow implementation if FFTW not available




//Wrap fftw library for multiple templated float types
template<typename mf_Float>
class FFTWwrapper{};

template<>
class FFTWwrapper<float>{
public:
  typedef fftwf_complex complexType;
  typedef fftwf_plan planType;
  
  static planType plan_many_dft(int rank, const int *n, int howmany,
			 complexType *in, const int *inembed,
			 int istride, int idist,
			 complexType *out, const int *onembed,
			 int ostride, int odist,
			 int sign, unsigned flags){
    return fftwf_plan_many_dft(rank,n,howmany,in,inembed,istride,idist,out,onembed,ostride,odist,sign,flags);
  }
  static complexType* alloc_complex(size_t n){ return fftwf_alloc_complex(n); }
  
  static void execute_dft(const planType p, complexType *in, complexType *out){
    fftwf_execute_dft(p,in,out);
  }
  static void destroy_plan(planType plan){ fftwf_destroy_plan(plan); }

  static void cleanup(){ fftwf_cleanup(); }
  static void free(void *p){ fftwf_free(p); }
};
template<>
class FFTWwrapper<double>{
public:
  typedef fftw_complex complexType;
  typedef fftw_plan planType;
  
  static planType plan_many_dft(int rank, const int *n, int howmany,
			 complexType *in, const int *inembed,
			 int istride, int idist,
			 complexType *out, const int *onembed,
			 int ostride, int odist,
			 int sign, unsigned flags){
    return fftw_plan_many_dft(rank,n,howmany,in,inembed,istride,idist,out,onembed,ostride,odist,sign,flags);
  }
  static complexType* alloc_complex(size_t n){ return fftw_alloc_complex(n); }
  
  static void execute_dft(const planType p, complexType *in, complexType *out){
    fftw_execute_dft(p,in,out);
  }
  static void destroy_plan(planType plan){ fftw_destroy_plan(plan); }

  static void cleanup(){ fftw_cleanup(); }
  static void free(void *p){ fftw_free(p); }
};

template<typename mf_Float>
class FFTplanContainer{
  typedef typename FFTWwrapper<mf_Float>::planType planType;
  typedef typename FFTWwrapper<mf_Float>::complexType complexType;
  planType plan;
  bool plan_allocd;
public:
  FFTplanContainer(): plan_allocd(false){}

  FFTplanContainer(int rank, const int *n, int howmany,
		   complexType *in, const int *inembed,
		   int istride, int idist,
		   complexType *out, const int *onembed,
		   int ostride, int odist,
		   int sign, unsigned flags): plan_allocd(false){
    setPlan(rank,n,howmany,in,inembed,istride,idist,out,onembed,ostride,odist,sign,flags);
  }
  
  void setPlan(int rank, const int *n, int howmany,
		  complexType *in, const int *inembed,
		  int istride, int idist,
		  complexType *out, const int *onembed,
		  int ostride, int odist,
		  int sign, unsigned flags){
    if(plan_allocd) FFTWwrapper<mf_Float>::destroy_plan(plan);

    plan = FFTWwrapper<mf_Float>::plan_many_dft(rank,n,howmany,in,inembed,istride,idist,out,onembed,ostride,odist,sign,flags);
    plan_allocd = true;
  }

  planType & getPlan(){ return plan; }
  const planType & getPlan() const{ return plan; }
  
  ~FFTplanContainer(){
    if(plan_allocd) FFTWwrapper<mf_Float>::destroy_plan(plan);
  } 
};



CPS_END_NAMESPACE 
#endif

