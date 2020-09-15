#ifndef CK_A2A_UTILS
#define CK_A2A_UTILS

#include <alg/alg_fix_gauge.h>
#include <alg/a2a/gsl_wrapper.h>
#include <alg/a2a/template_wizardry.h>
#include <util/spincolorflavormatrix.h>

CPS_START_NAMESPACE

//3x3 complex vector multiplication with different precision matrices and vectors
template<typename VecFloat, typename MatFloat>
void colorMatrixMultiplyVector(VecFloat* y, const MatFloat* u, const VecFloat* x){
	*y     =  *u      * *x     - *(u+1)  * *(x+1) + *(u+2)  * *(x+2)
		- *(u+3)  * *(x+3) + *(u+4)  * *(x+4) - *(u+5)  * *(x+5);
	*(y+1) =  *u      * *(x+1) + *(u+1)  * *x     + *(u+2)  * *(x+3)
		+ *(u+3)  * *(x+2) + *(u+4)  * *(x+5) + *(u+5)  * *(x+4);
	*(y+2) =  *(u+6)  * *x     - *(u+7)  * *(x+1) + *(u+8)  * *(x+2)
		- *(u+9)  * *(x+3) + *(u+10) * *(x+4) - *(u+11) * *(x+5);
	*(y+3) =  *(u+6)  * *(x+1) + *(u+7)  * *x     + *(u+8)  * *(x+3)
		+ *(u+9)  * *(x+2) + *(u+10) * *(x+5) + *(u+11) * *(x+4);
	*(y+4) =  *(u+12) * *x     - *(u+13) * *(x+1) + *(u+14) * *(x+2)
		- *(u+15) * *(x+3) + *(u+16) * *(x+4) - *(u+17) * *(x+5);
	*(y+5) =  *(u+12) * *(x+1) + *(u+13) * *x     + *(u+14) * *(x+3)
		+ *(u+15) * *(x+2) + *(u+16) * *(x+5) + *(u+17) * *(x+4);
}
//M^\dagger v

//0 ,1    2 ,3    4 ,5
//6 ,7    8 ,9    10,11
//12,13   14,15   16,17
//->
//0 ,-1   6 ,-7   12,-13
//2 ,-3   8 ,-9   14,-15 
//4 ,-5   10,-11  16,-17

template<typename VecFloat, typename MatFloat>
void colorMatrixDaggerMultiplyVector(VecFloat* y, const MatFloat* u, const VecFloat* x){
	*y     =  *u      * *x     + *(u+1)  * *(x+1) + *(u+6)  * *(x+2)	  
		+ *(u+7)  * *(x+3) + *(u+12)  * *(x+4) + *(u+13)  * *(x+5);
	*(y+1) =  *u      * *(x+1) - *(u+1)  * *x     + *(u+6)  * *(x+3)	  
		- *(u+7)  * *(x+2) + *(u+12)  * *(x+5) - *(u+13)  * *(x+4);	
	*(y+2) =  *(u+2)  * *x     + *(u+3)  * *(x+1) + *(u+8)  * *(x+2)	  
		+ *(u+9)  * *(x+3) + *(u+14) * *(x+4) + *(u+15) * *(x+5);	
	*(y+3) =  *(u+2)  * *(x+1) - *(u+3)  * *x     + *(u+8)  * *(x+3)	  
		- *(u+9)  * *(x+2) + *(u+14) * *(x+5) - *(u+15) * *(x+4);	
	*(y+4) =  *(u+4) * *x     + *(u+5) * *(x+1) + *(u+10) * *(x+2)	  
		+ *(u+11) * *(x+3) + *(u+16) * *(x+4) + *(u+17) * *(x+5);	
	*(y+5) =  *(u+4) * *(x+1) - *(u+5) * *x     + *(u+10) * *(x+3)	  
		- *(u+11) * *(x+2) + *(u+16) * *(x+5) - *(u+17) * *(x+4);
}

//Array *= with cps::Float(=double) input and arbitrary precision output
template<typename FloatOut,typename FloatIn>
void VecTimesEquFloat(FloatOut *out, FloatIn *in, const Float fac, const int len) 
{
#pragma omp parallel for
	for(int i = 0; i < len; i++) out[i] = in[i] * fac;
}

inline void getNodeWork(const int work, int &node_work, int &node_off, bool &do_work, const bool node_local = false){
  if(node_local){ node_work = work; node_off = 0; do_work = true; return; } //node does all the work

  int nodes = 1; for(int i=0;i<5;i++) nodes *= GJP.Nodes(i);
  int me = UniqueID();

  //Stolen from BFM :)
  int basework = work/nodes;
  int backfill = nodes-(work % nodes);
  node_work = (work+me)/nodes;
  node_off  = basework * me;
  if ( me > backfill ) 
    node_off+= (me-backfill);
  if(node_work > 0) do_work = true;
}

inline void thread_work(int &my_work, int &my_offset, const int total_work, const int me, const int team){
  my_work = total_work/team;
  my_offset = me * my_work;
  
  int rem = total_work - my_work * team;
  if(me < rem){
    ++my_work; //first rem threads mop up the remaining work
    my_offset += me; //each thread before me has gained one extra unit of work
  }else my_offset += rem; //after the first rem threads, the offset shift is uniform
}


inline void compute_overlap(std::vector<bool> &out, const std::vector<bool> &a, const std::vector<bool> &b){
  assert(a.size()==b.size());
  out.resize(a.size());
  for(int i=0;i<a.size();i++) out[i] = a[i] && b[i];
}

class NullObject
{
 public:
  NullObject(){}
};

//A class inheriting from this type must have template parameter T as a double or float
#define EXISTS_IF_DOUBLE_OR_FLOAT(T) public my_enable_if<is_double_or_float<mf_Float>::value,NullObject>::type

//Functions for performing global and timeslice sums of single or double precision quantities. Daiqian had to implement these himself as CPS can only do this with the Float=double type

// My global sum
template <typename T>
void QMP_sum_array(T *result, int len){
#ifdef USE_QMP
  if(sizeof(T) == sizeof(double)) {
    QMP_sum_double_array((double*)result, len);
  } else if(sizeof(T) == sizeof(float)) {
    QMP_sum_float_array((float*)result, len);
  } else {
    QMP_error("QMP_sum_array::data type not supported!\n");
  }
#else
  //CK: This only works for single-node code
  int nodes = 1; for(int i=0;i<4;i++) nodes *= cps::GJP.Nodes(i);
  if(nodes != 1){
    cps::ERR.General("","QMP_sum_array(T *result, int len)","Only implemented for QMP on parallel machines");
  }
  //do nothing!
#endif
}

#ifndef USE_QMP
  inline void QMP_sum_double_array(double *result, int len){
    //CK: This only works for single-node code
    int nodes = 1; for(int i=0;i<4;i++) nodes *= cps::GJP.Nodes(i);
    if(nodes != 1){
      cps::ERR.General("","QMP_sum_double_array fake definition","Not implemented on parallel machines: use QMP!");
    }
  }
  inline void QMP_sum_float_array(float *result, int len){
    //CK: This only works for single-node code
    int nodes = 1; for(int i=0;i<4;i++) nodes *= cps::GJP.Nodes(i);
    if(nodes != 1){
      cps::ERR.General("","QMP_sum_float_array fake definition","Not implemented on parallel machines: use QMP!");
    }
  }
#endif

//Look for contiguous blocks of indices in the idx_map, output a list of start,size pairs
inline void find_contiguous_blocks(std::vector<std::pair<int,int> > &blocks, const int idx_map[], int map_size){
  blocks.resize(0);
  std::pair<int,int> block(0,1); //start, size
  int prev = idx_map[0];
  for(int j_packed=1;j_packed<map_size;j_packed++){
    int j_unpacked = idx_map[j_packed];
    if(j_unpacked == prev+1){
      ++block.second;
    }else{
      blocks.push_back(block);
      block.first = j_packed;
      block.second = 1;      
    }
    prev = j_unpacked;
  }
  blocks.push_back(block);

  int sum = 0;
  for(int b=0;b<blocks.size();b++){
    //printf("Block %d, start %d, size %d\n",b,blocks[b].first,blocks[b].second);
    sum += blocks[b].second;
  }
  if(sum != map_size)
    ERR.General("find_contiguous_blocks","","Sum of block sizes %d, expect %d\n",sum,map_size);
}

template<typename T>
inline void resize_2d(std::vector<std::vector<T> > &v, const size_t i, const size_t j){
  v.resize(i);
  for(int a=0;a<i;a++) v[a].resize(j);
}
template<typename T>
inline void resize_3d(std::vector<std::vector<std::vector<T> > > &v, const size_t i, const size_t j, const size_t k){
  v.resize(i);
  for(int a=0;a<i;a++){
    v[a].resize(j);
    for(int b=0;b<j;b++)
      v[a][b].resize(k);
  }
}

inline std::complex<double> GSLtrace(const SpinColorFlavorMatrix& a, const SpinColorFlavorMatrix& b){

  const int scf_size = 24;
  std::complex<double> _a[scf_size][scf_size];
  std::complex<double> _bT[scf_size][scf_size];   //In-place transpose of b so rows are contiguous
  for(int i=0;i<scf_size;i++){
    int rem = i;
    int ci = rem % 3; rem /= 3;
    int si = rem % 4; rem /= 4;
    int fi = rem;
    
    for(int j=0;j<scf_size;j++){
      rem = j;
      int cj = rem % 3; rem /= 3;
      int sj = rem % 4; rem /= 4;
      int fj = rem;
      
      _bT[i][j] = b(sj,cj,fj, si,ci,fi);
      _a[i][j] = a(si,ci,fi, sj,cj,fj);
    }
  }

  double* ad = (double*)&_a[0][0];
  double* bd = (double*)&_bT[0][0];

  gsl_block_complex_struct ablock;
  ablock.size = 24*24;
  ablock.data = ad;

  gsl_vector_complex arow; //single row of a
  arow.block = &ablock;
  arow.owner = 0;
  arow.size = 24;
  arow.stride = 1;
  
  gsl_block_complex_struct bblock;
  bblock.size = 24*24;
  bblock.data = bd;

  gsl_vector_complex bcol; //single col of b
  bcol.block = &bblock;
  bcol.owner = 0;
  bcol.size = 24;
  bcol.stride = 1;

  //gsl_blas_zdotu (const gsl_vector_complex * x, const gsl_vector_complex * y, gsl_complex * dotu)
  //   //  a[0][0]*b[0][0] + a[0][1]*b[1][0] + a[0][2]*b[2][0] + ...
  //   //+ a[1][0]*b[0][1] + a[1][1]*b[1][1] + a[1][2]*b[2][1] + ....
  //   //...

  std::complex<double> out(0.0);
  gsl_complex tmp;
  for(int i=0;i<24;i++){
    arow.data = ad + 24*2*i; //i'th row offset
    bcol.data = bd + 24*2*i; //i'th col offset (remember we transposed it)

    gsl_blas_zdotu(&arow, &bcol, &tmp);
    reinterpret_cast<double(&)[2]>(out)[0] += GSL_REAL(tmp);
    reinterpret_cast<double(&)[2]>(out)[1] += GSL_IMAG(tmp);
  }
  return out;
}


//For a Nrows*Ncols matrix 'to' with elements in the standard order  idx=(Ncols*i + j), poke a submatrix into it with origin (i0,j0) and size (ni,nj)
template<typename T>
void pokeSubmatrix(T* to, const T* sub, const int Nrows, const int Ncols, const int i0, const int j0, const int ni, const int nj, const bool threaded = false){
  #define DOIT \
    for(int row = i0; row < i0+ni; row++){ \
      T* to_block = to + row*Ncols + j0;	  \
      const T* from_block = sub + (row-i0)*nj;	\
      memcpy(to_block,from_block,nj*sizeof(T));	\
    }
  if(threaded){
#pragma omp parallel for
    DOIT;
  }else{
    DOIT;
  }
  #undef DOIT
}
//For a Nrows*Ncols matrix 'from' with elements in the standard order  idx=(Ncols*i + j), get a submatrix with origin (i0,j0) and size (ni,nj) and store in sub
template<typename T>
void getSubmatrix(T* sub, const T* from, const int Nrows, const int Ncols, const int i0, const int j0, const int ni, const int nj, const bool threaded = false){
  #define DOIT \
    for(int row = i0; row < i0+ni; row++){		\
      const T* from_block = from + row*Ncols + j0;	\
      T* to_block = sub + (row-i0)*nj;			\
      memcpy(to_block,from_block,nj*sizeof(T));		\
    }
  if(threaded){
#pragma omp parallel for
    DOIT;
  }else{
    DOIT;
  }
  #undef DOIT
}


//Simple test allocator to find out when memory is allocated
template <typename T>
class mmap_allocator: public std::allocator<T>{
public:
  typedef size_t size_type;
  typedef T* pointer;
  typedef const T* const_pointer;

  template<typename _Tp1>
  struct rebind{
    typedef mmap_allocator<_Tp1> other;
  };

  pointer allocate(size_type n, const void *hint=0){
    fprintf(stderr, "Alloc %d bytes.\n", n*sizeof(T));
    return std::allocator<T>::allocate(n, hint);
  }

  void deallocate(pointer p, size_type n){
    fprintf(stderr, "Dealloc %d bytes (%p).\n", n*sizeof(T), p);
    return std::allocator<T>::deallocate(p, n);
  }

  mmap_allocator() throw(): std::allocator<T>() { fprintf(stderr, "Hello allocator!\n"); }
  mmap_allocator(const mmap_allocator &a) throw(): std::allocator<T>(a) { }
  template <class U>                    
  mmap_allocator(const mmap_allocator<U> &a) throw(): std::allocator<T>(a) { }
  ~mmap_allocator() throw() { }
};


CPS_END_NAMESPACE
#ifdef ARCH_BGQ
#include <spi/include/kernel/memory.h>
#else
#include <sys/sysinfo.h>
#endif
CPS_START_NAMESPACE

inline double byte_to_MB(const int b){
  return double(b)/1024./1024.;
}

inline void printMem(){
#ifdef ARCH_BGQ
  #warning "printMem using ARCH_BGQ"
  uint64_t shared, persist, heapavail, stackavail, stack, heap, guard, mmap;
  Kernel_GetMemorySize(KERNEL_MEMSIZE_SHARED, &shared);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_PERSIST, &persist);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAPAVAIL, &heapavail);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_STACKAVAIL, &stackavail);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_STACK, &stack);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAP, &heap);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_GUARD, &guard);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_MMAP, &mmap);

  if(!UniqueID()){
    printf("printMem: Allocated heap: %.2f MB, avail. heap: %.2f MB\n", (double)heap/(1024*1024),(double)heapavail/(1024*1024));
    printf("printMem: Allocated stack: %.2f MB, avail. stack: %.2f MB\n", (double)stack/(1024*1024), (double)stackavail/(1024*1024));
    printf("printMem: Memory: shared: %.2f MB, persist: %.2f MB, guard: %.2f MB, mmap: %.2f MB\n", (double)shared/(1024*1024), (double)persist/(1024*1024), (double)guard/(1024*1024), (double)mmap/(1024*1024));
  }
#else
#warning "printMem using NOARCH"
  /* unsigned long totalram;  /\* Total usable main memory size *\/ */
  /* unsigned long freeram;   /\* Available memory size *\/ */
  /* unsigned long sharedram; /\* Amount of shared memory *\/ */
  /* unsigned long bufferram; /\* Memory used by buffers *\/ */
  /* unsigned long totalswap; /\* Total swap space size *\/ */
  /* unsigned long freeswap;  /\* swap space still available *\/ */
  /* unsigned short procs;    /\* Number of current processes *\/ */
  /* unsigned long totalhigh; /\* Total high memory size *\/ */
  /* unsigned long freehigh;  /\* Available high memory size *\/ */
  /* unsigned int mem_unit;   /\* Memory unit size in bytes *\/ */

  struct sysinfo myinfo;
  sysinfo(&myinfo);
  double total_mem = myinfo.mem_unit * myinfo.totalram;
  total_mem /= (1024.*1024.);
  double free_mem = myinfo.mem_unit * myinfo.freeram;
  free_mem /= (1024.*1024.);
  
  if(!UniqueID()){
    printf("printMem: Memory: total: %.2f MB, avail: %.2f MB, used %.2f MB\n",total_mem, free_mem, total_mem-free_mem);
  }

  //# define PRINT_MALLOC_INFO    //Use of int means this is garbage for large memory systems
# ifdef PRINT_MALLOC_INFO
  struct mallinfo mi;
  mi = mallinfo();

  // int arena;     /* Non-mmapped space allocated (bytes) */
  // int ordblks;   /* Number of free chunks */
  // int smblks;    /* Number of free fastbin blocks */
  // int hblks;     /* Number of mmapped regions */
  // int hblkhd;    /* Space allocated in mmapped regions (bytes) */
  // int usmblks;   /* Maximum total allocated space (bytes) */
  // int fsmblks;   /* Space in freed fastbin blocks (bytes) */
  // int uordblks;  /* Total allocated space (bytes) */
  // int fordblks;  /* Total free space (bytes) */
  // int keepcost;  /* Top-most, releasable space (bytes) */

  if(!UniqueID()){
    printf("printMem: Malloc info: arena %f MB, ordblks %d, smblks %d, hblks %d, hblkhd %f MB, fsmblks %f MB, uordblks %f MB, fordblks %f MB, keepcost %f MB\n",
	   byte_to_MB(mi.arena), mi.ordblks, mi.smblks, mi.hblks, byte_to_MB(mi.hblkhd), byte_to_MB(mi.fsmblks), byte_to_MB(mi.uordblks), byte_to_MB(mi.fordblks), byte_to_MB(mi.keepcost) );
  }

# endif

  //# define PRINT_MALLOC_STATS  Also doesn't work well
# ifdef PRINT_MALLOC_STATS
  if(!UniqueID()) malloc_stats();
# endif
  
#endif
}

//Skip gauge fixing and set all gauge fixing matrices to unity
void gaugeFixUnity(Lattice &lat, const FixGaugeArg &fix_gauge_arg){
  FixGaugeType fix = fix_gauge_arg.fix_gauge_kind;
  int start = fix_gauge_arg.hyperplane_start;
  int step = fix_gauge_arg.hyperplane_step;
  int num = fix_gauge_arg.hyperplane_num;

  int h_planes[num];
  for(int i=0; i<num; i++) h_planes[i] = start + step * i;

  lat.FixGaugeAllocate(fix, num, h_planes);
  
#pragma omp parallel for
  for(int sf=0;sf<(GJP.Gparity()+1)*GJP.VolNodeSites();sf++){
    //s + vol*f
    int s = sf % GJP.VolNodeSites();
    int f = sf / GJP.VolNodeSites();
    
    const Matrix* mat = lat.FixGaugeMatrix(s,f);
    if(mat == NULL) continue;
    else{
      Matrix* mm = const_cast<Matrix*>(mat); //evil, I know, but it saves duplicating the accessor (which is overly complicated)
      mm->UnitMatrix();
    }
  }
}

//Set the complex number at pointer p to a random value of a chosen type
//Uses the current LRG for the given FermionFieldDimension. User should choose the range and the particular site-RNG themselves beforehand
template<typename mf_Float>
class RandomComplex{};

//Only for float and double, hence I have to control its access
template<typename mf_Float>
class RandomComplexBase{
 protected:
  template<typename T> friend class RandomComplex;
  
  static void rand(mf_Float *p, const RandomType type, const FermionFieldDimension frm_dim){
    static const Float PI = 3.14159265358979323846;
    Float theta = LRG.Urand(frm_dim);
  
    switch(type) {
    case UONE:
      p[0] = cos(2. * PI * theta);
      p[1] = sin(2. * PI * theta);
      break;
    case ZTWO:
      p[0] = theta > 0.5 ? 1 : -1;
      p[1] = 0;
      break;
    case ZFOUR:
      if(theta > 0.75) {
	p[0] = 1;
	p[1] = 0;
      }else if(theta > 0.5) {
	p[0] = -1;
	p[1] = 0;
      }else if(theta > 0.25) {
	p[0] = 0;
	p[1] = 1;
      }else {
	p[0] = 0;
	p[1] = -1;
      }
      break;
    default:
      ERR.NotImplemented("RandomComplexBase", "rand(...)");
    }
  }
};



template<typename T>
class RandomComplex<std::complex<T> > : public RandomComplexBase<T>{
public:
  static void rand(std::complex<T> *p, const RandomType &type, const FermionFieldDimension &frm_dim){
    RandomComplexBase<T>::rand( (T*)p, type, frm_dim);
  }
};

template<typename T, typename T_class>
struct _mult_sgn_times_i_impl{};

template<typename T>
struct _mult_sgn_times_i_impl<T,complex_double_or_float_mark>{
  inline static T doit(const int sgn, const T &val){
    return T( -sgn * val.imag(), sgn * val.real() ); // sign * i * val
  }
};

#ifdef USE_GRID
template<typename T>
struct _mult_sgn_times_i_impl<T,grid_vector_complex_mark>{
  inline static T doit(const int sgn, const T &val){
    return sgn == -1 ? timesMinusI(val) : timesI(val);
  }
};
#endif

// template<typename T, typename my_enable_if<is_complex_double_or_float<T>::value,int>::type = 0> //for standard complex types
// inline T multiplySignTimesI(const int sgn, const T &val){
//   return T( -sgn * val.imag(), sgn * val.real() ); // sign * i * val
// }

// #ifdef USE_GRID
// template<typename T, typename my_enable_if<is_grid_vector_complex<T>::value,int>::type = 0> //for Grid complex types
// inline T multiplySignTimesI(const int sgn, const T &val){
//   return sgn == -1 ? timesMinusI(val) : timesI(val);
// }
// #endif

template<typename T>
inline T multiplySignTimesI(const int sgn, const T &val){
  return _mult_sgn_times_i_impl<T,typename ComplexClassify<T>::type>::doit(sgn,val);
}

template<typename T, typename ComplexClass>
struct _cconj{};

template<typename T>
struct _cconj<T,complex_double_or_float_mark>{
  static inline T doit(const T &in){ return std::conj(in); }
};
#ifdef USE_GRID
template<typename T>
struct _cconj<T,grid_vector_complex_mark>{
  static inline T doit(const T &in){ return Grid::conjugate(in); }
};
#endif

template<typename T>
inline T cconj(const T& in){
  return _cconj<T,typename ComplexClassify<T>::type>::doit(in);
}

template<typename T>
std::complex<double> convertComplexD(const std::complex<T> &what){
  return what;
}
#ifdef USE_GRID
std::complex<double> convertComplexD(const Grid::vComplexD &what){  
  return Reduce(what);
}
std::complex<double> convertComplexD(const Grid::vComplexF &what){  
  return Reduce(what);
}
#endif


template<typename T>
void globalSumComplex(std::complex<T>* v, const int n){
  QMP_sum_array( (T*)v,2*n);
}
#ifdef USE_GRID
template<typename T>
struct _globalSumComplexGrid{
  static inline void doit(T *v, const int n){
    typedef typename T::scalar_type scalar_type; //an std::complex type
    typedef typename scalar_type::value_type floatType;    
    int vmult = sizeof(T)/sizeof(scalar_type);
    floatType * ptr = (floatType *)v; 
    QMP_sum_array(ptr,2*n*vmult);
  }
};

void globalSumComplex(Grid::vComplexD* v, const int n){
  _globalSumComplexGrid<Grid::vComplexD>::doit(v,n);
}
void globalSumComplex(Grid::vComplexF* v, const int n){
  _globalSumComplexGrid<Grid::vComplexF>::doit(v,n);
}
#endif

//The base G-parity momentum vector for quark fields with arbitrary sign
inline void GparityBaseMomentum(int p[3], const int sgn){
  for(int i=0;i<3;i++)
    if(GJP.Bc(i) == BND_CND_GPARITY)
      p[i] = sgn;
    else p[i] = 0;
}

#ifdef USE_MPI
//get MPI rank of this node
inline int getMyMPIrank(){
  int my_mpi_rank;
  int ret = MPI_Comm_rank(MPI_COMM_WORLD, &my_mpi_rank);
  if(ret != MPI_SUCCESS) ERR.General("A2AmesonField","read","Comm_rank failed\n");
  return my_mpi_rank;
}
//get the MPI rank of the node with UniqueID() == 0
inline int getHeadMPIrank(){
  int head_mpi_rank;
  int rank_tmp = UniqueID() == 0 ? getMyMPIrank() : 0;
  int ret = MPI_Allreduce(&rank_tmp,&head_mpi_rank, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD); //node is now the MPI rank corresponding to UniqueID == _node
  if(ret != MPI_SUCCESS) ERR.General("A2AmesonField","read","Reduce failed\n");
  return head_mpi_rank;
}

inline int node_lex(const int* coor, const int ndir){
  int out = 0;
  for(int i=ndir-1;i>=0;i--){
    out *= GJP.Nodes(i);
    out += coor[i];
  }
  return out;  
}

//Generate map to convert lexicographic node index from GJP to an MPI rank in MPI_COMM_WORLD
inline void getMPIrankMap(std::vector<int> &map){
  int nodes = 1;
  int my_node_coor[5];
  for(int i=0;i<5;i++){
    nodes*= GJP.Nodes(i);
    my_node_coor[i] = GJP.NodeCoor(i);
  }
  const int my_node_lex = node_lex( my_node_coor, 5 );
  const int my_mpi_rank = getMyMPIrank();

  int *node_map_send = (int*)malloc(nodes*sizeof(int));
  memset(node_map_send,0,nodes*sizeof(int));
  node_map_send[my_node_lex] = my_mpi_rank;

  map.resize(nodes);
  int ret = MPI_Allreduce(node_map_send, &map[0], nodes, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  assert(ret == MPI_SUCCESS);
  free(node_map_send);
}
#endif 
//Invert 3x3 complex matrix. Expect elements accessible as  row*3 + col
//0 1 2
//3 4 5
//6 7 8

//+ - +
//- + -
//+ - +

template<typename Zout, typename Zin>
void z3x3_invert(Zout* out, Zin const* in){
  out[0] = in[4]*in[8]-in[7]*in[5];
  out[1] = -in[3]*in[8]+in[6]*in[5];
  out[2] = in[3]*in[7]-in[6]*in[4];

  out[3] = -in[1]*in[8]+in[7]*in[2];
  out[4] = in[0]*in[8]-in[6]*in[2];
  out[5] = -in[0]*in[7]+in[6]*in[1];

  out[6] = in[1]*in[5]-in[4]*in[2];
  out[7] = -in[0]*in[5]+in[3]*in[2];
  out[8] = in[0]*in[4]-in[3]*in[1];
  
  Zout det = in[0]*out[0] + in[1]*out[1] + in[2]*out[2];

  out[0] /= det; out[1] /= det; out[2] /= det;
  out[3] /= det; out[4] /= det; out[5] /= det;
  out[6] /= det; out[7] /= det; out[8] /= det;
}

//A class that owns data via a pointer that has an assignment and copy constructor which does a deep copy.
template<typename T>
class PtrWrapper{
  T* t;
public:
  inline T& operator*(){ return *t; }
  inline T* operator->(){ return t; }
  inline T const& operator*() const{ return *t; }
  inline T const* operator->() const{ return t; }
  
  inline PtrWrapper(): t(NULL){};
  inline PtrWrapper(T* _t): t(_t){}
  inline ~PtrWrapper(){ if(t!=NULL) delete t; }

  inline const bool assigned() const{ return t != NULL; }
  
  inline void set(T* _t){
    if(t!=NULL) delete t;
    t = _t;
  }
  inline void free(){
    if(t!=NULL) delete t;
    t = NULL;
  }
  
  //Deep copies
  inline PtrWrapper(const PtrWrapper &r): t(NULL){
    if(r.t != NULL) t = new T(*r.t);
  }

  inline PtrWrapper & operator=(const PtrWrapper &r){
    if(t!=NULL){ delete t; t = NULL; }
    if(r.t!=NULL) t = new T(*r.t);
  } 
};



CPS_END_NAMESPACE

#endif
