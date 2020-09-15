//Implementation of basic matrix, vector types with underlying SIMD parallelization via Grid
//Note: I don't put them in CPS namespace to avoid conflicts.
#ifdef USE_GRID

#ifndef _GRID_SIMDLINALG_H
#define _GRID_SIMDLINALG_H

#include <Grid/Grid.h>

namespace Grid_A2A{
  
  //Basic vector dot product for testing
  template<typename T>
  T dot_basic(const std::vector<T> &A, const std::vector<T> &B){
    assert(A.size() == B.size());
    T out = 0.0;
    for(uint i=0;i<A.size();i++)
      out += A[i] * B[i];
    return out;
  }


  //Mapping of SIMD types
  template<typename T>
  struct SIMDvtype{};

  template<>
  struct SIMDvtype< std::complex<double> >{
    typedef Grid::vComplexD Vtype;
  };
  template<>
  struct SIMDvtype< std::complex<float> >{
    typedef Grid::vComplexF Vtype;
  };


  //SIMD paralellized vector
  template<typename T>
  class SIMDvector{
  public:
    typedef typename SIMDvtype<T>::Vtype Vtype;
  private:
    std::vector<Vtype, Grid::alignedAllocator<Vtype> > v;
  
    const int nsimd; //how many T can be fit into a packed Vtype

    uint _size_packed;
    uint _size;

    inline void alloc(uint n){
      //n is the unpacked size
      _size = n;
      bool overspill = (n % nsimd != 0); //size doesn't divide evenly over SIMD vectors, add a final SIMD vector padded with zeroes to mop up
      _size_packed = n/nsimd + (overspill ? 1:0);
      v.resize(_size_packed);
    }

  public:
    SIMDvector(): v(0), nsimd(Vtype::Nsimd()), _size(0){}
  
    SIMDvector(uint n): nsimd(Vtype::Nsimd()){
      alloc(n);
    }

    SIMDvector(SIMDvector<T> &&r): nsimd(Vtype::Nsimd()){
      _size = r._size; _size_packed = r._size_packed;
      v = std::move(r.v);
    }
    SIMDvector(const SIMDvector<T> &r){
      alloc(r._size);
      for(int i=0;i<_size_packed;i++) v[i] = r.v[i];
    }
  
    SIMDvector &operator=(const SIMDvector<T> &r){
      if(r._size != _size){
	alloc(r._size);
      }
      for(int i=0;i<_size_packed;i++) v[i] = r.v[i];
      return *this;
    }
  
    inline uint size() const{ return _size; }
    inline uint packedSize() const{ return _size_packed; }

    //sz is the unpacked size!
    inline void resize(uint sz){
      alloc(sz);
    }

    void import(const T* from, const size_t size){
      if(size!=_size){
	alloc(size);
      }
      int overspill = _size % nsimd;
      uint to = _size_packed - (overspill > 0 ? 1:0);
    
      T* d = const_cast<T*>(from);
      for(int i=0;i<to;i++){
	vset(v[i],d);
	d+=nsimd;
      }
    
      if(overspill){
	T block[nsimd];
	for(int i=0;i<overspill;i++) 
	  block[i] = d[i];
	for(int i=overspill; i<nsimd; i++)
	  block[i] = 0;
	vset(v[to],block);
      }
    }

    inline void import(const std::vector<T> &from){
      import(&from[0],from.size());
    }

    //Access to packed members
    inline Vtype & operator[](uint i){ return v[i]; }
    inline const Vtype & operator[](uint i) const{ return v[i]; }
  };

  //SIMD parallelized vector dot product
  template<typename T>
  T dot_simd(const SIMDvector<T> &A, const SIMDvector<T> &B){
    typedef typename SIMDvector<T>::Vtype Vtype;
    assert(A.size() == B.size());

    T out = 0.0;
    Vtype tmp;
    for(uint i=0;i<A.packedSize();i++){
      tmp = A[i] * B[i];
      out += Reduce(tmp);
    }
    return out;
  }

  template<typename T>
  T dot_simd(const std::vector<T> &A, const std::vector<T> &B){
    assert(A.size() == B.size());
    SIMDvector<T> Av; Av.import(A);
    SIMDvector<T> Bv; Bv.import(B);

    return dot_simd(Av,Bv);
  }


  //A basic matrix implementation for testing
  template<typename T>
  class Matrix{
    T* m;
    uint _rows;
    uint _cols;
    uint size;

    inline void alloc(uint r, uint c){
      _rows = r; _cols = c; size = r*c;
      m = (T*)malloc(size * sizeof(T));
    }
    inline void freeme(){
      if(m!=NULL) free(m);
    }
  public:
    Matrix(): m(NULL){}

    Matrix(uint r, uint c){
      alloc(r,c);
    }
    ~Matrix(){
      freeme();
    }
    Matrix(Matrix<T> &&r){
      _rows = r._rows;
      _cols = r._cols;
      size = r.size;
      m = r.m;
      r.m = NULL;
    }
    Matrix(const Matrix<T> &r){
      alloc(r._rows,r._cols);
      for(int i=0;i<size;i++) m[i] = r.m[i];
    }
  
    Matrix &operator=(const Matrix<T> &r){
      if(r._rows != _rows || r._cols != _cols){
	freeme();
	alloc(r._rows,r._cols);
      }
      for(int i=0;i<size;i++) m[i] = r.m[i];
      return *this;
    }
    inline const T & operator[](uint i) const{ return m[i]; }

    inline T & operator()(uint i, uint j){ return m[i*_cols + j]; }
    inline const T & operator()(uint i, uint j) const{ return m[i*_cols + j]; }
  
    inline uint rows() const{ return _rows; }
    inline uint cols() const{ return _cols; }  
  };

  //Comparison of matrices for testing
  inline bool equals(const Matrix<std::complex<double> > &A, const Matrix<std::complex<double> > &B, double tol=1e-9){
    if(A.rows()!=B.rows() || A.cols()!=B.cols()) return false;
    for(int i=0;i<A.rows()*A.cols();i++){
      double rd = fabs( std::real(A[i]) - std::real(B[i]) );
      if(rd > tol){ 
	printf("WARNING: Difference at elem %d real part: val A %g val B %g  diff %g\n",i, std::real(A[i]), std::real(B[i]), rd); 
	return false;
      }
      rd = fabs( std::imag(A[i]) - std::imag(B[i]) );
      if(rd > tol){ 
	printf("WARNING: Difference at elem %d imag part: val A %g val B %g  diff %g\n",i, std::imag(A[i]), std::imag(B[i]), rd); 
	return false;
      }
    }
    return true;
  }

  //Basic implementation of matrix multiplication
  template<typename T>
  Matrix<T> mult_basic(const Matrix<T> &A, const Matrix<T> &B){
    assert(A.cols() == B.rows());
    Matrix<T> out(A.rows(),B.cols());

    for(uint i=0;i<A.rows();i++){
      for(uint k=0;k<B.cols();k++){
	out(i,k) = 0.0;

	for(uint j=0;j<A.cols();j++)
	  out(i,k) += A(i,j)*B(j,k);
      }
    }
    return out;
  };


  //A matrix with SIMD vectors as rows (i.e. its element span an entire row)
  template<typename T>
  class SIMDrowMatrix{
    std::vector< SIMDvector<T> > _rows;
    typedef typename SIMDvector<T>::Vtype Vtype;
  public:
  
    void import(const Matrix<T> &from){
      _rows.resize(from.rows());
    
      int nsimd = Vtype::Nsimd(); //number of Ts in a V :)
      for(int r=0;r<_rows.size();r++){
	_rows[r].resize(from.cols());
      
	T* block_base = &from(r,0); //row major so column elements are consecutive
	for(int c_packed = 0; c_packed < _rows[r].packedSize(); c_packed++){
	  vset(_rows[r][c_packed], block_base);
	  block_base += nsimd;
	}
      }
    }

    inline const SIMDvector<T> & getRow(uint i) const{ return _rows[i]; }
    inline SIMDvector<T> & getRow(uint i){ return _rows[i]; }

    uint rows() const{ return _rows.size(); }
    uint cols() const{ return _rows[0].size(); }
  };

  //A matrix with SIMD vectors as columns
  template<typename T>
  class SIMDcolMatrix{
    std::vector< SIMDvector<T> > _cols;
    typedef typename SIMDvector<T>::Vtype Vtype;
  public:
  
    void import(const Matrix<T> &from){
      _cols.resize(from.cols());
    
      int nsimd = Vtype::Nsimd(); //number of Ts in a V :)
      T tmp[nsimd];
      int row_off = from.cols(); //row major

      for(int c=0;c<_cols.size();c++){
	_cols[c].resize(from.rows());

	T* col_p = &from(0,c);
	for(int r_packed = 0; r_packed < _cols[c].packedSize(); r_packed++){
	  for(int i=0;i<nsimd;i++){ //copy nsimd elements from the column to temp memory
	    tmp[i] = *col_p;
	    col_p += row_off;
	  }
	  vset(_cols[c][r_packed],tmp);
	}
      }
    }

    inline const SIMDvector<T> & getCol(uint i) const{ return _cols[i]; }
    inline SIMDvector<T> & getCol(uint i){ return _cols[i]; }

    uint rows() const{ return _cols[0].size(); }
    uint cols() const{ return _cols.size(); }
  };

  //SIMD implementation of matrix multiplication
  template<typename T>
  Matrix<T> mult_simd(const SIMDrowMatrix<T> &A, const SIMDcolMatrix<T> &B){
    Matrix<T> out(A.rows(),B.cols());

    for(uint i=0;i<A.rows();i++)
      for(uint k=0;k<B.cols();k++)
	out(i,k) = dot_simd(A.getRow(i),B.getCol(k));
    return out;
  };

  


  // inline std::complex<double> gridTrace(const SpinColorFlavorMatrix& a, const SpinColorFlavorMatrix& b){
  //   typedef typename SIMDvtype<std::complex<double> >::Vtype Vtype;
  //   const int scf_size = 24;
  //   const int nsimd = Vtype::Nsimd();
  //   assert( scf_size % nsimd == 0 ); //works for AVX1 and AVX2 (2 complex<double> per SIMD vector, nblocks = 12) and AVX512 (4 complex<double> per SIMD vector, nblocks = 6), will also work with AVX1024 when it comes around
  //   const int nblocks = scf_size / nsimd;

  //   //In-place transpose of b so rows are contiguous
  //   std::complex<double> bT[scf_size][scf_size];
  //   for(int i=0;i<scf_size;i++){
  //     int rem = i;
  //     int ci = rem % 3; rem /= 3;
  //     int si = rem % 4; rem /= 4;
  //     int fi = rem;

  //     for(int j=0;j<scf_size;j++){
  // 	rem = j;
  // 	int cj = rem % 3; rem /= 3;
  // 	int sj = rem % 4; rem /= 4;
  // 	int fj = rem;

  // 	bT[i][j] = b(sj,cj,fj, si,ci,fi);
  //     }
  //   }
    
  //   //  a[0][0]*b[0][0] + a[0][1]*b[1][0] + a[0][2]*b[2][0] + ...
  //   //+ a[1][0]*b[0][1] + a[1][1]*b[1][1] + a[1][2]*b[2][1] + ....
  //   //...

  //   //linearize index
  //   //  a[0]*b[0] + a[1]*b[24] + a[2]*b[48] + ....
  //   //+ a[24]*b[1] + a[25]*b[25] + a[26]*b[49] + ....
  //   //...

  //   Vtype res; zeroit(res);
  //   for(int i=0;i<scf_size;i++){
      



//     const int nfullblocks = scf_size / nsimd;
//     const int overspill = scf_size % nsimd;
//     const int nblocks = nfullblocks + (overspill ? 1:0);
    


//   Rcomplex trace(0.0,0.0);
//   for(int f2=0;f2<2;f2++){
//     for(int f1=0;f1<2;f1++){
//       for(int s2=0;s2<4;++s2){
// 	for(int c2=0;c2<3;++c2){
// 	  for(int s1=0;s1<4;++s1){
// 	    for(int c1=0;c1<3;++c1){
// 	      trace += a(s1,c1,f1,s2,c2,f2)*b(s2,c2,f2,s1,c1,f1);
// 	    }
// 	  }
// 	}
//       }
//     }
//   }
//   return trace;
// }




};





#endif

#endif
