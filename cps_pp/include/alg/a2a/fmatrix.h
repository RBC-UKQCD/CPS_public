#ifndef _FMATRIX_H
#define _FMATRIX_H

CPS_START_NAMESPACE

template<typename T>
class basicMatrix{
  T* tt;
  int rows, cols;
  int fsize; //number of elements
  
  void free(){
    if(tt!=NULL) sfree("basicMatrix","~basicMatrix","free",tt);
  }

  void alloc(const int &_rows, const int &_cols, T const* cp = NULL){
    if(_rows != rows || _cols != cols){
      free();
      rows = _rows; cols = _cols; fsize = rows*cols;
      tt = (T*)smalloc("basicMatrix", "basicMatrix", "alloc" , sizeof(T) * fsize);
    }
    if(cp != NULL) for(int i=0;i<fsize;i++) tt[i] = cp[i];
  }

public:
  basicMatrix(): rows(0),cols(0),fsize(0),tt(NULL){ }

  basicMatrix(const int &_rows, const int &_cols): tt(NULL){ 
    alloc(_rows,_cols);
  }
  basicMatrix(const basicMatrix<T> &r){
    alloc(r.rows,r.cols,r.tt);
  }
  
  T* ptr(){ return tt;}

  void resize(const int &_rows, const int &_cols){ alloc(_rows,_cols); }

  inline const T & operator()(const int &i, const int &j) const{ return tt[j + cols*i]; }
  inline T & operator()(const int &i, const int &j){ return tt[j + cols*i]; }

  inline const int &nRows() const{ return rows; }
  inline const int &nCols() const{ return cols; }

  ~basicMatrix(){
    free();
  }
};


//A matrix of complex numbers and some useful associated methods
template<typename mf_Complex>
class fMatrix{
  mf_Complex* tt;
  int rows, cols;
  int fsize; //number of elements
  
  void free_matrix(){
    if(tt!=NULL) sfree("fMatrix","~fMatrix","free",tt);
  }

  void alloc_matrix(const int _rows, const int _cols, mf_Complex const* cp = NULL){
    if(_rows != rows || _cols != cols){
      free_matrix();
      rows = _rows; cols = _cols; fsize = rows*cols;
      tt = (mf_Complex*)smalloc("fMatrix", "fMatrix", "alloc" , sizeof(mf_Complex) * fsize);
    }
    if(cp == NULL) zero();
    else for(int i=0;i<fsize;i++) tt[i] = cp[i];
  }

public:
  fMatrix(): rows(0),cols(0),fsize(0),tt(NULL){ }

  fMatrix(const int &_rows, const int &_cols): rows(0), cols(0), fsize(0),tt(NULL){ 
    alloc_matrix(_rows,_cols);
  }
  fMatrix(const fMatrix<mf_Complex> &r): rows(0), cols(0), fsize(0),tt(NULL){
    alloc_matrix(r.rows,r.cols,r.tt);
  }
  
  mf_Complex *ptr(){ return tt;}

  void resize(const int _rows, const int _cols){ alloc_matrix(_rows,_cols); }

  void zero(){ for(int i=0;i<fsize;i++) tt[i] = 0.0; }

  fMatrix & operator*=(const mf_Complex &r){ for(int i=0;i<fsize;i++) tt[i] *= r;  return *this; }
  fMatrix & operator*=(const typename mf_Complex::value_type &r){ for(int i=0;i<fsize*2;i++) ((typename mf_Complex::value_type*)tt)[i] *= r;  return *this; }

  fMatrix & operator+=(const fMatrix<mf_Complex> &r){ for(int i=0;i<fsize;i++) tt[i] += r.tt[i];  return *this; }
  
  inline const mf_Complex & operator()(const int i, const int j) const{ return tt[j + cols*i]; }
  inline mf_Complex & operator()(const int i, const int j){ return tt[j + cols*i]; }

  inline const int nRows() const{ return rows; }
  inline const int nCols() const{ return cols; }

  void nodeSum(){
    QMP_sum_array( (typename mf_Complex::value_type*)tt,2*fsize);
  }

  ~fMatrix(){
    free_matrix();
  }

  void write(const std::string &filename) const{
    FILE *p;
    if((p = Fopen(filename.c_str(),"w")) == NULL)
      ERR.FileA("fMatrix","write",filename.c_str());
    for(int r=0;r<rows;r++)
      for(int c=0;c<cols;c++)
	Fprintf(p,"%d %d %.16e %.16e\n",r,c, (*this)(r,c).real(), (*this)(r,c).imag());
    Fclose(p);
  }
};

//Rearrange an Lt*Lt matrix from ordering  tsnk, tsrc  to   tsrc,  tsep=tsnk-tsrc
template<typename mf_Complex>
void rearrangeTsrcTsep(fMatrix<mf_Complex> &m){
  int Lt = GJP.Tnodes()*GJP.TnodeSites();
  if(m.nRows()!=Lt || m.nCols()!=Lt) ERR.General("","rearrangeTsrcTsep(fMatrix<mf_Complex> &)","Expect an Lt*Lt matrix\n");

  fMatrix<mf_Complex> tmp(m);
  for(int tsnk=0;tsnk<Lt;tsnk++){
    for(int tsrc=0;tsrc<Lt;tsrc++){
      int tsep = (tsnk-tsrc+Lt) % Lt;
      m(tsrc,tsep) = tmp(tsnk,tsrc);
    }
  }
}



//A vector of complex numbers and some useful associated methods
template<typename mf_Complex>
class fVector{
  mf_Complex* tt;
  int fsize;
  
  void free_mem(){
    if(tt!=NULL) sfree("fVector","~fVector","free",tt);
  }

  void alloc_mem(const int _elems, mf_Complex const* cp = NULL){
    if(_elems != fsize){
      free_mem();
      fsize = _elems;
      tt = (mf_Complex*)smalloc("fVector", "fVector", "alloc" , sizeof(mf_Complex) * fsize);
    }
    if(cp == NULL) zero();
    else for(int i=0;i<fsize;i++) tt[i] = cp[i];
  }

public:
  fVector(): fsize(0),tt(NULL){ }

  fVector(const int _elems): fsize(0),tt(NULL){ 
    alloc_mem(_elems);
  }
  fVector(const fVector<mf_Complex> &r): fsize(0),tt(NULL){
    alloc_mem(r.fsize,r.tt);
  }
  
  mf_Complex *ptr(){ return tt;}

  void resize(const int _elems){ alloc_mem(_elems); }

  void zero(){ for(int i=0;i<fsize;i++) tt[i] = mf_Complex(0,0); }

  fVector & operator*=(const mf_Complex &r){ for(int i=0;i<fsize;i++) tt[i] *= r;  return *this; }
  fVector & operator*=(const typename mf_Complex::value_type &r){ for(int i=0;i<fsize*2;i++) ((typename mf_Complex::value_type*)tt)[i] *= r;  return *this; }

  inline const mf_Complex & operator()(const int &i) const{ return tt[i]; }
  inline mf_Complex & operator()(const int &i){ return tt[i]; }

  inline const int &size() const{ return fsize; }

  void nodeSum(){
    QMP_sum_array( (typename mf_Complex::value_type*)tt,2*fsize);
  }

  ~fVector(){
    free_mem();
  }

  void write(const std::string &filename) const{
    FILE *p;
    if((p = Fopen(filename.c_str(),"w")) == NULL)
      ERR.FileA("fVector","write",filename.c_str());
    for(int i=0;i<fsize;i++)
      Fprintf(p,"%d %.16e %.16e\n",i, tt[i].real(), tt[i].imag());
    Fclose(p);
  }
};




//Array of complex with optional threading
template<typename mf_Complex, typename AllocPolicy = StandardAllocPolicy>
class basicComplexArray: public AllocPolicy{
protected:
  int thread_size; //size of each thread unit
  int nthread;
  int size; //total size
  mf_Complex *con;
public:
  basicComplexArray(): size(0), con(NULL){}
  basicComplexArray(const int &_thread_size, const int &_nthread = 1): size(0), con(NULL){
    resize(_thread_size,_nthread);
  }
  void free_mem(){
    if(con != NULL){ AllocPolicy::_free(con); con = NULL; }
  }
  void resize(const int &_thread_size, const int &_nthread = 1){
    free_mem();
    thread_size = _thread_size; nthread = _nthread;
    size = _thread_size * _nthread;
    this->_alloc(&con, size*sizeof(mf_Complex));
    memset((void*)con, 0, size * sizeof(mf_Complex));
  }
  ~basicComplexArray(){
    free_mem();
  }
  inline const mf_Complex & operator[](const int i) const{ return con[i]; }
  inline mf_Complex & operator[](const int i){ return con[i]; }

  inline mf_Complex & operator()(const int i, const int thread){ return con[i + thread * thread_size]; }

  int nElementsTotal() const{
    return size;
  }
  int nElementsPerThread() const{
    return thread_size;
  }
  int nThreads() const{
    return nthread;
  }
    
  //Sum (reduce) over all threads
  void threadSum(){
    if(nthread == 1) return;
    basicComplexArray<mf_Complex,AllocPolicy> tmp(thread_size,1);
    
#pragma omp parallel for
    for(int i=0;i<thread_size;i++){
      for(int t=0;t<nthread;t++)
	tmp.con[i] += con[i + t*thread_size];
    }
    AllocPolicy::_free(con);
    con = tmp.con;
    nthread = 1;
    size = tmp.size;

    tmp.con = NULL;
  }
  void nodeSum(){
    globalSumComplex(this->con,size);
    //QMP_sum_array( (typename mf_Complex::value_type*)con,2*size);
  }
};


CPS_END_NAMESPACE

#endif
