#ifndef _SPARSE_MATRIX_H
#define _SPARSE_MATRIX_H

//The total number of high mode indices is  nhits * nspincolors * nflavors * Lt, and for each of those there are Lt 3d spin/color/flavor solution vectors
//Depending on the level of dilution the vast majority of entries are zero; in fact only the V vector has no non-zero entries
//It is therefore imperative to have sparse vectors and matrices to hold these


// inline void computeUnion(std::vector<bool> &into, const std::vector<bool> &a, const std::vector<bool> &b){
//   assert(a.size()==b.size());
//   into.resize(a.size());
//   for(int i=0;i<into.size();i++) into[i] = a[i] && b[i];
// }

// template<typename mf_Float>
// class A2AsparseVector{
//   std::vector<bool> zeroes;
//   std::vector<CPSfermion4D<mf_Float> *> elems;
//   CPSfermion4D<mf_Float> zero; //defaults to zero anyway
//   size_t len;
  
//   void clear(){
//     zeroes.resize(0);
//     for(int i=0;i<elems.size();i++) if(elems[i] != NULL) delete elems[i];
//     elems.resize(0);
//   }

// public:
//   //Setup a sparse vector of length _len. Reserves space for pointers but does not assign anything
//   SparseVector(const size_t &_len): len(_len), zeroes(_len,true), elems(_len,NULL){
//   }
//   SparseVector(): len(0){
//   }
//   const CPSfermion4D<mf_Float> & operator()(const size_t &idx){
//     if(zeroes[idx]) return zero;
//     else return *elems[idx];
//   }
//   void set(const size_t &idx, const CPSfermion4D<mf_Float> &to){
//     zeroes[idx] = false;
//     elems[idx] = new CPSfermion4D<mf_Float>(to);
//   }
//   const size_t &size() const{ return len; } //the supposed maximum length of the vector (not checked ofc!)
  
//   //Deletes everything already in vector
//   void resize(const size_t &to){ 
//     clear();
//     len = to; 
//     zeroes.resize(len,true);
//     elems.resize(len,NULL);
//   }

//   ~SparseVector(){ clear(); }
// };

// template<typename mf_Float>
// class A2AsparseMesonField{
//   std::vector<bool> row_zeroes;
//   std::vector<bool> col_zeroes;

//   std::vector<std::complex<mf_Float> > elems;
//   size_t side_length;
  
//   void clear(){
//     zeroes.resize(0);
//     elems.resize(0);
//   }

// public:
//   //Setup a sparse vector of length _len. Reserves space for pointers but does not assign anything
//   SparseVector(const size_t &_len): len(_len), zeroes(_len,true), elems(_len,std::complex<mf_Float>(0,0) ){
//   }
//   SparseVector(): len(0){
//   }
//   const CPSfermion4D<mf_Float> & operator()(const size_t &idx){
//     return 

//     if(zeroes[idx]) return zero;
//     else return *elems[idx];
//   }
//   void set(const size_t &idx, const CPSfermion4D<mf_Float> &to){
//     zeroes[idx] = false;
//     elems[idx] = new CPSfermion4D<mf_Float>(to);
//   }
//   const size_t &size() const{ return len; } //the supposed maximum length of the vector (not checked ofc!)
  
//   //Deletes everything already in vector
//   void resize(const size_t &to){ 
//     clear();
//     len = to; 
//     zeroes.resize(len,true);
//     elems.resize(len,NULL);
//   }

//   ~SparseVector(){ clear(); }
// };









#if 0
//The total number of high mode indices is  nhits * nspincolors * nflavors * Lt, and for each of those there are Lt 3d spin/color/flavor solution vectors
//Depending on the level of dilution the vast majority of entries are zero; in fact only the V vector has no non-zero entries
//It is therefore imperative to have sparse vectors and matrices to hold these
template<typename T>
class BasicZeroGenerator{
public:
  static void zero(T &z){ z = 0.0; }
};

template<typename ElementType, template <typename> class ZeroGenerator = BasicZeroGenerator>
class SparseVector{
  std::map<size_t,ElementType> sparse_map;
  const ElementType zero;
  size_t len;
public:
  SparseVector(const size_t &_len): len(_len){
    ZeroGenerator<ElementType>::zero(zero);
  }
  SparseVector(): len(0){
    ZeroGenerator<ElementType>::zero(zero);
  }
  const ElementType & operator()(const size_t &idx){
    typename std::map<size_t,ElementType>::const_iterator it = sparse_map.find(idx);
    return it == sparse_map.end() ? zero : *it;
  }
  void set(const size_t &idx, const ElementType &to){
    sparse_map[idx] = to;
  }
  const size_t &size() const{ return len; } //the supposed maximum length of the vector (not checked ofc!)
  
  void resize(const size_t &to){ len = to; }
};

template<typename ElementType, template <typename> class ZeroGenerator = BasicZeroGenerator>
class SparseMatrix{
  typedef SparseVector<ElementType,ZeroGenerator> Vtype;

  typename std::map<size_t,Vtype> sparse_col_map;
  const ElementType zero;

  size_t rows, cols;
public:
  SparseMatrix(const size_t &_rows, const size_t &_cols){
    ZeroGenerator<ElementType>::zero(zero);
  }
    
  const ElementType & operator()(const size_t &i1, const size_t &i2){
    typename std::map<size_t,Vtype>::const_iterator it = sparse_col_map.find(i1);
    return it == sparse_col_map.end() ? zero : it->operator()(i2);
  }
  void set(const size_t &i1, const size_t &i2, const ElementType &to){
    if(sparse_col_map.count(i1) == 0) .....
    sparse_col_map[i1](i2) = to;
  }

  const size_t &nRows() const{ return rows; }
  const size_t &nCols() const{ return cols; }
};
  
template<typename ElementType, template <typename> class ZeroGenerator = BasicZeroGenerator>
void dot(SparseVector<ElementType,ZeroGenerator> &out, const SparseMatrix<ElementType,ZeroGenerator> &M, const SparseVector<ElementType,ZeroGenerator> &v){
  
}
#endif


#endif



