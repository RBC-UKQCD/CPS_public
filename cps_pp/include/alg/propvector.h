CPS_START_NAMESPACE
#ifndef PROP_VECTOR_H
#define PROP_VECTOR_H
CPS_END_NAMESPACE

#include<config.h>
#include <alg/propagatorcontainer.h>
#include<vector>

CPS_START_NAMESPACE

//container for multiple props with destructor that deletes all props
template<typename T>
class PointerArray{
public:
  T & operator[](const int &idx);

  void set(const int &idx, T* to);

  T& append(T* v);

  void clear();
  int size() const;
  
  PointerArray();
  virtual ~PointerArray();

private:
  std::vector<T*> ptrs;
};


class PropVector: public PointerArray<PropagatorContainer>{
public:
  PropagatorContainer & addProp(PropagatorArg &arg);
  PropVector();
};

#ifdef USE_BFM
class LanczosVector: public PointerArray<LanczosContainer>{
public:
  LanczosContainer & add(LanczosContainerArg &arg);
  LanczosVector();
};
#endif

template<typename T>
PointerArray<T>::PointerArray(){ }
template<typename T>
PointerArray<T>::~PointerArray(){ for(int i=0;i<ptrs.size();i++) delete ptrs[i]; }

template<typename T>
T & PointerArray<T>::operator[](const int &idx){ return *ptrs[idx]; }

template<typename T>
void PointerArray<T>::clear(){
  for(int i=0;i<ptrs.size();i++) delete ptrs[i];
  ptrs.resize(0);
}
template<typename T>
int PointerArray<T>::size() const{ return ptrs.size(); }

template<typename T>
void PointerArray<T>::set(const int &idx, T* to){
  delete ptrs[idx]; ptrs[idx] = to;
}

template<typename T>
T& PointerArray<T>::append(T* v){
  ptrs.push_back(v);
  return *v;
}


#endif
CPS_END_NAMESPACE
