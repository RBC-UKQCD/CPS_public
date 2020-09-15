#ifndef _A2A_DILUTIONS_H
#define _A2A_DILUTIONS_H

#include<alg/a2a/fmatrix.h>

CPS_START_NAMESPACE

//NOTE: As it stands all dilutions do not contain any information other than the A2Aparams which is the same for all dilutions. This allows us to create instances of one dilution from another or from an A2Aparams alonee 
//without propagating extra information. If you change this you will probably have to do some extra work!
struct modeIndexSet{
  int hit;
  int spin_color;
  int flavor;
  int time;
  modeIndexSet(): hit(-1),spin_color(-1),flavor(-1),time(-1){}
};

//Various stages of high-mode dilution. Undiluted indices are referred to as 'packed'
class StandardIndexDilution: public A2Aparams{
public:
  enum { UndilutedIndices = 0 };

  StandardIndexDilution(): A2Aparams(){}
  StandardIndexDilution(const A2AArg &_args): A2Aparams(_args){}
  StandardIndexDilution(const A2Aparams &_p): A2Aparams(_p){}

  //Map dilution indices to a 'mode' index for the high-mode part of the propagator
  //This is a combination of the hit index, timeslice-block index and spin,color (and flavor) indices
  inline int indexMap(const int &hit, const int &tblock, const int &spin_color, const int &flavor = 0) const{
    return spin_color + nspincolor*( flavor + nflavors*( tblock + ntblocks * hit) );
  }
  inline void indexUnmap(int idx, int &hit, int &tblock, int &spin_color, int &flavor) const{
    spin_color = idx % nspincolor; idx/=nspincolor;
    flavor = idx % nflavors; idx/=nflavors;
    tblock = idx % ntblocks; idx/=ntblocks;
    hit = idx;
  }
  inline int indexMap(const modeIndexSet &mset) const{
    return indexMap(mset.hit, mset.time, mset.spin_color, mset.flavor);
  }
  inline void indexUnmap(int idx,modeIndexSet &mset) const{
    return indexUnmap(idx,mset.hit, mset.time, mset.spin_color, mset.flavor);
  }
  int getNmodes() const{ return nv; }

  //Compute the mapping between a full StandardIndex and the packed indices. Returns a vector<bool> indicating which elements are non-zero and a map of StandardIndex -> packed index
  //The field with which this index has a generalized coordinate (spincolor,flavor,time) which must be provided in 'field_coord'. 
  //Note: which subset of the generalized coordinate are needed depends on the mapping. It will throw a runtime error if the information is not available (could possibly be made compile time with some thought)
  void getIndexMapping(std::vector<int> &map, std::vector<bool> &non_zeroes, const modeIndexSet &field_coord) const{
    non_zeroes.resize(nv,true);
    map.resize(nv); for(int i=0;i<nv;i++) map[i] = i; //trivial mapping
  }

  static std::string name(){ return "StandardIndexDilution"; }
};

//An index containing only the hit, flavor and spin-color indices
class TimePackedIndexDilution: public A2Aparams{

public:
  enum { UndilutedIndices = 1 };

  TimePackedIndexDilution(): A2Aparams(){}
  TimePackedIndexDilution(const A2AArg &_args): A2Aparams(_args){}
  TimePackedIndexDilution(const A2Aparams &_p): A2Aparams(_p){}

  inline int indexMap(const int &hit, const int &spin_color, const int &flavor = 0) const{
    return spin_color + nspincolor * (flavor + nflavors*hit);
  }
  inline void indexUnmap(int idx, int &hit, int &spin_color, int &flavor) const{
    spin_color = idx % nspincolor; idx/=nspincolor;
    flavor = idx % nflavors; idx/=nflavors;
    hit = idx;
  }
  inline int indexMap(const modeIndexSet &mset) const{
    return indexMap(mset.hit, mset.spin_color, mset.flavor);
  }
  inline void indexUnmap(int idx,modeIndexSet &mset) const{
    return indexUnmap(idx,mset.hit, mset.spin_color, mset.flavor);
  }
  int getNmodes() const{ return nl + nhits*nspincolor*nflavors; }

  //Compute the mapping between a full StandardIndex and the packed indices. Returns a vector<bool> indicating which elements are non-zero and a map of StandardIndex -> packed index
  //The field with which this index has a generalized coordinate (spincolor,flavor,time) which must be provided in 'field_coord'. 
  //Note: which subset of the generalized coordinate are needed depends on the mapping. It will throw a runtime error if the information is not available (could possibly be made compile time with some thought)
  void getIndexMapping(std::vector<int> &map, std::vector<bool> &non_zeroes, const modeIndexSet &field_coord) const{
    assert(field_coord.time != -1);

    non_zeroes.resize(nv,false); map.resize(nv);
    for(int i=0;i<nl;i++){ non_zeroes[i] = true; map[i] = i; } //low modes mapping trivial
    
    const StandardIndexDilution &unpacked = static_cast< const StandardIndexDilution &>(*this);

    //Loop over packed modes and fill gaps
    for(int p = 0; p < nhits*nspincolor*nflavors; p++){
      modeIndexSet pparams;  indexUnmap(p, pparams); pparams.time = field_coord.time;
      int idx_packed = nl + p;
      int idx_unpacked = nl + unpacked.indexMap(pparams);
      non_zeroes[idx_unpacked] = true;
      map[idx_unpacked] = idx_packed;
    }
  }
  static std::string name(){ return "TimePackedIndexDilution"; }
};

//An index containing only the hit and spin-color indices
class TimeFlavorPackedIndexDilution: public A2Aparams{

public:
  enum { UndilutedIndices = 2 };

  TimeFlavorPackedIndexDilution(): A2Aparams(){}
  TimeFlavorPackedIndexDilution(const A2AArg &_args): A2Aparams(_args){}
  TimeFlavorPackedIndexDilution(const A2Aparams &_p): A2Aparams(_p){}

  //Mapping used by W_fftw for high 'modes', where only the spin/color index has been diluted out
  inline int indexMap(const int &hit, const int &spin_color) const{
    return spin_color + nspincolor * hit;
  }
  inline void indexUnmap(int idx, int &hit, int &spin_color) const{
    spin_color = idx % nspincolor; idx/=nspincolor;
    hit = idx;
  }
  inline int indexMap(const modeIndexSet &mset) const{
    return indexMap(mset.hit, mset.spin_color);
  }
  inline void indexUnmap(int idx,modeIndexSet &mset) const{
    return indexUnmap(idx,mset.hit, mset.spin_color);
  }
  int getNmodes() const{ return nl + nhits*nspincolor; }

  //Compute the mapping between a full StandardIndex and the packed indices. Returns a vector<bool> indicating which elements are non-zero and a map of StandardIndex -> packed index
  //The field with which this index has a generalized coordinate (spincolor,flavor,time) which must be provided in 'field_coord'. 
  //Note: which subset of the generalized coordinate are needed depends on the mapping. It will throw a runtime error if the information is not available (could possibly be made compile time with some thought)
  void getIndexMapping(std::vector<int> &map, std::vector<bool> &non_zeroes, const modeIndexSet &field_coord) const{
    assert(field_coord.time != -1);
    assert(field_coord.flavor != -1);

    non_zeroes.resize(nv,false); map.resize(nv);
    for(int i=0;i<nl;i++){ non_zeroes[i] = true; map[i] = i; } //low modes mapping trivial
    
    const StandardIndexDilution &unpacked = static_cast< const StandardIndexDilution &>(*this);

    //Loop over packed modes and fill gaps
    for(int p = 0; p < nhits*nspincolor; p++){
      modeIndexSet pparams;  indexUnmap(p, pparams); pparams.time = field_coord.time; pparams.flavor = field_coord.flavor; 
      int idx_packed = nl + p;
      int idx_unpacked = nl + unpacked.indexMap(pparams);
      non_zeroes[idx_unpacked] = true;
      map[idx_unpacked] = idx_packed;
    }
  }
  static std::string name(){ return "TimeFlavorPackedIndexDilution"; }
};


//An index containing only the hit index
class FullyPackedIndexDilution: public A2Aparams{
public:
  enum { UndilutedIndices = 3 };

  FullyPackedIndexDilution(): A2Aparams(){}
  FullyPackedIndexDilution(const A2AArg &_args): A2Aparams(_args){}
  FullyPackedIndexDilution(const A2Aparams &_p): A2Aparams(_p){}

  //Mapping used by W_fftw for high 'modes', where only the spin/color index has been diluted out
  inline int indexMap(const int &hit) const{
    return hit;
  }
  inline void indexUnmap(int idx, int &hit) const{
    hit = idx;
  }
  inline int indexMap(const modeIndexSet &mset) const{
    return indexMap(mset.hit);
  }
  inline void indexUnmap(int idx,modeIndexSet &mset) const{
    return indexUnmap(idx,mset.hit);
  }

  int getNmodes() const{ return nl + nhits*nspincolor; }

  //Compute the mapping between a full StandardIndex and the packed indices. Returns a vector<bool> indicating which elements are non-zero and a map of StandardIndex -> packed index
  //The field with which this index has a generalized coordinate (spincolor,flavor,time) which must be provided in 'field_coord'. 
  //Note: which subset of the generalized coordinate are needed depends on the mapping. It will throw a runtime error if the information is not available (could possibly be made compile time with some thought)
  void getIndexMapping(std::vector<int> &map, std::vector<bool> &non_zeroes, const modeIndexSet &field_coord) const{
    assert(field_coord.time != -1);
    assert(field_coord.flavor != -1);
    assert(field_coord.spin_color != -1);

    non_zeroes.resize(nv,false); map.resize(nv);
    for(int i=0;i<nl;i++){ non_zeroes[i] = true; map[i] = i; } //low modes mapping trivial
    
    const StandardIndexDilution &unpacked = static_cast< const StandardIndexDilution &>(*this);

    //Loop over packed modes and fill gaps
    for(int p = 0; p < nhits; p++){
      modeIndexSet pparams;  indexUnmap(p, pparams); pparams.time = field_coord.time; pparams.flavor = field_coord.flavor; pparams.spin_color = field_coord.spin_color;
      int idx_packed = nl + p;
      int idx_unpacked = nl + unpacked.indexMap(pparams);
      non_zeroes[idx_unpacked] = true;
      map[idx_unpacked] = idx_packed;
    }
  }
  static std::string name(){ return "FullyPackedIndexDilution"; }
};











CPS_END_NAMESPACE
#endif
