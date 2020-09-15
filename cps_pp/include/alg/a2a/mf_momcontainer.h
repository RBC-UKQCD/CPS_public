#ifndef _PIPI_MFCONTAINER_H
#define _PIPI_MFCONTAINER_H

#include<alg/a2a/threemomentum.h>
CPS_START_NAMESPACE

//We must construct meson fields with a number of different total momenta. This class holds the fields and allows access in a flexible and transparent manner
//The ThreeMomentum is the pion momentum
//The class owns the meson fields it stores, and they are deleted when it is destroyed
template<typename mf_Policies>
class MesonFieldMomentumContainer{
private:
  typedef std::map<ThreeMomentum, std::vector<A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorVfftw> >* > MapType; //vector is the time index of the meson field
  MapType mf; //store pointers so we don't have to copy
  
public:
  std::vector<A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorVfftw> > const* getPtr(const ThreeMomentum &p) const{
    typename MapType::const_iterator it = mf.find(p);
    if(it == mf.end()) return NULL;
    else return it->second;
  }
  const std::vector<A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorVfftw> >& get(const ThreeMomentum &p) const{
    typename MapType::const_iterator it = mf.find(p);
    if(it == mf.end()){
      std::cout << "MesonFieldMomentumContainer::get Cannot find meson field with ThreeMomentum " << p.str() << std::endl; std::cout.flush();
      exit(-1);
    }      
    else return *it->second;
  }
  std::vector<A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorVfftw> >& get(const ThreeMomentum &p){
    typename MapType::iterator it = mf.find(p);
    if(it == mf.end()){
      std::cout << "MesonFieldMomentumContainer::get Cannot find meson field with ThreeMomentum " << p.str() << std::endl; std::cout.flush();
      exit(-1);
    }
    else return *it->second;
  }

  void printMomenta(std::ostream &os) const{
    for(typename MapType::const_iterator it = mf.begin(); it != mf.end(); it++)
      os << it->first.str() << "\n";
  }

  std::vector<A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorVfftw> >& copyAdd(const ThreeMomentum &p, const std::vector<A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorVfftw> > &mfield){   
    mf[p] = new std::vector<A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorVfftw> >(mfield);
    return *mf[p];
  }
  std::vector<A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorVfftw> >& moveAdd(const ThreeMomentum &p, std::vector<A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorVfftw> > &mfield){
    mf[p] = new std::vector<A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorVfftw> >(mfield.size());
    for(int i=0;i<mfield.size();i++) mf[p]->operator[](i).move(mfield[i]);
    mfield.swap(std::vector<A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorVfftw> >());
    return *mf[p];
  }
  
  bool contains(const ThreeMomentum &p) const{ return mf.count(p) != 0; }

  ~MesonFieldMomentumContainer(){
    for(typename MapType::iterator it = mf.begin(); it != mf.end(); it++) delete it->second;
  }
};


CPS_END_NAMESPACE
#endif

