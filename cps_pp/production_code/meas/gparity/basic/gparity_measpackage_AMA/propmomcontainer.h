#ifndef _PROP_MOM_CONTAINER_H
#define _PROP_MOM_CONTAINER_H

#include <set>
#include "propwrapper.h"

CPS_START_NAMESPACE

class PropMomContainer{
  std::map<std::string, PropWrapper> props;
public:
  void insert(const PropWrapper &prop, const std::string &tag){
    if(props.count(tag) != 0) ERR.General("PropMomContainer","insert","Attempting to insert duplicate of prop with tag '%s'\n", tag.c_str());
    props[tag] = prop;
  }
  PropWrapper & get(const std::string &tag){
    if(props.count(tag) == 0) ERR.General("PropMomContainer","get","Could not find prop with tag '%s'\n", tag.c_str());
    return props[tag];
  } 
  const PropWrapper & get(const std::string &tag) const{
    std::map<std::string, PropWrapper>::const_iterator it = props.find(tag);
    if(it == props.end()) ERR.General("PropMomContainer","get","Could not find prop with tag '%s'\n", tag.c_str());
    return it->second;
  } 

  void clear(){
    //This container takes ownership of the memory associated with the propagators.
    //As a propagator may appear in multiple entries we must avoid double deletion
    std::set<QPropW*> deleted;
    for(std::map<std::string, PropWrapper>::iterator it = props.begin(); it != props.end(); it++){
      for(int f=0;f<1+GJP.Gparity();f++){
	QPropW* p = it->second.getPtr(f);
	if(!deleted.count(p)){ 
	  deleted.insert(p);
	  delete p;
	}
      }
    }
    props.clear();
  }
  void printAllTags() const{
    if(!UniqueID()){
      printf("Propagators stored:\n");
      for(std::map<std::string, PropWrapper>::const_iterator it = props.begin(); it != props.end(); it++)
	std::cout << it->first << '\n';
    }
  }

  //For debugging purposes print the time dependence of the 3d volume sum of the matrix norm2 of the prop, for each prop
  void writePropNormTdep(const std::string &results_dir){
    SpinColorFlavorMatrix tmp_scf[omp_get_max_threads()];
    WilsonMatrix tmp_sc[omp_get_max_threads()];
    int vol3d = GJP.VolNodeSites()/GJP.TnodeSites();

    std::string filename = results_dir + "/prop_norms.dat";

    FILE *p;
    if((p = Fopen(filename.c_str(),"w")) == NULL)
      ERR.FileA("PropMomContainer","writePropNormTdep",filename.c_str());

    for(std::map<std::string, PropWrapper>::const_iterator it = props.begin(); it != props.end(); it++){
      basicComplexArray<Rcomplex> pnorms(GJP.TnodeSites()*GJP.Tnodes(), omp_get_max_threads());
      const std::string &tag = it->first;
      const PropWrapper &prop = it->second;
      
      for(int t=0;t<GJP.TnodeSites();t++){
	int t_glb = t + GJP.TnodeCoor()*GJP.TnodeSites();
#pragma omp parallel for
	for(int x=0;x<vol3d;x++){
	  int me = omp_get_thread_num();
	  if(GJP.Gparity()){
	    prop.siteMatrix(tmp_scf[me],x + vol3d*t);
	    pnorms(t_glb,me) += tmp_scf[me].norm();
	  }else{
	    prop.siteMatrix(tmp_sc[me],x + vol3d*t);
	    pnorms(t_glb,me) += tmp_sc[me].norm();
	  }
	}      
      }
      pnorms.threadSum();
      pnorms.nodeSum();

      Fprintf(p,"%s",tag.c_str());
      for(int t=0;t<GJP.TnodeSites()*GJP.Tnodes();t++)
	Fprintf(p," %.16e",pnorms[t].real());
      Fprintf(p,"\n");
    }
    Fclose(p);
  }

  //Takes ownership and deletes props when destroyed
  ~PropMomContainer(){
    clear();
  }
};



CPS_END_NAMESPACE

#endif
