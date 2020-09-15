#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <util/qcdio.h>
#ifdef PARALLEL
#include <comms/sysfunc_cps.h>
#endif
#include <comms/scu.h>
#include <comms/glb.h>

#include <util/lattice.h>
#include <util/time_cps.h>
#include <util/smalloc.h>

#include <util/command_line.h>

#include<unistd.h>
#include<config.h>

#ifdef USE_BFM 

//CK: these are redefined by BFM (to the same values)
#undef ND
#undef SPINOR_SIZE
#undef HALF_SPINOR_SIZE
#undef GAUGE_SIZE
#endif

#ifdef USE_BFM
#include <alg/lanc_arg.h>
#include <alg/eigen/Krylov_5d.h>
#endif

#include <alg/propmanager.h>

CPS_START_NAMESPACE


PropVector PropManager::props;
#ifdef USE_BFM
LanczosVector PropManager::eig; 
#endif

PropagatorContainer & PropManager::getProp(const char *tag){
  for(int i=0;i<props.size();i++) if(props[i].tagEquals(tag)) return props[i];
  ERR.General("PropManager","getProp(const char *tag)","Prop '%s' does not exist!\n",tag);
};

PropagatorContainer & PropManager::addProp(PropagatorArg &arg){
  return props.addProp(arg);
}

#ifdef USE_BFM
LanczosContainer & PropManager::addLanczos(LanczosContainerArg &arg){
  return eig.add(arg);
}
#endif

void PropManager::setup(JobPropagatorArgs &prop_args){
#ifdef USE_BFM
  for(int i=0;i< prop_args.lanczos.lanczos_len; i++) addLanczos(prop_args.lanczos.lanczos_val[i]);
#endif
  for(int i=0;i< prop_args.props.props_len; i++) addProp(prop_args.props.props_val[i]);

  //after all props are loaded in, the attributes of props that combine other props are copied over
  for(int i=0;i<props.size();i++) if(props[i].type() == QPROPW_TYPE) props[i].convert<QPropWcontainer>().propCombSetupAttrib();
}

void PropManager::startNewTraj(){
  for(int i=0;i<props.size();i++) props[i].deleteProp();
#ifdef USE_BFM
  for(int i=0;i<eig.size();i++){
    eig[i].deleteEig();
    eig[i].reloadGauge(); //will re-import gauge into internal bfm object when calc is next called
  }
#endif
}

void PropManager::clear(){
  props.clear();
#ifdef USE_BFM
  eig.clear();
#endif
}

void PropManager::calcProps(Lattice &latt){
#ifdef USE_BFM
  for(int i=0;i<eig.size();i++) eig[i].calcEig(latt);
#endif

  for(int i=0;i<props.size();i++){
    props[i].readProp(latt);
    props[i].calcProp(latt); //will calculate if not read
  }
}

#ifdef USE_BFM
LanczosContainer& PropManager::getLanczos(const char *tag){
  for(int i=0;i<eig.size();i++) if(eig[i].tagEquals(tag)) return eig[i];
  ERR.General("PropManager","getLanczos(const char *tag)","Lanczos instance '%s' does not exist!\n",tag);
}
#endif

CPS_END_NAMESPACE


