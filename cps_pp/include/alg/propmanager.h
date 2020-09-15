CPS_START_NAMESPACE
#ifndef PROP_MANAGER_H
#define PROP_MANAGER_H
CPS_END_NAMESPACE

#include <config.h>
#include <alg/propvector.h>
#include <alg/prop_attribute_arg.h>

#ifdef USE_BFM
namespace BFM_Krylov{
template<class T> class Lanczos_5d;
}
#endif

CPS_START_NAMESPACE

//static class for managing propagators
class PropManager{
private:
  static PropVector props;
#ifdef USE_BFM
  static LanczosVector eig; //generate/store sets of eigenvectors, recallable by tag for use in A2A props, AMA etc
#endif
public:
  static PropagatorContainer & getProp(const char *tag); //get propagator indexed by tag
  static PropagatorContainer & addProp(PropagatorArg &arg); //this does not actually invert the propagator until getProp is called on the PropagatorContainer or using the above function on the vector

#ifdef USE_BFM
  static LanczosContainer & getLanczos(const char *tag);
  static LanczosContainer & addLanczos(LanczosContainerArg &arg);
#endif

  static void setup(JobPropagatorArgs &prop_args); //set up props from an args struct
  static void startNewTraj(); //deletes all existing props, but keeps the containers so they can be recreated on this new configuration
  static void calcProps(Lattice &latt); //read/invert all propagators
  static void clear(); //deletes all props and containers
};

#endif
CPS_END_NAMESPACE
