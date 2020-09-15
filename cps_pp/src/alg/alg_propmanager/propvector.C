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

#include <alg/qpropw.h>
#include <alg/propvector.h>

CPS_START_NAMESPACE

PropVector::PropVector(): PointerArray<PropagatorContainer>(){}

PropagatorContainer & PropVector::addProp(PropagatorArg &arg){
  PropagatorContainer* p = PropagatorContainer::create(arg.generics.type);
  p->setup(arg);
  return append(p);
}

#ifdef USE_BFM
LanczosVector::LanczosVector(): PointerArray<LanczosContainer>(){}

LanczosContainer & LanczosVector::add(LanczosContainerArg &arg){
  return append(new LanczosContainer(arg));
}
#endif

CPS_END_NAMESPACE
