#include "alg/prop_attribute_arg.h"
#include <string.h>
CPS_START_NAMESPACE

PropagatorArg::~PropagatorArg(){
  //if(attributes.attributes_val) delete attributes.attributes_val;
}
  
JobPropagatorArgs::~JobPropagatorArgs(){
  //if(props.props_val) delete props.props_val;
}
JobPropagatorArgs::JobPropagatorArgs(){
  props.props_len = 0; 
  props.props_val = NULL;
  lanczos.lanczos_len = 0; 
  lanczos.lanczos_val = NULL;
}

//type mapping
AttrType GenericPropAttrArg::getType (  ){ return GENERIC_PROP_ATTR; }
AttrType PointSourceAttrArg::getType (  ){ return POINT_SOURCE_ATTR; }
AttrType WallSourceAttrArg::getType (  ){ return WALL_SOURCE_ATTR; }
AttrType VolumeSourceAttrArg::getType (  ){ return VOLUME_SOURCE_ATTR; }
AttrType MomentumAttrArg::getType (  ){ return MOMENTUM_ATTR; }
AttrType PropIOAttrArg::getType (  ){ return PROP_IO_ATTR; }
AttrType GparityFlavorAttrArg::getType (  ){ return GPARITY_FLAVOR_ATTR; }
AttrType CGAttrArg::getType (  ){ return CG_ATTR; }
AttrType GaugeFixAttrArg::getType (  ){ return GAUGE_FIX_ATTR; }
AttrType MomCosAttrArg::getType (  ){ return MOM_COS_ATTR; }
AttrType PropCombinationAttrArg::getType (  ){ return PROP_COMBINATION_ATTR; }
AttrType GparityOtherFlavPropAttrArg::getType ( ){ return GPARITY_OTHER_FLAV_PROP_ATTR; }
AttrType GparityComplexConjSourcePartnerPropAttrArg::getType ( ){ return GPARITY_COMPLEX_CONJ_SOURCE_PARTNER_PROP_ATTR; }
AttrType TwistedBcAttrArg::getType( ){ return TWISTED_BC_ATTR; }
AttrType StoreMidpropAttrArg::getType( ){ return STORE_MIDPROP_ATTR; }
AttrType A2AAttrArg::getType( ){ return A2A_ATTR; }
AttrType DeflatedCGAttrArg::getType( ){ return DEFLATED_CG_ATTR; }
//As AttrArg types are used in a union, we cannot write copy constructors. Instead write a clone function to perform a deep copy,
//and the automatically generated trivial copy constructor is shallow.

template<typename StructType>
class clone_lambda{
public:
  StructType *to;

  void clone(StructType &into, StructType &from){
    to = &into;
    from.datamem_lambda(*this);
  }
  
  template<typename T>
  void doit(T StructType::* memptr, StructType &obj);
  template<typename T,int size>
  void doit_vector(T (StructType::* memptr)[size], StructType &obj);

  // template<typename T>
  // void clone_it(T &
};

template<typename StructType>
template<typename T>
void clone_lambda<StructType>::doit(T StructType::* memptr, StructType &obj){
  (*to).*memptr = obj.*memptr;
}
template<typename StructType>
template<typename T, int size>
void clone_lambda<StructType>::doit_vector(T (StructType::* memptr)[size], StructType &obj){
  for(int s=0;s<size;s++)
    ((*to).*memptr)[s] = (obj.*memptr)[s];
}



GenericPropAttrArg GenericPropAttrArg::clone(){
  GenericPropAttrArg out;
  //clone_lambda<GenericPropAttrArg> cl; cl.clone(out,*this);

  out.tag = new char[strlen(tag)+1]; strcpy(out.tag,tag);
  out.mass = mass;
  memcpy(out.bc,bc,4*sizeof(BndCndType));
  return out;
}
PointSourceAttrArg PointSourceAttrArg::clone(){
  PointSourceAttrArg out;
  memcpy(out.pos,pos,4*sizeof(int));
  return out;
}
WallSourceAttrArg WallSourceAttrArg::clone(){
  return *this;
}
VolumeSourceAttrArg VolumeSourceAttrArg::clone(){
  return *this;
}
MomentumAttrArg MomentumAttrArg::clone(){
  MomentumAttrArg out;
  memcpy(out.p,p,3*sizeof(int));
  return out;
}
PropIOAttrArg PropIOAttrArg::clone(){
  PropIOAttrArg out(*this);
  out.qio_filename = new char[strlen(qio_filename)+1]; strcpy(out.qio_filename,qio_filename);
  return out;
}
GparityFlavorAttrArg GparityFlavorAttrArg::clone(){
  return *this;
}
CGAttrArg CGAttrArg::clone(){
  return *this;
}
GaugeFixAttrArg GaugeFixAttrArg::clone(){
  return *this;
}
MomCosAttrArg MomCosAttrArg::clone(){
  return *this;
}
PropCombinationAttrArg PropCombinationAttrArg::clone(){
  PropCombinationAttrArg out;
  out.deep_copy(*this);
  return out;
}
GparityOtherFlavPropAttrArg GparityOtherFlavPropAttrArg::clone(){
  GparityOtherFlavPropAttrArg out;
  out.deep_copy(*this);
  return out;
}
TwistedBcAttrArg TwistedBcAttrArg::clone(){
  TwistedBcAttrArg out;
  out.deep_copy(*this);
  return out;
}
StoreMidpropAttrArg StoreMidpropAttrArg::clone(){
  return *this;
}
A2AAttrArg A2AAttrArg::clone(){
  A2AAttrArg out;
  out.deep_copy(*this);
  return out;
}

DeflatedCGAttrArg DeflatedCGAttrArg::clone(){
  DeflatedCGAttrArg out;
  out.deep_copy(*this);
  return out;
}


CPS_END_NAMESPACE
