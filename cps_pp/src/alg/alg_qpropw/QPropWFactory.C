#include <config.h>
#include <alg/qpropw.h>

CPS_START_NAMESPACE

QPropW* QPropWFactory::Create( 
  Lattice &lat, SourceType type, 
  QPropWArg *arg, CommonArg *c_arg, void *e_arg
){
  QPropW *qpropw=NULL;
  switch (type) {
    case WALL:
      qpropw = new QPropWWallSrc(lat, arg, c_arg);
      break;
    case MOM:
      qpropw = new QPropWMomSrc(lat, arg, (int *)e_arg, c_arg);
      break;
    case VOLUME:
      qpropw = new QPropWVolSrc(lat, arg, c_arg);
      break;
    case POINT:
      qpropw = new QPropWPointSrc(lat, arg, c_arg);
      break;
    case RANDWALL:
      qpropw = new QPropWRandWallSrc(lat, arg, (QPropWRandArg *)e_arg, c_arg);
      break;
    case RANDVOLUME:
      qpropw = new QPropWRandVolSrc(lat, arg, (QPropWRandArg *)e_arg, c_arg);
      break;
    case RANDSLAB:
      qpropw = new QPropWRandSlabSrc(lat, arg, (QPropWSlabArg *)e_arg, c_arg);
      break;
    case EXP:
      qpropw = new QPropWExpSrc(lat, arg, (QPropWExpArg *)e_arg, c_arg);
      break;
    default:
      ERR.General("QPropWFactory","Create","SourceType %d not recognized\n",type);
      break;
  }
  return qpropw;
}

void QPropWFactory::Destroy(QPropW *qp){
  delete qp;
}

CPS_END_NAMESPACE
