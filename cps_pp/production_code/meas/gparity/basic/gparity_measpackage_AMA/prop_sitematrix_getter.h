#ifndef _AMA_PROP_SITEMATRIX_GETTER_H
#define _AMA_PROP_SITEMATRIX_GETTER_H

#include "propwrapper.h"
#include "propmomcontainer.h"
#include "prop_tag.h"

CPS_START_NAMESPACE

//A class to manage the getting of the propagator for a particular local 3d coordinate and a *global* source->sink separation tdis_glb.  tdis_glb is linear, i.e. the user doesn't need to worry about the modulo Lt wrapping
//The class knows about the periodicity of the propagator
//Note, the resulting site must be available on this node - no comms are performed. This is not a huge restriction as the shifts we perform are in units of Lt
class PropSiteMatrixGetter{
public:
  virtual void siteMatrix(WilsonMatrix &into, const int x3d_lcl, const int tdis_glb, const PropSplane splane = SPLANE_BOUNDARY) const = 0;
  virtual void siteMatrix(SpinColorFlavorMatrix &into, const int x3d_lcl, const int tdis_glb, const PropSplane splane = SPLANE_BOUNDARY) const = 0;
  virtual ~PropSiteMatrixGetter(){}
};

//For propagators with modulo-Lt periodicity or antiperiodicity
class PropSiteMatrixStandard: public PropSiteMatrixGetter{
  PropWrapper prop;
  BndCndType bc;
  int tsrc;

  inline void get4dcoordAndSign(int &x4d_lcl, int &sgn, const int x3d_lcl, const int tdis_glb) const{
    int Lt = GJP.Tnodes()*GJP.TnodeSites();
    sgn = 1;
    int t = tdis_glb + tsrc;
    while(t<0 || t>=Lt){
      if(t<0) t+= Lt;
      else t-=Lt;
      sgn *= (BND_CND_APRD ? -1 : 1); //APRD   G(t-Lt) = -G(t)
    }
    int t_lcl = t - GJP.TnodeCoor() * GJP.TnodeSites();
    assert(t_lcl>=0 && t_lcl<GJP.TnodeSites());
    x4d_lcl = x3d_lcl + GJP.VolNodeSites()/GJP.TnodeSites()*t_lcl;
  }


public:
  PropSiteMatrixStandard(const PropWrapper &_prop, const BndCndType _bc, const int _tsrc): prop(_prop), bc(_bc), tsrc(_tsrc){}

  void siteMatrix(WilsonMatrix &into, const int x3d_lcl, const int tdis_glb, const PropSplane splane = SPLANE_BOUNDARY) const{
    int sgn, x4d_lcl;
    get4dcoordAndSign(x4d_lcl,sgn,x3d_lcl,tdis_glb);
    prop.siteMatrix(into,x4d_lcl,splane);
    if(!sgn) into *= Float(-1);
  }
  void siteMatrix(SpinColorFlavorMatrix &into, const int x3d_lcl, const int tdis_glb, const PropSplane splane = SPLANE_BOUNDARY) const{
    int sgn, x4d_lcl;
    get4dcoordAndSign(x4d_lcl,sgn,x3d_lcl,tdis_glb);
    prop.siteMatrix(into,x4d_lcl,splane);
    if(!sgn) into *= Float(-1);
  }    
};


//For F=P+A and B=P-A propagators with modulo-2Lt periodicity
//Utilize F(t+Lt) = B(t) and B(t+Lt) = F(t)
//Need both the F and B propagators. The 'base' propagator is either F or B, and the 'shift' propagator is B or F respectively
class PropSiteMatrixFB: public PropSiteMatrixGetter{
  PropWrapper prop_base;
  PropWrapper prop_shift;
  int tsrc;

  inline PropWrapper const* get4dcoordAndProp(int &x4d_lcl, const int x3d_lcl, const int tdis_glb) const{
    int Lt = GJP.Tnodes()*GJP.TnodeSites();
    int t = tdis_glb + tsrc;
    PropWrapper const* props[2] = { &prop_base, &prop_shift };
    int use = 0;
    while(t<0 || t>=Lt){
      if(t<0) t+= Lt;
      else t-=Lt;
      use = (use + 1) % 2;
    }
    int t_lcl = t - GJP.TnodeCoor() * GJP.TnodeSites();
    if(t_lcl < 0 || t_lcl >= GJP.TnodeSites()){
      printf("PropSiteMatrixFB::get4dcoordAndProp for Node %d with t range %d to %d, tdis_glb=%d and tsrc=%d, got lattice t=%d and prop %d. Lattice t is not on node!\n",
	     UniqueID(),
	     GJP.TnodeCoor()*GJP.TnodeSites(), (GJP.TnodeCoor()+1)*GJP.TnodeSites()-1,
	     tdis_glb,tsrc,t,use); 
      fflush(stdout);
      exit(-1);
    }
    x4d_lcl = x3d_lcl + GJP.VolNodeSites()/GJP.TnodeSites()*t_lcl;
    // if(!UniqueID()){ printf("PropSiteMatrixFB::get4dcoordAndProp for Node %d with t range %d to %d, tdis_glb=%d and tsrc=%d, got lattice t=%d and prop %d.\n",
    // 			    UniqueID(),
    // 			    GJP.TnodeCoor()*GJP.TnodeSites(), (GJP.TnodeCoor()+1)*GJP.TnodeSites()-1,
    // 			    tdis_glb,tsrc,t,use); fflush(stdout); }
    return props[use];
  }


public:
  PropSiteMatrixFB(const PropWrapper &_prop_base, const PropWrapper &_prop_shift, const int _tsrc): prop_base(_prop_base),prop_shift(_prop_shift),tsrc(_tsrc){}

  void siteMatrix(WilsonMatrix &into, const int x3d_lcl, const int tdis_glb, const PropSplane splane = SPLANE_BOUNDARY) const{
    int x4d_lcl;
    PropWrapper const* prop_use = get4dcoordAndProp(x4d_lcl,x3d_lcl,tdis_glb);
    prop_use->siteMatrix(into,x4d_lcl,splane);
  }
  void siteMatrix(SpinColorFlavorMatrix &into, const int x3d_lcl, const int tdis_glb, const PropSplane splane = SPLANE_BOUNDARY) const{
    int x4d_lcl;
    PropWrapper const* prop_use = get4dcoordAndProp(x4d_lcl,x3d_lcl,tdis_glb);
    prop_use->siteMatrix(into,x4d_lcl,splane);
  }    
};


class PropSiteMatrixGetterFactory{
public:
  static PropSiteMatrixGetter* get(const PropMomContainer &props, const QuarkType lh, const PropPrecision prec, const int tsrc, const ThreeMomentum &p, const TbcStatus &tbc){
    std::string tag = propTag(lh,prec,tsrc,p,tbc);
    if(!tbc.isCombinedType()){
      return new PropSiteMatrixStandard(props.get(tag), tbc.getSingleBc(), tsrc);
    }else{
      TbcCombination base_comb = tbc.getCombinedBc(); 
      assert(base_comb == CombinationF || base_comb == CombinationB);
      TbcCombination shift_comb = base_comb == CombinationF ? CombinationB : CombinationF;
      
      std::string shift_tag = propTag(lh,prec,tsrc,p,shift_comb);
      
      const PropWrapper &base_prop = props.get(tag);
      const PropWrapper &shift_prop = props.get(shift_tag);
      return new PropSiteMatrixFB(base_prop,shift_prop,tsrc);
    }
  }
};





template<typename MatrixType>
class WallSinkPropSiteMatrixGetter{
public:
  virtual void siteMatrix(MatrixType &into, const int tdis_glb) const = 0;
  virtual ~WallSinkPropSiteMatrixGetter(){}
};

template<typename MatrixType>
class WallSinkPropSiteMatrixStandard: public WallSinkPropSiteMatrixGetter<MatrixType>{
  WallSinkProp<MatrixType> prop;
  BndCndType bc;
  int tsrc;

  inline void getSinkTimeAndSign(int &t, int &sgn, const int tdis_glb) const{
    int Lt = GJP.Tnodes()*GJP.TnodeSites();
    sgn = 1;
    t = tdis_glb + tsrc;
    while(t<0 || t>=Lt){
      if(t<0) t+= Lt;
      else t-=Lt;
      sgn *= (BND_CND_APRD ? -1 : 1); //APRD   G(t-Lt) = -G(t)
    }
  }

public:
  WallSinkPropSiteMatrixStandard(const PropWrapper &prop_in, const BndCndType _bc, const int _tsrc, const ThreeMomentum &p_snk, Lattice &lat, const bool gauge_fix_sink = true): bc(_bc), tsrc(_tsrc), prop(gauge_fix_sink){
    prop.setProp(prop_in);
    prop.compute(lat,p_snk);
  }
  WallSinkPropSiteMatrixStandard(const PropWrapper &prop_in, const BndCndType _bc, const int _tsrc, const double *p_snk, Lattice &lat, const bool gauge_fix_sink = true): bc(_bc), tsrc(_tsrc), prop(gauge_fix_sink){
    prop.setProp(prop_in);
    prop.compute(lat,p_snk);
  }

  void siteMatrix(MatrixType &into, const int tdis_glb) const{
    int sgn, t_glb;
    getSinkTimeAndSign(t_glb,sgn,tdis_glb);
    into = prop(t_glb);
    if(!sgn) into *= Float(-1);
  }

  WallSinkProp<MatrixType> & getWallSinkProp(){ return prop; }
};

template<typename MatrixType>
class WallSinkPropSiteMatrixFB: public WallSinkPropSiteMatrixGetter<MatrixType>{
  WallSinkProp<MatrixType> prop_base;
  WallSinkProp<MatrixType> prop_shift;
  int tsrc;

  inline WallSinkProp<MatrixType> const* getSinkTimeAndProp(int &t, const int tdis_glb) const{
    int Lt = GJP.Tnodes()*GJP.TnodeSites();
    t = tdis_glb + tsrc;
    WallSinkProp<MatrixType> const* props[2] = { &prop_base, &prop_shift };
    int use = 0;
    while(t<0 || t>=Lt){
      if(t<0) t+= Lt;
      else t-=Lt;
      use = (use + 1) % 2;
    }
    return props[use];
  }

public:
  WallSinkPropSiteMatrixFB(const PropWrapper &prop_base_in, const PropWrapper &prop_shift_in, const int _tsrc, 
			   const ThreeMomentum &p_snk, Lattice &lat, const bool gauge_fix_sink = true): tsrc(_tsrc),prop_base(gauge_fix_sink), prop_shift(gauge_fix_sink){
    prop_base.setProp(prop_base_in);
    prop_base.compute(lat,p_snk);

    prop_shift.setProp(prop_shift_in);
    prop_shift.compute(lat,p_snk);
  }


  void siteMatrix(MatrixType &into, const int tdis_glb) const{
    int t_glb;
    WallSinkProp<MatrixType> const* prop_use = getSinkTimeAndProp(t_glb,tdis_glb);
    into = prop_use->operator()(t_glb);
  }
};


template<typename MatrixType>
class WallSinkPropSiteMatrixGetterFactory{
public:
  static WallSinkPropSiteMatrixGetter<MatrixType>* get(const PropMomContainer &props, const QuarkType lh, const PropPrecision prec, const int tsrc, 
						       const ThreeMomentum &p_src, const ThreeMomentum &p_snk, const TbcStatus &tbc, Lattice &lat, const bool gauge_fix_sink = true){
    std::string tag = propTag(lh,prec,tsrc,p_src,tbc);
    if(!tbc.isCombinedType()){
      return new WallSinkPropSiteMatrixStandard<MatrixType>(props.get(tag), tbc.getSingleBc(), tsrc, p_snk, lat, gauge_fix_sink);
    }else{
      TbcCombination base_comb = tbc.getCombinedBc(); 
      assert(base_comb == CombinationF || base_comb == CombinationB);
      TbcCombination shift_comb = base_comb == CombinationF ? CombinationB : CombinationF;
      
      std::string shift_tag = propTag(lh,prec,tsrc,p_src,shift_comb);
      
      const PropWrapper &base_prop = props.get(tag);
      const PropWrapper &shift_prop = props.get(shift_tag);
      return new WallSinkPropSiteMatrixFB<MatrixType>(base_prop,shift_prop,tsrc, p_snk, lat, gauge_fix_sink);
    }
  }
};



CPS_END_NAMESPACE

#endif
