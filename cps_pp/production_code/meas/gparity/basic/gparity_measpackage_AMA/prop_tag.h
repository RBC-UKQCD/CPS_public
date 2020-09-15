#ifndef _AMA_PROP_TAG_H
#define _AMA_PROP_TAG_H

#include "enums.h"

CPS_START_NAMESPACE

class TbcStatus{
  BndCndType basic_bc;
  TbcCombination bc_comb;
public:
  TbcStatus(const BndCndType bbc): basic_bc(bbc), bc_comb(Single){}
  TbcStatus(const TbcCombination bcomb): bc_comb(bcomb){}
  
  inline std::string getTag() const{
    if(bc_comb == Single){
      switch(basic_bc){
      case BND_CND_PRD:
	return std::string("P");
	break;
      case BND_CND_APRD:
	return std::string("A");
	break;
      default:
	ERR.General("TbcStatus","getTag","Unknown TBC\n");
	break;
      }
    }else return std::string(bc_comb == CombinationF ? "F" : "B");
  }
  inline bool isCombinedType() const{
    return bc_comb != Single;
  }
  inline BndCndType getSingleBc() const{
    assert(!isCombinedType());
    return basic_bc;
  }
  inline TbcCombination getCombinedBc() const{
    assert(isCombinedType());
    return bc_comb;
  }
  void swapTbcCombination(){
    assert(isCombinedType());
    TbcCombination bc_comb_other = bc_comb == CombinationF ? CombinationB : CombinationF;
    bc_comb = bc_comb_other;
  }


};

inline static std::string propTag(const QuarkType lh, const PropPrecision prec, const int tsrc, const ThreeMomentum &p, const TbcStatus &tbc){
  std::ostringstream tag;
  tag << (lh == Light ? "light" : "heavy");
  tag << (prec == Sloppy ? "_sloppy" : "_exact");
  tag << "_t" << tsrc;
  tag << "_mom" << p.file_str();
  tag << "_tbc_" << tbc.getTag();
  return tag.str();
}
inline static std::string propTag(const QuarkType lh, const PropPrecision prec, const int tsrc, const ThreeMomentum &p, const BndCndType tbc){
  return propTag(lh,prec,tsrc,p,TbcStatus(tbc));
}
inline static std::string propTag(const QuarkType lh, const PropPrecision prec, const int tsrc, const ThreeMomentum &p, const TbcCombination bc_comb){
  return propTag(lh,prec,tsrc,p,TbcStatus(bc_comb));
}

CPS_END_NAMESPACE

#endif
