#ifndef _REQUIRED_MOMENTA_H
#define _REQUIRED_MOMENTA_H

#include<alg/a2a/threemomentum.h>

CPS_START_NAMESPACE

//A class to compute and store the W^dagger V meson field momenta
template<typename MomComputePolicy>
class RequiredMomentum: public MomComputePolicy{
  std::vector<ThreeMomentum>  wdag_mom;
  std::vector<ThreeMomentum>  vmom;
  std::vector<ThreeMomentum>  wdag_mom_alt;
  std::vector<ThreeMomentum>  vmom_alt;
     
public:
  RequiredMomentum(){
    //We require meson fields with both p and -p for each non-identical momentum direction.     
    int ngp = 0; for(int i=0;i<3;i++){
      if(GJP.Bc(i) == BND_CND_GPARITY) ngp++;
      if(i>0 && GJP.Bc(i) == BND_CND_GPARITY && GJP.Bc(i-1) != BND_CND_GPARITY){ ERR.General("RequiredMomentum","RequiredMomentum","Expect G-parity directions to be consecutive\n"); }     //(as it is setup here we expect the G-parity directions to be consecutive, i.e.  x, or x and y, or x and y and z)
    }
    this->setupMomenta(ngp); //in MomComputePolicy
  }

  int nMom() const{ return vmom.size(); }
  int nAltMom() const{ return vmom_alt.size(); }
  
  //Get the twist momentum for the W field. Use alternative = true to access the second G-parity momentum averaged with the first to improve rotational symmetry
  ThreeMomentum getWmom(const int &i, const bool &alternative = false) const{ return ThreeMomentum::negative(alternative ? wdag_mom_alt[i] : wdag_mom[i]); } //negative because p(W) = -p(W^dag)
  //Get the twist momentum for the V field. Use alternative = true to access the second G-parity momentum averaged with the first to improve rotational symmetry
  ThreeMomentum getVmom(const int &i, const bool &alternative = false) const{ return alternative ? vmom_alt[i] : vmom[i]; }

  ThreeMomentum getMesonMomentum(const int &i) const{ return wdag_mom[i] + vmom[i]; }

  //Add a W^dag and V momentum (respectively) from a string in the form "(%d,%d,%d) + (%d%d,%d)", plus the negatives of these
  //Use alternative = true to specify second G-parity momentum averaged with the first to improve rotational symmetry
  void addPandMinusP(const std::string &p, const bool &alternative = false){
    std::pair<ThreeMomentum,ThreeMomentum> p2 = ThreeMomentum::parse_str_two_mom(p);
    std::vector<ThreeMomentum> &wto = alternative ? wdag_mom_alt : wdag_mom;
    std::vector<ThreeMomentum> &vto = alternative ? vmom_alt : vmom;

    wto.push_back(p2.first); vto.push_back(p2.second);
    wto.push_back(-p2.first); vto.push_back(-p2.second);

    if(alternative){
      if(vmom_alt.size() != vmom.size()){ ERR.General("RequiredMomentum","addPandMinusP","Alternative momentum combination must be added after standard combination\n"); }
      for(int i=vmom.size()-2; i<vmom.size(); i++){
	ThreeMomentum ptot_alt = vmom_alt[i] + wdag_mom_alt[i];
	ThreeMomentum ptot_orig = vmom[i] + wdag_mom[i];
	if(ptot_alt != ptot_orig){ ERR.General("RequiredMomentum","addPandMinusP","Alternative momentum combination must have same total momentum as standard combination\n"); } 
      }
    }
  }
  //Add a W^dag and V momentum (respectively) from a string in the form "(%d,%d,%d) + (%d%d,%d)"
  //Use alternative = true to specify second G-parity momentum averaged with the first to improve rotational symmetry
  void addP(const std::string &p, const bool &alternative = false){
    std::pair<ThreeMomentum,ThreeMomentum> p2 = ThreeMomentum::parse_str_two_mom(p);
    std::vector<ThreeMomentum> &wto = alternative ? wdag_mom_alt : wdag_mom;
    std::vector<ThreeMomentum> &vto = alternative ? vmom_alt : vmom;

    wto.push_back(p2.first); vto.push_back(p2.second);
    if(alternative){
      if(vmom_alt.size() != vmom.size()){ ERR.General("RequiredMomentum","addPandMinusP","Alternative momentum combination must be added after standard combination\n"); }
      int i = vmom.size()-1;
      ThreeMomentum ptot_alt = vmom_alt[i] + wdag_mom_alt[i];
      ThreeMomentum ptot_orig = vmom[i] + wdag_mom[i];
      if(ptot_alt != ptot_orig){ ERR.General("RequiredMomentum","addP","Alternative momentum combination must have same total momentum as standard combination\n"); } 
    }

  }
};



CPS_END_NAMESPACE

#endif
