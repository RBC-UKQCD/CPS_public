#ifndef _MESON_MOMENTA_H
#define _MESON_MOMENTA_H

#include <alg/a2a/threemomentum.h>

CPS_START_NAMESPACE

//Storage for unique quark momenta needed (for pairs +p and -p, we can trivially reconstruct the flipped momentum using the G-parity conjugate relation without having to do extra inversions)
class QuarkMomenta{
  std::vector<ThreeMomentum> mom;
public:
  inline void add(const ThreeMomentum &val){
    //if(!UniqueID()) std::cout << "QuarkMomenta adding mom " << val.str() << '\n';
    for(int i=0;i<mom.size();i++) 
      if(mom[i] == val || (GJP.Gparity() && mom[i] == -val)){
	//if(!UniqueID()) std::cout << "Skipping as it (or its negative) already exists: " << mom[i].str() << '\n';
	return;
      }

    mom.push_back(val);
  }
  inline int nMom() const{ return mom.size(); }
  inline const ThreeMomentum & getMom(const int pidx) const{ return mom[pidx]; }
};

//Two-point functions of the form 
//< \bar\psi_2(-p_psi,t) g5 A g5 \psi_1(-p_psibar,t) \bar\psi_1(p_psibar,0) g5 B g5 \psi_2(p_psi,0) >
//where momenta are associated with the quark field operators in the *source*

//These form correlators
//< [prop2(t,0;-p_psi)]^dag A prop1(t,0;p_psibar) B>

//There are 5 momenta values we might be interested in: cf MomentumOf enum

//DaggeredProp is the momentum of the propagator that will be daggered (prop2 above)

class MesonMomenta{
  typedef std::pair<QuarkType,ThreeMomentum> Mtype;
  std::vector<Mtype> p_psi;
  std::vector<Mtype> p_psibar;
  
public:
  inline int nMom() const{ return p_psi.size(); }
  
  inline ThreeMomentum getMomentum(const MomentumOf what, const int mom_idx) const{
    switch(what){
    case SrcPsi:
      return p_psi[mom_idx].second;
    case DaggeredProp:
      return -p_psi[mom_idx].second; //we require a propagator of the opposite momentum, to which we apply g5-hermiticity
    case SrcPsiBar:
    case UndaggeredProp:
      return p_psibar[mom_idx].second;
      break;
    case Total: 
      return p_psi[mom_idx].second + p_psibar[mom_idx].second;
    default:
      ERR.General("MesonMomenta","getMomentum","Unknown momentum label\n");
      break;
    }
  }

  inline QuarkType getQuarkType(const MomentumOf what, const int mom_idx) const{
    switch(what){
    case SrcPsi:
    case DaggeredProp:
      return p_psi[mom_idx].first;
      break;
    case SrcPsiBar:
    case UndaggeredProp:
      return p_psibar[mom_idx].first;
      break;
    default:
      ERR.General("MesonMomenta","getMomentum","Invalid choice\n");
      break;
    }
  }

  void printAllCombs(const std::string & descr = "") const{
    if(!UniqueID()){
      printf("(psibar,psi) Momentum combinations %s:\n",descr.c_str());
      for(int i=0;i<nMom();i++){
	std::cout << "(" << (p_psibar[i].first == Light ? "light" : "heavy") << ") " << p_psibar[i].second.str() << " + ";
	std::cout << "(" << (p_psi[i].first == Light ? "light" : "heavy") << ") " << p_psi[i].second.str() << '\n';
      }
    }
  }
    

  //Add a p(\bar\psi) and p(\psi) from a string in the form "(%d,%d,%d) + (%d%d,%d)"
  void addP(const std::string &p, const QuarkType qtype1,  const QuarkType qtype2){
    std::pair<ThreeMomentum,ThreeMomentum> p2 = ThreeMomentum::parse_str_two_mom(p);
    p_psibar.push_back(Mtype(qtype1,p2.first ));
    p_psi   .push_back(Mtype(qtype2,p2.second));
  }

  //Add to QuarkMomenta all the required *propagator* source momenta of a particular quark species (heavy/light)
  //Assumes gamma5-hermiticity is used for the backwards propagator (cf below)
  void appendQuarkMomenta(const QuarkType qtype,QuarkMomenta &qmom) const{
    //Psibar is associated with the undaggered prop
    //For the propagator   \sum_{x_src} e^{-ip_psibar x_src}<\psi(x_snk,t)\bar\psi(x_src,0)> 
    //the prop momentum is the same as the psibar momentum
    for(int i=0;i<p_psibar.size();i++)   
      if(p_psibar[i].first == qtype) qmom.add(p_psibar[i].second);

    //Psi is associated with the daggered prop. 
    //To get      \sum_{x_src} e^{-ip_psi x_src}<\psi(x_src,0)\bar\psi(x_snk,t)> 
    //we compute  \sum_{x_src} \gamma^5 [ e^{-i [-p_psi] x_snk}<\psi(x_snk,0)\bar\psi(x_src,0)> ]^\dagger \gamma^5
    //Where the prop computed has the *opposite* momentum, -p_psi
    for(int i=0;i<p_psi.size();i++) 
      if(p_psi[i].first == qtype) qmom.add(-p_psi[i].second);
  }

};



inline int getNgp(){
  int ngp = 0;
  for(int i=0;i<3;i++){
    if(GJP.Bc(i) == BND_CND_GPARITY) ngp++;
    if(i>0 && GJP.Bc(i) == BND_CND_GPARITY && GJP.Bc(i-1) != BND_CND_GPARITY){ ERR.General("","getNgp","Expect G-parity directions to be consecutive\n"); } //(as it is setup here we expect the G-parity directions to be consecutive, i.e.  x, or x and y, or x and y and z)
  }
  return ngp;
}



class PionMomenta{
public:
  static void setup(MesonMomenta &into, bool include_alternative_mom = true){
    int ngp = getNgp();
    
    if(ngp == 0){
      //p_pi = (0,0,0)
      into.addP("(0,0,0) + (0,0,0)",Light,Light);
    }else if(ngp == 1){
      //p_pi = (2,0,0)     (units of pi/2L)    
      into.addP("(1,0,0) + (1,0,0)",Light,Light); 
      if(include_alternative_mom) into.addP("(-1,0,0) + (3,0,0)",Light,Light); 
    }else if(ngp == 2){
      //Along G-parity direction:
      //p_pi = (2,2,0)     (units of pi/2L)  
      into.addP("(1,1,0) + (1,1,0)",Light,Light);
      if(include_alternative_mom) into.addP("(3,3,0) + (-1,-1,0)",Light,Light);

      //Along off-diagonal direction:      
      //p_pi = (-2,2,0)
      into.addP("(-3,1,0) + (1,1,0)",Light,Light); 
      if(include_alternative_mom) into.addP("(-1,3,0) + (-1,-1,0)",Light,Light);
    }else if(ngp == 3){
      //p_pi = (2,2,2)     (units of pi/2L)  
      into.addP("(1,1,1) + (1,1,1)",Light,Light); 
      if(include_alternative_mom) into.addP("(3,3,3) + (-1,-1,-1)",Light,Light);

      //p_pi = (-2,2,2) //only do one off-diagonal as we have symmetry around that axis
      into.addP("(-3,1,1) + (1,1,1)",Light,Light);
      if(include_alternative_mom) into.addP("(-1,3,3) + (-1,-1,-1)",Light,Light);
    }else{
      ERR.General("PionMomenta","setup","ngp cannot be >3\n");
    }
  }
};


//This is the stationary \bar\psi \gamma^5 \psi pseudoscalar flavor-singlet 
class LightFlavorSingletMomenta{
public:
  static void setup(MesonMomenta &into){
    int ngp = getNgp();
    
    if(ngp == 0){
      //p_pi = (0,0,0)
      into.addP("(0,0,0) + (0,0,0)",Light,Light);
    }else if(ngp == 1){
      //p_pi = (0,0,0)     (units of pi/2L)    
      into.addP("(1,0,0) + (-1,0,0)",Light,Light); 
    }else if(ngp == 2){
      //Along G-parity direction:
      //p_pi = (0,0,0)     (units of pi/2L)  
      into.addP("(1,1,0) + (-1,-1,0)",Light,Light);
    }else if(ngp == 3){
      //p_pi = (0,0,0)     (units of pi/2L)  
      into.addP("(1,1,1) + (-1,-1,-1)",Light,Light); 
    }else{
      ERR.General("LightFlavorSingletMomenta","setup","ngp cannot be >3\n");
    }
  }
};

//Note the heavy quark propagator (the one that is daggered - cf Eq 191 of the paper) is assigned the - momentum
class KaonMomenta{
public:
  static void setup(MesonMomenta &into){
    int ngp = getNgp();
    
    if(ngp == 0){
      //p_pi = (0,0,0)
      into.addP("(0,0,0) + (0,0,0)",Light,Heavy);
    }else if(ngp == 1){
      //p_pi = (0,0,0)     (units of pi/2L)    
      into.addP("(1,0,0) + (-1,0,0)",Light,Heavy); 
    }else if(ngp == 2){
      //Along G-parity direction:
      //p_pi = (0,0,0)     (units of pi/2L)  
      into.addP("(1,1,0) + (-1,-1,0)",Light,Heavy);
    }else if(ngp == 3){
      //p_pi = (0,0,0)     (units of pi/2L)  
      into.addP("(1,1,1) + (-1,-1,-1)",Light,Heavy); 
    }else{
      ERR.General("KaonMomenta","setup","ngp cannot be >3\n");
    }
  }
};

CPS_END_NAMESPACE

#endif
