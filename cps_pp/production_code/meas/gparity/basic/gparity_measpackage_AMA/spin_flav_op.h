#ifndef _SPIN_FLAV_H
#define _SPIN_FLAV_H

CPS_START_NAMESPACE

//Assumes momenta are in units of \pi/2L, and must be *odd integer* (checked)
inline static int getProjSign(const int p[3]){
  if(!GJP.Gparity()){ ERR.General("","getProjSign","Requires GPBC in at least one direction\n"); }

  //Sign is exp(i\pi n_p)
  //where n_p is the solution to  p_j = \pi/2L( 1 + 2n_p )
  //Must be consistent for all j

  int np;
  for(int j=0;j<3;j++){
    if(GJP.Bc(j)!=BND_CND_GPARITY) continue;

    if(abs(p[j]) %2 != 1){ ERR.General("","getProjSign","Component %d of G-parity momentum (%d,%d,%d) is invalid as it is not an odd integer!\n",j,p[0],p[1],p[2]); }
    int npj = (p[j] - 1)/2;
    if(j == 0) np = npj;
    else if(abs(npj)%2 != abs(np)%2){ 
      ERR.General("","getProjSign","Momentum component %d of G-parity momentum (%d,%d,%d) is invalid because it doesn't differ from component 0 by multiple of 2pi (4 in these units). Got np(0)=%d, np(j)=%d\n",j,p[0],p[1],p[2],np,npj); 
    }
  }
  int sgn = (abs(np) % 2 == 0 ? 1 : -1); //exp(i\pi n_p) = exp(-i\pi n_p)  for all integer n_p
  //if(!UniqueID()){ printf("getProjSign got sign %d (np = %d) for p=(%d,%d,%d)pi/2L\n",sgn,np,p[0],p[1],p[2]); fflush(stdout); }

  return sgn;
}

inline FlavorMatrix getProjector(const ThreeMomentum &p){
  const int proj_sign = getProjSign(p.ptr());
  //1/2(1 +/- sigma2)
  FlavorMatrix proj(0.0); 
  proj(0,0) = proj(1,1) = 0.5;
  proj(0,1) = Complex(0,-0.5*proj_sign);
  proj(1,0) = Complex(0,0.5*proj_sign);
  return proj;
}

template<typename MatrixType>
class SrcSnkOp{
public:
  virtual void rightMultiply(MatrixType &prop) const = 0;
};

class BasicOp : public SrcSnkOp<WilsonMatrix>{
  SpinMatrixType smat;
public:
  BasicOp(const SpinMatrixType _smat): smat(_smat){}
  void rightMultiply(WilsonMatrix &prop) const{
    switch(smat){
    case(gamma1):
      prop.gr(0); break;
    case(gamma2):
      prop.gr(1); break;
    case(gamma3):
      prop.gr(2); break;
    case(gamma4):
      prop.gr(3); break;
    case(gamma5):
      prop.gr(-5); break;
    case(spin_unit):
    default:
      break;
    }
  }
};



class BasicGparityOp : public SrcSnkOp<SpinColorFlavorMatrix>{
  SpinMatrixType smat;
  FlavorMatrixType fmat;
public:
  BasicGparityOp(const SpinMatrixType _smat, const FlavorMatrixType _fmat): smat(_smat), fmat(_fmat){}
  void rightMultiply(SpinColorFlavorMatrix &prop) const{
    switch(smat){
    case(gamma1):
      prop.gr(0); break;
    case(gamma2):
      prop.gr(1); break;
    case(gamma3):
      prop.gr(2); break;
    case(gamma4):
      prop.gr(3); break;
    case(gamma5):
      prop.gr(-5); break;
    case(spin_unit):
    default:
      break;
    }
    prop.pr(fmat);
  }
};

class GparityOpWithFlavorProject : public BasicGparityOp{
  ThreeMomentum p_psi;
public:
  GparityOpWithFlavorProject(const SpinMatrixType _smat, const FlavorMatrixType _fmat, const ThreeMomentum &_p_psi, bool _use_wrong_proj_sign = false): BasicGparityOp(_smat,_fmat), p_psi(_p_psi){}

  std::string printProj() const{
    const int proj_sign = getProjSign(p_psi.ptr());
    std::ostringstream os; os << "1/2(1" << (proj_sign == 1 ? '+' : '-') << "s2)";
    return os.str();
  }

  void rightMultiply(SpinColorFlavorMatrix &prop) const{
    this->BasicGparityOp::rightMultiply(prop);
    FlavorMatrix proj = getProjector(p_psi);    
    prop *= proj;
  }
};



CPS_END_NAMESPACE

#endif
