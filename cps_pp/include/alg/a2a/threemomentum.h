#ifndef _THREEMOMENTUM_H
#define _THREEMOMENTUM_H

CPS_START_NAMESPACE

class ThreeMomentum{
  int p[3];
public:
  ThreeMomentum(){ p[0]=p[1]=p[2]=0; }
  explicit ThreeMomentum(const int all){ p[0]=p[1]=p[2]=all; }
  ThreeMomentum(const int px, const int py, const int pz){ p[0]=px; p[1]=py; p[2]=pz; }
  ThreeMomentum(const int _p[3]){ memcpy(p,_p,3*sizeof(int)); }
  
  int &operator()(const int i){ return p[i]; }
  const int operator()(const int i) const{ return p[i]; }

  bool operator<(const ThreeMomentum &r) const{
    for(int i=0;i<3;i++){
      if(p[i] < r.p[i]) return true;
      else if(p[i] > r.p[i]) return false;
    }//lexicographic comparison
    return false; //all equal
  }
  bool operator>(const ThreeMomentum &r) const{
    return (*this)!=r && ! ( (*this) < r );
  }
  
  inline static ThreeMomentum negative(const ThreeMomentum &p){
    return ThreeMomentum( -p(0), -p(1), -p(2) );
  }
  ThreeMomentum operator-() const{ return ThreeMomentum( -p[0], -p[1], -p[2] ); }

  bool operator==(const ThreeMomentum &r) const{ return p[0]==r.p[0] && p[1]==r.p[1] && p[2]==r.p[2]; }
  bool operator!=(const ThreeMomentum &r) const{ return p[0]!=r.p[0] || p[1]!=r.p[1] || p[2]!=r.p[2]; }

  //Parse a string "(%d,%d,%d)" into a ThreeMomentum  
  static ThreeMomentum parse_str(const std::string &p){
    ThreeMomentum out;
    int ret = sscanf(p.c_str(),"(%d,%d,%d)",&out.p[0],&out.p[1],&out.p[2]);
    if(ret != 3){ ERR.General("ThreeMomentum","parse_str","Could not parse '%s'\n",p.c_str()); }
    return out;
  }
  //Parse a string "(%d,%d,%d) + (%d,%d,%d)"  (spaces included) into a pair of ThreeMomentum
  static std::pair< ThreeMomentum, ThreeMomentum > parse_str_two_mom(const std::string &pp){
    std::pair< ThreeMomentum, ThreeMomentum  > out;
    
    int ret = sscanf(pp.c_str(),"(%d,%d,%d) + (%d,%d,%d)",
		     &out.first.p[0],&out.first.p[1],&out.first.p[2],
		     &out.second.p[0],&out.second.p[1],&out.second.p[2]
		     );
    if(ret != 6){ ERR.General("ThreeMomentum","parse_str_two_mom","Could not parse '%s'\n",pp.c_str()); }
    return out;
  }
  //Get a pointer to the underlying array
  int const* ptr() const{ return &p[0]; }
  int * ptr(){ return &p[0]; }

  //Return momentum as string "(%d,%d,%d)"
  std::string str() const{
    std::ostringstream os; os << "(" << p[0] << "," << p[1] << "," << p[2] << ")";
    return os.str();
  }
  //Return momentum as string suitable for filenames. Negative momenta are prefixed with '_' so
  //e.g.  p=(1,-2,3) becomes 1_23
  //Daiqian writes meson momentum in units of pi/L, whereas the natural G-parity units for ThreeMomentum are pi/2L.
  //Here we allow for a numerical factor to divide by.
  std::string file_str(const int divide_by = 1) const{
    std::ostringstream os; 
    for(int i=0;i<3;i++){
      if(p[i] < 0) os << '_';
      os << abs(p[i]/divide_by);
    }
    return os.str();
  }


  ThreeMomentum &operator+=(const ThreeMomentum &r){ for(int i=0;i<3;i++) p[i] += r.p[i]; return *this; }
  ThreeMomentum operator+(const ThreeMomentum &r) const{ ThreeMomentum out(*this); out += r; return out; }

  ThreeMomentum &operator-=(const ThreeMomentum &r){ for(int i=0;i<3;i++) p[i] -= r.p[i]; return *this; }
  ThreeMomentum operator-(const ThreeMomentum &r) const{ ThreeMomentum out(*this); out -= r; return out; }

  ThreeMomentum &operator*=(const int r){ for(int i=0;i<3;i++) p[i] *= r; return *this; }
  ThreeMomentum operator*(const int r) const{ ThreeMomentum out(*this); out *= r; return out; }
  
  //Convert the momenta into lattice units
  void latticeUnits(Float p_latt[3]) const{
    for(int i=0;i<3;i++){
      double L = GJP.NodeSites(i)*GJP.Nodes(i);
      double unit;
      switch(GJP.Bc(i)){
      case BND_CND_PRD:
	unit = 2.*M_PI/L; break;
      case BND_CND_APRD:
	unit = M_PI/L; break;
      case BND_CND_GPARITY:
	unit = M_PI/2./L; break;
      default:
	ERR.General("ThreeMomentum","latticeUnits","Unknown BC\n");
      }
      p_latt[i] = p[i] * unit;
    }
  }

};

//exp(-i p.x)
//if xlocal == false, it assumes the coordinate is global rather than node-local
inline std::complex<double> mom_phase(const ThreeMomentum &p, const int x[3], const bool xlocal = true){
  Float p_latt[3];
  p.latticeUnits(p_latt);

  double pdx = 0.;
  for(int i=0;i<3;i++){
    int xi = x[i];
    if(xlocal) xi += GJP.NodeCoor(i)*GJP.NodeSites(i);
    pdx += p_latt[i] * xi;
  }
  return std::complex<double>( cos(pdx), -sin(pdx) );
}


CPS_END_NAMESPACE
#endif
