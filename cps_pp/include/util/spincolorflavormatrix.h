#ifndef SPIN_COLOR_FLAVOR_MATRIX_H
#define SPIN_COLOR_FLAVOR_MATRIX_H

#include<config.h>
#include <alg/alg_base.h>
#include <alg/qpropw.h>
#include <alg/prop_attribute_arg.h>
#include <alg/propagatorcontainer.h>
#include <alg/propmanager.h>
#include <alg/wilson_matrix.h>
#include <util/spinflavormatrix.h>
#include <util/flavormatrix.h>

CPS_START_NAMESPACE

enum SpinMatrixType { gamma1, gamma2, gamma3, gamma4, gamma5, spin_unit };
enum PropSplane { SPLANE_BOUNDARY, SPLANE_MIDPOINT }; //Usual boundary 5d propagator or the midpoint propagator

class SpinColorFlavorMatrix{
protected:
  WilsonMatrix wmat[2][2];
  const char *cname;

  inline static WilsonMatrix & getSite(const int site, const int flav, QPropW &from, const PropSplane splane){
    return splane == SPLANE_BOUNDARY ? from.SiteMatrix(site,flav) : from.MidPlaneSiteMatrix(site,flav);
  }

  inline static WilsonMatrix & getSite(const int site, const int flav, QPropWcontainer &from,  Lattice &lattice, const PropSplane splane){
    return splane == SPLANE_BOUNDARY ? from.getProp(lattice).SiteMatrix(site,flav) : from.getProp(lattice).MidPlaneSiteMatrix(site,flav);
  }
public:

  SpinColorFlavorMatrix(QPropWcontainer &from, Lattice &lattice, const int &site, const PropSplane &splane = SPLANE_BOUNDARY): cname("SpinColorFlavorMatrix"){
    generate(from,lattice,site,splane);
  }
  SpinColorFlavorMatrix():cname("SpinColorFlavorMatrix"){
  }
  SpinColorFlavorMatrix(const SpinColorFlavorMatrix &from): cname("SpinColorFlavorMatrix"){
    wmat[0][0] = from.wmat[0][0];
    wmat[0][1] = from.wmat[0][1];
    wmat[1][0] = from.wmat[1][0];
    wmat[1][1] = from.wmat[1][1];
  }
  SpinColorFlavorMatrix(const Float &rhs): cname("SpinColorFlavorMatrix"){
    for(int i=0;i<2;i++)
      for(int j=0;j<2;j++)
	wmat[i][j] = rhs;
  }

  SpinColorFlavorMatrix(const SpinMatrixType &g, const FlavorMatrixType &f): cname("SpinColorFlavorMatrix"){
    static const int spin_map[] = { 0,1,2,3,-5 };
    Unit();
    if(g != spin_unit) gr(spin_map[ (int)g - (int)gamma1 ]);
    if(f != sigma0) pr(f);
  }

  void generate(QPropW &from_f0, QPropW &from_f1, const int site, const PropSplane splane = SPLANE_BOUNDARY);
  void generate(QPropWcontainer &from_f0, QPropWcontainer &from_f1, Lattice &lattice, const int site, const PropSplane splane = SPLANE_BOUNDARY);


  //Use the propagator conjugate relation to compute the 2x2 flavor matrix propagator using a  propagator from a single flavor and one with the complex conjugate of the source used for the first
  void generate_from_cconj_pair(QPropW &from, QPropW &from_conj, const int from_flav, const int site, const PropSplane splane = SPLANE_BOUNDARY);
  void generate_from_cconj_pair(QPropWcontainer &from, QPropWcontainer &from_conj, Lattice &lattice, const int site, const PropSplane splane = SPLANE_BOUNDARY);

  //Use the prop conj relation for a real source to compute the full 2x2 propagator  
  void generate_from_real_source(QPropWcontainer &from, Lattice &lattice, const int &site, const PropSplane &splane = SPLANE_BOUNDARY);

  void generate(QPropWcontainer &from, Lattice &lattice, const int &site, const PropSplane &splane = SPLANE_BOUNDARY);

  //We can use the G-parity prop-conj relation to flip the source momentum a posteriori providing the source has the form  e^{-ipx} \eta    where \eta obeys  \sigma_2 C\gamma^5 \eta^* = \eta\sigma_2 C\gamma^5
  //This applies to pretty much all standard source types
  void flipSourceMomentum();

  
  //multiply on left by a flavor matrix
  SpinColorFlavorMatrix & pl(const FlavorMatrixType &type){
    if(type == F0){
      wmat[1][0] = 0.0;
      wmat[1][1] = 0.0;
      return *this;
    }else if(type == F1){
      wmat[0][0] = 0.0;
      wmat[0][1] = 0.0;
      return *this;
    }else if(type == Fud){
      WilsonMatrix _00(wmat[0][0]);
      WilsonMatrix _01(wmat[0][1]);
      wmat[0][0] = wmat[1][0];
      wmat[0][1] = wmat[1][1];
      wmat[1][0] = _00;
      wmat[1][1] = _01;
      return *this;
    }else if(type == sigma0){
      return *this;
    }else if(type == sigma1){
      return pl(Fud);
    }else if(type == sigma2){
      WilsonMatrix _i00(wmat[0][0]); _i00*=Complex(0.0,1.0);
      WilsonMatrix _i01(wmat[0][1]); _i01*=Complex(0.0,1.0);
      wmat[0][0] = wmat[1][0]; wmat[0][0]*=Complex(0.0,-1.0);
      wmat[0][1] = wmat[1][1]; wmat[0][1]*=Complex(0.0,-1.0);
      wmat[1][0] = _i00;
      wmat[1][1] = _i01;
      return *this;
    }else if(type == sigma3){
      wmat[1][0]*=-1.0;
      wmat[1][1]*=-1.0;
      return *this;
    }
    ERR.General(cname,"pl(const FlavorMatrixType &type)","Unknown FlavorMatrixType");
  }
  //multiply on right by a flavor matrix
  SpinColorFlavorMatrix & pr(const FlavorMatrixType &type){
    if(type == F0){
      wmat[0][1] = 0.0;
      wmat[1][1] = 0.0;
      return *this;
    }else if(type == F1){
      wmat[0][0] = 0.0;
      wmat[1][0] = 0.0;
      return *this;
    }else if(type == Fud){
      WilsonMatrix _00(wmat[0][0]);
      WilsonMatrix _10(wmat[1][0]);
      wmat[0][0] = wmat[0][1];
      wmat[1][0] = wmat[1][1];
      wmat[0][1] = _00;
      wmat[1][1] = _10;
      return *this;
    }else if(type == sigma0){
      return *this;
    }else if(type == sigma1){
      return pr(Fud);
    }else if(type == sigma2){
      WilsonMatrix _mi00(wmat[0][0]); _mi00 *= Complex(0.0,-1.0);
      WilsonMatrix _mi10(wmat[1][0]); _mi10 *= Complex(0.0,-1.0);
      wmat[0][0] = wmat[0][1]; wmat[0][0] *= Complex(0.0,1.0); 
      wmat[1][0] = wmat[1][1]; wmat[1][0] *= Complex(0.0,1.0);
      wmat[0][1] = _mi00;
      wmat[1][1] = _mi10;
      return *this;
    }else if(type == sigma3){
      wmat[0][1]*=-1.0;
      wmat[1][1]*=-1.0;
      return *this;
    }
    ERR.General(cname,"pr(const FlavorMatrixType &type)","Unknown FlavorMatrixType");
  }
  WilsonMatrix FlavorTrace() const{
    WilsonMatrix out(wmat[0][0]);
    out+=wmat[1][1];
    return out;
  }
  Complex Trace() const{
    Complex out(0);
    for(int f=0;f<2;f++) out += wmat[f][f].Trace();
    return out;
  }
  Matrix SpinFlavorTrace() const{
    return SpinTrace(FlavorTrace());
  }
  inline FlavorSpinMatrix ColorTrace() const{ 
    FlavorSpinMatrix out;
    for(int f1 = 0; f1 < 2; f1++) 
      for(int f2 = 0; f2 < 2; f2++) 
	wmat[f1][f2].ColorTrace(out(f1,f2));

	//out(f2,f1) = cps::ColorTrace(wmat[f2][f1]);
    return out;
  }  


  SpinColorFlavorMatrix operator*(const SpinColorFlavorMatrix& rhs) const{
    SpinColorFlavorMatrix out;
    out.wmat[0][0] = wmat[0][0]*rhs.wmat[0][0] + wmat[0][1]*rhs.wmat[1][0];
    out.wmat[1][0] = wmat[1][0]*rhs.wmat[0][0] + wmat[1][1]*rhs.wmat[1][0];
    out.wmat[0][1] = wmat[0][0]*rhs.wmat[0][1] + wmat[0][1]*rhs.wmat[1][1];
    out.wmat[1][1] = wmat[1][0]*rhs.wmat[0][1] + wmat[1][1]*rhs.wmat[1][1];
    return out;
  }

  SpinColorFlavorMatrix& operator*=(const SpinColorFlavorMatrix& rhs){
    SpinColorFlavorMatrix cp(*this);
    wmat[0][0] = cp.wmat[0][0]*rhs.wmat[0][0] + cp.wmat[0][1]*rhs.wmat[1][0];
    wmat[1][0] = cp.wmat[1][0]*rhs.wmat[0][0] + cp.wmat[1][1]*rhs.wmat[1][0];
    wmat[0][1] = cp.wmat[0][0]*rhs.wmat[0][1] + cp.wmat[0][1]*rhs.wmat[1][1];
    wmat[1][1] = cp.wmat[1][0]*rhs.wmat[0][1] + cp.wmat[1][1]*rhs.wmat[1][1];
    return *this;
  }
  SpinColorFlavorMatrix& operator*=(const FlavorMatrix& rhs){
    SpinColorFlavorMatrix cp(*this);
    wmat[0][0] = cp.wmat[0][0]*rhs(0,0) + cp.wmat[0][1]*rhs(1,0);
    wmat[1][0] = cp.wmat[1][0]*rhs(0,0) + cp.wmat[1][1]*rhs(1,0);
    wmat[0][1] = cp.wmat[0][0]*rhs(0,1) + cp.wmat[0][1]*rhs(1,1);
    wmat[1][1] = cp.wmat[1][0]*rhs(0,1) + cp.wmat[1][1]*rhs(1,1);
    return *this;
  }

  SpinColorFlavorMatrix& operator*=(const Float& rhs){
    for(int i=0;i<2;i++)
      for(int j=0;j<2;j++)
	wmat[i][j]*=rhs;
    return *this;
  }
  SpinColorFlavorMatrix operator*(const Float& rhs) const{
    SpinColorFlavorMatrix out(*this);
    out*=rhs;
    return out;
  }
  SpinColorFlavorMatrix& operator*=(const Complex& rhs){
    for(int i=0;i<2;i++)
      for(int j=0;j<2;j++)
	wmat[i][j]*=rhs;
    return *this;
  }
  SpinColorFlavorMatrix operator*(const Complex& rhs) const{
    SpinColorFlavorMatrix out(*this);
    out*=rhs;
    return out;
  }
  SpinColorFlavorMatrix& operator*=(const WilsonMatrix& rhs){
    for(int i=0;i<2;i++)
      for(int j=0;j<2;j++)
	wmat[i][j]*=rhs;
    return *this;
  }
  SpinColorFlavorMatrix operator*(const WilsonMatrix& rhs) const{
    SpinColorFlavorMatrix out(*this);
    out*=rhs;
    return out;
  }

  SpinColorFlavorMatrix& operator+=(const SpinColorFlavorMatrix& rhs){
    wmat[0][0]+= rhs.wmat[0][0];
    wmat[0][1]+= rhs.wmat[0][1];
    wmat[1][0]+= rhs.wmat[1][0];
    wmat[1][1]+= rhs.wmat[1][1];
    return *this;
  }
  SpinColorFlavorMatrix operator+(const SpinColorFlavorMatrix& rhs) const{
    SpinColorFlavorMatrix out(*this);
    out+=rhs;
    return out;
  }
  SpinColorFlavorMatrix& operator-=(const SpinColorFlavorMatrix& rhs){
    wmat[0][0]-= rhs.wmat[0][0];
    wmat[0][1]-= rhs.wmat[0][1];
    wmat[1][0]-= rhs.wmat[1][0];
    wmat[1][1]-= rhs.wmat[1][1];
    return *this;
  }
  SpinColorFlavorMatrix operator-(const SpinColorFlavorMatrix& rhs) const{
    SpinColorFlavorMatrix out(*this);
    out-=rhs;
    return out;
  }
    
  SpinColorFlavorMatrix& operator= (const Float& rhs){
    for(int i=0;i<2;i++)
      for(int j=0;j<2;j++)
	wmat[i][j] = rhs;
    return *this;
  }
  SpinColorFlavorMatrix& operator= (const Complex& rhs){
    for(int i=0;i<2;i++)
      for(int j=0;j<2;j++)
	wmat[i][j] = WilsonMatrix(rhs);
    return *this;
  }
  SpinColorFlavorMatrix& operator=(const SpinColorFlavorMatrix &from){
    wmat[0][0] = from.wmat[0][0];
    wmat[0][1] = from.wmat[0][1];
    wmat[1][0] = from.wmat[1][0];
    wmat[1][1] = from.wmat[1][1];
    return *this;
  }
  
  SpinColorFlavorMatrix& gl(int dir){
    for(int i=0;i<2;i++)
      for(int j=0;j<2;j++)
	wmat[i][j].gl(dir);
    return *this;
  }
  SpinColorFlavorMatrix& gr(int dir){
    for(int i=0;i<2;i++)
      for(int j=0;j<2;j++)
	wmat[i][j].gr(dir);
    return *this;
  }

  //Left-multiply this by \gamma^5\gamma^\mu
  inline SpinColorFlavorMatrix& glAx(int mu){
    for(int i=0;i<2;i++)
      for(int j=0;j<2;j++)
	wmat[i][j].glAx(mu);
    return *this;
  }


  SpinColorFlavorMatrix glL(int dir){
		SpinColorFlavorMatrix out;
    for(int i=0;i<2;i++)
      for(int j=0;j<2;j++)
				out(i,j) = this->wmat[i][j].glL(dir);
    return out;
  }
  SpinColorFlavorMatrix glR(int dir){
		SpinColorFlavorMatrix out;
    for(int i=0;i<2;i++)
      for(int j=0;j<2;j++)
				out(i,j) = this->wmat[i][j].glR(dir);
    return out;
  }

  SpinColorFlavorMatrix& ccl(int dir){
    for(int i=0;i<2;i++)
      for(int j=0;j<2;j++)
	wmat[i][j].ccl(dir);
    return *this;
  }
  SpinColorFlavorMatrix& ccr(int dir){
    for(int i=0;i<2;i++)
      for(int j=0;j<2;j++)
	wmat[i][j].ccr(dir);
    return *this;
  }
 //! left multiply by 1/2(1-gamma_5)
  SpinColorFlavorMatrix& glPL(){
    for(int i=0;i<2;i++)
      for(int j=0;j<2;j++)
	wmat[i][j].glPL();
    return *this;
  }
  //! left multiply by 1/2(1+gamma_5)
  SpinColorFlavorMatrix& glPR(){
    for(int i=0;i<2;i++)
      for(int j=0;j<2;j++)
	wmat[i][j].glPR();
    return *this;
  }

  //! right multiply by 1/2(1-gamma_5)
  SpinColorFlavorMatrix& grPL(){
    for(int i=0;i<2;i++)
      for(int j=0;j<2;j++)
	wmat[i][j].grPL();
    return *this;
  }

  //! right multiply by 1/2(1+gamma_5)
  SpinColorFlavorMatrix& grPR(){
    for(int i=0;i<2;i++)
      for(int j=0;j<2;j++)
	wmat[i][j].grPR();
    return *this;
  }
  IFloat norm() const{
    IFloat out(0.0);
    for(int i=0;i<2;i++)
      for(int j=0;j<2;j++)
	out += wmat[i][j].norm();
    return out;
  }
  //! hermitean conjugate
  SpinColorFlavorMatrix & hconj(){
    WilsonMatrix hc01 = wmat[0][1]; hc01.hconj();
    wmat[0][0].hconj();
    wmat[1][1].hconj();
    wmat[0][1] = wmat[1][0];
    wmat[0][1].hconj();
    wmat[1][0] = hc01;
    return *this;
  }
 
  //! complex conjugate
  SpinColorFlavorMatrix & cconj(){
    for(int i=0;i<2;i++)
      for(int j=0;j<2;j++)
	wmat[i][j].cconj();
    return *this;
  }

  //spin color and flavor transpose
  SpinColorFlavorMatrix & transpose(){
    WilsonMatrix _01tr = wmat[0][1]; _01tr.transpose();
    wmat[0][0].transpose();
    wmat[1][1].transpose();
    wmat[0][1] = wmat[1][0];
    wmat[0][1].transpose();
    wmat[1][0] = _01tr;
    return *this;
  }
  //just flavor
  SpinColorFlavorMatrix & tranpose_flavor(){
    WilsonMatrix _01 = wmat[0][1];
    wmat[0][1] = wmat[1][0];
    wmat[1][0] = _01;
    return *this;
  }
  SpinColorFlavorMatrix & transpose_color(){
    wmat[0][0].transpose_color();
    wmat[0][1].transpose_color();
    wmat[1][0].transpose_color();
    wmat[1][1].transpose_color();
	return *this;
  }


  inline Complex& operator()(int s1, int c1, int f1, int s2, int c2, int f2){
    return wmat[f1][f2](s1,c1,s2,c2);
  }  
  inline const Complex operator()(int s1, int c1, int f1, int s2, int c2, int f2) const{
    return wmat[f1][f2](s1,c1,s2,c2);
  }
  WilsonMatrix &operator()(int f1,int f2){
    return wmat[f1][f2];
  }
  const WilsonMatrix &operator()(int f1,int f2) const{
    return wmat[f1][f2];
  }
  SpinColorFlavorMatrix& LeftTimesEqual(const SpinColorFlavorMatrix& lhs){
    SpinColorFlavorMatrix tmp(*this);
    for(int f1=0;f1<2;f1++) for(int f2=0;f2<2;f2++){
	wmat[f1][f2] = 0.0;
	for(int f3=0;f3<2;f3++) wmat[f1][f2] += lhs.wmat[f1][f3]*tmp.wmat[f3][f2];
    }
    return *this;
  }
  SpinColorFlavorMatrix& LeftTimesEqual(const WilsonMatrix& lhs){
    for(int f1=0;f1<2;f1++) for(int f2=0;f2<2;f2++) wmat[f1][f2].LeftTimesEqual(lhs);
    return *this;
  }
  SpinColorFlavorMatrix& LeftTimesEqual(const Matrix& lhs){
    for(int f1=0;f1<2;f1++) for(int f2=0;f2<2;f2++) wmat[f1][f2].LeftTimesEqual(lhs);
    return *this;
  }

  void Unit(){
    wmat[0][0].Unit();
    wmat[0][1] = 0;
    wmat[1][0] = 0;
    wmat[1][1].Unit();	
  }
};

inline Complex Trace(const SpinColorFlavorMatrix& a){ return a.Trace(); }
Complex Trace(const SpinColorFlavorMatrix& a, const SpinColorFlavorMatrix& b);
Matrix SpinFlavorTrace(const SpinColorFlavorMatrix& a, const SpinColorFlavorMatrix& b);
SpinFlavorMatrix ColorTrace(const SpinColorFlavorMatrix& a, const SpinColorFlavorMatrix& b);
FlavorSpinMatrix ColorTraceFS(const SpinColorFlavorMatrix& a, const SpinColorFlavorMatrix& b);



typedef struct { wilson_vector f[2]; } flav_spin_color_vector;
typedef struct { flav_spin_color_vector c[3]; } color_flav_spin_color_vector;
typedef struct { color_flav_spin_color_vector d[4]; } spin_color_flav_spin_color_vector;
typedef struct { spin_color_flav_spin_color_vector f[2]; } spin_color_flavor_matrix;

class FlavorSpinColorMatrix{
public:
  enum PropSplane { SPLANE_BOUNDARY, SPLANE_MIDPOINT }; //Usual boundary 5d propagator or the midpoint propagator
protected:
  spin_color_flavor_matrix p;

  inline static WilsonMatrix & getSite(const int &site, const int &flav, QPropWcontainer &from,  Lattice &lattice, const PropSplane &splane){
    return splane == SPLANE_BOUNDARY ? from.getProp(lattice).SiteMatrix(site,flav) : from.getProp(lattice).MidPlaneSiteMatrix(site,flav);
  }
public:
  Rcomplex* ptr() { return reinterpret_cast<Rcomplex*>(&p); }
  const Rcomplex* ptr() const { return reinterpret_cast<const Rcomplex*>(&p); }

  //Complex pointer offset to element
  inline static int wilson_matrix_map(const int &rs,const int&rc, const int &cs, const int &cc){
    return (cc+3*(cs+4*(rc+3*rs)));
  }
  inline static void wilson_matrix_invmap(int idx, int &rs, int&rc, int &cs, int &cc){
    cc = idx%3; idx/=3;
    cs = idx%4; idx/=4;
    rc = idx%3; idx/=3;
    rs = idx;
  }
  inline static int fsc_matrix_map(const int &rf, const int &rs,const int&rc, const int &cf, const int &cs, const int &cc){
    return (cc+3*(cs+4*(cf+2*(rc+3*(rs+4*rf)))));
  }
  inline static int fsc_matrix_invmap(int idx, int &rf, int &rs, int&rc, int &cf, int &cs, int &cc){
    cc = idx%3; idx/=3;
    cs = idx%4; idx/=4;
    cf = idx%2; idx/=2;
    rc = idx%3; idx/=3;
    rs = idx%4; idx/=4;
    rf = idx;
  }

  FlavorSpinColorMatrix(QPropWcontainer &from, Lattice &lattice, const int &site, const PropSplane &splane = SPLANE_BOUNDARY){
    generate(from,lattice,site,splane);
  }
  FlavorSpinColorMatrix(){}

  FlavorSpinColorMatrix(const FlavorSpinColorMatrix &from){
    Float* top = (Float*)&p;
    Float* fromp = (Float*)&from.p;
    memcpy((void*)top,(void*)fromp,1152*sizeof(Float));
  }
  FlavorSpinColorMatrix(const Float &rhs){
    Float* to = (Float*)&p;
    for(int i=0;i<1152;++i) to[i] = rhs;
  }
  Complex& operator()(int f1, int s1, int c1, int f2, int s2, int c2){
    return p.f[f1].d[s1].c[c1].f[f2].d[s2].c[c2];
  }  
  const Complex& operator()(int f1, int s1, int c1, int f2, int s2, int c2) const{
    return p.f[f1].d[s1].c[c1].f[f2].d[s2].c[c2];
  }
  WilsonMatrix operator()(int f1,int f2){
    WilsonMatrix out;
    for(int s1=0;s1<4;++s1) for(int c1=0;c1<3;++c1)
    for(int s2=0;s2<4;++s2) for(int c2=0;c2<3;++c2)
       out(s1,c1,s2,c2) = p.f[f1].d[s1].c[c1].f[f2].d[s2].c[c2];
  }

  void generate(QPropWcontainer &from_f0, QPropWcontainer &from_f1, Lattice &lattice, const int &site, const PropSplane &splane = SPLANE_BOUNDARY){
    WilsonMatrix* wmat[2][2] = { { &getSite(site,0,from_f0,lattice,splane), &getSite(site,0,from_f1,lattice,splane) },
				 { &getSite(site,1,from_f0,lattice,splane), &getSite(site,1,from_f1,lattice,splane) } };

    for(int f1=0;f1<2;++f1) for(int s1=0;s1<4;++s1) for(int c1=0;c1<3;++c1)
    for(int f2=0;f2<2;++f2) for(int s2=0;s2<4;++s2) for(int c2=0;c2<3;++c2)
       p.f[f1].d[s1].c[c1].f[f2].d[s2].c[c2] = (*wmat[f1][f2])(s1,c1,s2,c2);
  }

  void generate(QPropWcontainer &from, Lattice &lattice, const int &site, const PropSplane &splane = SPLANE_BOUNDARY);

  //multiply on left by a flavor matrix
  FlavorSpinColorMatrix & pl(const FlavorMatrixType &type);
  //multiply on right by a flavor matrix
  FlavorSpinColorMatrix & pr(const FlavorMatrixType &type);

  WilsonMatrix FlavourTrace(){
    WilsonMatrix out(0);
    for(int f=0;f<2;++f)
    for(int s1=0;s1<4;++s1) for(int c1=0;c1<3;++c1)
    for(int s2=0;s2<4;++s2) for(int c2=0;c2<3;++c2)
       out(s1,c1,s2,c2) += p.f[f].d[s1].c[c1].f[f].d[s2].c[c2];

    return out;
  }
  Rcomplex Trace(){
    Rcomplex out(0);
    for(int f=0;f<2;++f) for(int s=0;s<4;++s) for(int c=0;c<3;++c)
	out += p.f[f].d[s].c[c].f[f].d[s].c[c];					      
    return out;
  }

  FlavorSpinColorMatrix operator*(const FlavorSpinColorMatrix& rhs){
    FlavorSpinColorMatrix out;
    for(int f1=0;f1<2;++f1) for(int s1=0;s1<4;++s1) for(int c1=0;c1<3;++c1)
    for(int f2=0;f2<2;++f2) for(int s2=0;s2<4;++s2) for(int c2=0;c2<3;++c2){
      out(f1,s1,c1,f2,s2,c2) = 0.0;
      for(int f3=0;f3<2;++f3) for(int s3=0;s3<4;++s3) for(int c3=0;c3<3;++c3) out(f1,s1,c1,f2,s2,c2) +=  (*this)(f1,s1,c1,f3,s3,c3) * rhs(f3,s3,c3,f2,s2,c2);
    }
    return out;
  }
  
  FlavorSpinColorMatrix& operator*=(const FlavorSpinColorMatrix& rhs){
    FlavorSpinColorMatrix cp(*this);
    for(int f1=0;f1<2;++f1) for(int s1=0;s1<4;++s1) for(int c1=0;c1<3;++c1)
    for(int f2=0;f2<2;++f2) for(int s2=0;s2<4;++s2) for(int c2=0;c2<3;++c2){
      p.f[f1].d[s1].c[c1].f[f2].d[s2].c[c2] = 0.0;
      for(int f3=0;f3<2;++f3) for(int s3=0;s3<4;++s3) for(int c3=0;c3<3;++c3) p.f[f1].d[s1].c[c1].f[f2].d[s2].c[c2] +=  cp(f1,s1,c1,f3,s3,c3) * rhs(f3,s3,c3,f2,s2,c2);
    }
    return *this;
  }

  FlavorSpinColorMatrix& operator*=(const Float& rhs){
    Float *to = (Float*)&p;
    for(int i=0;i<1152;++i) to[i]*=rhs;
    return *this;
  }
  FlavorSpinColorMatrix operator*(const Float& rhs){
    FlavorSpinColorMatrix out(*this);
    out*=rhs;
    return out;
  }
  FlavorSpinColorMatrix& operator*=(const Rcomplex& rhs){
    Complex *to = (Complex*)&p;
    for(int i=0;i<576;++i) to[i]*=rhs;
    return *this;
  }
  FlavorSpinColorMatrix operator*(const Rcomplex& rhs){
    FlavorSpinColorMatrix out(*this);
    out*=rhs;
    return out;
  }
  FlavorSpinColorMatrix& operator*=(const WilsonMatrix& rhs){
    FlavorSpinColorMatrix cp(*this);
    for(int f1=0;f1<2;++f1) for(int s1=0;s1<4;++s1) for(int c1=0;c1<3;++c1)
    for(int f2=0;f2<2;++f2) for(int s2=0;s2<4;++s2) for(int c2=0;c2<3;++c2){
      p.f[f1].d[s1].c[c1].f[f2].d[s2].c[c2] = 0.0;
      for(int s3=0;s3<4;++s3) for(int c3=0;c3<3;++c3) p.f[f1].d[s1].c[c1].f[f2].d[s2].c[c2] +=  cp(f1,s1,c1,f2,s3,c3) * rhs(s3,c3,s2,c2);
    }
    return *this;
  }
  FlavorSpinColorMatrix operator*(const WilsonMatrix& rhs){
    FlavorSpinColorMatrix out(*this);
    out*=rhs;
    return out;
  }

  FlavorSpinColorMatrix& operator+=(const FlavorSpinColorMatrix& rhs){
    Float* top = (Float*)&p;
    Float* addp = (Float*)&rhs.p;
    for(int i=0;i<1152;++i) top[i] += addp[i];
    return *this;
  }
  FlavorSpinColorMatrix operator+(const FlavorSpinColorMatrix& rhs){
    FlavorSpinColorMatrix out(*this);
    out+=rhs;
    return out;
  }
  FlavorSpinColorMatrix& operator-=(const FlavorSpinColorMatrix& rhs){
    Float* top = (Float*)&p;
    Float* addp = (Float*)&rhs.p;
    for(int i=0;i<1152;++i) top[i] -= addp[i];
    return *this;
  }
  FlavorSpinColorMatrix operator-(const FlavorSpinColorMatrix& rhs){
    FlavorSpinColorMatrix out(*this);
    out-=rhs;
    return out;
  }
    
  FlavorSpinColorMatrix& operator= (const Float& rhs){
    Float* top = (Float*)&p;
    for(int i=0;i<1152;++i) top[i] = rhs;
    return *this;
  }
  FlavorSpinColorMatrix& operator= (const Rcomplex& rhs){
    Complex *to = (Complex*)&p;
    for(int i=0;i<576;++i) to[i]=rhs;
    return *this;
  }
  FlavorSpinColorMatrix& operator=(const FlavorSpinColorMatrix &from){
    Float* top = (Float*)&p;
    Float* fromp = (Float*)&from.p;
    memcpy((void*)top,(void*)fromp,1152*sizeof(Float));
    return *this;
  }
  
  FlavorSpinColorMatrix& gl(int dir);
  FlavorSpinColorMatrix& gr(int dir);

  FlavorSpinColorMatrix& ccl(int dir);//Note, this has the same incorrect ordering as for WilsonMatrix:   ccl(1)  =  C^{-1}M   ccl(-1) = CM
  FlavorSpinColorMatrix& ccr(int dir);

 //! left multiply by 1/2(1-gamma_5)
  FlavorSpinColorMatrix& glPL();

  //! left multiply by 1/2(1+gamma_5)
  FlavorSpinColorMatrix& glPR();

  //! right multiply by 1/2(1-gamma_5)
  FlavorSpinColorMatrix& grPL();

  //! right multiply by 1/2(1+gamma_5)
  FlavorSpinColorMatrix& grPR();

  IFloat norm() const{
    Float out(0.0);
    Complex *c = (Complex*)&p;
    for(int i=0;i<576;++i) out += std::norm(c[i]);
    return out;
  }
  //! hermitean conjugate
  void hconj(){
    FlavorSpinColorMatrix tmp(*this);   
    for(int f1=0;f1<2;++f1) for(int s1=0;s1<4;++s1) for(int c1=0;c1<3;++c1)
    for(int f2=0;f2<2;++f2) for(int s2=0;s2<4;++s2) for(int c2=0;c2<3;++c2){
	  p.f[f1].d[s1].c[c1].f[f2].d[s2].c[c2] = Complex(tmp.p.f[f2].d[s2].c[c2].f[f1].d[s1].c[c1].real(), -tmp.p.f[f2].d[s2].c[c2].f[f1].d[s1].c[c1].imag());
    }
  }
 
  //! complex conjugate
  void cconj(){
    Float *f = (Float*)&p;
    for(int i=1;i<1152;i+=2) f[i]*=-1;
  }

  //spin color and flavor transpose
  void transpose(){
    FlavorSpinColorMatrix tmp(*this);   
    for(int f1=0;f1<2;++f1) for(int s1=0;s1<4;++s1) for(int c1=0;c1<3;++c1)
    for(int f2=0;f2<2;++f2) for(int s2=0;s2<4;++s2) for(int c2=0;c2<3;++c2)
      p.f[f1].d[s1].c[c1].f[f2].d[s2].c[c2] = tmp.p.f[f2].d[s2].c[c2].f[f1].d[s1].c[c1];
  }
  //just flavor
  void tranpose_flavor(){
    FlavorSpinColorMatrix tmp(*this); 
    for(int f1=0;f1<2;++f1) for(int s1=0;s1<4;++s1) for(int c1=0;c1<3;++c1)
    for(int f2=0;f2<2;++f2) for(int s2=0;s2<4;++s2) for(int c2=0;c2<3;++c2)
      p.f[f1].d[s1].c[c1].f[f2].d[s2].c[c2] = tmp.p.f[f2].d[s1].c[c1].f[f1].d[s2].c[c2];
  }

  FlavorSpinColorMatrix& LeftTimesEqual(const FlavorSpinColorMatrix& lhs){
    FlavorSpinColorMatrix cp(*this);
    for(int f1=0;f1<2;++f1) for(int s1=0;s1<4;++s1) for(int c1=0;c1<3;++c1)
    for(int f2=0;f2<2;++f2) for(int s2=0;s2<4;++s2) for(int c2=0;c2<3;++c2){
      p.f[f1].d[s1].c[c1].f[f2].d[s2].c[c2] = 0.0;
      for(int f3=0;f3<2;++f3) for(int s3=0;s3<4;++s3) for(int c3=0;c3<3;++c3) p.f[f1].d[s1].c[c1].f[f2].d[s2].c[c2] +=  lhs(f1,s1,c1,f3,s3,c3) * cp(f3,s3,c3,f2,s2,c2);
    }
    return *this;
  }

  FlavorSpinColorMatrix& LeftTimesEqual(const WilsonMatrix& lhs){
    FlavorSpinColorMatrix cp(*this);
    for(int f1=0;f1<2;++f1) for(int s1=0;s1<4;++s1) for(int c1=0;c1<3;++c1)
    for(int f2=0;f2<2;++f2) for(int s2=0;s2<4;++s2) for(int c2=0;c2<3;++c2){
      p.f[f1].d[s1].c[c1].f[f2].d[s2].c[c2] = 0.0;
      for(int s3=0;s3<4;++s3) for(int c3=0;c3<3;++c3) p.f[f1].d[s1].c[c1].f[f2].d[s2].c[c2] +=  lhs(s1,c1,s3,c3) * cp(f1,s3,c3,f2,s2,c2);
    }
    return *this;
  }
  FlavorSpinColorMatrix& LeftTimesEqual(const Matrix& lhs){
    FlavorSpinColorMatrix cp(*this);
    for(int f1=0;f1<2;++f1) for(int s1=0;s1<4;++s1) for(int c1=0;c1<3;++c1)
    for(int f2=0;f2<2;++f2) for(int s2=0;s2<4;++s2) for(int c2=0;c2<3;++c2){
      p.f[f1].d[s1].c[c1].f[f2].d[s2].c[c2] = 0.0;
      for(int c3=0;c3<3;++c3) p.f[f1].d[s1].c[c1].f[f2].d[s2].c[c2] +=  lhs(c1,c3) * cp(f1,s1,c3,f2,s2,c2);
    }
    return *this;
  }

  void unit(){
    *this = 0.0;
    for(int f=0;f<2;f++) for(int s=0;s<4;s++) for(int c=0;c<3;c++) p.f[f].d[s].c[c].f[f].d[s].c[c] = 1.0;
  }
};


inline static SpinColorFlavorMatrix ColorTranspose(const SpinColorFlavorMatrix &m){
  SpinColorFlavorMatrix out;
  for(int f1=0;f1<2;++f1) for(int s1=0;s1<4;++s1) for(int c1=0;c1<3;++c1)
  for(int f2=0;f2<2;++f2) for(int s2=0;s2<4;++s2) for(int c2=0;c2<3;++c2)
    out(s1,c1,f1,s2,c2,f2) = m(s1,c2,f1,s2,c1,f2);
  return out;
}




















CPS_END_NAMESPACE

#endif

