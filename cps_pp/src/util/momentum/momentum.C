#include <config.h>
#include <util/momentum.h>
#include <math.h>

CPS_START_NAMESPACE

//#define TWO_PI 6.283185308
#define TWO_PI (M_PI*2)
void ThreeMom::CalcLatMom(void){
  pp[0] = TWO_PI/(GJP.XnodeSites()*GJP.Xnodes()) ;
  pp[1] = TWO_PI/(GJP.YnodeSites()*GJP.Ynodes()) ;
  pp[2] = TWO_PI/(GJP.ZnodeSites()*GJP.Znodes()) ;
}
#undef TWO_PI

void ThreeMom::setZeroMomFlag(void){
  ZeroMom = true ;
  for(int i(0);i<3;i++)
    if(p[i]!=0)
      { 
	ZeroMom = false ;
	break ;
      }
}

// exp(-i p*x) 
Complex  ThreeMom::Fact(Site& s){  
  if(ZeroMom)
    {
      Complex F(1.0,0.0) ;
      return F ;
    }
  
  Float xp(0.0) ;
  
  for(int i = 0 ; i< 3; i++)
    if(p[i]!=0)
      {
	int xx ;
	xx = p[i]*(s.physCoor(i)) ;
	xp += pp[i]*xx ;
      }
  
  Complex F(cos(xp), -sin(xp) );
  
  return F ;
}

// cos(p1*x1)*cos(p2*x2)*cos(p3*x3)
Complex  ThreeMom::FactCos(Site& s){  
  if(ZeroMom)
    {
      Complex F(1.0,0.0) ;
      return F ;
    }
  
  Float ccc(1.0) ;
  
  for(int i = 0 ; i< 3; i++)
    if(p[i]!=0)
      {
	int xx ;
	xx = p[i]*(s.physCoor(i)) ;
	ccc *= cos(pp[i]*xx) ;
      }
  
  Complex F(ccc, 0 );
  
  return F ;
}

Complex ThreeMom::Fact(Site& s, int *sx) // exp(-i p*(x+sx)) 
{
  // WARNING we have to take care of the boundary condition factor
  // As things are now it only work with periodic BC in spatial
  // directions
  if(ZeroMom)
    {
      Complex F(1.0,0.0) ;
      return F ;
    }
  
  Float xp(0.0) ;
  for(int i = 0 ; i< 3; i++)
    if(p[i]!=0)
      {
	int xx ;
	xx = p[i]*(s.physCoor(i) + sx[i]) ;
	xp += pp[i]*xx ;
      }
  
  Complex F(cos(xp), -sin(xp) );
  
  return F ;
}

ThreeMom::ThreeMom() 
{ 
  CalcLatMom() ; 
  p[0]=0 ;
  p[1]=0 ;
  p[2]=0 ;
  
  setZeroMomFlag() ;
}

ThreeMom::ThreeMom(const int *q)
{
  CalcLatMom() ;
  p[0]=q[0] ;
  p[1]=q[1] ;
  p[2]=q[2] ;
  
  setZeroMomFlag() ;
}

ThreeMom::ThreeMom(int q0,int q1,int q2)
{
  CalcLatMom() ;
  p[0]=q0 ;
  p[1]=q1 ;
  p[2]=q2 ;
  
  setZeroMomFlag() ;
}

ThreeMom::ThreeMom(const ThreeMom& rhs)
{
  p[0] = rhs.p[0] ;
  p[1] = rhs.p[1] ;
  p[2] = rhs.p[2] ;
  pp[0] = rhs.pp[0] ;
  pp[1] = rhs.pp[1] ;
  pp[2] = rhs.pp[2] ;
  ZeroMom = rhs.ZeroMom ;
}

ThreeMomTwist::ThreeMomTwist() : ThreeMom()
{
  for (int i=0 ; i<3 ; i++)
    pp[i]*=0.5;
}

ThreeMomTwist::ThreeMomTwist(const int *q) : ThreeMom(q)
{
  for (int i=0 ; i<3 ; i++)
    pp[i]*=0.5;
}

ThreeMomTwist::ThreeMomTwist(int q0,int q1,int q2) : ThreeMom(q0,q1,q2)
{
  for (int i=0 ; i<3 ; i++)
    pp[i]*=0.5;
}

ThreeMomTwist::ThreeMomTwist(const ThreeMomTwist& rhs) : ThreeMom(rhs)
{ }


CPS_END_NAMESPACE
