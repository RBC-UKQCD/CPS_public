#include <util/momentum.h>

#define TWO_PI 6.283185308
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

ThreeMom::ThreeMom(int *q)
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
