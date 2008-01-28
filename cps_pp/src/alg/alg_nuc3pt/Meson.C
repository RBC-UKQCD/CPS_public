//------------------------------------------------------------------
//
// Meson.C
//
// Implementation of class Meson.
// Can be used to constuct mesons with any algorith that uses QPropW
//
// February 2002
//
// Kostas Orginos
//
//------------------------------------------------------------------
#include <alg/meson.h>

CPS_START_NAMESPACE

Meson::Meson(int mu, int nu,char *src_i)
{
  char *cname = "Meson";
  char *fname = "Meson()";
  VRB.Func(cname,fname);

  gamma[0] = mu ;
  gamma[1] = nu ;

  m[0] = m[1] = -1111; 

  strcpy(src,src_i);
  setGamma(mu,nu) ;
}

Meson::Meson(int mu, char *src_i)
{
  char *cname = "Meson";
  char *fname = "Meson()";
  VRB.Func(cname,fname);

  gamma[0] = mu ;
  gamma[1] = 1969 ;

  m[0] = m[1] = -1111; 

  setGamma(mu) ;
  strcpy(src,src_i);
}

Meson::Meson(char *src_i)
{
  char *cname = "Meson";
  char *fname = "Meson()";
  VRB.Func(cname,fname);

  gamma[0] = 1969 ;
  gamma[1] = 1969 ;

  m[0] = m[1] = -1111; 

  setGamma();
  strcpy(src,src_i);
}




void Meson::calcMeson(QPropW& q1, QPropW& q2)//adds to the correlation function
{
  char *cname = "Meson";
  char *fname = "calcMeson()";
  VRB.Func(cname,fname);

  // m[i] should have been set by setMass!
  if((m[0]!=q1.Mass())||(m[1]!=q2.Mass())) 
  ERR.General(cname,fname,"OOPS: Masses changed!\n");
  
  int t_src(q1.SourceTime());
  if(t_src!=q2.SourceTime())
    ERR.General(cname,fname,"Source times not equal!\n");
  

  int t ;
  int first_index(0) ;
  int second_index(0) ;

  WilsonMatrix temp;

  if( (gamma[0]==-5) || ((gamma[0]>-1) && (gamma[0]<4)) )
    first_index = 1 ;

  if( (gamma[1]==-5) || ((gamma[1]>-1) && (gamma[1]<4)) )
    second_index = 1 ;
  // If both indices are out of range then the Scalar meson is constructed

  int Nt(func.TimeSize());

  for(int i=0; i< GJP.VolNodeSites(); i++)
    {
      t=i/(GJP.VolNodeSites()/GJP.TnodeSites()); // t coordinate in node 
      t += GJP.TnodeCoor()*GJP.TnodeSites(); // plus t node offset physical t
      
      temp = q2[i] ;
      temp.hconj() ;
      temp.gr(-5).gl(-5);
      //temp is the antiquark

      if(first_index)  temp.gl(gamma[0]).gr(gamma[0]) ;
      if(second_index) temp.gl(gamma[1]).gr(gamma[1]) ;
      
      //Shift the time to make source be at t=0
      func[(t + Nt - t_src)%Nt] += Trace(q1[i],temp) ; 
      // meson correlator point sink
    }

  // Global sum the correlators
  func.GlobalSum() ;

}

void Meson::calcMidPointPion(QPropW& q1, QPropW& q2)//adds to the correlation function
{
  char *cname = "Meson";
  char *fname = "calcMidPointPion()";
  VRB.Func(cname,fname);

  // m[i] should have been set by setMass!
  if((m[0]!=q1.Mass())||(m[1]!=q2.Mass())) 
  ERR.General(cname,fname,"OOPS: Masses changed!\n");
  
  int t_src(q1.SourceTime());
  if(t_src!=q2.SourceTime())
    ERR.General(cname,fname,"Source times not equal!\n");
  

  int t ;

  WilsonMatrix temp;

  int Nt(func.TimeSize());

  for(int i=0; i< GJP.VolNodeSites(); i++)
    {
      t=i/(GJP.VolNodeSites()/GJP.TnodeSites()); // t coordinate in node 
      t += GJP.TnodeCoor()*GJP.TnodeSites(); // plus t node offset physical t
      
      temp = q2(i) ; // the "()" returns the midpoint prop.
      temp.hconj() ;
      //temp is the antiquark

      //Shift the time to make source be at t=0
      func[(t + Nt - t_src)%Nt] += Trace(q1(i),temp) ; 
      // meson correlator point sink
    }

  // Global sum the correlators
  func.GlobalSum() ;

}

void Meson::setSrc(char *src_i)
{
  if(strlen(src_i)<9)
    strcpy(src,src_i); 
  else
    {
      char *cname="Meson()";
      char *fname="setSrc" ;
      VRB.Warn(cname,fname,"SRC tag too long!\n");
      strncpy(src,src_i,8);
      src[8] ='\0' ;
    }  
}

CPS_END_NAMESPACE
