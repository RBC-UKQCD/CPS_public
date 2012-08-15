#include <config.h>
#include  <alg/diquark.h>
// This makes sure that the static memeber initialization
// only gets compiled once
#define COMPILE_THE_EPSILON
#include "epsilon.h"
#undef COMPILE_THE_EPSILON

CPS_START_NAMESPACE

void Diquark::Project(WilsonVector& res, ProjectType P)
{
  char *fname = "Project(FermVec&,WilsV&)";
  char *cname = "Diquark";

  VRB.Func(cname, fname);

  res.Zero() ;
  //  Shoichi's comments
  // Projection of baryon propagator or vector and scalar 3pt_Function
  //  Tr{ (1 + gamma4) x Baryon_Propagator } for
  //    dd = 0 i.e. decay direction along t
  //  Tr{ (1 - gamma4) x Baryon_Propagator } for
  //    dd = 1 i.e. decay direction along t
  // Projection for axial and tensor 3pt_function
  //  Tr{ (1 + gamma4) x (i*gamma5*gammaZ) x 3pt_Function } for
  //    dd = 2 i.e. decay direction along t
  //  Tr{ (1 - gamma4) x (i*gamma5*gammaZ) x 3pt_Function } for
  //    dd = 3 i.e. decay direction along t
  // Projection for axial and tensor 3pt_function
  //  Tr{ (1 + gamma4) x (gamma5*gammaY) x 3pt_Function } for
  //    dd = 4 i.e. decay direction along t
  //  Tr{ (1 - gamma4) x (gamma5*gammaY) x 3pt_Function } for
  //    dd = 5 i.e. decay direction along t
  // Projection for axial and tensor 3pt_function
  //  Tr{ (1 + gamma4) x (i*gamma5*gammaX) x 3pt_Function } for
  //    dd = 6 i.e. decay direction along t
  //  Tr{ (1 - gamma4) x (i*gamma5*gammaX) x 3pt_Function } for
  //    dd = 7 i.e. decay direction along t
  // Projection for magnetic momment 3pt_function
  //  Tr{ (1 + gamma4) x (i*gammaX*gammaY) x 3pt_Function } for
  //    dd = 8 i.e. decay direction along t
  //  Tr{ (1 - gamma4) x (i*gammaX*gammaY) x 3pt_Function } for
  //    dd = 9 i.e. decay direction along t

  //Comments bellow here have been checked against the chiral
  //basis used in the code. Kostas Orginos
  switch (P) {
  case PPAR: //1+gamma4
    res += q[0][0] ;
    res += q[1][1] ;
    res += q[2][2] ;
    res += q[3][3] ;
    res += q[3][1] ;
    res += q[2][0] ;
    res += q[1][3] ;
    res += q[0][2] ;
    break;
  case NPAR: //1-gamma4
    res +=q[0][0];
    res +=q[1][1];
    res +=q[2][2];
    res +=q[3][3];
    res -=q[3][1];
    res -=q[2][0];
    res -=q[1][3];
    res -=q[0][2];
    break;
  case PPAR_5Z: //(1+gamma4) i gamma5 gammaZ
    res -=q[0][2];
    res +=q[1][3];
    res -=q[2][0];
    res +=q[3][1];
    res -=q[0][0];
    res +=q[1][1];
    res -=q[2][2];
    res +=q[3][3];
    break;
  case NPAR_5Z: //(1-gamma4) i gamma5 gammaZ
    res -=q[0][2];
    res +=q[1][3];
    res -=q[2][0];
    res +=q[3][1];
    res +=q[0][0];
    res -=q[1][1];
    res +=q[2][2];
    res -=q[3][3];
    break;
  case PPAR_5Y: //(1+gamma4)  gamma5 (-gammaY) !! This is what Shoichi did! 
    res +=q[0][3];
    res -=q[1][2];
    res +=q[2][1];
    res -=q[3][0];
    res +=q[0][1];
    res -=q[1][0];
    res +=q[2][3];
    res -=q[3][2];
    break;
  case NPAR_5Y: //(1-gamma4) gamma5 (-gammaY) !! This is what Shoichi did! 
    res +=q[0][3];
    res -=q[1][2];
    res +=q[2][1];
    res -=q[3][0];
    res -=q[0][1];
    res +=q[1][0];
    res -=q[2][3];
    res +=q[3][2];
    break;
  case PPAR_5X: //(1+gamma4) i gamma5 gammaX
    res -=q[0][3];
    res -=q[1][2];
    res -=q[2][1];
    res -=q[3][0];
    res -=q[0][1];
    res -=q[1][0];
    res -=q[2][3];
    res -=q[3][2];
    break;
  case NPAR_5X://(1-gamma4) i gamma5 gammaX
    res -=q[0][3];
    res -=q[1][2];
    res -=q[2][1];
    res -=q[3][0];
    res +=q[0][1];
    res +=q[1][0];
    res +=q[2][3];
    res +=q[3][2];
    break;
  case PPAR_XY: //(1+gamma4) i * gammaX gammaY
    res +=q[0][0];
    res -=q[1][1];
    res +=q[2][2];
    res -=q[3][3];
    res +=q[0][2];
    res +=q[2][0];
    res -=q[1][3];
    res -=q[3][1];
    break;
  case NPAR_XY://(1-gamma4) i * gammaX gammaY
    res +=q[0][0];
    res -=q[1][1];
    res +=q[2][2];
    res -=q[3][3];
    res -=q[0][2];
    res -=q[2][0];
    res +=q[1][3];
    res +=q[3][1];
    break;
  default: 
    ERR.General(cname, fname, "invalid projection");
  }
}



/*! 
   For a given spin and color at the source this calculates:
  \f[
  Diq(\gamma',\gamma)[\delta',a'] = 
  \epsilon^{a',b',c'}\epsilon_{color,b,c}Q_1(\delta',b';spin,b)
                                         Q_2(\gamma',c';\gamma,c) -
  \epsilon^{a',b',c'}\epsilon_{color,b,c}Q_3(\delta',b';\gamma,c)
                                         Q_4(\gamma',c';spin,c)
  \f]
  
  \f$Q_1\f$ is asumed to be the quark propagator multiplied by
  \f$C\gamma_5\f$ left and right

  \f[
  Q_1(\delta',\delta) = (C\gamma_5)_{\delta',\beta'} Q(\beta',\beta)
  (C\gamma_5)_{\beta,\delta}
  \f]

  \f$ Q_2 \f$ is a quark propagator
  \f[
  Q_3(\delta',\beta) =  (C\gamma_5)_{\delta',\beta'} Q(\beta',\beta) 
  Q_4(\beta',\delta) =  Q(\beta',\beta)(C\gamma_5)_{\beta,\delta}  
  \f]
  I use an extra minus sing from that in the notes since the
  ccr (or ccl) has a minus sign wrong in their definition
**/
void Diquark::D_diquark(WilsonMatrix& Q1, WilsonMatrix& Q2, 
			WilsonMatrix& Q3, WilsonMatrix& Q4, 
			int spin, int color)
{

  Epsilon eps1 ; //source
  Epsilon eps2 ; //sink


  Float sign ;
  Complex pp1 ;
  Complex pp2 ;
  
  for(int gamma_p=0;gamma_p<4;gamma_p++)
    for(int gamma=0;gamma<4;gamma++)
      {
	q[gamma_p][gamma].Zero() ;
	// source color contruction (d     in notes)
	for(eps1.Begin(color);eps1.Contracting();eps1.Next(color)) 
	  // sink color contruction (d'     in notes) = eps2.a()
		for(eps2.Begin();eps2.Contracting();eps2.Next()) 
		  {
		    sign = eps1.sign()*eps2.sign() ;
		    // (delta' in notes) 
		    for (int delta_p=0;delta_p<4;delta_p++)
		      {
			//First Term
			pp1=Q1.wmat().d[delta_p].c[eps2.b()].d[spin ].c[eps1.b()] ;
			pp2=Q2.wmat().d[gamma_p].c[eps2.c()].d[gamma].c[eps1.c()] ;
			q[gamma_p][gamma].d[delta_p].c[eps2.a()] += sign*pp1*pp2 ;

			//Second Term
			pp1=Q3.wmat().d[delta_p].c[eps2.b()].d[gamma].c[eps1.c()] ;
			pp2=Q4.wmat().d[gamma_p].c[eps2.c()].d[spin ].c[eps1.b()] ;
			q[gamma_p][gamma].d[delta_p].c[eps2.a()] -= sign*pp1*pp2 ;
		      }
		  }
      }
}



/*! 
  For a given spin and color at the source this calculates:
  \f[
  Diq(\gamma',\gamma)[\delta',a'] = 
  \epsilon^{a',b',c'}\epsilon_{color,b,c}Q_1(\delta',b';spin,b)
                                         Q_2(\gamma',c';\gamma,c) +
  \epsilon^{a',b',c'}\epsilon_{color,b,c}Q_1(\beta',b';\beta,b)
                                         Q_2(\beta',c';\beta,c)
					 \delta^{\delta',\gamma'}
					 \delta^{spin,\gamma} -
  \f]
  \f[
 -\epsilon^{a',b',c'}\epsilon_{color,b,c}Q_1(\beta',c';spin,b)
                                         Q_2(\beta',b';\beta,c)
					 \delta^{\delta',\gamma'} -
  \epsilon^{a',b',c'}\epsilon_{color,b,c}Q_1(\delta',b';\beta,c)
                                         Q_2(\gamma',c';\beta,b)
					 \delta^{spin,\gamma}

  \f]

  \f$Q_1\f$ is asumed to be the quark propagator multiplied by
  \f$C\gamma_5\f$ left and right

  \f[
  Q_1(\delta',\delta) = (C\gamma_5)_{\delta',\beta'} Q(\beta',\beta)
  (C\gamma_5)_{\beta,\delta}
  \f]

  \f$ Q_2 \f$ is a quark propagator

  I use an extra minus sing from that in the notes since the
  ccr (or ccl) has a minus sign wrong in their definition
**/
void Diquark::U_diquark(WilsonMatrix& Q1, WilsonMatrix& Q2, int spin,int color)
{

  Epsilon eps1 ; //source
  Epsilon eps2 ; //sink


  Float sign ;
  Complex pp1 ;
  Complex pp2 ;
  
  for(int gamma_p=0;gamma_p<4;gamma_p++)
    for(int gamma=0;gamma<4;gamma++)
      {
	q[gamma_p][gamma].Zero() ;
	// source color contruction (d     in notes)
	for(eps1.Begin(color);eps1.Contracting();eps1.Next(color)) 
	  // sink color contruction (d'     in notes) = eps2.a()
		for(eps2.Begin();eps2.Contracting();eps2.Next()) 
		  {
		    sign = eps1.sign()*eps2.sign() ;
		    // (delta' in notes) 
		    for (int delta_p=0;delta_p<4;delta_p++)
		      {
			//First Term
			pp1=Q1.wmat().d[delta_p].c[eps2.b()].d[spin ].c[eps1.b()] ;
			pp2=Q2.wmat().d[gamma_p].c[eps2.c()].d[gamma].c[eps1.c()] ;
			q[gamma_p][gamma].d[delta_p].c[eps2.a()] += sign*pp1*pp2 ;

			//Second Term
			if((gamma_p == delta_p)&&(gamma == spin))
			  {
			    for(int beta_p=0;beta_p<4;beta_p++)
			      for(int beta=0;beta<4;beta++)
				{
				  pp1=Q1.wmat().d[beta_p].c[eps2.b()].d[beta].c[eps1.b()] ;
				  pp2=Q2.wmat().d[beta_p].c[eps2.c()].d[beta].c[eps1.c()] ;
				  q[gamma_p][gamma].d[delta_p].c[eps2.a()] += sign*pp1*pp2 ;
				}
			  }

			//Third term
			if(gamma_p == delta_p)
			  {
			    for(int beta_p=0;beta_p<4;beta_p++)
			      {
				pp1=Q1.wmat().d[beta_p].c[eps2.c()].d[spin ].c[eps1.b()] ;
				pp2=Q2.wmat().d[beta_p].c[eps2.b()].d[gamma].c[eps1.c()] ;
				q[gamma_p][gamma].d[delta_p].c[eps2.a()] -= sign*pp1*pp2 ;
			      }
			  }

			//Forth term
			if(gamma == spin)
			  {
			    for(int beta=0;beta<4;beta++)
			      {
				pp1=Q1.wmat().d[delta_p].c[eps2.b()].d[beta].c[eps1.c()] ;
				pp2=Q2.wmat().d[gamma_p].c[eps2.c()].d[beta].c[eps1.b()] ;
				q[gamma_p][gamma].d[delta_p].c[eps2.a()] -= sign*pp1*pp2 ;
			      }
			  }

		      }
		  }
      }
}

CPS_END_NAMESPACE
