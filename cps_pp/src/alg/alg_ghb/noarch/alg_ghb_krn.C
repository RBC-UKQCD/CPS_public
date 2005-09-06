#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
/*!\file
  \brief Routines used by the AlgGheatBath class methods:

  $Id: alg_ghb_krn.C,v 1.9 2005-09-06 20:33:36 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2005-09-06 20:33:36 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_ghb/noarch/alg_ghb_krn.C,v 1.9 2005-09-06 20:33:36 chulwoo Exp $
//  $Id: alg_ghb_krn.C,v 1.9 2005-09-06 20:33:36 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: alg_ghb_krn.C,v $
//  $Revision: 1.9 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_ghb/noarch/alg_ghb_krn.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// alg_ghb.C
//
// AlgGheatBath is derived from Alg and is relevant to the 
// gauge heat bath algorithm. The type of gauge action is
// determined by the argument to the constructor.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <stdlib.h>	// exit()
#include <util/qcdio.h>
#include <math.h>
#include<alg/alg_ghb.h>
#include<alg/common_arg.h>
#include<alg/ghb_arg.h>
#include<util/lattice.h>
#include<util/gjp.h>
#include<util/random.h>
#include<util/smalloc.h>
#include<util/vector.h>
#include<util/verbose.h>
#include<util/error.h>
CPS_START_NAMESPACE

void m_multiply2r( Float *AB, Float *B );
void m_multiply2l( Float *AB, Float *A );
void m_multiply3( Float *AB, Float *A, Float *B );
void m_add( Float *AplusB, Float *A, Float *B );
void m_equal( Float *A, Float *B );
void m_identity( Float *A );
void m_conjugate( Float *A );
void m_invert( Float *matrix );
void m_rand( Float* eps, Float squeeze );
Float m_determinantR( Float* mtx );
Float m_determinantI( Float* mtx );
Float absR( Float x );
void m_subtract( Float* AminusB, Float* A, Float* B );
void m_zero( Float* x );
void m_print( Float *A );


static const int MATRIXSIZE = 18;

static int	NHITS = 10;
static Float	SMALL = 0.2;
static Float	EPSILON[MATRIXSIZE];
static Float	MTEMP1[MATRIXSIZE];
static Float	MTEMP2[MATRIXSIZE];

//*******************************************************
//     Here is a Metropolis update kernel...  		*
//*******************************************************
void metropolis_kernel(  Float *sigma, Float *U)
{
//  LRG.SetInterval(1,-1);

VRB.Debug(" Beta = %f\n", GJP.Beta() );
  Float fGamma = (2./3.) * GJP.Beta();
 
  // These variables are associated with the actual update. 
  Float  old_action;
  Float  new_action;
  Float  accept_probability;
  Float  die_roll;
 
VRB.Debug("Entering update core metro new with nhits=%d and small=%f.\n",
 NHITS, SMALL );

  for(int HIT_CNT=0; HIT_CNT<NHITS; HIT_CNT++)
  {

    {
      Float a,b,c,d;
      Float norm;
      Float select;
      Float* v;

      v = EPSILON;
    
      a = SMALL * LRG.Urand(); if(a>0)a*=-1; a+=1.0;
      b = LRG.Urand();
      c = LRG.Urand();
      d = LRG.Urand();
    

VRB.Debug("unnormalized: a=%+6.5f b=%+6.5f c=%+6.5f d=%+6.5f\n",a,b,c,d);

      norm = sqrt( (1.0 - a*a) / (b*b + c*c + d*d) );
      b *= norm; c *= norm; d *= norm;

VRB.Debug("normalized:   a=%+6.5f b=%+6.5f c=%+6.5f d=%+6.5f\n",a,b,c,d);

    
      select = LRG.Urand();
    
      if ( select < 0.33333333 ) {
        *(v++)=a  ; *(v++)=c  ; *(v++)=0.0;
        *(v++)=-c ; *(v++)=a  ; *(v++)=0.0;
        *(v++)=0.0; *(v++)=0.0; *(v++)=1.0;

        *(v++)=b  ; *(v++)=d  ; *(v++)=0.0;
        *(v++)=d  ; *(v++)=-b ; *(v++)=0.0;
        *(v++)=0.0; *(v++)=0.0; *(v++)=0.0;
      }
      else
        if ( select > 0.66666666 ) {
    
          *(v++)=1.0; *(v++)=0.0; *(v++)=0.0;
          *(v++)=0.0; *(v++)=a  ; *(v++)=c  ;
          *(v++)=0.0; *(v++)=-c ; *(v++)=a  ;
 
          *(v++)=0.0; *(v++)=0.0; *(v++)=0.0;
          *(v++)=0.0; *(v++)=b  ; *(v++)=d  ;
          *(v++)=0.0; *(v++)=d  ; *(v++)=-b ;
        }
        else {
          *(v++)=a  ; *(v++)=0.0; *(v++)=c  ;
          *(v++)=0.0; *(v++)=1.0; *(v++)=0.0;
          *(v++)=-c ; *(v++)=0.0; *(v++)=a  ;

          *(v++)=b  ; *(v++)=0.0; *(v++)=d  ;
          *(v++)=0.0; *(v++)=0.0; *(v++)=0.0;
          *(v++)=d  ; *(v++)=0.0; *(v++)=-b ;
        }

    }

VRB.Debug("epsilon matrix:\n"); m_print( EPSILON );

      
    // Now begin the update process. 

    m_multiply3( MTEMP1, sigma, U );
    old_action = fGamma * (9. - .5 * (MTEMP1[0]+MTEMP1[4]+MTEMP1[8]) );

VRB.Debug("Old action %e\n", old_action );


    m_multiply3( MTEMP1, EPSILON, U );
    m_multiply3( MTEMP2, sigma, MTEMP1 );
    new_action = fGamma * (9. - .5 * (MTEMP2[0]+MTEMP2[4]+MTEMP2[8]) );

VRB.Debug("New action %e\n", new_action );


    accept_probability = exp( old_action - new_action );

VRB.Debug("Accept chance %e\n", accept_probability );


    die_roll = LRG.Urand();
    if( die_roll < 0 ) die_roll *= -1;

VRB.Debug("die roll %e\n", die_roll );


    if( die_roll < accept_probability ) {
VRB.Debug("ACCEPT\n\n");
      m_multiply3( MTEMP1, EPSILON, U );
      m_equal( U, MTEMP1 );
    }
else
VRB.Debug("REJECT\n\n");

  }             /* End HIT_CNT loop for multiple hits.           */
 
}

//*******************************************************
//       Here is a Heat Bath update kernel...   	
// This routine is based on the algorithm presented by Cabbibo and
// Marrinari for updating a SU(n) link element by seperately
// updating SU(2) subgroups of that matrix using for each the
// SU(2) heat bath algorithm of creutz. 
//*******************************************************
void  cmhb_kernel( Float *sigma, Float *u)
{
//    LRG.SetInterval(1, -1);

#if 1    
    // This tells which subblock is being updated: 
  int	subblock=0;

    // Here is the product of the link to be updated and the environment,
    // which is the sum of the six plquettes which contain u. 
  Float   u_sigma[MATRIXSIZE];

    // Choosing an SU(2) subblock of u_sigma, it can be expressed
    // in terms of four real numbers r[0] through r[3], where 
    // R = 1*r[0] + I*r_vector DOT sigma_vector,
    // where sigma_vector refers to the three pauli matrices. 
  Float  r[4];

    // The subblock of u_sigma will not in general be SU(2), it will
    // have to be seperated into a real SU(2) part and an imaginary
    // SU(2), both needed to be rescaled to be SU(2).  k is the
    // rescaling parameter used there. 
  Float  k;
  Float  alp;
  Float	R, R_, R__, R___, X, X_;
  Float  C, A, delta;

    // A boltzman distributed vector which represents an SU(2) matrix. 
  Float  bd[4];

    // Variables required to get a random 3-vector on a sphere. 
  Float  phi, sin_theta, cos_theta, sin_phi;

    // and a normalization for those vectors 
  Float  scale;

    // Once a boltzman distributed SU(2) matrix has been constructed, 
    // it is necessary to imbed it into an SU(3) matrix (unit on the last 
    // diagonal)
    // so that it can be multiplied by the origional link. 
  Float  alpha[MATRIXSIZE];
 
  Float  mtx_r[MATRIXSIZE];
  Float  mtx_bd[MATRIXSIZE];
 
  for(subblock=0; subblock<=1; subblock++)
  {
      // Construct the sum of plaquettes containing the link "u". 
    m_multiply3( u_sigma, u, sigma );

      // Sort out the first SU2 subblock of that sum matrix and 
      // taking the upper 2x2 subblock of u_sigma, it can be expressed
      // as:    k_1 ( r[0] + r_vector DOT sigma_vector )
      //      + k_2 ( s[0] + s_vector DOT sigma_vector )
      // extract the SU2 part called r and the coefficiant k. 
    if( subblock==0 )
    {
      r[0] = ( u_sigma[0]  + u_sigma[4]  ) / 2.;
      r[1] = ( u_sigma[10] + u_sigma[12] ) / 2.;
      r[2] = ( u_sigma[1]  - u_sigma[3]  ) / 2.;
      r[3] = ( u_sigma[9]  - u_sigma[13] ) / 2.;
    }
    else 
    {
      r[0] = ( u_sigma[4]  + u_sigma[8]  ) / 2.;
      r[1] = ( u_sigma[14] + u_sigma[16] ) / 2.;
      r[2] = ( u_sigma[5]  - u_sigma[7]  ) / 2.;
      r[3] = ( u_sigma[13] - u_sigma[17] ) / 2.;
    }

    k = sqrt( r[0]*r[0] + r[1]*r[1] + r[2]*r[2] + r[3]*r[3] );

    r[0] = r[0]/k;
    r[1] = r[1]/k;
    r[2] = r[2]/k;
    r[3] = r[3]/k;
 
//    Float lambda = (2. * GJP.Beta() * k) / 3.;
    // The 3 is the three of SU(3). 
   
    //  Generate the zero-th component of the boltzman distributed
    //  using the Kennedy Pendleton algorithm: Phys Lett 156B (393).  

    int good_z = 0;
    while( good_z == 0 )
    {
      Float	my_pi = 3.14159265358979323846;

      alp  = (2./3.)* GJP.Beta() * k;

      R__  =  absR(LRG.Urand() );
      R_   =  absR(LRG.Urand() );
      R    =  absR(LRG.Urand() );
      R___ =  absR(LRG.Urand() );

      X  = -(log(R ))/alp ;
      X_ = -(log(R_))/alp ;
      C  = cos( 2*my_pi*R__ )*cos( 2*my_pi*R__ );
      A  = X * C;
      delta = X_ + A;
      if( (R___*R___) <= (1.-delta/2.) ) good_z = 1;
    }
/* ---------  For Debugging purposes
    int good_z = 0;
    int counter = 0;
    while( good_z == 0 )
    {
      Float     my_pi = 3.14159265358979323846;

      alp  = (2./3.)* GJP.Beta() * k;

      R__  =  absR(0.553422 );
      R_   =  absR(-0.098343 );
      R    =  absR(-0.873420 );
      R___ =  absR(0.204563 );
      X  = -(log(R ))/alp ;
      X_ = -(log(R_))/alp ;
      C  = cos( 2*my_pi*R__ )*cos( 2*my_pi*R__ );
      A  = X * C;
      delta = X_ + A;

      if(counter++ == 3) good_z = 1;
    }
// -----------------------------------
*/

    bd[0] = 1 - delta;

    //  Generate points on the surface of a sphere by directly
    // generating theta and phi with the apropriate weighting. 

    scale = sqrt( 1. - bd[0]*bd[0] );
 
    // Generate the three components of bd on a unit sphere. 

    cos_theta = LRG.Urand();
    Float	my_pi = 3.14159265358979323846;
    phi = LRG.Urand()*my_pi;   // phi e [-pi,pi) 

/* 
//--------------  For Debugging Purposes
    cos_theta = -0.912342;
    Float       my_pi = 3.14159265358979323846;
    phi = 0.243134*my_pi;   // phi e [-pi,pi)
// ---------------------------------------------
*/ 

    sin_theta = sqrt( 1. - cos_theta*cos_theta );  // sin(Theta) > 0 
    sin_phi =   sqrt( 1. - cos(phi)*cos(phi) );
    if( phi < 0 ) sin_phi *= -1;

    // x = sin(Theta) x cos(Phi) 
    bd[1] = sin_theta * cos(phi);

    // y = sin(Theta) x sin(Phi) 
    bd[2] = sin_theta * sin_phi;

    // z = cos(Theta) 
    bd[3] = cos_theta;
 
    bd[1]*=scale;
    bd[2]*=scale;
    bd[3]*=scale;

    m_identity( mtx_r );
    if(subblock==0)
    {
      *(mtx_r + 0)  = r[0]; mtx_r[9]  = r[3];
      *(mtx_r + 1)  = r[2]; mtx_r[10] = r[1];
      *(mtx_r + 3)  =-r[2]; mtx_r[12] = r[1];
      *(mtx_r + 4)  = r[0]; mtx_r[13] =-r[3];
    }
    else
    {
      *(mtx_r + 4)  = r[0]; mtx_r[13] = r[3];
      *(mtx_r + 5)  = r[2]; mtx_r[14] = r[1];
      *(mtx_r + 7)  =-r[2]; mtx_r[16] = r[1];
      *(mtx_r + 8)  = r[0]; mtx_r[17] =-r[3];
    }
    m_invert( mtx_r );

    m_identity( mtx_bd );

    if(subblock==0)
    {
      *(mtx_bd + 0)  = bd[0]; mtx_bd[9]  = bd[3];
      *(mtx_bd + 1)  = bd[2]; mtx_bd[10] = bd[1];
      *(mtx_bd + 3)  =-bd[2]; mtx_bd[12] = bd[1];
      *(mtx_bd + 4)  = bd[0]; mtx_bd[13] =-bd[3];
    }
    else
    {
      *(mtx_bd + 4)  = bd[0]; mtx_bd[13] = bd[3];
      *(mtx_bd + 5)  = bd[2]; mtx_bd[14] = bd[1];
      *(mtx_bd + 7)  =-bd[2]; mtx_bd[16] = bd[1];
      *(mtx_bd + 8)  = bd[0]; mtx_bd[17] =-bd[3];
    }

    m_multiply3( alpha, mtx_bd, mtx_r );

    // Now that the update factor is computed multiply U by that
    // factor.  
    m_multiply2l( u, alpha );

  }
#endif  
}

CPS_END_NAMESPACE
