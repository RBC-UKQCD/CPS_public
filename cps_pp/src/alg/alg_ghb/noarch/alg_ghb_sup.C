#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
/*!\file
  \brief Routines used by the AlgGheatBath class methods:

  $Id: alg_ghb_sup.C,v 1.6 2004-08-18 11:57:38 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-08-18 11:57:38 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_ghb/noarch/alg_ghb_sup.C,v 1.6 2004-08-18 11:57:38 zs Exp $
//  $Id: alg_ghb_sup.C,v 1.6 2004-08-18 11:57:38 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: alg_ghb_sup.C,v $
//  $Revision: 1.6 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_ghb/noarch/alg_ghb_sup.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// alg_ghb_support.C
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



void m_multiply2r( Float *AB, Float *B )
/* Two operand multiply with implicit right multiplication,
** analogous to AB *= B, meaning AB = AB * B. */
{
  Float	scratch1[18];
  m_equal( scratch1, AB );
  m_multiply3( AB, scratch1, B );
}

void m_multiply2l( Float *AB, Float *A )
/* Two operand multiply with implicit left multiplication,
** analogous to AB *= A, meaning AB = A * AB. */
{
  Float	scratch1[18];
  m_equal( scratch1, AB );
  m_multiply3( AB, A, scratch1 );
}

void m_multiply3( Float *AB, Float *A, Float *B )
/* Three operand multiply for matrices, analogous to C=A*B */
{
  /* First multiply the components. */
  *(AB + 0)  = (*(A + 0) * *(B + 0) ) + (*(A + 1) * *(B + 3) )
             + (*(A + 2) * *(B + 6) ) - (*(A + 9) * *(B + 9) )
             - (*(A + 10)* *(B + 12)) - (*(A + 11)* *(B + 15));
  *(AB + 1)  = (*(A + 0) * *(B + 1) ) + (*(A + 1) * *(B + 4) )
             + (*(A + 2) * *(B + 7) ) - (*(A + 9) * *(B + 10))
             - (*(A + 10)* *(B + 13)) - (*(A + 11)* *(B + 16));
  *(AB + 2)  = (*(A + 0) * *(B + 2) ) + (*(A + 1) * *(B + 5) )
             + (*(A + 2) * *(B + 8) ) - (*(A + 9) * *(B + 11)) 
             - (*(A + 10)* *(B + 14)) - (*(A + 11)* *(B + 17));
  *(AB + 3)  = (*(A + 3) * *(B + 0) ) + (*(A + 4) * *(B + 3) )
             + (*(A + 5) * *(B + 6) ) - (*(A + 12)* *(B + 9) )
             - (*(A + 13)* *(B + 12)) - (*(A + 14)* *(B + 15));
  *(AB + 4)  = (*(A + 3) * *(B + 1) ) + (*(A + 4) * *(B + 4) )
             + (*(A + 5) * *(B + 7) ) - (*(A + 12)* *(B + 10))
             - (*(A + 13)* *(B + 13)) - (*(A + 14)* *(B + 16));
  *(AB + 5)  = (*(A + 3) * *(B + 2) ) + (*(A + 4) * *(B + 5) )
             + (*(A + 5) * *(B + 8) ) - (*(A + 12)* *(B + 11))
             - (*(A + 13)* *(B + 14)) - (*(A + 14)* *(B + 17));
  *(AB + 6)  = (*(A + 6) * *(B + 0) ) + (*(A + 7) * *(B + 3) )
             + (*(A + 8) * *(B + 6) ) - (*(A + 15)* *(B + 9) )
             - (*(A + 16)* *(B + 12)) - (*(A + 17)* *(B + 15));
  *(AB + 7)  = (*(A + 6) * *(B + 1) ) + (*(A + 7) * *(B + 4) )
             + (*(A + 8) * *(B + 7) ) - (*(A + 15)* *(B + 10))
             - (*(A + 16)* *(B + 13)) - (*(A + 17)* *(B + 16));
  *(AB + 8)  = (*(A + 6) * *(B + 2) ) + (*(A + 7) * *(B + 5) )
             + (*(A + 8) * *(B + 8) ) - (*(A + 15)* *(B + 11))
             - (*(A + 16)* *(B + 14)) - (*(A + 17)* *(B + 17));
 
  *(AB + 9)  = (*(A + 0) * *(B + 9) ) + (*(A + 1) * *(B + 12))
             + (*(A + 2) * *(B + 15)) + (*(A + 9) * *(B + 0) )
             + (*(A + 10)* *(B + 3) ) + (*(A + 11)* *(B + 6) );
  *(AB + 10) = (*(A + 0) * *(B + 10)) + (*(A + 1) * *(B + 13))
             + (*(A + 2) * *(B + 16)) + (*(A + 9) * *(B + 1) )
             + (*(A + 10)* *(B + 4) ) + (*(A + 11)* *(B + 7) );
  *(AB + 11) = (*(A + 0) * *(B + 11)) + (*(A + 1) * *(B + 14))
             + (*(A + 2) * *(B + 17)) + (*(A + 9) * *(B + 2) )
             + (*(A + 10)* *(B + 5) ) + (*(A + 11)* *(B + 8) );
  *(AB + 12) = (*(A + 3) * *(B + 9) ) + (*(A + 4) * *(B + 12))
             + (*(A + 5) * *(B + 15)) + (*(A + 12)* *(B + 0) )
             + (*(A + 13)* *(B + 3) ) + (*(A + 14)* *(B + 6) );
  *(AB + 13) = (*(A + 3) * *(B + 10)) + (*(A + 4) * *(B + 13))
             + (*(A + 5) * *(B + 16)) + (*(A + 12)* *(B + 1) )
             + (*(A + 13)* *(B + 4) ) + (*(A + 14)* *(B + 7) );
  *(AB + 14) = (*(A + 3) * *(B + 11)) + (*(A + 4) * *(B + 14))
             + (*(A + 5) * *(B + 17)) + (*(A + 12)* *(B + 2) )
             + (*(A + 13)* *(B + 5) ) + (*(A + 14)* *(B + 8) );
  *(AB + 15) = (*(A + 6) * *(B + 9) ) + (*(A + 7) * *(B + 12))
             + (*(A + 8) * *(B + 15)) + (*(A + 15)* *(B + 0) )
             + (*(A + 16)* *(B + 3) ) + (*(A + 17)* *(B + 6) );
  *(AB + 16) = (*(A + 6) * *(B + 10)) + (*(A + 7) * *(B + 13))
             + (*(A + 8) * *(B + 16)) + (*(A + 15)* *(B + 1) )
             + (*(A + 16)* *(B + 4) ) + (*(A + 17)* *(B + 7) );
  *(AB + 17) = (*(A + 6) * *(B + 11)) + (*(A + 7) * *(B + 14))
             + (*(A + 8) * *(B + 17)) + (*(A + 15)* *(B + 2) )
             + (*(A + 16)* *(B + 5) ) + (*(A + 17)* *(B + 8) );
}
  
void  m_add( Float *AplusB, Float *A, Float *B )
{

  *(AplusB + 0)  =  *(A + 0)  + *(B + 0) ;
  *(AplusB + 9)  =  *(A + 9)  + *(B + 9) ;
  *(AplusB + 1)  =  *(A + 1)  + *(B + 1) ;
  *(AplusB + 10) =  *(A + 10) + *(B + 10);
  *(AplusB + 2)  =  *(A + 2)  + *(B + 2) ;
  *(AplusB + 11) =  *(A + 11) + *(B + 11);
  *(AplusB + 3)  =  *(A + 3)  + *(B + 3) ;
  *(AplusB + 12) =  *(A + 12) + *(B + 12);
  *(AplusB + 4)  =  *(A + 4)  + *(B + 4) ;
  *(AplusB + 13) =  *(A + 13) + *(B + 13);
  *(AplusB + 5)  =  *(A + 5)  + *(B + 5) ;
  *(AplusB + 14) =  *(A + 14) + *(B + 14);
  *(AplusB + 6)  =  *(A + 6)  + *(B + 6) ;
  *(AplusB + 15) =  *(A + 15) + *(B + 15);
  *(AplusB + 7)  =  *(A + 7)  + *(B + 7) ;
  *(AplusB + 16) =  *(A + 16) + *(B + 16);
  *(AplusB + 8)  =  *(A + 8)  + *(B + 8) ;
  *(AplusB + 17) =  *(A + 17) + *(B + 17);

}

void  m_equal( Float *A, Float *B )
{
  *(A + 0)  = *(B + 0) ;
  *(A + 9)  = *(B + 9) ;
  *(A + 1)  = *(B + 1) ;
  *(A + 10) = *(B + 10);
  *(A + 2)  = *(B + 2) ;
  *(A + 11) = *(B + 11);
  *(A + 3)  = *(B + 3) ;
  *(A + 12) = *(B + 12);
  *(A + 4)  = *(B + 4) ;
  *(A + 13) = *(B + 13);
  *(A + 5)  = *(B + 5) ;
  *(A + 14) = *(B + 14);
  *(A + 6)  = *(B + 6) ;
  *(A + 15) = *(B + 15);
  *(A + 7)  = *(B + 7) ;
  *(A + 16) = *(B + 16);
  *(A + 8)  = *(B + 8) ;
  *(A + 17) = *(B + 17);
}

void m_identity( Float *A )
{
  *(A + 0)  = 1.; *(A + 9)  = 0.;
  *(A + 1)  = 0.; *(A + 10) = 0.;
  *(A + 2)  = 0.; *(A + 11) = 0.;
  *(A + 3)  = 0.; *(A + 12) = 0.;
  *(A + 4)  = 1.; *(A + 13) = 0.;
  *(A + 5)  = 0.; *(A + 14) = 0.;
  *(A + 6)  = 0.; *(A + 15) = 0.;
  *(A + 7)  = 0.; *(A + 16) = 0.;
  *(A + 8)  = 1.; *(A + 17) = 0.;
}

void  m_conjugate( Float *A )
{
  Float	scratch1[18];
  m_equal( scratch1, A );

  *(A + 0)  = + *(scratch1 + 0) ;
  *(A + 9)  = - *(scratch1 + 9) ;
  *(A + 1)  = + *(scratch1 + 3) ;
  *(A + 10) = - *(scratch1 + 12);
  *(A + 2)  = + *(scratch1 + 6) ;
  *(A + 11) = - *(scratch1 + 15);
  *(A + 3)  = + *(scratch1 + 1) ;
  *(A + 12) = - *(scratch1 + 10);
  *(A + 4)  = + *(scratch1 + 4) ;
  *(A + 13) = - *(scratch1 + 13);
  *(A + 5)  = + *(scratch1 + 7) ;
  *(A + 14) = - *(scratch1 + 16);
  *(A + 6)  = + *(scratch1 + 2) ;
  *(A + 15) = - *(scratch1 + 11);
  *(A + 7)  = + *(scratch1 + 5) ;
  *(A + 16) = - *(scratch1 + 14);
  *(A + 8)  = + *(scratch1 + 8) ;
  *(A + 17) = - *(scratch1 + 17);
}

void m_invert( Float *matrix )
{
  Float	scratch1[18];
  Float	scratch2[18];

  /* Define a cofactor matrix. */
  m_identity( scratch1 );

  /* Define a temp to use for the inverse. */
  m_identity( scratch2 );
 
  /* Create a scratch1actor matrix. */
  *(scratch1 + 0)  = (*(matrix + 4) * *(matrix + 8) -
		    *(matrix + 13)* *(matrix + 17))
                   - (*(matrix + 5) * *(matrix + 7) -
		    *(matrix + 14)* *(matrix + 16));
  *(scratch1 + 9)  = (*(matrix + 4) * *(matrix + 17)+
		    *(matrix + 13)* *(matrix + 8) )
                   - (*(matrix + 5) * *(matrix + 16)+
		    *(matrix + 14)* *(matrix + 7) );
  *(scratch1 + 1)  = (*(matrix + 3) * *(matrix + 8) -
		    *(matrix + 12)* *(matrix + 17))
                   - (*(matrix + 5) * *(matrix + 6) -
		    *(matrix + 14)* *(matrix + 15));
  *(scratch1 + 10) = (*(matrix + 3) * *(matrix + 17)+
		    *(matrix + 12)* *(matrix + 8) )
                   - (*(matrix + 5) * *(matrix + 15)+
		    *(matrix + 14)* *(matrix + 6) );
  *(scratch1 + 2)  = (*(matrix + 3) * *(matrix + 7) -
		    *(matrix + 12)* *(matrix + 16))
                   - (*(matrix + 4) * *(matrix + 6) -
		    *(matrix + 13)* *(matrix + 15));
  *(scratch1 + 11) = (*(matrix + 3) * *(matrix + 16)+
		    *(matrix + 12)* *(matrix + 7) )
                   - (*(matrix + 4) * *(matrix + 15)+
		    *(matrix + 13)* *(matrix + 6) );
  *(scratch1 + 3)  = (*(matrix + 1) * *(matrix + 8) -
		    *(matrix + 10)* *(matrix + 17))
                   - (*(matrix + 2) * *(matrix + 7) -
		    *(matrix + 11)* *(matrix + 16));
  *(scratch1 + 12) = (*(matrix + 1) * *(matrix + 17)+
		    *(matrix + 10)* *(matrix + 8) )
                   - (*(matrix + 2) * *(matrix + 16)+
		    *(matrix + 11)* *(matrix + 7) );
  *(scratch1 + 4)  = (*(matrix + 0) * *(matrix + 8) -
		    *(matrix + 9) * *(matrix + 17))
                   - (*(matrix + 2) * *(matrix + 6) -
		    *(matrix + 11)* *(matrix + 15));
  *(scratch1 + 13) = (*(matrix + 0) * *(matrix + 17)+
		    *(matrix + 9) * *(matrix + 8) )
                   - (*(matrix + 2) * *(matrix + 15)+
		    *(matrix + 11)* *(matrix + 6) );
  *(scratch1 + 5)  = (*(matrix + 0) * *(matrix + 7) -
		    *(matrix + 9) * *(matrix + 16))
                   - (*(matrix + 1) * *(matrix + 6) -
		    *(matrix + 10)* *(matrix + 15));
  *(scratch1 + 14) = (*(matrix + 0) * *(matrix + 16)+
		    *(matrix + 9) * *(matrix + 7) )
                   - (*(matrix + 1) * *(matrix + 15)+
		    *(matrix + 10)* *(matrix + 6) );
  *(scratch1 + 6)  = (*(matrix + 1) * *(matrix + 5) -
		    *(matrix + 10)* *(matrix + 14))
                   - (*(matrix + 2) * *(matrix + 4) -
		    *(matrix + 11)* *(matrix + 13));
  *(scratch1 + 15) = (*(matrix + 1) * *(matrix + 14)+
		    *(matrix + 10)* *(matrix + 5) )
                   - (*(matrix + 2) * *(matrix + 13)+
		    *(matrix + 11)* *(matrix + 4) );
  *(scratch1 + 7)  = (*(matrix + 0) * *(matrix + 5) -
		    *(matrix + 9) * *(matrix + 14))
                   - (*(matrix + 2) * *(matrix + 3) -
		    *(matrix + 11)* *(matrix + 12));
  *(scratch1 + 16) = (*(matrix + 0) * *(matrix + 14)+
		    *(matrix + 9) * *(matrix + 5) )
                   - (*(matrix + 2) * *(matrix + 12)+
		    *(matrix + 11)* *(matrix + 3) );
  *(scratch1 + 8)  = (*(matrix + 0) * *(matrix + 4) -
		    *(matrix + 9) * *(matrix + 13))
                   - (*(matrix + 1) * *(matrix + 3) -
		    *(matrix + 10)* *(matrix + 12));
  *(scratch1 + 17) = (*(matrix + 0) * *(matrix + 13)+
		    *(matrix + 9) * *(matrix + 4) )
                   - (*(matrix + 1) * *(matrix + 12)+
		    *(matrix + 10)* *(matrix + 3) );
 
  *(scratch1 + 1)  = - *(scratch1 + 1) ;
  *(scratch1 + 10) = - *(scratch1 + 10);
  *(scratch1 + 3)  = - *(scratch1 + 3) ;
  *(scratch1 + 12) = - *(scratch1 + 12);
  *(scratch1 + 5)  = - *(scratch1 + 5) ;
  *(scratch1 + 14) = - *(scratch1 + 14);
  *(scratch1 + 7)  = - *(scratch1 + 7) ;
  *(scratch1 + 16) = - *(scratch1 + 16);
  
  *(scratch2 + 0)  = *(scratch1 + 0) ;
  *(scratch2 + 9)  = *(scratch1 + 9) ;  
  *(scratch2 + 1)  = *(scratch1 + 3) ;
  *(scratch2 + 10) = *(scratch1 + 12);  
  *(scratch2 + 2)  = *(scratch1 + 6) ;
  *(scratch2 + 11) = *(scratch1 + 15);  
  *(scratch2 + 3)  = *(scratch1 + 1) ;
  *(scratch2 + 12) = *(scratch1 + 10);  
  *(scratch2 + 4)  = *(scratch1 + 4) ;
  *(scratch2 + 13) = *(scratch1 + 13);  
  *(scratch2 + 5)  = *(scratch1 + 7) ;
  *(scratch2 + 14) = *(scratch1 + 16);  
  *(scratch2 + 6)  = *(scratch1 + 2) ;
  *(scratch2 + 15) = *(scratch1 + 11);  
  *(scratch2 + 7)  = *(scratch1 + 5) ;
  *(scratch2 + 16) = *(scratch1 + 14);  
  *(scratch2 + 8)  = *(scratch1 + 8) ;
  *(scratch2 + 17) = *(scratch1 + 17);  
 
  m_equal( matrix, scratch2 );
}

void  m_rand( Float* eps, Float squeeze)
{
#define	R	0
#define	I	1
  LRG.SetInterval(-1, 1);

  Float  mag_squared;
  Float  AdotBprime[2];
  Float  aA[3][2], A[3][2], bB[3][2], B[3][2], Bprime[3][2], C[3][2];

  m_identity( eps );

  /* Generate A near i-hat. */
  aA[0][R] = 1-squeeze*(1.+LRG.Urand());
  aA[0][I] = squeeze*LRG.Urand();
  aA[1][R] = squeeze*LRG.Urand();
  aA[1][I] = squeeze*LRG.Urand();
  aA[2][R] = squeeze*LRG.Urand();
  aA[2][I] = squeeze*LRG.Urand();
 
  /* Normalize A. */
  mag_squared = aA[0][R]*aA[0][R] + aA[0][I]*aA[0][I]
              + aA[1][R]*aA[1][R] + aA[1][I]*aA[1][I]
              + aA[2][R]*aA[2][R] + aA[2][I]*aA[2][I];
  A[0][R] = aA[0][R] / sqrt( mag_squared );
  A[0][I] = aA[0][I] / sqrt( mag_squared );
  A[1][R] = aA[1][R] / sqrt( mag_squared );
  A[1][I] = aA[1][I] / sqrt( mag_squared );
  A[2][R] = aA[2][R] / sqrt( mag_squared );
  A[2][I] = aA[2][I] / sqrt( mag_squared );

  /* Generate Bprime near j-hat. */
  Bprime[0][R] = squeeze*LRG.Urand();
  Bprime[0][I] = squeeze*LRG.Urand();
  Bprime[1][R] = 1-squeeze*(1.+LRG.Urand());
  Bprime[1][I] = squeeze*LRG.Urand();
  Bprime[2][R] = squeeze*LRG.Urand();
  Bprime[2][I] = squeeze*LRG.Urand();

  /* Now calculate B which is defined by 
  ** Bprime = aA + bB by taking A.Bprime.
  ** This gives:
  ** AdotBprime = a(AdotA) +b(AdotB) = a,
  ** So:
  ** bB = Bprime - (AdotBprime)A      */
  AdotBprime[R] = 
    A[0][R]*Bprime[0][R] - (- A[0][I])*Bprime[0][I]
  + A[1][R]*Bprime[1][R] - (- A[1][I])*Bprime[1][I]
  + A[2][R]*Bprime[2][R] - (- A[2][I])*Bprime[2][I];
  AdotBprime[I] = 
    (- A[0][I])*Bprime[0][R] + A[0][R]*Bprime[0][I]
  + (- A[1][I])*Bprime[1][R] + A[1][R]*Bprime[1][I]
  + (- A[2][I])*Bprime[2][R] + A[2][R]*Bprime[2][I];

  bB[0][R] = Bprime[0][R]-(AdotBprime[R]*A[0][R]-AdotBprime[I]*A[0][I]);
  bB[0][I] = Bprime[0][I]-(AdotBprime[I]*A[0][R]+AdotBprime[R]*A[0][I]);
  bB[1][R] = Bprime[1][R]-(AdotBprime[R]*A[1][R]-AdotBprime[I]*A[1][I]);
  bB[1][I] = Bprime[1][I]-(AdotBprime[I]*A[1][R]+AdotBprime[R]*A[1][I]);
  bB[2][R] = Bprime[2][R]-(AdotBprime[R]*A[2][R]-AdotBprime[I]*A[2][I]);
  bB[2][I] = Bprime[2][I]-(AdotBprime[I]*A[2][R]+AdotBprime[R]*A[2][I]);

  // Now that we have bB, I should calculate b and devide it out.
  //  To this end, calculate BdotB which is |b|^2 
  mag_squared = bB[0][R]*bB[0][R] + bB[0][I]*bB[0][I]
              + bB[1][R]*bB[1][R] + bB[1][I]*bB[1][I]
              + bB[2][R]*bB[2][R] + bB[2][I]*bB[2][I];
  B[0][R] = bB[0][R] / sqrt( mag_squared );
  B[0][I] = bB[0][I] / sqrt( mag_squared );
  B[1][R] = bB[1][R] / sqrt( mag_squared );
  B[1][I] = bB[1][I] / sqrt( mag_squared );
  B[2][R] = bB[2][R] / sqrt( mag_squared );
  B[2][I] = bB[2][I] / sqrt( mag_squared );

  C[0][R] =      (A[1][R]*  B[2][R]  - (-A[1][I])*(-B[2][I]))
               - (A[2][R]*  B[1][R]  - (-A[2][I])*(-B[1][I]));
  C[0][I] =      (A[1][R]*(-B[2][I]) + (-A[1][I])*  B[2][R] )
               - (A[2][R]*(-B[1][I]) + (-A[2][I])*  B[1][R] );

  C[1][R] = - (  (A[0][R]*  B[2][R]  - (-A[0][I])*(-B[2][I]))
               - (A[2][R]*  B[0][R]  - (-A[2][I])*(-B[0][I]))  );
  C[1][I] = - (  (A[0][R]*(-B[2][I]) + (-A[0][I])*  B[2][R] )
               - (A[2][R]*(-B[0][I]) + (-A[2][I])*  B[0][R] )  );

  C[2][R] =      (A[0][R]*  B[1][R]  - (-A[0][I])*(-B[1][I]))
               - (A[1][R]*  B[0][R]  - (-A[1][I])*(-B[0][I]));
  C[2][I] =      (A[0][R]*(-B[1][I]) + (-A[0][I])*  B[1][R] )
               - (A[1][R]*(-B[0][I]) + (-A[1][I])*  B[0][R] );
 
  *(eps + 0)  = A[0][R]; *(eps + 9)  = A[0][I];
  *(eps + 1)  = A[1][R]; *(eps + 10) = A[1][I];
  *(eps + 2)  = A[2][R]; *(eps + 11) = A[2][I];
 
  *(eps + 3)  = B[0][R]; *(eps + 12) = B[0][I];
  *(eps + 4)  = B[1][R]; *(eps + 13) = B[1][I];
  *(eps + 5)  = B[2][R]; *(eps + 14) = B[2][I];

  *(eps + 6)  = C[0][R]; *(eps + 15) = C[0][I];
  *(eps + 7)  = C[1][R]; *(eps + 16) = C[1][I];
  *(eps + 8)  = C[2][R]; *(eps + 17) = C[2][I];
 
  /* In general, detailed balance requires that any matrix which can
  ** be generated should have equal probability with it's own inverse.
  ** It should not be obvious that this is the case at this time, so
  ** at this point I flip a coin, heads I continue, tails I invert 
  ** the matrix eps.  */
  if(  absR( LRG.Urand() ) < .5 )
    m_invert( eps );

}

#define	R00	(*(mtx+ 0))
#define	R01	(*(mtx+ 1))
#define	R02	(*(mtx+ 2))
#define	R10	(*(mtx+ 3))
#define	R11	(*(mtx+ 4))
#define	R12	(*(mtx+ 5))
#define	R20	(*(mtx+ 6))
#define	R21	(*(mtx+ 7))
#define	R22	(*(mtx+ 8))
#define	I00	(*(mtx+ 9))
#define	I01	(*(mtx+10))
#define	I02	(*(mtx+11))
#define	I10	(*(mtx+12))
#define	I11	(*(mtx+13))
#define	I12	(*(mtx+14))
#define	I20	(*(mtx+15))
#define	I21	(*(mtx+16))
#define	I22	(*(mtx+17))

Float  m_determinantR( Float* mtx )
{
  Float	result;

  /* Calculate the real part */
  result = 
   + R00*R11*R22 - R00*I11*I22 - I00*R11*I22 - I00*I11*R22
   - R00*R12*R21 + R00*I12*I21 + I00*R12*I21 + I00*I12*R21 
   - R01*R10*R22 + R01*I10*I22 + I01*R10*I22 + I01*I10*R22
   + R01*R12*R20 - R01*I12*I20 - I01*R12*I20 - I01*I12*R20 
   + R02*R10*R21 - R02*I10*I21 - I02*R10*I21 - I02*I10*R21
   - R02*R11*R20 + R02*I11*I20 + I02*R11*I20 + I02*I11*R20; 

  return( result );
}
    
Float  m_determinantI( Float* mtx )
{
  Float	result;

    /* Calculate the imaginary part */
  result = 
   - I00*I11*I22 + I00*R11*R22 + R00*I11*R22 + R00*R11*I22
   + I00*I12*I21 - I00*R12*R21 - R00*I12*R21 - R00*R12*I21 
   + I01*I10*I22 - I01*R10*R22 - R01*I10*R22 - R01*R10*I22
   - I01*I12*I20 + I01*R12*R20 + R01*I12*R20 + R01*R12*I20
   - I02*I10*I21 + I02*R10*R21 + R02*I10*R21 + R02*R10*I21
   + I02*I11*I20 - I02*R11*R20 - R02*I11*R20 - R02*R11*I20;

  return( result );
}

Float  absR( Float x )
{
  if( x>=0 ) return( x );
  else return( -x );
}

void  m_subtract( Float* AminusB, Float* A, Float* B )
{

  *(AminusB + 0)  =  *(A + 0)  - *(B + 0) ;
  *(AminusB + 9)  =  *(A + 9)  - *(B + 9) ;
  *(AminusB + 1)  =  *(A + 1)  - *(B + 1) ;
  *(AminusB + 10) =  *(A + 10) - *(B + 10);
  *(AminusB + 2)  =  *(A + 2)  - *(B + 2) ;
  *(AminusB + 11) =  *(A + 11) - *(B + 11);
  *(AminusB + 3)  =  *(A + 3)  - *(B + 3) ;
  *(AminusB + 12) =  *(A + 12) - *(B + 12);
  *(AminusB + 4)  =  *(A + 4)  - *(B + 4) ;
  *(AminusB + 13) =  *(A + 13) - *(B + 13);
  *(AminusB + 5)  =  *(A + 5)  - *(B + 5) ;
  *(AminusB + 14) =  *(A + 14) - *(B + 14);
  *(AminusB + 6)  =  *(A + 6)  - *(B + 6) ;
  *(AminusB + 15) =  *(A + 15) - *(B + 15);
  *(AminusB + 7)  =  *(A + 7)  - *(B + 7) ;
  *(AminusB + 16) =  *(A + 16) - *(B + 16);
  *(AminusB + 8)  =  *(A + 8)  - *(B + 8) ;
  *(AminusB + 17) =  *(A + 17) - *(B + 17);

}

void m_zero( Float* x ) { for( int i=0; i<18; i++ ) *(x++)=0.; }

void m_print( Float *A )
{

  VRB.Debug("# 0x%10x:\n", A );
  VRB.Debug("# | %+9.6f%+9.6fi ", *(A + 0) , *(A + 9)  );
  VRB.Debug("| %+9.6f%+9.6fi ",   *(A + 1) , *(A + 10) );
  VRB.Debug("| %+9.6f%+9.6fi  |", *(A + 2) , *(A + 11) );
  VRB.Debug("\n");
  VRB.Debug("# | %+9.6f%+9.6fi ", *(A + 3) , *(A + 12) );
  VRB.Debug("| %+9.6f%+9.6fi ",   *(A + 4) , *(A + 13) );
  VRB.Debug("| %+9.6f%+9.6fi  |", *(A + 5) , *(A + 14) );
  VRB.Debug("\n");
  VRB.Debug("# | %+9.6f%+9.6fi ", *(A + 6) , *(A + 15) );
  VRB.Debug("| %+9.6f%+9.6fi ",   *(A + 7) , *(A + 16) );
  VRB.Debug("| %+9.6f%+9.6fi  |", *(A + 8) , *(A + 17) );
  VRB.Debug("\n");
  VRB.Debug("\n");

}

CPS_END_NAMESPACE
