#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-08-18 11:57:39 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_inst/alg_inst.C,v 1.8 2004-08-18 11:57:39 zs Exp $
//  $Id: alg_inst.C,v 1.8 2004-08-18 11:57:39 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: alg_inst.C,v $
//  $Revision: 1.8 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_inst/alg_inst.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// alg_inst.h
//
// Header file for the AlgInst class.
//
// AlgInst is derived from Alg and is relevant to the 
// generation of an instanton configuration.
//
//----- Begin Adrian
// ???
//
//      This operator creates a single instanton at
//      the center of the lattice.  It is aware of the
//      overall size of the lattice, and uses this, to
//      place the instanton at the actual center of the
//      whole lattice, regardless of the number of processors.
//        The instanton requires two arguments.
//              0       'Instanton Size'
//              1       average value of delta-nu.  If this
//                       argument is +1 or -1 then the move
//                       is forced, if it is 0, then the
//                       move is in a random direction and
//                       there is an accept reject step.
//              2       For squashed instantons, this is
//                        the value of "r_max"
//      The return value is the change in the action of the 
//      configuration as a result of the operation.
//
//        The actual 'instanton' configuration is the solution
//      outlined in Coleman(p297):
//              g1 = (1/r) * (x_4 + i x_dot_sigma )
//              f(r) = r^2/(r^2/rho^2)
//              A_mu = f(r) g1 d_mu [g1]^-1
//
//----- End Adrian
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <stdlib.h>
#include <util/qcdio.h>
#include <math.h>
#include <alg/alg_inst.h>
#include <alg/common_arg.h>
#include <alg/inst_arg.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/smalloc.h>
#include <util/vector.h>
#include <util/verbose.h>
#include <util/error.h>
CPS_START_NAMESPACE


static void AddInstanton( Matrix*, Float[4], InstMethod );
static int Eta( int i, int mu, int nu );
static int EtaBar( int i, int mu, int nu );
static Float Rho( Float, Float, Float ); // function gives the radial
					 // dependance of rho, to
					 // "squash" the instanton.


//------------------------------------------------------------------
// Constructor 
//------------------------------------------------------------------
AlgInst::AlgInst(Lattice& latt, 
 	         CommonArg *c_arg,
	         InstArg *arg) : 
	         Alg(latt, c_arg) 
{
  cname = "AlgInst";
  char *fname = "AlgInst(L&,CommonArg*,InstArg*)";
  VRB.Func(cname,fname);

  // Initialize the argument pointer
  //----------------------------------------------------------------
  if(arg == 0)
    ERR.Pointer(cname,fname, "arg");
  alg_inst_arg = arg;

}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
AlgInst::~AlgInst() {
  char *fname = "~AlgInst()";
  VRB.Func(cname,fname);

}


//------------------------------------------------------------------
// run() runs the code
//------------------------------------------------------------------
void AlgInst::run()
{
  char *fname = "run()";
  VRB.Func(cname,fname);

  //----------------------------------------------------------------
  // Check if anisotropy is present and exit 
  //----------------------------------------------------------------
  if(GJP.XiBare() != 1 ||
     GJP.XiV()    != 1 ||
     GJP.XiVXi()  != 1   ){
    ERR.General(cname,fname,
    "XiBare=%g, XiV=%g, XiVXi=%g : Not implemented for anisotropy\n",
                GJP.XiBare(), GJP.XiV(), GJP.XiVXi());
  }

  // Set the Lattice pointer and inst_arg
  //----------------------------------------------------------------
  Lattice& lat = AlgLattice();

  // Set the parameters using the argument structure
  //----------------------------------------------------------------
  InstType eIType = alg_inst_arg->inst_kind;
  Float rRhoCutoff = alg_inst_arg->rho_cutoff;
  Float rCharge = alg_inst_arg->charge;
  Float rRho = alg_inst_arg->rho;
  int n1 = alg_inst_arg->n1;
  int n2 = alg_inst_arg->n2;
  int nCharge = int(rCharge);

  //----------------------------------------------------------------
  // Define a temporary matrix to be used as buffer
  //----------------------------------------------------------------
  Matrix mTmpMat;

  //----------------------------------------------------------------
  // Adrian's code
  //----------------------------------------------------------------


  // Set the desired instanton size.
  Float	rRho2 = rRho*rRho;

  // Whipe the lattice clean before trying to lay down an instanton.
  //???
  //  void  order( int job_id );
  //  order( job_id );
  lat.SetGfieldOrd();

  if( eIType == REGULAR_SQUASHED 
     || eIType == REGULAR_SQUASHED_TRANSFORMED )
  VRB.Flow(cname,fname,"Beginning to calculate instanton %4.3f,%d,%f)\n", 
	   IFloat(rRho), nCharge, IFloat(rRhoCutoff) );
  else
  VRB.Flow(cname,fname,"Beginning to calculate instanton %4.3f,%d)\n", 
	   IFloat(rRho), nCharge);

VRB.Debug("---\n\n  eIType = %d\n\n---\n", (int)eIType );
VRB.Debug("---\n\n  n1,n2 = %d,%d\n\n---\n", n1,n2 );

  Float arInstCenter[4];
  arInstCenter[0] = (GJP.TnodeSites()*GJP.Tnodes() - 1)/2.;
  arInstCenter[1] = (GJP.XnodeSites()*GJP.Xnodes() - 1)/2.;
  arInstCenter[2] = (GJP.YnodeSites()*GJP.Ynodes() - 1)/2.;
  arInstCenter[3] = (GJP.ZnodeSites()*GJP.Znodes() - 1)/2.;

//  unsigned auX[4];
  int auX[4];
  for( auX[0]=0; auX[0]<GJP.TnodeSites(); auX[0]++ )
  for( auX[1]=0; auX[1]<GJP.XnodeSites(); auX[1]++ )
  for( auX[2]=0; auX[2]<GJP.YnodeSites(); auX[2]++ )
  for( auX[3]=0; auX[3]<GJP.ZnodeSites(); auX[3]++ ) {

    // Define y to be the point x on this node, located with
    // reference to the instanton's center.
    Float arY[4] = {
      GJP.TnodeCoor()*GJP.TnodeSites() + auX[0] - arInstCenter[0],
      GJP.XnodeCoor()*GJP.XnodeSites() + auX[1] - arInstCenter[1],
      GJP.YnodeCoor()*GJP.YnodeSites() + auX[2] - arInstCenter[2],
      GJP.ZnodeCoor()*GJP.ZnodeSites() + auX[3] - arInstCenter[3] };

    Float r2=0.0;
    unsigned uMu;
    for( uMu=0; uMu<4; uMu++ )
      r2 += arY[uMu]*arY[uMu];

    Float  arTheta[4];
    int auX_can[4];
    int uMu_can[4];

    //???
    auX_can[0] = auX[1];
    auX_can[1] = auX[2];
    auX_can[2] = auX[3];
    auX_can[3] = auX[0];
    uMu_can[0] = 3; 
    uMu_can[1] = 0;
    uMu_can[2] = 1;
    uMu_can[3] = 2;


    for( uMu=0; uMu<4; uMu++ ) {

      //???
      Matrix* pmThisLink = lat.GaugeField();
      pmThisLink += uMu_can[uMu] + lat.GsiteOffset(auX_can);


      switch( eIType ) {
        case REGULAR:
        case REGULAR_TRANSFORMED: {

          Float t1 = sqrt( rRho2 + r2 - arY[uMu]*arY[uMu] );
          Float t2 = rRho2 + r2 + arY[uMu];

          Float rPhi = - atan(t1/t2)/t1;

          // The link U_mu is exp( i Theta dot Sigma ), Theta
          // has only three components, 1,2,3 but Theta[0] is
          // used for the magnitude of theta.  The components
          // are then normalized.

          for( int j=0; j<4; j++ ) arTheta[j]=0.0;

          // Gor each generator of SU(2)...
          for( int ii = 1; ii<4; ii++ ) {

            // For each direction...
            for( int uNu = 0; uNu<4; uNu++ ) 
              if(nCharge==+1)
                arTheta[ii] += EtaBar(ii,uMu,uNu) * arY[uNu] * rPhi;
              else
                arTheta[ii] += Eta(ii,uMu,uNu) * arY[uNu] * rPhi;

            arTheta[0] += arTheta[ii] * arTheta[ii];
          }
          arTheta[0] = sqrt( (double)arTheta[0] );
          arTheta[1] /= arTheta[0];
          arTheta[2] /= arTheta[0];
          arTheta[3] /= arTheta[0];

        }
	AddInstanton( pmThisLink, arTheta, alg_inst_arg->inst_method );
        break;

        case REGULAR_SQUASHED:
        case REGULAR_SQUASHED_TRANSFORMED: {

          Float arZ[4] = { arY[0], arY[1], arY[2], arY[3] };
	  //??????
          arZ[uMu] = arZ[uMu] + 1.0;

          Float z2 = arZ[0]*arZ[0] + arZ[1]*arZ[1]
                  + arZ[2]*arZ[2] + arZ[3]*arZ[3];

          Float z=sqrt((double)z2);
          Float r=sqrt((double)r2);


          Float rRho_x = Rho(r, rRho, rRhoCutoff );
          Float rRho_z = Rho(z, rRho, rRhoCutoff );

          Float rPhi;

          rPhi  = - 1./(  r2 + rRho_x*rRho_x );
          rPhi += - 1./(  z2 + rRho_z*rRho_z );
          rPhi /= 2.;

          for( int j=0; j<4; j++ ) arTheta[j]=0.;

          for( int ii = 1; ii<4; ii++ ) {

            for( int uNu = 0; uNu<4; uNu++ ) {
              if(nCharge==+1)
                arTheta[ii] += EtaBar(ii,uMu,uNu) * arY[uNu] * rPhi;
              else
                arTheta[ii] += Eta(ii,uMu,uNu) * arY[uNu] * rPhi;
            }

            arTheta[0] += arTheta[ii] * arTheta[ii];
          }
          arTheta[0] = sqrt( (double)arTheta[0] );
          arTheta[1] /= arTheta[0];
          arTheta[2] /= arTheta[0];
          arTheta[3] /= arTheta[0];
        }
	AddInstanton( pmThisLink, arTheta, alg_inst_arg->inst_method );
        break;

        case SINGULAR: {

          Float b1 = sqrt( (double)(rRho2 + r2 - arY[uMu]*arY[uMu]) );
          Float b2 = rRho2 + r2 + arY[uMu];

          Float a1 = sqrt( (double)(r2 - arY[uMu]*arY[uMu]) );
          Float a2 = r2 + arY[uMu];
         
          Float rPhi = atan(b1/b2)/b1 - atan(a1/a2)/a1;

          // The link U_mu is exp( i Theta dot Sigma ), Theta
          // has only three components, 1,2,3 but Theta[0] is
          // used for the magnitude of theta.  The components
          // are then normalized.

          for( int j=0; j<4; j++ ) arTheta[j]=0.;

          for( int ii = 1; ii<=3; ii++ ) {

            for( int uNu = 0; uNu<=3; uNu++ ) {
              if(nCharge==+1)
                arTheta[ii] += Eta(ii,uMu,uNu)*arY[uNu]*rPhi;
              else
                arTheta[ii] += EtaBar(ii,uMu,uNu)*arY[uNu]*rPhi;
            }

            arTheta[0] += arTheta[ii]*arTheta[ii];
          }
          arTheta[0] = sqrt( (double)arTheta[0] );
          arTheta[1] /= arTheta[0];
          arTheta[2] /= arTheta[0];
          arTheta[3] /= arTheta[0];

        }
	AddInstanton( pmThisLink, arTheta, alg_inst_arg->inst_method );
        break;

        case CONSTANT_FIELD: {
   
          Float My_PI = 3.14159265358979323846;

          int anY[4] = {
            GJP.TnodeCoor()*GJP.TnodeSites() + auX[0],
            GJP.XnodeCoor()*GJP.XnodeSites() + auX[1],
            GJP.YnodeCoor()*GJP.YnodeSites() + auX[2],
            GJP.ZnodeCoor()*GJP.ZnodeSites() + auX[3] };

          int nL0 = GJP.TnodeSites()*GJP.Tnodes();
          int nL1 = GJP.XnodeSites()*GJP.Xnodes();
          int nL2 = GJP.YnodeSites()*GJP.Ynodes();
          int nL3 = GJP.ZnodeSites()*GJP.Znodes();

          Float rOmega1 = (2.*My_PI*n1)/(Float)(nL1*nL2);
          Float rOmega2 = (2.*My_PI*n2)/(Float)(nL3*nL0);
   
          (*pmThisLink).UnitMatrix();
  
          Float* arLink = (Float*) pmThisLink;
  
          switch( uMu ) {
     
            case 1:
            arLink[0] = cos( + rOmega1 * anY[2] );
            arLink[1] = sin( + rOmega1 * anY[2] );
            arLink[8] = + arLink[0];
            arLink[9] = - arLink[1];
            break;
     
            case 2:
            if( anY[2] == nL2-1 ) {
              arLink[0] = cos( - rOmega1 * nL2 * anY[1] );
              arLink[1] = sin( - rOmega1 * nL2 * anY[1] );
              arLink[8] = + arLink[0];
              arLink[9] = - arLink[1];
            }
            break;
     
            case 3:
            arLink[0] = cos( - rOmega2 * anY[0] );
            arLink[1] = sin( - rOmega2 * anY[0] );
            arLink[8] = + arLink[0];
            arLink[9] = - arLink[1];
            break;
     
            case 0: // case 4: in the case of S&V's Euclidian notation
            if( anY[0] == nL0-1 ) {
              arLink[0] = cos( + rOmega2 * nL0 * anY[3] );
              arLink[1] = sin( + rOmega2 * nL0 * anY[3] );
              arLink[8] = + arLink[0];
              arLink[9] = - arLink[1];
            }
            break;
     
            default:
            exit(-1);
          }

//char ccc;
//switch(uMu){case 0:ccc='t';break;case 1:ccc='x';break;case 2:ccc='y';break;case 3:ccc='z';break;}
//VRB.Debug("%d %d %d %d %c ", anY[0],anY[1],anY[2],anY[3], ccc );
//VRB.Debug("%+.20e %+.20e ", arLink[0], arLink[1] );
//VRB.Debug("%+.20e %+.20e ", arLink[2], arLink[3] );
//VRB.Debug("%+.20e %+.20e ", arLink[4], arLink[5] );
//VRB.Debug("%+.20e %+.20e ", arLink[6], arLink[7] );
//VRB.Debug("%+.20e %+.20e ", arLink[8], arLink[9] );
//VRB.Debug("%+.20e %+.20e ", arLink[10], arLink[11] );
//VRB.Debug("%+.20e %+.20e ", arLink[12], arLink[13] );
//VRB.Debug("%+.20e %+.20e ", arLink[14], arLink[15] );
//VRB.Debug("%+.20e %+.20e\n", arLink[16], arLink[17] );
        }

        break;
      }



      // Now for the types which require a gauge transformation,
      // do that. 
      if( eIType == REGULAR_TRANSFORMED 
       || eIType == REGULAR_SQUASHED_TRANSFORMED ) {

        // Here we have an instanton, and we would
        // like to gauge transform it here on the lattice with
        // a (now manifestly nonsingular) gauge transformation, to
        // supress the boundary effects.

        // It is convenient to define a point z which is at the
        // other end of the link U_mu(y), z = y + mu_hat. Then
        // the link U_mu(y) becomes, under the transformation
        // g which takes a fermion psi(y) to g(y).psi(y), into 
	// U_mu(y) -> g(y) U_mu(y) g(z)^-1.

        Float arZ[4] = { arY[0], arY[1], arY[2], arY[3] };
        arZ[uMu] += 1.0;

        Matrix	mG;		// matrix G(y)
        Matrix	mG_Inv;		// matrix G(z)^-1
	Float *prG = (Float *) &mG;
	Float *prG_Inv = (Float *) &mG_Inv;

	//??????
	mG.UnitMatrix();
	mG_Inv.UnitMatrix();

        Float z2 = arZ[0]*arZ[0] + arZ[1]*arZ[1]
                + arZ[2]*arZ[2] + arZ[3]*arZ[3];

        Float z=sqrt((double)z2);
        Float r=sqrt((double)r2);

        prG[0] =  arZ[0]/z;			// G(z)
        prG[1] = -nCharge * arZ[3]/z;
        prG[2] = -nCharge * arZ[2]/z;
        prG[3] = -nCharge * arZ[1]/z;
        prG[6] =  nCharge * arZ[2]/z;
        prG[7] = -nCharge * arZ[1]/z;
        prG[8] =  arZ[0]/z;
        prG[9] =  nCharge * arZ[3]/z;

        prG_Inv[0] =  arY[0]/r;		// G(y)^-1
        prG_Inv[1] =  nCharge * arY[3]/r;
        prG_Inv[2] =  nCharge * arY[2]/r;
        prG_Inv[3] =  nCharge * arY[1]/r;
        prG_Inv[6] = -nCharge * arY[2]/r;
        prG_Inv[7] =  nCharge * arY[1]/r;
        prG_Inv[8] =  arY[0]/r;
        //prG_Inv[9] = -nCharge * arY[3]/r;

	mTmpMat.DotMEqual(*pmThisLink, mG);
	pmThisLink->DotMEqual(mG_Inv, mTmpMat);


      }


    }
  }

}

static int Eta( int i, int mu, int nu )
// This is the SO(4)-SU(2) mixing tensor
// defined by Rajaraman
// which naturally arises in the definition of the 
// regular (non-singular) gauge instanton.
{
  int eta=0;

  if( nu == 0 )
    if( i == mu ) eta = +1;
  if( mu == 0 )
    if( i == nu ) eta = -1;
        
  if( i == 1 && mu==2 && nu==3 ) eta = +1;
  if( i == 1 && mu==3 && nu==2 ) eta = -1;
  if( i == 2 && mu==3 && nu==1 ) eta = +1;
  if( i == 2 && mu==1 && nu==3 ) eta = -1;
  if( i == 3 && mu==1 && nu==2 ) eta = +1;
  if( i == 3 && mu==2 && nu==1 ) eta = -1;

  return eta;
}

static int EtaBar( int i, int mu, int nu )
// This is the space-time-group mixing tensor
// origionally defined by t'Hooft (read from Rajaraman)
// which naturally arises in the definition of the 
// singular gauge instanton.
{
  int eta=0;

  if( nu == 0 && mu != 0 )
    if( i == mu ) eta = -1;
  if( mu == 0 && nu != 0 )
    if( i == nu ) eta = +1;
        
  if( i == 1 && mu==2 && nu==3 ) eta = +1;
  if( i == 1 && mu==3 && nu==2 ) eta = -1;
  if( i == 2 && mu==3 && nu==1 ) eta = +1;
  if( i == 2 && mu==1 && nu==3 ) eta = -1;
  if( i == 3 && mu==1 && nu==2 ) eta = +1;
  if( i == 3 && mu==2 && nu==1 ) eta = -1;

  return eta;
}


/* This is the function which squashes "squashed-instantons" */
static Float Rho( Float r, Float rho, Float cut ) 
{
  Float  rRho_x;

  if( r < cut ) {
    rRho_x = rho * ( 1. - ( r / cut )  );
  }
  else {
    rRho_x = 0.;
  }

  return rRho_x;

}


static void AddInstanton( Matrix* pmThisLink, Float arTheta[4],
  InstMethod IM ) 
{
      Matrix mPhase;
      mPhase.UnitMatrix();
      Float *prPhase = (Float *) &mPhase;

      // Ok, now I have an exponential of an SU(2) matrix,
      // exp( i theta_0 * SUM[j] theta_j sigma_j )
      //  = cos( theta_0 ) + i SUM[j] theta_j sigma_j sin( theta_0 )
      //
      // Transcribe this into the corner of an SU(3) matrix:

      prPhase[0] =                cos( arTheta[0] );
      prPhase[1] =   arTheta[3] * sin( arTheta[0] );

      prPhase[2] =   arTheta[2] * sin( arTheta[0] );
      prPhase[3] =   arTheta[1] * sin( arTheta[0] );

      prPhase[6] = - arTheta[2] * sin( arTheta[0] );
      prPhase[7] =   arTheta[1] * sin( arTheta[0] );

      prPhase[8] =                cos( arTheta[0] );
      prPhase[9] = - arTheta[3] * sin( arTheta[0] );

      switch( IM ) {
        case DESTROY:
          *pmThisLink = mPhase;
          break;
        case ADD: {
          Matrix mTmpMat;
          mTmpMat.DotMEqual(*pmThisLink, mPhase);
          *pmThisLink = mTmpMat;
        } break;
      }
}

CPS_END_NAMESPACE
