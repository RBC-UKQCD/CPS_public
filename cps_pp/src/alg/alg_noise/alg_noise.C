#include<config.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------
/*!\file
  \brief Definitions of the AlgNoise class methods.

  $Id: alg_noise.C,v 1.8 2004-09-02 16:56:37 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-09-02 16:56:37 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_noise/alg_noise.C,v 1.8 2004-09-02 16:56:37 zs Exp $
//  $Id: alg_noise.C,v 1.8 2004-09-02 16:56:37 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: alg_noise.C,v $
//  $Revision: 1.8 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_noise/alg_noise.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <stdlib.h>
#include <util/qcdio.h>
#include <math.h>
#include <alg/alg_noise.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/smalloc.h>
#include <util/vector.h>
#include <util/random.h>
#include <util/verbose.h>
#include <util/error.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
/*!
  \param latt The lattice on which to compute the new gauge configuration.
  \param c_arg The common argument structure for all algorithms.
  \param arg  Arguments specific to this algorithm.
 */
//------------------------------------------------------------------
AlgNoise::AlgNoise(Lattice& latt,
		   CommonArg *c_arg,
		   NoiseArg *arg):    Alg(latt, c_arg)
{
  cname = "AlgNoise";
  const char *fname = "AlgNoise";
  VRB.Func(cname,fname);

  // Initialize the argument pointer
  //----------------------------------------------------------------
  if(arg == 0)
    ERR.Pointer(cname,fname, "arg");
  alg_noise_arg = arg;

}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
AlgNoise::~AlgNoise() {
  const char *fname = "~AlgNoise";
  VRB.Func(cname,fname);

}



//------------------------------------------------------------------
//! Performs the computation.
/*!
  \post The results are written to the file specified in the common_arg
  structure.
*/
//------------------------------------------------------------------
void AlgNoise::run()
{
  const char *fname = "run";
  VRB.Func(cname,fname);

  // Set the Lattice pointer and noise_arg
  //----------------------------------------------------------------
  Lattice& lat = AlgLattice();

  // Set the parameters using the argument structure
  //----------------------------------------------------------------
  NoiseType NKind = alg_noise_arg->noise_kind;
  Float size = alg_noise_arg->size;

  VRB.Flow(cname,fname,"Beginning to add noise %4.3f\n", IFloat(size));

  int X[4];
  for( X[0]=0; X[0]<GJP.TnodeSites(); X[0]++ )
  for( X[1]=0; X[1]<GJP.XnodeSites(); X[1]++ )
  for( X[2]=0; X[2]<GJP.YnodeSites(); X[2]++ )
  for( X[3]=0; X[3]<GJP.ZnodeSites(); X[3]++ ) {

      LRG.AssignGenerator(X);

      for( int Mu=0; Mu<4; Mu++ ) {

	  // The coefficients of the Gell-Mann matrices, from
	  // which m in g will be constructed (via M in Lg ).
	  IFloat  gm_coeff[8] = {0,0,0,0,0,0,0,0};
	  Matrix mNoise_g;	// noisy matrix in the group  SU(3)
	  Matrix mNoise_Lg;	// noisy matrix in the Lie algebra of su(3)

	  Matrix* pmThisLink = lat.GaugeField();
	  pmThisLink += Mu + lat.GsiteOffset(X);

	  // Generate 8 random numbers, the coefficients of the
	  // Gell-Mann matrices, from which M in L(g) will be constructed.
	  switch( NKind ) {

	  case FLAT: 
	      LRG.SetInterval(-size,size);
	      for( int i=0; i<8; i++ ) gm_coeff[i] = LRG.Urand();
	      break;

	  case GAUSSIAN: 
	      LRG.SetSigma( size*size );
	      for( int i=0; i<8; i++ ) gm_coeff[i] = LRG.Grand();
	      break;

	  }

	  // Construct M in L(g)
	  mNoise_Lg.AntiHermMatrix( gm_coeff );

	  mNoise_g = Exponentiate_Matrix( mNoise_Lg, 8 );

	  Matrix mTmpMat;
	  mTmpMat.DotMEqual( mNoise_g, *pmThisLink );
	  *pmThisLink = mTmpMat;

      }
  }

}

Matrix AlgNoise::Exponentiate_Matrix( Matrix iA, int order )
{

  Matrix mg;
  mg.ZeroMatrix();

  Matrix iAn;	// iA to the power n.
  iAn.UnitMatrix();
  Float nFactorial = 1;

  for( int n=0; n<=order; n++ ) {
    Matrix u = iAn;
    vecTimesEquFloat( (IFloat*)&u, 1./nFactorial, 18 );
    vecAddEquVec( (IFloat*)&mg, (const IFloat*)&u, 18 );

    mDotMEqual( (IFloat*)&u, (const IFloat*)&iAn, (const IFloat*)&iA);
    iAn = u;
    nFactorial *= (n+1);
  }
  
  mg.Unitarize();
  return mg;
}

CPS_END_NAMESPACE
