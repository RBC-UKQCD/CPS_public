//------------------------------------------------------------------
//
// fermion_vector.h
//
// Header file for the FermionVectorTp  class for  Wilson-like quarks.
//
//
//------------------------------------------------------------------
#ifndef INCLUDED_FERMIONVEC_H
#define INCLUDED_FERMIONVEC_H

#include <math.h>
#include <string.h>
#include <util/gjp.h>
#include <util/lattice.h>
#include <util/random.h>
#include <util/rcomplex.h>
#include <util/vector.h>
#include <util/wilson.h>
#include <alg/alg_base.h>
#include <alg/common_arg.h>
#include <util/site.h>
#include <util/momentum.h>
#include "wilson_matrix.h"

CPS_START_NAMESPACE

class FermionVectorTp {

private:
  Float *fv;

  char *cname ;

public:
  // CREATORS
  FermionVectorTp();
 ~FermionVectorTp();

  // ACCESSORS
  void   print(void) const;
  Float *data (void) const { return fv; }

  // MANIPULATORS

  void SetPointSource( int color, int spin, int x, int y, int z, int t);

  /*! Point source, fixed into Landau gauge */
  void SetGFPointSource( Lattice& lat, int colour, int spin,
                         int x, int y, int z, int t );

  void SetWallSource ( int color, int spin, int time_slice);
  void SetWallSource ( int color, int spin, int time_slice, Float* src);
  void SetBoxSource  ( int color, int spin, int bstart, int bend, 
		    int time_slice);

  /*! Gauge fix sink - Coulomb gauge only */
  void GaugeFixSink       ( Lattice& lat, int dir);

  /*! Gauge fix sink - Landau gauge */
  void LandauGaugeFixSink ( Lattice& lat );
  void SetLandauGaugeMomentaSource ( Lattice& lat , 
                                     int colour   , 
                                     int spin     ,
                                     int p[]     );


  void SetLandauWallSource( Lattice& lat, int spin, int where );
  void GFWallSource       ( Lattice& lat, int spin, int dir, int source_time );

  void ZeroSource();
  void SetVolSource(int color, int spin);
  void SetVolSource(int color, int spin, Float* src) ;

  ///
  //  Tom's momentum stuff
  //  void setMomentum(int* mom, int color, int spin);
  //  void FermionVectorTp::setMomSource(int color, int spin, int source_time,
  //                                   int* mom);

  void SetMomSource(int color, int spin, int source_time, ThreeMom& mom) ;

  void ChiralToDirac();
  void DiracToChiral();

  const Float& operator[](int i);
  /*! Copies a WilsonVector into the site i of the FermionVectorTp */
  void CopyWilsonVec(int i, WilsonVector& WV) ;

   /*! Copies the sink  of a WilsonMatrix  into 
     the site i of the FermionVectorTp                             */
  void CopyWilsonMatSink(int i, int spin, int color, WilsonMatrix& WM) ;

  /* source for the NPR of derivative operators */
  void SetGFDerivativeSource(Lattice& lat,int color, int spin,
			     int d, int x, int y, int z, int t) ;
};

CPS_END_NAMESPACE

#endif

