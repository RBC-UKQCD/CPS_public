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
#include <vector>
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
#include <alg/wilson_matrix.h>

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

  Float Norm();

  // MANIPULATORS

  void SetPointSource( int color, int spin, int x, int y, int z, int t);

  /*! Point source, fixed into Landau gauge */
  void SetGFPointSource( Lattice& lat, int colour, int spin,
                         int x, int y, int z, int t );

  void SetWallSource ( int color, int spin, int time_slice);
  void SetWallSource ( int color, int spin, int time_slice, Float* src);

  // if src_offset is non NULL use the array as offsets for the starts and ends of the box soruce
  void SetBoxSource  ( int color, int spin, int bstart, int bend, 
		       int time_slice,
		       int* src_offset=0);

  // Sets a 4D box source. If you want to set a 3D xyz box, set
  // size[3] == 1 and glb_x[3] to the global time slice you want.
  //
  // With this function in disposal, why do we need separate functions
  // for point/wall/box sources?
  void Set4DBoxSource(int color,
                      int spin,
                      const int start[4], // global starting location in x, y, z and t directions
                      const int size[4], // global size in x, y, z and t directions
                      const Float mom[4]); // momentum

  void SetZ3BWall(int color, int spin, int t, const int size[3],
                  const std::vector<Rcomplex> &rand_num);

  /*! Gauge fix sink - Coulomb gauge only */
  void GaugeFixSink       ( Lattice& lat, int dir, int unfix=0);

  /*! Gauge fix sink - Landau gauge */
  void LandauGaugeFixSink ( Lattice& lat );
  void SetLandauGaugeMomentaSource ( Lattice& lat , 
                                     int colour   , 
                                     int spin     ,
                                     int p[]     );


  void SetLandauWallSource( Lattice& lat, int spin, int where );
  void GFWallSource       ( Lattice& lat, int spin, int dir, int source_time );

  void ZeroSource();
//  void SetVolSourceEqualZero();
  void SetVolSource(int color, int spin);
  void SetVolSource(int color, int spin, Float* src) ;

  //Peter and HueyWen for general function source
  void SetGFLfuncSource(Lattice& lat, int color, int spin,
                Float (*func)(int gx,int gy,int gz,int gt));

  // === for exponential smeared source ================================
  void SetExpSource( int color, int spin, int x, int y, int z, int t,
                     Float A, Float B, Float C) ;
  void dumpqrkvec( int t);
  // ===================================================================

  // Wuppertal-like smearing on fv
  void GaussianSmearVector(Lattice& lat,
                           int spin,
                           int iter,
                           Float omega,
                           int source_time);
  // Wuppertal-like smearing on fv
  void GaussianSmearVector(Lattice& lat,
                           int spin,
                           int iter,
                           Float omega);

  ///
  //  Tom's momentum stuff
  //  void setMomentum(int* mom, int color, int spin);
  //  void FermionVectorTp::setMomSource(int color, int spin, int source_time,
  //                                   int* mom);

  void SetMomSource(int color, int spin, int source_time, ThreeMom& mom) ;
  void SetMomCosSource(int color, int spin, int source_time, ThreeMom& mom) ;
  void SetMomCosTwistSource(int color, int spin, int source_time, ThreeMomTwist& mom); // Use this if you are using twisted boundary conditions

  void ChiralToDirac();
  void DiracToChiral();

  const Float& operator[](int i);
  FermionVectorTp& operator*=(Float f) ;

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

