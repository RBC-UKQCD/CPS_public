#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-08-18 11:57:36 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/alg/qpropw.h,v 1.4 2004-08-18 11:57:36 zs Exp $
//  $Id: qpropw.h,v 1.4 2004-08-18 11:57:36 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.4 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/alg/qpropw.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// propw.h
//
// Header file for the QPropW and WilsonMatrix classes for 
// Wilson-like quarks.
//
// For now this is specific to three colors. The constructor
// will exit if the number of colors is not equal to three.
//
// This file contains the declarations of the QPropW class.
//
//------------------------------------------------------------------


#ifndef INCLUDED_PROPW_H
#define INCLUDED_PROPW_H

CPS_END_NAMESPACE
#include <math.h>
#include <string.h>
#include <util/gjp.h>
#include <util/lattice.h>
#include <util/random.h>
#include <util/rcomplex.h>
#include <alg/spin_matrix.h>
#include <util/vector.h>
#include <util/wilson.h>
#include <alg/wilson_matrix.h>
#include <alg/alg_base.h>
#include <alg/common_arg.h>
CPS_START_NAMESPACE

const int WALL = 11;
const int RWALL = 12;

class FermionVectorTp {

private:
  Float *fv_;

  FermionVectorTp(FermionVectorTp&);
  FermionVectorTp& operator=(FermionVectorTp&);

public:
  // CREATORS
  FermionVectorTp();
 ~FermionVectorTp();

  // ACCESSORS
  void print(void) const;
  Float * data(void) const { return fv_; }

  // MANIPULATORS
  void setPointSource(int color, int spin, int x, int y, int z, int t);
  void setWallSource(int color, int spin, int time_slice);
  void setWallSource(int color, int spin, int time_slice, Float* src);
  void FermionVectorTp::GaugeFixSink(Lattice& lat, int dir);
  void FermionVectorTp::LandauGaugeFixSink(Lattice& lat, int dir);
  void FermionVectorTp::setLandauWallSource(Lattice& lat, int color, int spin,
                int source_time);
  void FermionVectorTp::setGFWallSource(Lattice& lat, int color, int spin,
                int dir, int source_time);
  void setVolSource();
  void setVolSource(int color, int spin);

  const Float& operator[](int i);
};

// The QPropW class.
//------------------------------------------------------------------
class QPropW
{

  protected:

    WilsonMatrix* prop;
    
  public:

    // CONSTRUCTORS
    QPropW();
    
    void CG(Lattice& lat, CgArg* cg, FermionVectorTp&, 
		FermionVectorTp&, int&, Float&);
    void CGDwf(Lattice& lat, CgArg* cg, FermionVectorTp&, Vector*);

    // operator functions
    QPropW& operator=(const QPropW& rhs);
    WilsonMatrix& operator[](int i);

};

class QPropWRandWallSrc: public QPropW
{

    Float* rsrc;

  public:

    // CONSTRUCTORS
    QPropWRandWallSrc();
    QPropWRandWallSrc(Lattice& lat, CgArg* cg, int source_time, 
			int seed, CommonArg*);
    QPropWRandWallSrc::QPropWRandWallSrc(Lattice& lat, CgArg *arg,
                QPropWRandWallSrc* prop1, QPropWRandWallSrc* prop2 );
    
    //MANIPULATORS
    const Rcomplex& rand_src(int i);

    //operator functions

    // DESTRUCTOR
    void Delete();
    ~QPropWRandWallSrc();


};

class QPropWWallSrc: public QPropW
{

  public:

    // CONSTRUCTORS
    QPropWWallSrc();
    QPropWWallSrc(Lattice& lat, CgArg* cg, int source_time,
			CommonArg* c_arg);
    QPropWWallSrc::QPropWWallSrc(Lattice& lat, CgArg *arg,
                	QPropWWallSrc& prop1, QPropWWallSrc& prop2 );
    QPropWWallSrc::QPropWWallSrc(Lattice& lat, CgArg *arg,
			QPropWWallSrc* prop1, QPropWWallSrc* prop2 );

    // DESTRUCTOR
    ~QPropWWallSrc();

};

class QPropWVolSrc: public QPropW
{

  public:

    // CONSTRUCTORS
    QPropWVolSrc(Lattice& lat, CgArg* cg);

    // DESTRUCTOR
    ~QPropWVolSrc();

};


#endif


CPS_END_NAMESPACE
