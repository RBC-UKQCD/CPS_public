#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:52 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/alg/qpropw.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
//  $Id: qpropw.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.2  2001/06/19 18:11:32  anj
//  Serious ANSIfication.  Plus, degenerate double64.h files removed.
//  Next version will contain the new nga/include/double64.h.  Also,
//  Makefile.gnutests has been modified to work properly, propagating the
//  choice of C++ compiler and flags all the way down the directory tree.
//  The mpi_scu code has been added under phys/nga, and partially
//  plumbed in.
//
//  Everything has newer dates, due to the way in which this first alteration was handled.
//
//  Anj.
//
//  Revision 1.2  2001/05/25 06:16:00  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: qpropw.h,v $
//  $Revision: 1.1.1.1 $
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
#include<util/gjp.h>
#include<util/lattice.h>
#include<util/random.h>
#include<util/rcomplex.h>
#include<alg/spin_matrix.h>
#include<util/vector.h>
#include<util/wilson.h>
#include<alg/wilson_matrix.h>
#include<alg/alg_base.h>
#include<alg/common_arg.h>
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
