#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2003-07-24 16:53:53 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/gw_hb.h,v 1.2 2003-07-24 16:53:53 zs Exp $
//  $Id: gw_hb.h,v 1.2 2003-07-24 16:53:53 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.2  2001/06/19 18:13:17  anj
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
//  Revision 1.2  2001/05/25 06:16:09  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: gw_hb.h,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/gw_hb.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// gw_hb.h
//
// Header file for the GwHb base class. There are no derived
// classes. GwHb is the front-end for the Wilson gauge action
// heat bath code.
//
//------------------------------------------------------------------


#ifndef INCLUDED_GW_HB_H
#define INCLUDED_GW_HB_H

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/vector.h>
#include <alg/ghb_arg.h>
CPS_START_NAMESPACE

//------------------------------------------------------------------
//
// GwHb base class.
//
//------------------------------------------------------------------
class GwHb
{
 private:
  char *cname;                     // Class name.

 protected:
  Lattice& lat;                    // Lattice object.
  Matrix *gauge_field;             // pointer to the gauge field

  
 public:
  GwHb(Lattice& latt);             // Lattice object.

  virtual ~GwHb();

  void HeatBath(GhbArg *ghb_arg);
      // The heat bath algorithm

};

#endif

CPS_END_NAMESPACE
