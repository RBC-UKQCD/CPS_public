#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-08-18 11:57:37 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/gw_hb.h,v 1.4 2004-08-18 11:57:37 zs Exp $
//  $Id: gw_hb.h,v 1.4 2004-08-18 11:57:37 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: gw_hb.h,v $
//  $Revision: 1.4 $
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
