#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-06-04 21:13:59 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_s_spect/aots_s.C,v 1.5 2004-06-04 21:13:59 chulwoo Exp $
//  $Id: aots_s.C,v 1.5 2004-06-04 21:13:59 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: aots_s.C,v $
//  $Revision: 1.5 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_s_spect/aots_s.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
// aots_s.C
CPS_END_NAMESPACE
#include <alg/aots_s.h>
#include <util/gjp.h>
#include <util/verbose.h>
CPS_START_NAMESPACE

//-----------------------------------------------------------------
// static member initialization
//-----------------------------------------------------------------
char Aots::cname[] = "Aots";

//-----------------------------------------------------------------
// CTOR
//-----------------------------------------------------------------
Aots::Aots(int s, int e, int stride)
: start(s), end(e), step(stride), current (s)
{
  char *fname = "Aots(int, int, int)";
  VRB.Func(cname, fname);

  if ( (step <= 0) || (start < 0) || (end < 0) || (start > end) || 
       (start >= GJP.TnodeSites() * GJP.Tnodes()) || 
       (end >= GJP.TnodeSites() * GJP.Tnodes()))
    ERR.General(cname, fname, "Invalid Aots Parameters\n");
}

CPS_END_NAMESPACE
