#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of Fdwf class.

  $Id: f_dwf.C,v 1.4 2004-04-27 03:51:20 cwj Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: cwj $
//  $Date: 2004-04-27 03:51:20 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/f_dwf/f_dwf.C,v 1.4 2004-04-27 03:51:20 cwj Exp $
//  $Id: f_dwf.C,v 1.4 2004-04-27 03:51:20 cwj Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: f_dwf.C,v $
//  $Revision: 1.4 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/f_dwf/f_dwf.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// f_dwf.C
//
// Fdwf is derived from FwilsonTypes and is relevant to
// domain wall fermions
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <stdio.h>
#include <math.h>
#include <util/lattice.h>
#include <util/dirac_op.h>
#include <util/dwf.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/vector.h>
#include <util/random.h>
#include <util/error.h>
#include <comms/scu.h>
#include <comms/glb.h>
USING_NAMESPACE_CPS

Fdwf::Fdwf(){
}

Fdwf::~Fdwf(){
}
