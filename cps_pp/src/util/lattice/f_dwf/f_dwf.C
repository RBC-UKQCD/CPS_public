#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of Fdwf class.

  $Id: f_dwf.C,v 1.5 2004-06-04 21:14:12 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-06-04 21:14:12 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/f_dwf/f_dwf.C,v 1.5 2004-06-04 21:14:12 chulwoo Exp $
//  $Id: f_dwf.C,v 1.5 2004-06-04 21:14:12 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: f_dwf.C,v $
//  $Revision: 1.5 $
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
