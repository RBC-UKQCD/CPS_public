#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of FdwfBase class.

  $Id: f_dwf_base_force.C,v 1.3 2008-02-08 18:35:08 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/f_dwf_base/qmp/f_dwf_base_force.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// f_dwf_base.C
//
// FdwfBase is derived from FwilsonTypes and is relevant to
// domain wall fermions, moved from old Fdwf class
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/qcdio.h>
#include <math.h>
#include <util/lattice.h>
#include <util/dirac_op.h>
#include <util/dwf.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/vector.h>
#include <util/random.h>
#include <util/error.h>
#include <util/time_cps.h>
#include <comms/scu.h> // GRF
#include <comms/glb.h>
CPS_START_NAMESPACE



#undef PROFILE
// CJ: change start
//------------------------------------------------------------------
// EvolveMomFforce(Matrix *mom, Vector *chi, Float mass, 
//                 Float dt):
// It evolves the canonical momentum mom by dt
// using the fermion force.
//------------------------------------------------------------------
ForceArg FdwfBase::EvolveMomFforce(Matrix *mom, Vector *chi, 
			   Float mass, Float dt){
  ERR.NotImplemented("FdwfBase", "EvolveMomFforce");
}
CPS_END_NAMESPACE
