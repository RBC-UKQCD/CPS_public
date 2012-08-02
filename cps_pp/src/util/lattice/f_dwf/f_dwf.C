#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of Fdwf class.

*/
//--------------------------------------------------------------------
//  CVS keywords
//
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
#include <config.h>
#include <util/lattice.h>
#include <util/time_cps.h>
USING_NAMESPACE_CPS

Fdwf::Fdwf() : FdwfBase(){
  cname = "Fdwf";
}

Fdwf::~Fdwf(){
}

#undef PROFILE
ForceArg Fdwf::EvolveMomFforce(Matrix *mom, Vector *chi,
                           Float mass, Float step_size){
  char *fname = "EvolveMomFforce()";
#ifdef PROFILE
  Float time = -dclock();
  ForceFlops=0;
#endif
  ForceArg Fdt = FdwfBase::EvolveMomFforce(mom, chi, mass, step_size);
#ifdef PROFILE
  time += dclock();
  print_flops(cname,fname,ForceFlops,time);
#endif
     
  return Fdt;
}
