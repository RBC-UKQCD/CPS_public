#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief Alg class methods.
  
  $Id: alg_base.C,v 1.6 2004/08/18 11:57:38 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004/08/18 11:57:38 $
//  $Header: /space/cvs/cps/cps++/src/alg/alg_base/alg_base.C,v 1.6 2004/08/18 11:57:38 zs Exp $
//  $Id: alg_base.C,v 1.6 2004/08/18 11:57:38 zs Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $RCSfile: alg_base.C,v $
//  $Revision: 1.6 $
//  $Source: /space/cvs/cps/cps++/src/alg/alg_base/alg_base.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// alg_base.C
//
// Alg is the base abstract class
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <stdlib.h>	// exit()
#include <alg/alg_base.h>
#include <alg/common_arg.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/error.h>
CPS_START_NAMESPACE

/*!
  \param latt The lattice and actions with which the algorithm is to act.
  \param c_arg The common argument structure for all algorithms, apparently.
*/
Alg::Alg(Lattice & latt, 
	 CommonArg *c_arg) :
	 alg_lattice(latt)
{
  cname = "Alg";
  char *fname = "Alg(L&,CommonArg*)";
  VRB.Func(cname,fname);
  
  // Set the common argument pointer
  //----------------------------------------------------------------
  if(c_arg == 0)
    ERR.Pointer(cname,fname, "common_arg");
  common_arg = c_arg;
}


Alg::~Alg() {
  char *fname = "~Alg()";
  VRB.Func(cname,fname);

  //???
}


Lattice& Alg::AlgLattice()
{
  return alg_lattice;
}




CPS_END_NAMESPACE
