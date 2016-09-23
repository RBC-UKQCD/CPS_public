#include <config.h>
#ifdef USE_SSE
#include "../sse/wilson_end.C"
#else
CPS_START_NAMESPACE
/*! \file
  \brief  Routine used internally in the DiracOpWilson class.

  $Id: wilson_end.C,v 1.5 2011/03/04 11:25:28 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2011/03/04 11:25:28 $
//  $Header: /space/cvs/cps/cps++/src/util/dirac_op/d_op_wilson/noarch/wilson_end.C,v 1.5 2011/03/04 11:25:28 chulwoo Exp $
//  $Id: wilson_end.C,v 1.5 2011/03/04 11:25:28 chulwoo Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $Revision: 1.5 $
//  $Source: /space/cvs/cps/cps++/src/util/dirac_op/d_op_wilson/noarch/wilson_end.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/wilson.h>
#include <util/smalloc.h>
#include <util/verbose.h>
CPS_START_NAMESPACE

  
void wilson_end( Wilson *wilson_p)
{
  char *cname = " ";
  char *fname = "wilson_end(Wilson*)";
  VRB.Func(cname,fname);

  VRB.Sfree(cname,fname, "af[1]", wilson_p->af[1]);
  sfree(wilson_p->af[1]);

  VRB.Sfree(cname,fname, "af[0]", wilson_p->af[0]);
  sfree(wilson_p->af[0]);

  VRB.Sfree(cname,fname, "ptr", wilson_p->ptr);
  sfree(wilson_p->ptr);

}

CPS_END_NAMESPACE
#endif
