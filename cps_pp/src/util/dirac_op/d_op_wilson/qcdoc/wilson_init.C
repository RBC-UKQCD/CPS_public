/*!\file
  Wilson Dirac operator code for QCDOC

  $Id: wilson_init.C,v 1.7 2004-08-18 11:57:56 zs Exp $
*/

#include <util/wilson.h>
#include <util/gjp.h>
#include <util/smalloc.h>
#include <util/verbose.h>
#include <util/error.h>
#include <util/wfm.h>

CPS_START_NAMESPACE

/*----------------------------------------------------------------------*/
/* Pointer to Wilson type structure is ignored here                     */
/*----------------------------------------------------------------------*/
static WilsonArg wil;

void wilson_end(Wilson *wilson_p)
{
  VRB.Func("","wilson_end()");
  wilson_compat_end(wilson_p);
  wfm_end(&wil); 
}

void wilson_init(Wilson *wilson_p)  
{
  char *cname = "wfm";
  char *fname = "wilson_init(Wilson *)";


/*----------------------------------------------------------------------*/
/* Set sublattice direction sizes                                       */
/*----------------------------------------------------------------------*/

  wil.local_latt[0] = GJP.XnodeSites();
  wil.local_latt[1] = GJP.YnodeSites();
  wil.local_latt[2] = GJP.ZnodeSites();
  wil.local_latt[3] = GJP.TnodeSites();

/*----------------------------------------------------------------------*/
/* Set whether comms are necessary                                      */
/*----------------------------------------------------------------------*/

  if ( GJP.Xnodes() > 1 )  wil.local_comm[0] = 0; 
  else  wil.local_comm[0] = 1; 

  if ( GJP.Ynodes() > 1 )  wil.local_comm[1] = 0; 
  else  wil.local_comm[1] = 1; 

  if ( GJP.Znodes() > 1 )  wil.local_comm[2] = 0; 
  else  wil.local_comm[2] = 1; 

  if ( GJP.Tnodes() > 1 )  wil.local_comm[3] = 0; 
  else  wil.local_comm[3] = 1; 

  //Temporary hack to force non-local comms
  wil.local_comm[0] = 0;
  wil.local_comm[1] = 0;
  wil.local_comm[2] = 0;
  wil.local_comm[3] = 0;

  wfm_init (&wil);
  wilson_compat_init(wilson_p);

}

/*Wrapper used by the clover code*/
void wilson_dslash(Float *chi, 
		   Float *u, 
		   Float *psi, 
		   int cb,
		   int dag,
		   Wilson *wp)
{
  wfm_dslash(chi,u,psi,cb,dag);
}

CPS_END_NAMESPACE
