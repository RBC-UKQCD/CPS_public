#include <config.h>
/*!\file
  Wilson Dirac operator code for QCDOC

  $Id: wilson_init.C,v 1.4 2009/03/23 19:13:32 chulwoo Exp $
*/

#include <util/wilson.h>
#include <util/gjp.h>
#include <util/smalloc.h>
#include <util/verbose.h>
#include <util/error.h>
#include "wfm.h"

CPS_START_NAMESPACE

/* twisted mass fermion modifications                                   */
/*                                                                      */
/* - WilsonArg struct:  add ptr to spinor_tmp in                        */
/*                                                                      */
/* - wfm::init:   set *spinor_tmp to address of temporary spinor        */
/*                                                                      */
/* - wilson_compat_init.C:  arguments are prt to the WilsonArg struct   */
/* and ptr to the Wilson struct; its uses *spinor_tmp from WilsonArg    */
/* to define temporary checkerboards af[0] and af[1] in Wilson struct;  */
/* af[0] and af[1] are used in d_op_wilson and d_op_wilsonTm            */
/*                                                                      */
/*----------------------------------------------------------------------*/
/* Pointer to Wilson type structure                                     */
/*----------------------------------------------------------------------*/
static WilsonArg wil;

void wilson_set_sloppy(bool sloppy){
  wil.SloppyPrecision = sloppy;
}

void wilson_end(Wilson *wilson_p)
{
  VRB.Func("","wilson_end()");
  wilson_compat_end(wilson_p);
  wfm_vec_end(&wil); 
}

void wilson_init(Wilson *wilson_p)  
{
  char *cname = "wfm";
  char *fname = "wilson_init(Wilson *)";


  wil.CoreCount=1;
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

#if 0
  //Temporary hack to force non-local comms
  wil.local_comm[0] = 0;
  wil.local_comm[1] = 0;
  wil.local_comm[2] = 0;
  wil.local_comm[3] = 0;
#endif

  wfm_vec_init (&wil);
//printf("wilson_init::spinor_tmp %p\n", wil.spinor_tmp);
//~~ twisted mass fermions:  added second argument WilsonArg *wil
//~~ for transfer of *spinor_tmp from WilsonArg to Wilson
  wilson_compat_init(wilson_p, &wil);

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

//void wfm_dslash_two( Float *chi0, Float *chi1,
//                     Float *u,
//                     Float *psi0, Float *psi1,
//                     int cb0, int cb1, int dag)
/*Wrapper used by the DWF code*/
void wilson_dslash_two(Float *chi0, Float *chi1,
		   Float *u, 
		   Float *psi0, Float *psi1, 
		   int cb0, int cb1,
		   int dag,
		   Wilson *wp)
{
  wfm_dslash_two(chi0,chi1,u,psi0,psi1,cb0,cb1,dag);
}

CPS_END_NAMESPACE
