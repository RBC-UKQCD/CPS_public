#include <util/wilson.h>
#include <util/gjp.h>
#include <util/smalloc.h>
#include <util/verbose.h>
#include <util/error.h>
#include <util/wfm.h>

/*----------------------------------------------------------------------*/
/* Pointer to Wilson type structure is ignored here                     */
/*----------------------------------------------------------------------*/
static WilsonArg wil;
void wilson_end(Wilson *wilson_p)
{
  wfm_end(&wil); 
}

void wilson_init(Wilson *wilson_p)  
{

  char *cname = "wfm";
  char *fname = "wilson_init(Wilson)";

  VRB.Func(cname,fname);

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

  wil.local_comm[0] = 0;
  wil.local_comm[1] = 0;
  wil.local_comm[2] = 0;
  wil.local_comm[3] = 0;

  wfm_init (&wil);

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
