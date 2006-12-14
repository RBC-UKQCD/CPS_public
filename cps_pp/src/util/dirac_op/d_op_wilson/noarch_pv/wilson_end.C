#include<config.h>
CPS_START_NAMESPACE
/****************************************************************************/
/* 10/6/96                                                                  */
/*                                                                          */
/* wilson_end:                                                              */
/*                                                                          */
/* This routine frees any memory that was allocated by wilson_init          */
/*                                                                          */
/****************************************************************************/

/*--------------------------------------------------------------------------*/
/* Include header files                                                     */
/*--------------------------------------------------------------------------*/
CPS_END_NAMESPACE
#include<util/wilson.h>
#include<util/smalloc.h>
#include<util/verbose.h>
CPS_START_NAMESPACE


void wilson_end( Wilson *wilson_p)
{
  int i;
  char *cname = " ";
  char *fname = "wilson_end(Wilson*)";
  VRB.Func(cname,fname);

  for(i=0; i<4; i++){
    VRB.Sfree(cname,fname, "af[i]", wilson_p->af[i]);
    sfree(wilson_p->af[i]);
    VRB.Sfree(cname,fname, "ab[i]", wilson_p->ab[i]);
    sfree(wilson_p->ab[i]);
  }

  VRB.Sfree(cname,fname, "spinor_tmp", wilson_p->spinor_tmp);
  sfree(wilson_p->spinor_tmp);

  VRB.Sfree(cname,fname, "ptr", wilson_p->ptr);
  sfree(wilson_p->ptr);

}
CPS_END_NAMESPACE
