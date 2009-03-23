/*! \file

$Id: cps_compat.C,v 1.3 2009-03-23 19:13:32 chulwoo Exp $
*/
	
  
#include <util/wilson.h>
#include <util/gjp.h>
#include <util/smalloc.h>
#include <util/verbose.h>
#include <util/error.h>
#include "wfm.h"
CPS_START_NAMESPACE

/*----------------------------------------------------------------------*/
/* for backward dwf and clover compatibility in CPS                     */
/*----------------------------------------------------------------------*/
void wilson_compat_init(Wilson *wilson_p);
void wilson_compat_end(Wilson *wilson_p);

void wilson_compat_init(Wilson *wilson_p, WilsonArg *wil)
{
  int size;
  char *cname = "wfm";
  char *fname = "wilson_compat_init(Wilson *)";
  int spinor_words;

        /*--------------------------------------------------------------------------*/
        /* Reserve memory for the node sublattice sizes                             */
        /*--------------------------------------------------------------------------*/
          size = 4*sizeof(int);
          wilson_p->ptr = (int *) smalloc(size);
          if( wilson_p->ptr == 0)
            ERR.Pointer(cname,fname, "ptr");
          VRB.Smalloc(cname,fname,
                      "ptr", wilson_p->ptr, size);

        /*--------------------------------------------------------------------------*/
        /* Set the node sublattice sizes                                            */
        /*--------------------------------------------------------------------------*/
          wilson_p->ptr[0] = GJP.XnodeSites();
          wilson_p->ptr[1] = GJP.YnodeSites();
          wilson_p->ptr[2] = GJP.ZnodeSites();
          wilson_p->ptr[3] = GJP.TnodeSites();
          wilson_p->vol[0] = wilson_p->ptr[0] * wilson_p->ptr[1] *
                             wilson_p->ptr[2] * wilson_p->ptr[3] / 2;
          wilson_p->vol[1] = wilson_p->vol[0];
// TB modified:
//	  wilson_p->af[0] = NULL;
//	  wilson_p->af[1] = NULL;
//~~ 
//~~ twisted mass fermions:  defines half-fields Wilson.af[0] & Wilson.af[1]
//~~ in temporary spinor temporary full field 
//~~ af[0], af[1] used in d_op_wilson & d_op_wilsonTm members 
//~~ MatPc, MatPcDag, and MatPcDagMatPc
//~~ 
          spinor_words = SPINOR_SIZE * wilson_p->vol[0];
          wilson_p->af[0] = wil->spinor_tmp;
          wilson_p->af[1] = wil->spinor_tmp + spinor_words;
//printf("set af[] temps. size= %d\n",spinor_words);
//printf("cps_compat_init::af[0]=%p\n",wilson_p->af[0]);
//printf("cps_compat_init::af[1]=%p\n",wilson_p->af[1]);
//~~
// End TB
/*----------------------------------------------------------------------*/
/* end of for backward and dwf compatibility                            */
/*----------------------------------------------------------------------*/

}
void wilson_compat_end(Wilson *wilson_p)
{
/*----------------------------------------------------------------------*/
/* for backward and dwf compatibility                                   */
/*----------------------------------------------------------------------*/
  sfree(wilson_p->ptr);
}

CPS_END_NAMESPACE
