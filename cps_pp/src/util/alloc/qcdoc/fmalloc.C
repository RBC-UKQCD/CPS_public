#include<config.h>
/*!\file
  \brief  Implementation of dynamic memory management routines.	

  $Id: fmalloc.C,v 1.8 2004-10-27 15:16:38 zs Exp $
*/

#include <util/error.h>
#include <util/verbose.h>
#include <qcdoc_align.h>
#include <qalloc.h>
#include <util/smalloc.h>

CPS_START_NAMESPACE

void* fmalloc(size_t request,
	      const char *vname, const char *fname, const char *cname){
    void *p = qalloc(QFAST,request);
    if (!p) return smalloc(request, cname, fname, vname);
    VRB.Smalloc(cname, fname, vname, p, request);
    return p;
}

void ffree(void* p,
	   const char *vname, const char *fname, const char *cname){
    sfree(p, cname, fname, vname);
}



CPS_END_NAMESPACE
