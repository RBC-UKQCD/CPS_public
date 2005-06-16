#include<config.h>
/*!\file
  \brief  Implementation of dynamic memory management routines.	

  $Id: smalloc.C,v 1.12 2005-06-16 07:23:12 chulwoo Exp $
*/

#include <util/error.h>
#include <util/verbose.h>
#include <qcdoc_align.h>
#include <qalloc.h>

CPS_START_NAMESPACE

void* smalloc(size_t request,
	      const char *vname, const char *fname, const char *cname){
    void *p = qalloc(QCOMMS,request);
    if (!p) ERR.Pointer(cname, fname, vname);
    VRB.Smalloc(cname, fname, vname, p, request);
    return p;
}

void sfree(void* p, const char *vname, const char *fname, const char *cname){
    VRB.Sfree(cname, fname, vname, p);
    qfree(p);
}

void sclear(void){};





CPS_END_NAMESPACE
