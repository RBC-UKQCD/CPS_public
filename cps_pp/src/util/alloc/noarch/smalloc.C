#include<config.h>
/*!\file
  \brief  Implementation of dynamic memory management routines.	

  $Id: smalloc.C,v 1.4 2004-10-27 14:24:31 zs Exp $
*/

#include <util/error.h>
#include <util/verbose.h>
#include <stdlib.h>

CPS_START_NAMESPACE

void* smalloc(size_t request,
	      const char *vname="", const char *fname="smalloc", const char *cname=""){
    void *p = malloc(request);
    if(!p) ERR.Pointer(cname, fname, vname);
    VRB.Smalloc(cname, fname, vname, p, request);
    return p;
}

void sfree(void* p, const char *vname, const char *fname, const char *cname){
    VRB.Sfree(cname, fname, vname, p);
    free(p);
}



void sclear(void){};





CPS_END_NAMESPACE
