#include<config.h>
/*!\file
  \brief  Implementation of dynamic memory management routines.	

  $Id: pmalloc.C,v 1.3 2004-09-02 16:58:08 zs Exp $
*/

#include <util/verbose.h>
#include <util/error.h>
#include <stdlib.h>

CPS_START_NAMESPACE

void* pmalloc(int request){
    void *p = malloc(request);
    if(!p) ERR.Pointer("","pmalloc","");
    VRB.Pmalloc("","pmalloc","", p, request);
    return p;
}

void pfree(void* p){
    VRB.Pfree("","pfree","",p);
    free(p);
}

void pclear(void){};





CPS_END_NAMESPACE
