#include<config.h>
/*!\file
  \brief  Implementation of dynamic memory management routines.	

  $Id: pmalloc.C,v 1.4 2004-10-27 14:24:31 zs Exp $
*/

#include <util/verbose.h>
#include <util/error.h>
#include <stdlib.h>
#include <util/smalloc.h>

CPS_START_NAMESPACE

void* pmalloc(size_t request){
    return smalloc(request, "", "pmalloc");
}

void pfree(void* p){
    VRB.Pfree("","pfree","",p);
    free(p);
}

void pclear(void){};





CPS_END_NAMESPACE
