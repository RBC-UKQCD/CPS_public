#include<config.h>
/*!\file
  \brief  Implementation of dynamic memory management routines.	

  $Id: pmalloc.C,v 1.2 2006/12/14 17:53:51 chulwoo Exp $
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
