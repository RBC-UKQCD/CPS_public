#include<config.h>
/*!\file
  \brief  Implementation of dynamic memory management routines.	

  $Id: pmalloc.C,v 1.7 2004-10-27 14:24:31 zs Exp $
*/

#include <util/error.h>
#include <util/verbose.h>
#include <qcdoc_align.h>
#include <qalloc.h>
#include <util/smalloc.h>

CPS_START_NAMESPACE

void* pmalloc(size_t request){
    return smalloc(request, "", "pmalloc");
}

void pfree(void* p){
    VRB.Pfree("","pfree","",p);
    qfree(p);
}

void pclear(void){};





CPS_END_NAMESPACE
