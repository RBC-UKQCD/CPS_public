#include<config.h>
/*!\file
  \brief  Implementation of dynamic memory management routines.	

  $Id: pmalloc.C,v 1.6 2004-09-02 16:57:07 zs Exp $
*/

#include <util/error.h>
#include <util/verbose.h>
#include <qcdoc_align.h>
#include <qalloc.h>

CPS_START_NAMESPACE

void* pmalloc(int request){
    void* p = qalloc(QCOMMS,request);
    VRB.Pmalloc("","pmalloc","", p, request);
    if(!p) ERR.Pointer("","pmalloc","");
    return p;
}

void pfree(void* p){
    VRB.Pfree("","pfree","",p);
    qfree(p);
}

void pclear(void){};





CPS_END_NAMESPACE
