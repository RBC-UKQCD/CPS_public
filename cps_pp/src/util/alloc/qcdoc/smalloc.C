#include<config.h>
/*!\file
  \brief  Implementation of dynamic memory management routines.	

  $Id: smalloc.C,v 1.8 2004-09-07 05:21:48 chulwoo Exp $
*/

#include <util/error.h>
#include <util/verbose.h>
#include <qcdoc_align.h>
#include <qalloc.h>

CPS_START_NAMESPACE

void* smalloc(int request){
    void *p = qalloc(QCOMMS,request);
#if 0
    if (!p) ERR.Pointer("","smalloc","");
    VRB.Smalloc("","smalloc","", p, request);
#endif
    return p;
}

void sfree(void* p){
#if 0
    VRB.Sfree("","sfree","",p);
#endif
    qfree(p);
}

void sclear(void){};





CPS_END_NAMESPACE
