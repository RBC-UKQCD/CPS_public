#include<config.h>
/*!\file
  \brief  Implementation of dynamic memory management routines.	

  $Id: smalloc.C,v 1.7 2004-09-02 16:57:09 zs Exp $
*/

#include <util/error.h>
#include <util/verbose.h>
#include <qcdoc_align.h>
#include <qalloc.h>

CPS_START_NAMESPACE

void* smalloc(int request){
    void *p = qalloc(QCOMMS,request);
    if (!p) ERR.Pointer("","smalloc","");
    VRB.Smalloc("","smalloc","", p, request);
    return p;
}

void sfree(void* p){
    VRB.Sfree("","sfree","",p);
    qfree(p);
}

void sclear(void){};





CPS_END_NAMESPACE
