#include<config.h>
/*!\file
  \brief  Implementation of dynamic memory management routines.	

  $Id: fmalloc.C,v 1.5 2004-09-07 05:21:48 chulwoo Exp $
*/

#include <util/error.h>
#include <util/verbose.h>
#include <qcdoc_align.h>
#include <qalloc.h>

CPS_START_NAMESPACE

void* fmalloc(int request){
    void *p = qalloc(QFAST,request);
    if (!p) p = qalloc(QCOMMS,request);
#if 0
    if (!p) ERR.Pointer("","fmalloc","");
    VRB.Smalloc("","fmalloc","", p, request);
#endif 
    return p;
}

void ffree(void* p){
//  VRB.Sfree("","ffree","",p);
  qfree(p);
}



CPS_END_NAMESPACE
