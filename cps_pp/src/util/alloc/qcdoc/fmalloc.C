#include<config.h>
/*!\file
  \brief  Implementation of dynamic memory management routines.	

  $Id: fmalloc.C,v 1.6 2004-10-14 22:09:59 chulwoo Exp $
*/

#include <util/error.h>
#include <util/verbose.h>
#include <qcdoc_align.h>
#include <qalloc.h>

CPS_START_NAMESPACE

void* fmalloc(int request){
    void *p = qalloc(QFAST,request);
    if (!p) p = qalloc(QCOMMS,request);
    return p;
}

void ffree(void* p){
  qfree(p);
}



CPS_END_NAMESPACE
