#include<config.h>
/*!\file
  \brief  Implementation of dynamic memory management routines.	

  $Id: fmalloc.C,v 1.4 2004-09-07 05:21:48 chulwoo Exp $
*/

#include <stdlib.h>
#include <util/verbose.h>
#include <util/error.h>

CPS_START_NAMESPACE

void* fmalloc(int request){
    void *p =  malloc(request);
//    if(!p) ERR.Pointer("","fmalloc","");
//    VRB.Smalloc("","fmalloc","", p, request);
    return p;
}

void ffree(void* p){
//    VRB.Sfree("","ffree","",p);
    free(p);
}


CPS_END_NAMESPACE
