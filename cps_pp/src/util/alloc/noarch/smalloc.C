#include<config.h>
/*!\file
  \brief  Implementation of dynamic memory management routines.	

  $Id: smalloc.C,v 1.3 2004-09-02 16:58:10 zs Exp $
*/

#include <util/error.h>
#include <util/verbose.h>
#include <stdlib.h>

CPS_START_NAMESPACE

void* smalloc(int request){
    void *p = malloc(request);
    if(!p) ERR.Pointer("","smalloc","");
    VRB.Smalloc("","smalloc","", p, request);
    return p;
}

void sfree(void* p){
    VRB.Sfree("","sfree","",p);
    free(p);
}



void sclear(void){};





CPS_END_NAMESPACE
