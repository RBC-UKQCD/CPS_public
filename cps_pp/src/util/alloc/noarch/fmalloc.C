#include<config.h>
/*!\file
  \brief  Implementation of dynamic memory management routines.	

  $Id: fmalloc.C,v 1.8 2012-03-26 13:50:11 chulwoo Exp $
*/

#include <stdlib.h>
#include <util/verbose.h>
#include <util/error.h>
#include <util/smalloc.h>
#include <stdlib.h> //  posix_memalign

CPS_START_NAMESPACE

#define ALLOC_MEMALIGN_NUM 4096
void* fmalloc(size_t request,
	      const char vname[], const char fname[], const char cname[]){

    return smalloc(request, vname, fname, cname);
    
}

void* fmalloc(size_t request){
    return smalloc(request);
}

void ffree(void* p, const char vname[], const char fname[], const char cname[]){
    sfree(p, cname, fname, vname);
}
void ffree(void* p){
    sfree(p);
}


CPS_END_NAMESPACE
