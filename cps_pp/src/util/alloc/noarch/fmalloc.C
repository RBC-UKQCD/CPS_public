#include<config.h>
/*!\file
  \brief  Implementation of dynamic memory management routines.	

  $Id: fmalloc.C,v 1.8 2012/03/26 13:50:11 chulwoo Exp $
*/

#include <stdlib.h>
#include <util/verbose.h>
#include <util/error.h>
#include <util/smalloc.h>
#include <stdlib.h> //  posix_memalign

CPS_START_NAMESPACE

#define ALLOC_MEMALIGN_NUM 4096
void* fmalloc( const char cname[], const char fname[], const char vname[], size_t request){

    return smalloc(cname,fname,vname,request);
    
}

#if 0
void* fmalloc(size_t request){
    return smalloc(request);
}
#endif

void ffree( const char cname[], const char fname[], const char vname[], void *p){
    sfree(cname, fname, vname,p);
}
#if 0
void ffree(void* p){
    sfree(p);
}
#endif


CPS_END_NAMESPACE
