#include<conf.h>
#include<config.h>
/*!\file
  \brief  Implementation of dynamic memory management routines.	

*/

#include <util/error.h>
#include <util/verbose.h>
#include <stdlib.h>
#include <util/smalloc.h>


CPS_START_NAMESPACE

void* smalloc( const char cname[], const char fname[], const char vname[], size_t request){

  void *p;
    if (request<=0)
	ERR.General(cname,fname,"smalloc requested with size %d!\n",request);
//#ifdef HAVE_POSIX_MEMALIGN
#if 0
#define ALLOC_MEMALIGN_NUM 512
  if( posix_memalign((void**)&p, ALLOC_MEMALIGN_NUM, request) ) ERR.Pointer(cname, fname, vname);
#else 
    p = malloc(request);
    if(!p) ERR.Pointer(cname, fname, vname);
#endif
    VRB.Smalloc(cname, fname, vname, p, request);
    return p;
}

void sfree( const char cname[], const char fname[], const char vname[], void* p ){
    VRB.Sfree(cname, fname, vname, p);
    free(p);
}

void sclear(void){};

CPS_END_NAMESPACE
