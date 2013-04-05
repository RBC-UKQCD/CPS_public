#include<conf.h>
#include<config.h>
/*!\file
  \brief  Implementation of dynamic memory management routines.	

  $Id: smalloc.C,v 1.11 2013-04-05 17:46:30 chulwoo Exp $
*/

#include <util/error.h>
#include <util/verbose.h>
#include <stdlib.h>
#include <util/smalloc.h>


CPS_START_NAMESPACE

void* smalloc(size_t request,
	      const char vname[], const char fname[], const char cname[]){

#ifdef HAVE_POSIX_MEMALIGN
#define ALLOC_MEMALIGN_NUM 4096
  void *p;
  if( posix_memalign((void**)&p, ALLOC_MEMALIGN_NUM, request) ) ERR.Pointer(cname, fname, vname);
    VRB.Smalloc(cname, fname, vname, p, request);
#else 
    if (request<=0)
	ERR.General(cname,fname,"smalloc requested with size %d!\n",request);
    void *p = malloc(request);
    if(!p) ERR.Pointer(cname, fname, vname);
#endif
    return p;
}

void* smalloc(size_t request){
    if (request<=0)
	ERR.General("","","smalloc requested with size %d!\n",request);
    void *p = malloc(request);
    if(!p) ERR.Pointer("", "", "");
    VRB.Smalloc("", "", "", p, request);
    return p;
}

void sfree(void* p, const char vname[], const char fname[], const char cname[]){
    VRB.Sfree(cname, fname, vname, p);
    free(p);
}

void sfree(void* p){
    VRB.Sfree("", "", "", p);
    free(p);
}

void sclear(void){};

CPS_END_NAMESPACE
