#include<config.h>
#include<string.h>
CPS_START_NAMESPACE

CPS_END_NAMESPACE
#include <util/pmalloc.h>
#include <util/error.h>
#include <util/verbose.h>
#include <stdlib.h>
#include <emalloc.h>
CPS_START_NAMESPACE

static unsigned long extra_heap[0*1048576] LOCATE(".dheap");

//! Allocate memory
/*!
  \param request The amount of memory (in bytes) to allocate
  \return A pointer to the allocated memory
*/
void* pmalloc(int request){
  void* ptr;
  ptr = malloc(request);
  VRB.Pmalloc("","pmalloc(i)","", ptr, request);
  if(ptr==NULL)
	ERR.Pointer("","pmalloc(i)","");
  bzero((char *)ptr,request);
  return ptr;
}

//! Free allocate memory
/*!
  \param p Pointer to the memory to be freed.
*/
void pfree(void* p){
  VRB.Pfree("","pfree(v*)","",p);
  free((char*) p);
}

void pclear(void){};





CPS_END_NAMESPACE
