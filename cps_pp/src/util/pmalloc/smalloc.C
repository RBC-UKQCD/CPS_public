#include<config.h>
CPS_START_NAMESPACE

CPS_END_NAMESPACE
#include <util/smalloc.h>
#include <util/verbose.h>
#include <stdlib.h>
#include <emalloc.h>
CPS_START_NAMESPACE

static unsigned long extra_heap[1048576] LOCATE(".dheap");

//! Allocate memory
/*!
  \param request The amount of memory (in bytes) to allocate
  \return A pointer to the allocated memory
*/
void* smalloc(int request){
  void* ptr;
  ptr = malloc(request);
  VRB.Smalloc("","smalloc(i)","", ptr, request);
  return ptr;
}

//! Free allocate memory
/*!
  \param p Pointer to the memory to be freed.
*/
void sfree(void* p){
  VRB.Sfree("","sfree(v*)","",p);
  free((char*) p);
}

void sclear(void){};





CPS_END_NAMESPACE
