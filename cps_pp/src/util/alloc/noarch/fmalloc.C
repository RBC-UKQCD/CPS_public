#include<config.h>
CPS_START_NAMESPACE

CPS_END_NAMESPACE
#include <util/smalloc.h>
#include <util/verbose.h>
#include <stdlib.h>
CPS_START_NAMESPACE

//! Allocate memory
/*!
  \param request The amount of memory (in bytes) to allocate
  \return A pointer to the allocated memory
*/
void* fmalloc(int request){
  void* ptr;
  ptr = malloc(request);
  return ptr;
}

//! Free allocate memory
/*!
  \param p Pointer to the memory to be freed.
*/
void ffree(void* p){
  free((char*) p);
}






CPS_END_NAMESPACE
