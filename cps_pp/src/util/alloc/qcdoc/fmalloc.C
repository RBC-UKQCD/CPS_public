#include<config.h>
#include<string.h>
CPS_START_NAMESPACE

CPS_END_NAMESPACE
#include <util/smalloc.h>
#include <util/error.h>
#include <util/verbose.h>
#include <stdlib.h>
#include <util/qcdio.h>
#include <qcdoc_align.h>
#include <qalloc.h>
CPS_START_NAMESPACE

//! Allocate memory
/*!
  \param request The amount of memory (in bytes) to allocate
  \return A pointer to the allocated memory
*/
void* fmalloc(int request){
  void* ptr;
  ptr = qalloc(QFAST,request);
  if (ptr == NULL)
  ptr = qalloc(QCOMMS,request);
  if (ptr == NULL){
  ERR.Pointer("","smalloc(i)","");
  exit(1);
  }
  VRB.Smalloc("","fmalloc(i)","", ptr, request);
  return ptr;
}

//! Free allocate memory
/*!
  \param p Pointer to the memory to be freed.
*/

void ffree(void* p){
  VRB.Sfree("","sfree(v*)","",p);
  qfree((char*) p);
}



CPS_END_NAMESPACE
