#include<config.h>
#include<string.h>
CPS_START_NAMESPACE

CPS_END_NAMESPACE
#include <util/smalloc.h>
#include <util/error.h>
#include <util/verbose.h>
#include <stdlib.h>
#include <stdio.h>
#include <qcdoc_align.h>
//#include <qalloc.h>
CPS_START_NAMESPACE

//static unsigned long extra_heap[1024*1024] LOCATE(".ddr_t_heap");
static unsigned long extra_heap[1024*1024] LOCATE(".dheap");

//! Allocate memory
/*!
  \param request The amount of memory (in bytes) to allocate
  \return A pointer to the allocated memory
*/
void* smalloc(int request){
  void* ptr;
//  ptr = qalloc(QCOMMS,request);
  ptr = malloc(request);
//  printf("request =%d ptr=%p\n",request,ptr);
  if (ptr == NULL){
  ERR.Pointer("","smalloc(i)","");
  exit(1);
  }
  VRB.Smalloc("","smalloc(i)","", ptr, request);
  bzero((char *)ptr,request);
  return ptr;
}

//! Free allocate memory
/*!
  \param p Pointer to the memory to be freed.
*/
void sfree(void* p){
  VRB.Sfree("","sfree(v*)","",p);
//  qfree((char*) p);
  free((char*) p);
}

void sclear(void){};





CPS_END_NAMESPACE
