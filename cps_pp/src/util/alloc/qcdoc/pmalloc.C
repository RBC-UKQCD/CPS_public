#include<config.h>
#include<string.h>
CPS_START_NAMESPACE

CPS_END_NAMESPACE
#include <util/pmalloc.h>
#include <util/error.h>
#include <util/verbose.h>
#include <stdlib.h>
#include <qcdoc_align.h>
#include <qalloc.h>
CPS_START_NAMESPACE

//static unsigned long extra_heap[0*1048576] LOCATE(".ddr_t_heap");

//! Allocate memory
/*!
  \param request The amount of memory (in bytes) to allocate
  \return A pointer to the allocated memory
*/
void* pmalloc(int request){
  void* ptr;
  ptr = qalloc(QCOMMS,request);
  VRB.Pmalloc("","pmalloc(i)","", ptr, request);
  if(ptr==NULL)
	ERR.Pointer("","pmalloc(i)","");
  return ptr;
}

//! Free allocate memory
/*!
  \param p Pointer to the memory to be freed.
*/
void pfree(void* p){
  VRB.Pfree("","pfree(v*)","",p);
  qfree((char*) p);
}

void pclear(void){};





CPS_END_NAMESPACE
