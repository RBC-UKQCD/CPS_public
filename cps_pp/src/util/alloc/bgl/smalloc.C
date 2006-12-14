#include<config.h>
CPS_START_NAMESPACE
/*--------------------------------------------------------------------*/
/*!\file
  \brief  Dynamic memory management routines.	
*/
/*--------------------------------------------------------------------*/

CPS_END_NAMESPACE
#include <util/smalloc.h>
#include <util/verbose.h>
#include <util/error.h>
#include <stdlib.h>
CPS_START_NAMESPACE


//--------------------------------------------------------------------
// If BGL use MEM_ALIGN byte alligned memory allocation
//--------------------------------------------------------------------

// Define the allignment boundary in bytes
#define MEM_ALIGN 32

// Define the maximum number of pointers that can be allocated at
// any single time
#define MAX_PTR_NUM 1000

// Define the list tag that identifies a free place in the array
#define FREE_INDEX 0xFFFFFFFF

// Define a first time lock that will initialize the list arrays
int smalloc_first_time = 1;

// This array hold all the pointers allocated by malloc
unsigned smalloc_addr_list[MAX_PTR_NUM];

// This array hold all the pointers allocated by malloc after 
// they are alligned backwards on an MEM_ALIGN byte boundary
unsigned ma_smalloc_addr_list[MAX_PTR_NUM];


//! Allocate memory
/*!
  \param request The amount of memory (in bytes) to allocate
  \return A pointer to the allocated memory
*/

//--------------------------------------------------------------------
// If BGL use MEM_ALIGN byte alligned memory allocation
//--------------------------------------------------------------------

//-------------------------------------------------------------------------
// qmalloc (MEM_ALIGN byte alligned)
//-------------------------------------------------------------------------

//! Allocate memory
/*!
  \param request The amount of memory (in bytes) to allocate
  \return A pointer to the allocated memory
*/
void* smalloc(int request){
  int index;
  unsigned addr;
  unsigned ma_addr;
  void* ptr;
  void* ma_ptr;

  // If this is the first call the initialize the list arrays
  // and set smalloc_first_time=0
  if(smalloc_first_time == 1){
    smalloc_first_time = 0;
    for(index =0; index < MAX_PTR_NUM; index++){
      smalloc_addr_list[index] = FREE_INDEX;
      ma_smalloc_addr_list[index] = FREE_INDEX;
    }
  }

  // Allocate requested memory plus an extra MEM_ALIGN Bytes
  ptr =  malloc(request+MEM_ALIGN);
  
  // Calculate a pointer that is alligned in the previous
  // MEM_ALIGN Byte boundary
  addr     = (unsigned) ptr;
  ma_addr = addr + MEM_ALIGN - (addr%MEM_ALIGN);
  ma_ptr  = (void *) ma_addr;

  // Log ptr and ma_ptr in the same index of the two arrays.
  for(index =0; index < MAX_PTR_NUM; index++){
    // Find the first free slot
    if(smalloc_addr_list[index] == FREE_INDEX){
      smalloc_addr_list[index] = addr;
      ma_smalloc_addr_list[index] = ma_addr;
      break;
    }
  }

  if(index == MAX_PTR_NUM){
    ERR.General("smalloc","smalloc","No more room in the MEM_ALIGN word allign list\n");
  }

  VRB.Smalloc("","smalloc(i)","", ma_ptr, request);
  return ma_ptr;
}

//! Free allocate memory
/*!
  \param p Pointer to the memory to be freed.
*/
void sfree(void* p){
  int index;
  unsigned ma_addr;
  unsigned addr;
  void *ptr;

  VRB.Sfree("","sfree(v*)","",p);

  ma_addr = (unsigned) p;

  // Find the original (non-MEM_ALIGN byte-alligned pointer from the list
  for(index =0; index < MAX_PTR_NUM; index++){
    if(ma_smalloc_addr_list[index] == ma_addr){
      addr = smalloc_addr_list[index];
      ptr = (void *) addr;
      ma_smalloc_addr_list[index] = FREE_INDEX;
      smalloc_addr_list[index]     = FREE_INDEX;
      break;
    }
  }

  if(index == MAX_PTR_NUM){
    ERR.General("sfree","sfree","pointer not found in the MEM_ALIGN byte allign list\n");
  }


  free((char*) ptr);
}


void* smalloc(size_t request,
	      const char *vname, const char *fname, const char *cname){
    void *p = malloc(request);
    if(!p) ERR.Pointer(cname, fname, vname);
    VRB.Smalloc(cname, fname, vname, p, request);
    return p;
}


void sfree(void* p, const char *vname, const char *fname, const char *cname){
    VRB.Sfree(cname, fname, vname, p);
    free(p);
}



void sclear(void){};





CPS_END_NAMESPACE
