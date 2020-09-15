#ifndef _smalloc_h
#define _smalloc_h                //!< Prevent multiple inclusion

#include<config.h>
#include <stdlib.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Declaration of dynamic memory management routines.	

*/


/*! \addtogroup mem_alloc Memory allocation
  @{
*/


//! Allocate memory
/*!
  With verbose reporting of allocation details, if the appropriate Verbose
  level is enabled.
  \param request The amount of memory (in bytes) to allocate
  \param vname The name of the variable pointing to the allocated memory.
  \param fname The name of the calling function or method
  \param cname The name of the calling class
  \return A pointer to the allocated memory
  \post Program exits with the appropriate Error code if allocation fails.  
*/

//! Allocate memory
/*!
  With verbose reporting of allocation details, if the appropriate Verbose
  level is enabled.
  \param cname The name of the calling class
  \param fname The name of the calling function or method
  \param vname The name of the variable pointing to the allocated memory.
  \param request The amount of memory (in bytes) to allocate
  \return A pointer to the allocated memory
  \post  Exits with the appropriate Error code if allocation fails.
*/

void* smalloc(const char cname[], const char fname[], const char vname[], size_t request);
inline void* smalloc(size_t request, const char vname[], const char fname[]="smalloc", const char cname[]=""){
return smalloc(cname,fname,vname,request);
}
inline void* smalloc(size_t request){ return smalloc("","","",request);}

//! Free allocate memory
/*!
  With verbose reporting of allocation details, if the appropriate Verbose
  level is enabled.
  \param p Pointer to the memory to be freed.
  \param vname The name of the variable pointing to the allocated memory.
  \param fname The name of the calling function or method
  \param cname The name of the calling class
  \post  Exits with the appropriate Error code if allocation fails.  
*/

//! Free memory
/*!
  With verbose reporting of details, if the appropriate Verbose
  level is enabled.
  \param cname The name of the calling class
  \param fname The name of the calling function or method
  \param vname The name of the variable pointing to the allocated memory.
  \param p Pointer to the memory to be freed.
*/
void sfree(const char cname[], const char fname[], const char vname[], void *p);
inline void sfree(void* p, const char vname[], const char fname[],const char cname[]){
   sfree(cname,fname,vname,p);
}
inline void sfree(void* p){ sfree ("","","",p);}

//! Doesn't appear to do anything.
void sclear();


//! Allocate memory
/*!
  Allocates in fast memory (EDRAM) on the QCDOC.
  If this fails, then allocation of transient DDR memory is attempted.
  With verbose reporting of allocation details, if the appropriate Verbose
  level is enabled.
  \param request The amount of memory (in bytes) to allocate
  \param vname The name of the variable pointing to the allocated memory.
  \param fname The name of the calling function or method
  \param cname The name of the calling class
  \return A pointer to the allocated memory
  \post  Exits with the appropriate Error code if allocation fails.  
*/
void* fmalloc(const char cname[], const char fname[], const char vname[], size_t request);
inline void* fmalloc(size_t request,
	      const char vname[], const char fname[]="fmalloc", const char cname[]=""){
      return fmalloc(cname,fname,vname,request);
}
inline void* fmalloc(size_t request){
      return fmalloc("","","",request);
}

//! Free allocate memory
/*!
  With verbose reporting of allocation details, if the appropriate Verbose
  level is enabled.  
  \param p Pointer to the memory to be freed.
  \param vname The name of the variable pointing to the allocated memory.
  \param fname The name of the calling function or method
  \param cname The name of the calling class
*/
void ffree(const char cname[], const char fname[], const char vname[], void *p);
inline void ffree(void* p, const char vname[], const char fname[]="ffree", const char cname[]=""){
     ffree(cname,fname,vname,p);
}
inline void ffree(void* p){
     ffree("","","",p);
}

//! Allocate memory
/*!
  With verbose reporting of allocation details, if the appropriate Verbose
  level is enabled.
  \param cname The name of the calling class
  \param fname The name of the calling function or method
  \param vname The name of the variable pointing to the allocated memory.
  \param request The amount of memory (in bytes) to allocate
  \return A pointer to the allocated memory
  \post  Exits with the appropriate Error code if allocation fails.
*/

//! Free memory
/*!
  With verbose reporting of details, if the appropriate Verbose
  level is enabled.
  \param cname The name of the calling class
  \param fname The name of the calling function or method
  \param vname The name of the variable pointing to the allocated memory.
  \param p Pointer to the memory to be freed.
*/

/*! @} */


CPS_END_NAMESPACE

#endif /* !_smalloc_h */
