#include <config.h>
/*!
  \file
  \brief Functions for dynamical array allocation.

  $Id: amalloc.C,v 1.4 2004-09-02 16:59:38 zs Exp $
*/
#include <util/smalloc.h>
#include <stdarg.h>
#include <stddef.h>
#include <util/error.h>

CPS_START_NAMESPACE

// ------------------------------------------------------------ 
/*!
  This is a utility routine for amalloc().
  \param size   Size (in bytes) of objects in the array.
  \param n_dim  Dimension of this subarray.
  \param n_ptr  Number of pointers we need to set.
  \param prev  Location of the pointers.
  \param start Location of this subarray.
  \param dimension The array dimensions.
  \param index  Recursion depth - the array index of this subarray
*/
//  ------------------------------------------------------------

void subarray(size_t size, int n_dim, int n_ptr,
	      void ***prev, void **start, int *dimension, int index){

    int i, dim = dimension[index];

    if(n_dim>0){    // Set up pointers to pointers  

	for(i=0; i<n_ptr; i++) prev[i] = start + i*dim;

	subarray(size, n_dim-1, n_ptr*dim, 
		 (void***)start,	// Pointers to be initialised 
		 start+n_ptr*dim,	// Pointing to locations starting here 
		 dimension, index+1);

    }else{		// Last recursion; set up pointers to data   

	for(i=0; i<n_ptr; i++) (char*)prev[i] = (char*)start+i*dim*size;
	
    }

}

            


/*!
  Dynamically allocates memory for a (pseudo) array.

  The array consists of a single contiguous block of data, so
  a single 'free' call will release it.  

  The data is stored in the usual C order.

  \param size The size (in bytes) of the elements of the array.
  \param n_dim The dimensionality of the array.
  \param ... The dimensions of the array.
  \return A pointer to the allocated memory
*/

void *amalloc(size_t size, int n_dim, ...){

/*
  The first bit of data in the array are pointers which point to the data
  in the next dimension, which are pointers to data in the next dimension,
  and so on, until the penultimate dimension, where the pointers point to
  the actual data.
*/

    int *dimension = (int*)smalloc(n_dim*sizeof(int));
    if(!dimension) ERR.Pointer("", "amalloc", "dimension");

    // Count the pointers and data elements in the array. 

    va_list ap;
    va_start(ap, n_dim); 

    int d, n_ptr = 0, n_data = 1;
    for(d = 0; d<n_dim; d++){
	dimension[d] = va_arg(ap, int); 
	n_data *= dimension[d];
	if(d<n_dim-1) n_ptr  += n_data;
    }

    va_end(ap);   

    // Allocate the memory for the data and the pointers. 

    void **start = (void**)smalloc((size_t)(n_data*size)+n_ptr*sizeof(void*));
    if(!start) return 0;//ERR.Pointer("", "amalloc", "start");

    // Set up the pointers.
		   
    void **p;
    subarray(size, n_dim-1, 1, &p, start, dimension, 0);

    sfree(dimension);
    return (void*)start;

}


CPS_END_NAMESPACE
