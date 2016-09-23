/*!
  \file
  \brief Functions for dynamical array allocation.

  $Id: pt_amalloc.C,v 1.3 2009/08/28 17:20:49 chulwoo Exp $
*/
#include <stdarg.h>
#include <stddef.h>
#include <stdlib.h>



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


static const size_t alignment = 128;

static void subarray(size_t size, int n_dim, int n_ptr,
	      void ***prev, void **start, int *dimension, int index){

    int i, dim = dimension[index];

    if(n_dim>0){    // Set up pointers to pointers  

	for(i=0; i<n_ptr; i++) prev[i] = start + i*dim;

	subarray(size, n_dim-1, n_ptr*dim, 
		 (void***)start,	// Pointers to be initialised 
		 start+n_ptr*dim,	// Pointing to locations starting here 
		 dimension, index+1);

    }else{		// Last recursion; set up pointers to data   

	// The data should lie on a 'align'-byte boundary, so, if necessary,
	// we move 'start' to the next 'align'-byte boundary.
	// We allocated additional space for this using 'align_pad' in amalloc.
	
	if(0!=(long)start%alignment)
	    start = (void**)(((long)start/alignment+1)*alignment);

  	for(i=0; i<n_ptr; i++){ 
  	    char **previ = (char**)&prev[i]; 
  	    *previ = (char*)start+i*dim*size; 
  	} 
	
    }

}

            


/*!
  Dynamically allocates memory for a (pseudo) array.

  The array consists of a single contiguous block of data, so
  a single 'free' call will release it.  

  The data is stored in the usual C order.

  \param allocator Pointer to the function allocating the memory <em>viz.</em>
  ::smalloc or ::fmalloc.
  \param size The size (in bytes) of the elements of the array.
  \param n_dim The dimensionality of the array.
  \param ... The dimensions of the array.
  \return A pointer to the allocated memory
*/

void *pt_amalloc(void*  (*allocator)(size_t, const char *vname,
			  const char *fname, const char *cname),
	      size_t size, int n_dim, ...){

/*
  The first bit of data in the array are pointers which point to the data
  in the next dimension, which are pointers to data in the next dimension,
  and so on, until the penultimate dimension, where the pointers point to
  the actual data.
*/

    int *dimension = (int*)malloc(n_dim*sizeof(int));

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

    // The allocator will align the pointers correctly. We should
    // align the data on a boundary which is a multiple of 128 bytes so let's
    // pad the array so there is space to realign the data if necessary.

  
    // if the size of all the pointers happens to be a multiple of the
    //alignment length  then no padding is necessary here, right?
    size_t align;   
    if(n_ptr*sizeof(void*)%alignment==0) align = 0;
    else align = alignment;

    void **start = (void**)allocator(n_data*size+align+n_ptr*sizeof(void*),"","amalloc","");


    // Set up the pointers.
		   
    void **p;
    subarray(size, n_dim-1, 1, &p, start, dimension, 0);

//    VRB.Smalloc("","amalloc","", start, n_data*size+align+n_ptr*sizeof(void*));
    
    free(dimension);
    return (void*)start;

}

