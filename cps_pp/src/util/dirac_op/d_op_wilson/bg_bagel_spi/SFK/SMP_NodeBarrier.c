/******************************************************************************/
/*                                                                            */
/*                             SMP MODE BARRIER                               */
/*                                                         30.06.07 S.Krieg   */
/******************************************************************************/

/* This code is largely due to P. Heidelberger, IBM Watson Research, IBM USA  */
/* Adapted (fixed to 4 threads) for the BGP SMP mode.                                  */

#include <pthread.h> 
#include <stdio.h>
#include <stdlib.h>
#include "SMP_NodeBarrier.h"

/* int pthread_barrier_init(pthread_barrier_t * barrier,  */
/* 			 const pthread_barrierattr_t * attr, unsigned count);  */

/* The pthread_barrier_init() function allocates the resources needed for the */ 
/* specified barrier and initializes that barrier with the attributes pointed */ 
/* to by the attr argument. When attr is NULL, the specified barrier is       */
/* initialized with the default attributes. You should always initialize a    */
/* barrier with this function before making any other reference to that       */
/* barrier.                                                                   */
/*                                                                            */
/* Parameters:                                                                */
/* attr                                                                       */
/*          Points to the attributes to be used when initializing a barrier   */
/*          with pthread_barrier_init().                                      */
/* barrier                                                                    */
/*          Points to the barrier to be destroyed by pthread_barrier_destroy()*/
/*          or initialized by pthread_barrier_init().                         */
/* count                                                                      */
/*          Indicates how many threads (minimum of one) must call the         */
/*          pthread_barrier_wait() function before any of those threads       */
/*          successfully return from the call.                                */


/* int pthread_barrier_destroy(pthread_barrier_t *barrier); */

/* Init the barrier structure contained in the var barrier. When the          */
/* pthread_barrier_init() function fails, the specified barrier is not        */
/* initialized.                                                               */
/*                                                                            */
/* Parameters:                                                                */
/* barrier                                                                    */
/*          Points to the barrier to be destroyed by pthread_barrier_destroy()*/
/*          or initialized by pthread_barrier_init().                         */


static pthread_barrier_t SMP_int_barrier; 
static pthread_barrierattr_t SMP_int_barrier_attr; 
static unsigned int not_initialized=1;

static int SMP_InitBarrier(void) 
{ 
    int rc=0; 
    unsigned int numThreads=4; 
    if (not_initialized){
	rc = pthread_barrier_init( &(SMP_int_barrier), 
				   &(SMP_int_barrier_attr), 
				   numThreads);
	not_initialized=0;
    }
    return rc; 

} 

void SMP_NodeBarrier(void)   
{ 
    int rc=0;
    if (not_initialized){
	printf("SMP_Barrier: Uninitialized call to SMP_Barrier. "
	       "Initializing.\n");
	rc = SMP_InitBarrier(); 
	if (rc!=0)
	{
	    printf("SMP_Barrier: FAILURE to initialize barrier.\n");
	    exit(rc);
	}
	not_initialized=0;
    }
    pthread_barrier_wait(&(SMP_int_barrier));   
    return; 
}

