#include<config.h>
#include<util/qcdio.h>
//These are routines not provided for by 
#include<mpi.h>
//-------------------------------------------------------------------
/*!\file
  \brief Definition of glb_sum routine.
*/
#include<comms/glb.h>
#include<comms/scu.h>
#include<util/gjp.h>
#include<comms/double64.h>
#include <comms/sysfunc_cps.h>
#include <comms/glb_sum_internal.h>
CPS_START_NAMESPACE
#ifndef USE_QMP
#define USE_QMP
#endif

    int glb_sum (long *send, const long n_elem)
    {
      long recv[n_elem];
      int ret= MPI_Allreduce ((long *) send, recv, n_elem, MPI_LONG, MPI_SUM,
			    MPI_COMM_WORLD);
	memcpy(send,recv,n_elem*sizeof(long));
	return ret;
    }

    int glb_sum (uint32_t * send,
			 const long n_elem)
    {
      uint32_t recv[n_elem];
      int ret= MPI_Allreduce ((uint32_t *) send, recv, n_elem, MPI_UNSIGNED,
			    MPI_SUM, MPI_COMM_WORLD);
	memcpy(send,recv,n_elem*sizeof(uint32_t));
	return ret;
    }

#if 1
    int glb_sum (int *send, const long n_elem)
    {
      int recv[n_elem];
      int ret= MPI_Allreduce ((int *) send, recv, n_elem, MPI_INT, MPI_SUM,
			    MPI_COMM_WORLD);
	memcpy(send,recv,n_elem*sizeof(int));
	return ret;
    }

    int glb_sum (double *send, const long n_elem)
    {
#if 1
      double recv[n_elem];
      int ret= MPI_Allreduce ((double *) send, recv, n_elem, MPI_DOUBLE, MPI_SUM,
			    MPI_COMM_WORLD);
	memcpy(send,recv,n_elem*sizeof(double));
#else
	QMP_sum_double_array(recv,n_elem);
#endif
        return 1;
    }
#endif
CPS_END_NAMESPACE
