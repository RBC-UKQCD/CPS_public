#include<config.h>
CPS_START_NAMESPACE
/*----------------------------------------------------------*/
/*! The MPI comms request handle manager: mpi_requests.h

  This handles the array of MPI resuests.  Currently, it is very dumb
  and uses a fixed size (MPI_REQ_BASE_SIZE) array to store the
  handles.  If this is exceeded, it just crashes out.

  This should: 
  - Use a linked list.  
  - Store the direction as well as the handle, in case we wish to
  implement the form of SCUTransComplete that only waits for the
  transfers in a specified direction to cease.

  A.N.Jackson: ajackson@epcc.ed.ac.uk                      
  -----------------------------------------------------------
  CVS keywords
 
  $Author: mcneile $
  $Date: 2003-06-22 13:34:52 $
  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/comms/mpi_requests.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
  $Id: mpi_requests.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
  $Name: not supported by cvs2svn $
  $Locker:  $
  $RCSfile: mpi_requests.h,v $
  $Revision: 1.1.1.1 $
  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/comms/mpi_requests.h,v $
  $State: Exp $  */
/*----------------------------------------------------------*/

CPS_END_NAMESPACE
#include<config.h>
CPS_START_NAMESPACE

/* Allow the MPI stuff to be switched out, thus avoiding compiler
   errors (for the time being). */
#ifdef INCLUDE_MPI_SCU

#ifndef INCLUDED_MPI_REQ_MAN
#define INCLUDED_MPI_REQ_MAN

CPS_END_NAMESPACE
#include<comms/sysfunc.h>
CPS_START_NAMESPACE

#define MPI_REQ_BASE_SIZE 100

class MPIRequestManager {
 private:
    MPI_Request *mpi_req;
    int num_req, max_req;

 public:
    //! Default constructor: Initialises the request array:
    MPIRequestManager() {
	max_req = MPI_REQ_BASE_SIZE;
	mpi_req = new MPI_Request[max_req];
	num_req = 0;
	for( int i=0; i < max_req; i++ ) mpi_req[i] = MPI_REQUEST_NULL;
    }
    //! Default destructor: Deletes the storage associated with the request array:
    ~MPIRequestManager() {
	delete[] mpi_req;
    }

    //! Add a request handle to the array:
    void AddRequest( MPI_Request req ) {
	mpi_req[num_req] = req;
	num_req++;
	if( num_req >= max_req ) {
	    // This should do something more intelligent than this:
	    printf("ERROR: Maximum request-handle limit has been exceeded!\n");
	    exit(EXIT_FAILURE);
	}
    }

    //! Get the pointer to the request array:
    MPI_Request* ReqArray() {
	return mpi_req;
    }

    //! Get the number of requests:
    int NumReq() {
	return num_req;
    }

    //! Empty the array:
    void Clear() {
	for( int i=0; i < max_req; i++ ) mpi_req[i] = MPI_REQUEST_NULL;
	num_req = 0;
    }

};

#endif

#endif /* INCLUDE_MPI_SCU */

CPS_END_NAMESPACE
