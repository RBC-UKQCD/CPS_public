#include<config.h>
CPS_START_NAMESPACE
/*----------------------------------------------------------*/
/*!\file
  \brief  Declaration and definition of the MPIRequestManager class

  $Id: mpi_requests.h,v 1.3 2003-10-27 16:45:38 zs Exp $
*/
/*----------------------------------------------------------*/
/* The MPI comms request handle manager: mpi_requests.h

  A.N.Jackson: ajackson@epcc.ed.ac.uk                      
  -----------------------------------------------------------
  CVS keywords
 
  $Author: zs $
  $Date: 2003-10-27 16:45:38 $
  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/comms/mpi_requests.h,v 1.3 2003-10-27 16:45:38 zs Exp $
  $Id: mpi_requests.h,v 1.3 2003-10-27 16:45:38 zs Exp $
  $Name: not supported by cvs2svn $
  $Locker:  $
  $RCSfile: mpi_requests.h,v $
  $Revision: 1.3 $
  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/comms/mpi_requests.h,v $
  $State: Exp $  */
/*----------------------------------------------------------*/

CPS_END_NAMESPACE
CPS_START_NAMESPACE


#ifndef INCLUDED_MPI_REQ_MAN
#define INCLUDED_MPI_REQ_MAN

CPS_END_NAMESPACE
#include <comms/sysfunc.h>
CPS_START_NAMESPACE

//! Maximum number of concurrent MPI requests.
#define MPI_REQ_BASE_SIZE 100 

/*!
  \defgroup mpicomms Code specific to the MPI interface
  \ingroup comms
*/

//! This handles the MPI requests.
/*!
  The MPI requests are used as identifying handles in non-blocking
  communications; see MPI literature for more more details, <e>e.g.</e>
  http://www-unix.mcs.anl.gov/mpi/mpi-standard/mpi-report-1.1/node44.htm#Node44
  
  Currently, it is very dumb and uses a fixed size (#MPI_REQ_BASE_SIZE)
  array to store the handles.  If this is exceeded, it just crashes out.

  \todo This should  
- Use a linked list.  
- Store the direction as well as the handle, in case we wish to
  implement the form of SCUTransComplete that only waits for the
  transfers in a specified direction to cease.
  
  \ingroup comms mpicomms
*/

class MPIRequestManager {
 private:
    MPI_Request *mpi_req;
    int num_req, max_req;

 public:

    //! Default constructor:
    /*!
      Initialises the list of requests.
      \post All MPI requests are initially set to MPI_REQUEST_NULL.
    */
    MPIRequestManager() {
	max_req = MPI_REQ_BASE_SIZE;
	mpi_req = new MPI_Request[max_req];
	num_req = 0;
	for( int i=0; i < max_req; i++ ) mpi_req[i] = MPI_REQUEST_NULL;
    }
    //! Default destructor.
    /*! Deletes the storage associated with the request array. */
    ~MPIRequestManager() {
	delete[] mpi_req;
    }

    //! Store a request handle.
    /*!
      \param req The request handle.
      \post The request is stored.
    */
    void AddRequest( MPI_Request req ) {
	mpi_req[num_req] = req;
	num_req++;
	if( num_req >= max_req ) {
	    // This should do something more intelligent than this:
	    printf("ERROR: Maximum request-handle limit has been exceeded!\n");
	    exit(EXIT_FAILURE);
	}
    }

    //! Get the list of requests.
    /*!
      \return A pointer to the first request.
    */
    MPI_Request* ReqArray() {
	return mpi_req;
    }

    //! Get the number of requests currently stored.
    /*!
      \return The number of stored handles.
    */
    int NumReq() {
	return num_req;
    }

    //! Empty the list.
    void Clear() {
	for( int i=0; i < max_req; i++ ) mpi_req[i] = MPI_REQUEST_NULL;
	num_req = 0;
    }

};

#endif

CPS_END_NAMESPACE
