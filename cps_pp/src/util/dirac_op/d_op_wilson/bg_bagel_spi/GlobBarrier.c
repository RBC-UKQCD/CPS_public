/******************************************************************************/
/*                                                                            */
/*                         Global Barriers using SPI                          */
/*                                                         30.06.07 S.Krieg   */
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <spi/bgp_SPI.h>
#include "NodeBarrier.h"
#include "GlobBarrier.h"

/** \brief Barrier using a given channel that calls a waitproc while idle
*
* Upon entry the given channel must expected be raised.  If not,
* the barrier will wait for it to raise.  It is assumed a low channel is
* from a previous barrier which is returning to idle state.
*
* This is the general case barrier.  While the barrier is waiting for a state
* change in the given global interrupt channel, it repeatedly calls the
* optional waitproc(waitarg).
*
* The waitproc must return an integer value.  A return of zero means the
* barrier should continue.  A non-zero value will break the barrier and the
* barrier function itself will return this same non-zero value (application
* defined).  This is useful for processing timeouts, etc.
*
* A typical use with no wait processing would be:
*      GlobInt_Barrier(0, NULL, NULL);
*
* \verbatim A typical use with a timeout counting iterations would be:
      int timeoutproc(void *arg)
      {
          int *countp = (int *)arg;
          if (++(*countp) > TIMEOUT)
              return 1;
          return 0;
      }
      int count = 0;
      int ret = GlobInt_Barrier(0, timeoutproc, &count);
      if (ret == 1)
              printf("barrier timeout!\n");
\endverbatim
* \param channel       Channel to use (e.g. 0..3)
* \param waitproc      Function to call while waiting (may be NULL)
* \param waitarg       Argument to pass to waitproc
*
* \returns 0 on success, else value of waitproc which breaks the barrier.
*/
/* __INLINE__ uint32_t GlobInt_Barrier(uint32_t channel, GlobInt_WaitProc_t waitproc, void* waitarg) */


/******************************************************************************/
/* Barrier using GLOBINT_CHANNEL_APP_BARRIER hardcoded.                       */
/* No waitproc (see above) is used.                                           */
/******************************************************************************/

int GlobBarrier(void){
    int rc=0;
    int pid = Kernel_PhysicalProcessorID();

    NodeBarrier();
    if (!pid){
	rc = GlobInt_Barrier(GLOBINT_CHANNEL_APP_BARRIER, NULL, NULL);
    }
    NodeBarrier();
    return rc;
}
