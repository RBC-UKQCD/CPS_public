/* begin_generated_IBM_copyright_prolog                             */
/*                                                                  */
/* This is an automatically generated copyright prolog.             */
/* After initializing,  DO NOT MODIFY OR MOVE                       */
/*  --------------------------------------------------------------- */
/*                                                                  */
/* (C) Copyright IBM Corp.  2007, 2007                              */
/* IBM CPL License                                                  */
/*                                                                  */
/*  --------------------------------------------------------------- */
/*                                                                  */
/* end_generated_IBM_copyright_prolog                               */
/**
 * \file spi/GlobInt.h
 */


#include <common/namespace.h>

__BEGIN_DECLS
#include <common/linkage.h>
#include <bpcore/bgp_types.h>
#include <bpcore/ppc450_inlines.h>
#include <bpcore/bgp_global_ints.h>

#include <stdio.h>

/*! \brief There are 4 Global Interrupt channels, commonly used as follows:
 *           0: Application Compute Nodes Barrier
 *           1: Kernel All Nodes Alert for Power Throttling
 *           2: Kernel Compute Nodes Barrier
 *           3: Kernel All Nodes (Compute + I/O) Barrier
 */
#define GLOBINT_CHANNEL_APP_BARRIER                  (0) // User-space access
#define GLOBINT_CHANNEL_KERNEL_THROTTLE              (1) // Privileged access, Global OR Indicates Over-Power
#define GLOBINT_CHANNEL_KERNEL_COMPUTE_NODES_BARRIER (2) // Priviliged access
#define GLOBINT_CHANNEL_KERNEL_ALL_NODES_BARRIER     (3) // Priviliged access


/*
 * Note:  the hardware does not seem to allow mtdcrux/mfdcrux when in supervisor
 *        mode.  Kernel usage of this SPI requires the use of the mtdcrx/mfdcrx
 *        instructions.
 */

#ifdef __KERNEL__
#define WriteGlobIntDCR _bgp_mtdcrx
#define ReadGlobIntDCR  _bgp_mfdcrx
#else
#define WriteGlobIntDCR _bgp_mtdcrux
#define ReadGlobIntDCR  _bgp_mfdcrux
#endif


/*!\brief Information returned by status reads.  Note that userEnables returns the uptree state in diagnostic read mode.
*
*/
typedef union {
              unsigned gi_status_word;         /*!< Anonymous union to access as a word */
              struct {
                     unsigned mode:4;          /*!< Mode 0-3 (1=sticky, 0=direct) */
                     unsigned stickyInt:4;     /*!< Sticky mode interrupt 0-3 (1=active) */
                     unsigned stickyActive:4;  /*!< Sticky mode in progress 0-3 (1=true) */
                     unsigned stickyType:4;    /*!< Sticky mode type 0-3 (0=OR, 1=AND) */
                     unsigned reserved0:1;
                     unsigned accessviolation:1;/*!< Access violation (1=true) */
                     unsigned reserved1:1;
                     unsigned parityErr:1;     /*!< Parity error (1=true) */
                     unsigned giRecvState:4;   /*!< Receive state 0-3 ("downtree" state) */
                     unsigned userEnables:4;   /*!< User-level enables 0-3 (1=enabled). */
                                               /* Actual "uptree" state when readMode=1. */
                     unsigned giSendState:4;   /*!< Send state 0-3 ("uptree" state) */
                     };
               } GlobInt_Status_t;

/*!\brief Function prototype needed for GlobInt_Barrier
 * \see GlobInt_Barrier
 */
typedef int (*GlobInt_WaitProc_t)(void *arg);

/*!\brief Global interrupt usage modes.  Needed for GlobInt_Initialize
 * \see GlobInt_Initialize
 */
typedef enum
{
	GLOBINT_MODE_DIRECT=0,
	GLOBINT_MODE_STICKYOR,
	GLOBINT_MODE_STICKYAND
} GLOBINT_MODE_t;


/*!\brief Initialize the global interrupts to a specific mode
 *\see GLOBINT_MODE_t
 *\return Success/Failure
 *\retval 0 Success
 *\retval -1 Invalid mode
*/
__INLINE__ int GlobInt_Initialize(uint32_t channel, GLOBINT_MODE_t mode)
{
	switch(mode)
    {
		case GLOBINT_MODE_DIRECT: // direct
			WriteGlobIntDCR(_BGP_DCR_GLOBINT_ASSERT_CH(channel), 0x0);
			break;
		case GLOBINT_MODE_STICKYOR:
			WriteGlobIntDCR(_BGP_DCR_GLOBINT_ASSERT_CH(channel), 0x0);
			WriteGlobIntDCR(_BGP_DCR_GLOBINT_SET_CH(channel), _BGP_DCR_GLOBINT_SET_ARM_OR);
			_bgp_msync();
			_bgp_Delay(50);
			break;
		case GLOBINT_MODE_STICKYAND:
			WriteGlobIntDCR(_BGP_DCR_GLOBINT_ASSERT_CH(channel), 0x1);
			WriteGlobIntDCR(_BGP_DCR_GLOBINT_SET_CH(channel), _BGP_DCR_GLOBINT_SET_ARM_AND);
			_bgp_msync();
			_bgp_Delay(50);
			break;
		default:
			return -1;
    }
	return 0;
}

/** \brief Set state of a global interrupt channel
* \note This also unsets "sticky" mode for the channel.
*
* \param channel       Channel to use (e.g. 0..3)
* \param state         State to set (0 or 1)
*/
__INLINE__ void GlobInt_SetState(uint32_t channel, uint32_t state)
{
	WriteGlobIntDCR(_BGP_DCR_GLOBINT_ASSERT_CH(channel), state);
}

/** \brief Get status of global interrupts
*
* \return Returns the value of the global interrupts
* \see GlobInt_Status_t
*
*/
__INLINE__ GlobInt_Status_t GlobInt_GetStatus(void)
{
#if 1
    union
    {
	uint32_t uint32;
	GlobInt_Status_t globint;
    } status;
    status.uint32 = ReadGlobIntDCR(_BGP_DCR_GLOBINT_STATUS);
    return status.globint;
#else
	return (GlobInt_Status_t)ReadGlobIntDCR(_BGP_DCR_GLOBINT_STATUS);
#endif
}

/** \brief Set sticky arm type for a given channel.
*
* \param channel       Channel to use (e.g. 0..3)
* \param armtype       Type to set (0=OR, 1=AND)
*/
__INLINE__ void GlobInt_SetARMType(uint32_t channel, uint32_t armtype)
{
	WriteGlobIntDCR(_BGP_DCR_GLOBINT_SET_CH(channel), armtype );
}

/** \brief Return true (non-zero) if channel's receive state is raised
*
* \param channel       Channel to test (e.g. 0..3)
* \return Indication of whether the global interrupt is raised or not on the specified channel
* \retval 0 Not-Raised
* \retval 1 Raised
*/
__INLINE__ uint32_t GlobInt_IsRaised(uint32_t channel)
{
    // NOTE: Status bits 20-23 contain the downtree status for channels 0-3 respectively.
    return (ReadGlobIntDCR(_BGP_DCR_GLOBINT_STATUS) & _BN(20+channel)) ? 1 : 0;
}

/** \brief Disarm a sticky global interrupt
*
* \param channel       Channel to disarm (e.g. 0..3)
*/
__INLINE__ void GlobInt_Disarm(int channel)
{
	WriteGlobIntDCR(_BGP_DCR_GLOBINT_DRIVE_CH(channel), 0);
}

/** \brief Fire a sticky global interrupt
*
* \param channel       Channel to fire (e.g. 0..3)
*/
__INLINE__ void GlobInt_Fire(int channel)
{
	WriteGlobIntDCR(_BGP_DCR_GLOBINT_DRIVE_CH(channel), 1);
}

/** \brief Return true (non-zero) if channel's sticky state is raised
*
* \param channel       Channel to test (e.g. 0..3)
* \return Indication of whether sticky bit was raised
* \retval 0 Not-raised
* \retval 1 Raised
*/
__INLINE__ uint32_t GlobInt_IsStickyRaised(uint32_t channel)
{
    // NOTE: Status bits 4-7 contain the stick bits for channels 0-3, respectively.
    return (ReadGlobIntDCR(_BGP_DCR_GLOBINT_STATUS) & _BN(4+channel)) ? 1 : 0;
}



__INLINE__ uint32_t GlobInt_InitBarrier(uint32_t channel)
{
  /* Wait for idle state. */
  while(! GlobInt_IsRaised(channel))
    {
    }
  /* Fire global interrupt to show we are entering the barrier. */
  GlobInt_Disarm(channel);
  GlobInt_SetARMType(channel, _BGP_DCR_GLOBINT_SET_ARM_AND);
  GlobInt_Fire(channel);
  return 0;
}

__INLINE__ uint32_t GlobInt_QueryDone(uint32_t channel)

{
  if (GlobInt_IsStickyRaised(channel))
    return 1;
  else
    return 0;
}

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
__INLINE__ uint32_t GlobInt_Barrier(uint32_t channel, GlobInt_WaitProc_t waitproc, void* waitarg)
{
  GlobInt_InitBarrier(channel);

  for (;;)
    {
      if(GlobInt_QueryDone(channel))
	break;

      if (waitproc) {
	/* The waitproc may do something that allows other nodes to enter
	 * the barrier (e.g. MPI advance)
	 */
	int ret = waitproc(waitarg);
	if (ret)
	  return ret;
      }
    }
  return 0;
}

__END_DECLS
