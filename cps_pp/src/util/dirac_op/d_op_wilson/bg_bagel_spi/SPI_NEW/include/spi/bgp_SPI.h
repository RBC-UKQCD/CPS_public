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
#ifndef _BGP_SPI_H_ // Prevent multiple inclusion
#define _BGP_SPI_H_



/*!
 * \file spi/bgp_SPI.h
 *
 * \brief BG/P Master Include File for all System Programming Intefaces (SPIs)
 */

#include <common/namespace.h>

__BEGIN_DECLS

#include <common/linkage.h>
#include <bpcore/bgp_types.h>
#include <bpcore/ppc450_inlines.h>
#include <common/bgp_personality.h>

#ifndef __INLINE__
#define __INLINE__ extern inline
#endif
#define SPI_DEPRECATED 1

#include <spi/kernel_interface.h>        // functions for Virtual to/from Physical Translation
#include <spi/lockbox_interface.h>       // Lockbox Functions, eg BGP_LockBox_Barrier_Group()
#include <spi/DMA.h>                     // DMA definitions
#include <spi/GlobInt.h>
#include <spi/UPC.h>                     // UPC functions
#include <spi/UPC_Events.h>              // UPC events

__END_DECLS



#endif // Add nothing below this line
