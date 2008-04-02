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
#ifndef _BGP_CLOCKSTOP_SPI_H_ // Prevent multiple inclusion  
#define _BGP_CLOCKSTOP_SPI_H_

#include <common/namespace.h>

__BEGIN_DECLS



#include <bpcore/bgp_dcrmap.h>
#include <common/bgp_bitnumbers.h>


#define _BGP_DCR_CLOCKSTOP_STATUS0          (_BGP_DCR_CLOCKSTOP + 0x04)
         // 32 LSb of 43-bit current counter


#define _BGP_DCR_CLOCKSTOP_STATUS1          (_BGP_DCR_CLOCKSTOP + 0x05)
#define   _CLOCKSTOP_STATUS1_OUTPUT           _BN(0)     // output value of "precise clock stop"
#define   _CLOCKSTOP_STATUS1_CURRENT(x)       _B11(31,x) // 11 MSb of 43-bit current counter val

__END_DECLS



#endif // Add nothing below this line

