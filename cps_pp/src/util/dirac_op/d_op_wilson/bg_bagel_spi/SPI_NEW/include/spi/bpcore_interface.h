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
 * \file spi/bpcore_interface.h
 */
#ifndef _BGP_BPCORE_INT_H_ // Prevent multiple inclusion
#define _BGP_BPCORE_INT_H_

#define _BGP_UA_SCRATCH      (0x4)          // eDRAM Scratch: 0 to 8MB
#define _BGP_PA_SCRATCH      (0x00000000)
#define _BGP_PS_SCRATCH      (8 * 1024 * 1024) 
#define _BGP_PM_SCRATCH      (0x007FFFFF)

/* ************************************************************************* */
/* DMA Non-Fatal Interrupt Request: Group 3 bits 00:31                       */
/* ************************************************************************* */

#define _BGP_IC_DMA_NFT_G3_HIER_POS   3
#define _BGP_IC_DMA_NFT_G3_UNIT_NUM   3
#define _BGP_IC_DMA_NFT_G3_UNIT_POS   0
#define _BGP_IC_DMA_NFT_G3_UNIT_SIZE  32
#define _BGP_IC_DMA_NFT_G3_UNIT_MASK  0xffffffff


#endif // Add nothing below this line

