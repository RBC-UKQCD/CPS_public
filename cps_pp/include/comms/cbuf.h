#include<config.h>

/*!\file
  \brief  Functions used by the data layout conversion functions.

  $Id: common.C,v 1.5 2004/09/02 16:59:26 zs Exp $
*/
#include <util/data_types.h>

CPS_START_NAMESPACE

#if 0
#ifdef __cplusplus
extern "C" {
#endif

void negate_link(unsigned link_size, IFloat *link) ;
void site2cram(IFloat *src, IFloat *dst, unsigned site_size) ;
void site2dram(IFloat *src, IFloat *dst, unsigned *link_tbl, unsigned site_size) ;

#ifdef __cplusplus
}
#endif
#endif

//! Negate a floating point array.
/*!
  \param link_size The length of the array.
  \param link The array
  \post All elements of \a link are negated.
*/
  
inline void negate_link(unsigned link_size, IFloat *link)
{
	unsigned idx ;

	for (idx=0; idx<link_size; idx++) *(link+idx) = - *(link+idx) ;
}

//! Copy an array.
/*!
  \param src The array to copy from.
  \param dst The array to copy to.
  \param site_size The length of the array.
*/

inline void site2cram(IFloat *src, IFloat *dst, unsigned site_size)
{
	unsigned offset ;

	for (offset=0; offset<site_size; offset++)
		*(dst+offset) = *(src+offset) ;
}

//! Copy and rearrange an array.
/*!
  The copy is of the form <em>\n
  dst[link_tbl[i]] = src[i]
  \n</em>
  \param src The array to copy from.
  \param dst The array to copy to.
  \param link_tbl The look-up table for the \a dst indices.
  \param site_size The length of the \a src array.
*/

inline void site2dram(IFloat *src, IFloat *dst, unsigned *link_tbl, unsigned site_size)
{
	unsigned offset ;

	for (offset=0; offset<site_size; offset++)
		*(dst+*(link_tbl+offset)) = *(src+offset) ;
}

#if 0
//------------------------------------------------------------------
//
// cbuf.h
//
// Various routines that relate to the circular buffer.
//
//------------------------------------------------------------------

#ifndef INCLUDED_CBUF_H
#define INCLUDED_CBUF_H


//------------------------------------------------------------------
// saveCbufCntrlReg:
// Saves the contents of the circular buffer control registers
// in the cbuf_CntrlReg_reg static array
//------------------------------------------------------------------
void saveCbufCntrlReg(void);


//------------------------------------------------------------------
// restoreCbufCntrlReg:
// Restores the contents of the circular buffer control registers
// from the cbuf_CntrlReg_reg static array
//------------------------------------------------------------------
void restoreCbufCntrlReg(void);


//------------------------------------------------------------------
// setCbufCntrlReg(reg_no, value):
// Sets the contents of the circular buffer control register 
// reg_no to value.
//------------------------------------------------------------------
void setCbufCntrlReg(int reg_no, unsigned int value);


//------------------------------------------------------------------
// printCbufCntrlReg:
// Prints the contents of the circular buffer control registers
// under the control of VRB.Flow
//------------------------------------------------------------------
void printCbufCntrlReg(void);

#endif

#endif



CPS_END_NAMESPACE
