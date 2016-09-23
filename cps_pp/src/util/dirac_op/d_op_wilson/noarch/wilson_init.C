#include <config.h>
#ifdef USE_SSE
#include "../sse/wilson_init.C"
#else
CPS_START_NAMESPACE
/*! \file
  \brief  Routine used internally in the DiracOpWilson class.
  
  $Id: wilson_init.C,v 1.6 2012/03/26 13:50:12 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2012/03/26 13:50:12 $
//  $Header: /space/cvs/cps/cps++/src/util/dirac_op/d_op_wilson/noarch/wilson_init.C,v 1.6 2012/03/26 13:50:12 chulwoo Exp $
//  $Id: wilson_init.C,v 1.6 2012/03/26 13:50:12 chulwoo Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $Revision: 1.6 $
//  $Source: /space/cvs/cps/cps++/src/util/dirac_op/d_op_wilson/noarch/wilson_init.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
/****************************************************************************/
/* 10/16/97                                                                 */
/*                                                                          */
/* wilson_int:                                                              */
/*                                                                          */
/* This routine performs all initializations needed before wilson funcs     */
/* are called. It sets the addressing related arrays and reserves memory    */
/* for the needed temporary buffers. It only needs to be called             */
/* once at the begining of the program (or after a wilson_end call)         */
/* before any number of calls to wilson funcs are made.                     */
/*                                                                          */
/* WARNING:                                                                 */
/*                                                                          */
/* This set of routines will work only if the node sublattices have         */
/* even number of sites in each direction.                                  */
/*                                                                          */
/****************************************************************************/

/*--------------------------------------------------------------------------*/
/* Include header files                                                     */
/*--------------------------------------------------------------------------*/
CPS_END_NAMESPACE
#include <util/data_types.h>
#include <util/wilson.h>
#include <util/gjp.h>
#include <util/smalloc.h>
#include <util/verbose.h>
#include <util/error.h>

#if TARGET == BGP
void spi_init(){}
#endif



CPS_START_NAMESPACE
#if TARGET == BGP
void wilson_set_sloppy(bool sloppy) {}
#endif

void wilson_init(Wilson *wilson_p)  /* pointer to Wilson type structure    */
{
  char *cname = " ";
  char *fname = "wilson_init(Wilson*)";
  VRB.Func(cname,fname);

  int spinor_words;             /* size of the spinor field on the         */
				/* sublattice checkerboard                 */
  int size;

/*--------------------------------------------------------------------------*/
/* Reserve memory for the node sublattice sizes                             */
/*--------------------------------------------------------------------------*/
  size = 4*sizeof(int);
  wilson_p->ptr = (int *) smalloc(size);
  if( wilson_p->ptr == 0)
    ERR.Pointer(cname,fname, "ptr");
  VRB.Smalloc(cname,fname,
	      "ptr", wilson_p->ptr, size);

/*--------------------------------------------------------------------------*/
/* Set the node sublattice sizes                                            */
/*--------------------------------------------------------------------------*/
  wilson_p->ptr[0] = GJP.XnodeSites();
  wilson_p->ptr[1] = GJP.YnodeSites();
  wilson_p->ptr[2] = GJP.ZnodeSites();
  wilson_p->ptr[3] = GJP.TnodeSites();
  wilson_p->vol[0] = wilson_p->ptr[0] * wilson_p->ptr[1] *
                     wilson_p->ptr[2] * wilson_p->ptr[3] / 2;
  wilson_p->vol[1] = wilson_p->vol[0];

/*--------------------------------------------------------------------------*/
/* Reserve memory for 2  temporary spinors                                  */
/* Use the af[] array to pass them (for this routines af are not            */
/* spin projected half spinors but instead they are full 4 component        */
/* spinors.)                                                                */
/*--------------------------------------------------------------------------*/
  spinor_words = SPINOR_SIZE * wilson_p->vol[0];

  wilson_p->af[0] = (IFloat *) smalloc(spinor_words*sizeof(Float));
  if(wilson_p->af[0] == 0)
    ERR.Pointer(cname,fname, "af[0]");
  VRB.Smalloc(cname,fname,
	      "af[0]", wilson_p->af[0], spinor_words*sizeof(Float));

  wilson_p->af[1] = (IFloat *) smalloc(spinor_words*sizeof(Float));
  if(wilson_p->af[1] == 0)
    ERR.Pointer(cname,fname, "af[1]");
  VRB.Smalloc(cname,fname,
	      "af[1]", wilson_p->af[1], spinor_words*sizeof(Float));

  VRB.Debug("x = %d\n", wilson_p->ptr[0]);
  VRB.Debug("y = %d\n", wilson_p->ptr[1]);
  VRB.Debug("z = %d\n", wilson_p->ptr[2]);
  VRB.Debug("t = %d\n", wilson_p->ptr[3]);
  VRB.Debug("vol0 = %d\n", wilson_p->vol[0]);
  VRB.Debug("vol1 = %d\n", wilson_p->vol[1]);
  VRB.Debug("af0 = %x\n", wilson_p->af[0]);
  VRB.Debug("af1 = %x\n", wilson_p->af[1]);
}


CPS_END_NAMESPACE
#endif
