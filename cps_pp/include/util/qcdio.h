#include<config.h>
CPS_START_NAMESPACE
/*----------------------------------------------------------*/
/*!\file
  \brief Prototypes of gauge configuration IO functions.

  $Id: qcdio.h,v 1.2 2003-07-24 16:53:53 zs Exp $
*/
/*2  A.N.Jackson: ajackson@epcc.ed.ac.uk                      
  -----------------------------------------------------------
   CVS keywords
 
   $Author: zs $
   $Date: 2003-07-24 16:53:53 $
   $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/qcdio.h,v 1.2 2003-07-24 16:53:53 zs Exp $
   $Id: qcdio.h,v 1.2 2003-07-24 16:53:53 zs Exp $
   $Name: not supported by cvs2svn $
   $Locker:  $
   $RCSfile: qcdio.h,v $
   $Revision: 1.2 $
   $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/qcdio.h,v $
   $State: Exp $  */
/*----------------------------------------------------------*/

#ifndef INCLUDED_QCDIO_H_
#define INCLUDED_QCDIO_H_

CPS_END_NAMESPACE
#include <stdio.h>
#include <util/data_types.h>
#include <util/lattice.h>
#include <util/qcdio_qprintf.h>
CPS_START_NAMESPACE

#ifndef GAUGE_CONF_PREC
/*! The default precision at which gauge configurations are stored in files. */
#define GAUGE_CONF_PREC (sizeof(float))
#endif

#ifndef SWAP_BYTE_ORDER
/*! When loading gauge configurations, should the byte order be swapped before
 * putting the results into memory?  To swap the bytes by default, set this
 to 1, or set it to 0 to leave the byte order alone by default. */
#define SWAP_BYTE_ORDER 0
#endif

#ifndef TRANSPOSE_THE_MATRICES
/*! When the SU(3) matrices have been loaded, the row vs. column convention
  must be the same as for the code which generated the configuration.  To
  transpose the matrix by default, set this to be 1, otherwise use 0. */
#define TRANSPOSE_THE_MATRICES 1
#endif

  /* ------------------------------------------------------
  	Gauge configuration I/O:
    ------------------------------------------------------ */
  
  //! Routine to define whether the loaded data should be normalized.
  void qcdio_set_normalize( int );

  //! Routine for loading a gauge configuration
  void qload_gauge( char* fprefix, Lattice& lat,
             int prec = GAUGE_CONF_PREC,
	     int swap = SWAP_BYTE_ORDER,
	     int transp = TRANSPOSE_THE_MATRICES );

  //! Routine for saving the current gauge configuration
  void qsave_gauge( char* fprefix, Lattice& lat,
             int prec = GAUGE_CONF_PREC,
             int swap = SWAP_BYTE_ORDER,
	     int transp = TRANSPOSE_THE_MATRICES );

#endif


CPS_END_NAMESPACE
