#include<config.h>
CPS_START_NAMESPACE
/*----------------------------------------------------------*/
/*! The QCD I/O Interface:

  $Id: qcdio.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $

  A.N.Jackson: ajackson@epcc.ed.ac.uk                      
  -----------------------------------------------------------
   CVS keywords
 
   $Author: mcneile $
   $Date: 2003-06-22 13:34:52 $
   $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/qcdio.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
   $Id: qcdio.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
   $Name: not supported by cvs2svn $
   $Locker:  $
   $RCSfile: qcdio.h,v $
   $Revision: 1.1.1.1 $
   $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/qcdio.h,v $
   $State: Exp $  */
/*----------------------------------------------------------*/

#ifndef INCLUDED_QCDIO_H_
#define INCLUDED_QCDIO_H_

CPS_END_NAMESPACE
#include <stdio.h>
#include<config.h>
#include<util/data_types.h>
#include<util/lattice.h>
#include<util/qcdio_qprintf.h>
CPS_START_NAMESPACE

#ifndef GAUGE_CONF_PREC
/*! When loading a guage configurations, the default precision is: */
#define GAUGE_CONF_PREC (sizeof(float))
#endif

#ifndef SWAP_BYTE_ORDER
/*! When loading gauge configurations, should the byte der be swapped before
 * putting the results into memory?  1 = yes, 0 = no. */
#define SWAP_BYTE_ORDER 0
#endif

#ifndef TRANSPOSE_THE_MATRICIES
/*! When the SU(3) matricies have been loaded, the row v. col convention must
 * be the same as for the code which generated the configuration.  To
 * transpose the matrix, set this to be 1, else use 0. */
#define TRANSPOSE_THE_MATRICIES 1
#endif

  /* ------------------------------------------------------
  	Gauge configuration I/O:
    ------------------------------------------------------ */
  
  //! Routine to define whether the loaded data should be normalized:
  void qcdio_set_normalize( int );

  //! Routine for loading a gauge configuration
  void qload_gauge( char* fprefix, Lattice& lat,
             int prec = GAUGE_CONF_PREC,
	     int swap = SWAP_BYTE_ORDER,
	     int transp = TRANSPOSE_THE_MATRICIES );

  //! Routine for saving the current gauge configuration
  void qsave_gauge( char* fprefix, Lattice& lat,
             int prec = GAUGE_CONF_PREC,
             int swap = SWAP_BYTE_ORDER,
	     int transp = TRANSPOSE_THE_MATRICIES );

#endif

CPS_END_NAMESPACE
