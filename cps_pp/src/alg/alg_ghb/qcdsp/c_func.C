#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-06-02 09:36:38 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_ghb/qcdsp/c_func.C,v 1.4 2004-06-02 09:36:38 zs Exp $
//  $Id: c_func.C,v 1.4 2004-06-02 09:36:38 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: c_func.C,v $
//  $Revision: 1.4 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_ghb/qcdsp/c_func.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------
//
//	These functions are C++ support routines for the
//  assembly code.  bad_sqrt() is an error call out of the
//  assembly optimized square root function, the others
//  are simple printing routines designed to be called
//  from assembly code in debugging situations.  In the
//  production code all comments to those functions are
//  commented out, and only bad_sqrt() is required.
//------------------------------------------------------------
#include	<stdlib.h>
#include	<stdio.h>
CPS_END_NAMESPACE
#include <util/data_types.h>
#include <util/verbose.h>
#include <util/error.h>
CPS_START_NAMESPACE

#if( _TARTAN )
extern "C" {
  void  bad_sqrt( void )
  {
    printf("# SQUARE ROOT FAILURE, PASSED NEGATIVE VALUE!\n");
    exit( -1 );
  }
 
  void  nullprint( void ) { printf("\n"); }
  void  dotprint()   { printf("."); }
  void  iprint( int x )   { printf("# %d\n", x); }
  void  fprint( IFloat x ) { printf("# %+f\n", x); }
  void  fprint2( IFloat x ) { printf("> %f\n", x); }
  void  xprint( unsigned x ) { printf("# %08x\n", x); }
}
#endif


CPS_END_NAMESPACE
