#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-01-13 20:38:59 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_ghb/qcdsp/c_func.C,v 1.3 2004-01-13 20:38:59 chulwoo Exp $
//  $Id: c_func.C,v 1.3 2004-01-13 20:38:59 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.2.10.1  2003/11/06 00:10:38  cwj
//  *** empty log message ***
//
//  Revision 1.1.1.1  2003/11/04 05:04:56  chulwoo
//
//  starting again
//
//
//  Revision 1.2  2003/07/24 16:53:53  zs
//  Addition of documentation via doxygen:
//  doxygen-parsable comment blocks added to many source files;
//  New target in makefile and consequent alterations to configure.in;
//  New directories and files under the doc directory.
//
//  Revision 1.6  2001/08/16 14:29:16  anj
//  Minor changes to allow compilation on all platforms, plus a better
//  regression testing script and more regressions output (for the serial
//  version of 4_1_0. Anj
//
//  Revision 1.5  2001/08/16 10:49:38  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.2  2001/06/19 18:11:23  anj
//  Serious ANSIfication.  Plus, degenerate double64.h files removed.
//  Next version will contain the new nga/include/double64.h.  Also,
//  Makefile.gnutests has been modified to work properly, propagating the
//  choice of C++ compiler and flags all the way down the directory tree.
//  The mpi_scu code has been added under phys/nga, and partially
//  plumbed in.
//
//  Everything has newer dates, due to the way in which this first alteration was handled.
//
//  Anj.
//
//  Revision 1.2  2001/05/25 06:15:59  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: c_func.C,v $
//  $Revision: 1.3 $
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
