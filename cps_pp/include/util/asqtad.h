#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-01-13 20:38:56 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/asqtad.h,v 1.2 2004-01-13 20:38:56 chulwoo Exp $
//  $Id: asqtad.h,v 1.2 2004-01-13 20:38:56 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.1.2.1  2003/11/05 16:13:21  mike
//  Initial attempt at producing working branch
//
//  Revision 1.1.1.1  2003/09/18 22:30:58  chulwoo
//  Mike's files for single node QCDOC + Parallel transport
//  I added some hacks for PARALLEL without MPI_SCU
//  PARALLEL=2 set PARALLEL without MPI_SCU
//
//
//  Revision 1.2  2003/07/24 16:53:53  zs
//  Addition of documentation via doxygen:
//  doxygen-parsable comment blocks added to many source files;
//  New target in makefile and consequent alterations to configure.in;
//  New directories and files under the doc directory.
//
//  Revision 1.4  2001/08/16 10:50:31  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.2  2001/06/19 18:13:19  anj
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
//  Revision 1.2  2001/05/25 06:16:09  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: asqtad.h,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/asqtad.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
/****************************************************************************/
/*                                                                          */
/* stag.h                                                                   */
/*                                                                          */
/* C header file for the staggered library.                                 */
/*                                                                          */
/****************************************************************************/

#ifndef INCLUDED_STAG_LIB_H
#define INCLUDED_STAG_LIB_H

//------------------------------------------------------------------
// The staggered dslash operator
// dir_flag is flag which takes value 0 when all direction contribute to D,
// 1 - when only the special anisotropic direction contributes to D,
// 2 - when all  except the special anisotropic direction. 
//------------------------------------------------------------------
extern "C" void asqtad_dirac(IFloat *f_out, IFloat *f_in, int cb, int dag , 
		      int dir_flag = 0);

//------------------------------------------------------------------
// Initialize all global variables and address tables needed by
// staggered dirac() function.
//------------------------------------------------------------------
extern "C" void asqtad_dirac_init(const void *gauge_field_addr);
extern "C" void asqtad_dirac_init_g ();

//------------------------------------------------------------------
// Destroy all address tables needed by staggered dirac() function
//------------------------------------------------------------------
extern "C" void asqtad_destroy_dirac_buf();
extern "C" void asqtad_destroy_dirac_buf_g();

#endif

CPS_END_NAMESPACE
