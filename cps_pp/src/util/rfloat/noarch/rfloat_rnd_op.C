#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Definition of rfloat operator methods.

  $Id: rfloat_rnd_op.C,v 1.2 2003-07-24 16:53:54 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2003-07-24 16:53:54 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/rfloat/noarch/rfloat_rnd_op.C,v 1.2 2003-07-24 16:53:54 zs Exp $
//  $Id: rfloat_rnd_op.C,v 1.2 2003-07-24 16:53:54 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.3  2001/08/16 10:50:39  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.2  2001/06/19 18:13:38  anj
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
//  Revision 1.2  2001/05/25 06:16:11  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: rfloat_rnd_op.C,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/rfloat/noarch/rfloat_rnd_op.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// rfloat_rnd_op.C
//
// These overloaded operators perform the actual rounding.
// This file should not be compiled for use on the DSP.
// This file is provided for consistency when compiling
// for a machine that does rounding. If one needs to compile
// for the DSP the rfload_rnd_op.asm file should be used
// instead. The code in the .asm file does the actual rounding.
//
// More overloaded operators have been directly defined
// in rfloat.h but they all use the operators below.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/rfloat.h>
CPS_START_NAMESPACE


rfloat& rfloat::operator+=(IFloat a)
{x += a;  return *this;}

rfloat& rfloat::operator-=(IFloat a)
{x -= a;  return *this;}

rfloat& rfloat::operator*=(IFloat a)
{x *= a;  return *this;}

rfloat& rfloat::operator/=(IFloat a)
{x /= a;  return *this;}

rfloat operator-(const rfloat& a)
{return -a.x;}


CPS_END_NAMESPACE
