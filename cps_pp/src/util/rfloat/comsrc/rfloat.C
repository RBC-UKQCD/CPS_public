#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Definition of rfloat methods.

  $Id: rfloat.C,v 1.2 2003-07-24 16:53:54 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2003-07-24 16:53:54 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/rfloat/comsrc/rfloat.C,v 1.2 2003-07-24 16:53:54 zs Exp $
//  $Id: rfloat.C,v 1.2 2003-07-24 16:53:54 zs Exp $
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
//  Revision 1.2  2001/06/19 18:13:37  anj
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
//  $RCSfile: rfloat.C,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/rfloat/comsrc/rfloat.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//  rfloat.C


CPS_END_NAMESPACE
#include <util/rfloat.h>
CPS_START_NAMESPACE


rfloat::rfloat(IFloat a) : x(a) {}

rfloat::rfloat(const rfloat& a) : x(a.x) {}

rfloat::~rfloat() {}


//-----------------------------------------------------------------
//  free friend functions
//-----------------------------------------------------------------

rfloat operator+(const rfloat& a, const rfloat& b)
	{ return rfloat(a.x) += b.x; }
rfloat operator+(double a, const rfloat& b)
	{ return rfloat(a) += b.x; }
rfloat operator+(const rfloat& a, double b)
	{ return rfloat(a.x) += b; }

rfloat operator-(const rfloat& a, const rfloat& b)
	{ return rfloat(a.x) -= b.x; }
rfloat operator-(double a, const rfloat& b)
	{ return rfloat(a) -= b.x; }
rfloat operator-(const rfloat& a, double b)
	{ return rfloat(a.x) -= b; }

rfloat operator*(const rfloat& a, const rfloat& b)
	{ return rfloat(a.x) *= b.x; }
rfloat operator*(double a, const rfloat& b)
	{ return rfloat(a) *= b.x; }
rfloat operator*(const rfloat& a, double b)
	{ return rfloat(a.x) *= b; }

rfloat operator/(const rfloat& a, const rfloat& b)
	{ return rfloat(a.x) /= b.x; }
rfloat operator/(double a, const rfloat& b)
	{ return rfloat(a) /= b.x; }
rfloat operator/(const rfloat& a, double b)
	{ return rfloat(a.x) /= b; }


CPS_END_NAMESPACE
