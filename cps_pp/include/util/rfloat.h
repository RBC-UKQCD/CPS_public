#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:52 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/rfloat.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
//  $Id: rfloat.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.5  2001/08/16 10:50:31  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.4  2001/07/19 13:02:24  anj
//  Removed the overloaded = operator, which appears to cause some problems on QCDSP.
//
//  Revision 1.3  2001/06/28 14:34:22  anj
//
//  The core ANSIfication should now be complete.  There are a few
//  remaining issues, but this version should compile anywhere and be
//  backward compatable with QCDSP (although this requires the top source
//  directory (.../phys/ to be added to the include path).
//
//  The serial GCC version has also been tested, and all test programs
//  appear to behave as they should (not to imply that they all work, but
//  I believe those that should work are ok).  There are minor differences
//  in the results due to rounding, (see example pbp_gccsun.dat files),
//  but that is all.
//
//  Anj.
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
//  $RCSfile: rfloat.h,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/rfloat.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
// rfloat.h
//

#ifndef INCLUDED_RFLOAT_H
#define INCLUDED_RFLOAT_H

class rfloat;
CPS_END_NAMESPACE
#include<util/data_types.h>
CPS_START_NAMESPACE

class rfloat {
    IFloat x;
public:
    rfloat(IFloat a = 0);
    rfloat(const rfloat& a);
    ~rfloat();

    rfloat& operator=(IFloat a)
    	{ x = a; return *this; }

    operator IFloat() const { return x; }
    // conversion. Even though the user-defined conversion is a bad
    	// practice in many situation, it's OK here, because the
	// behavior of rfloat is exactly the same as IFloat except
	// the overloaded operators.

    // Overloaded = operator to allow printf to print rfloats as IFloats:
    //IFloat operator=( const rfloat& a ) { return a.x; }

    //---------------------------------------------------------
    //  overloaded operators
    //---------------------------------------------------------
    	// binary plus
    friend rfloat operator+(const rfloat&, const rfloat&);
    friend rfloat operator+(double, const rfloat&);
    friend rfloat operator+(const rfloat&, double);

    	// binary minus
    friend rfloat operator-(const rfloat&, const rfloat&);
    friend rfloat operator-(double, const rfloat&);
    friend rfloat operator-(const rfloat&, double);

    	// binary multiply
    friend rfloat operator*(const rfloat&, const rfloat&);
    friend rfloat operator*(double, const rfloat&);
    friend rfloat operator*(const rfloat&, double);

    	// binary division
    friend rfloat operator/(const rfloat&, const rfloat&);
    friend rfloat operator/(double, const rfloat&);
    friend rfloat operator/(const rfloat&, double);


    	// prefix unary minus
    friend rfloat operator-(const rfloat&);


    rfloat& operator+=(IFloat a);
    rfloat& operator+=(const rfloat& a)
    	{ *this += a.x;  return *this; }

    rfloat& operator-=(IFloat a);
    rfloat& operator-=(const rfloat& a)
    	{ *this -= a.x;  return *this; }

    rfloat& operator*=(IFloat a);
    rfloat& operator*=(const rfloat& a)
    	{ *this *= a.x;  return *this; }

    rfloat& operator/=(IFloat a);
    rfloat& operator/=(const rfloat& a)
    	{ *this /= a.x;  return *this; }
};


//-----------------------------------------------------------------
// Free freind functions
//-----------------------------------------------------------------
rfloat operator-(const rfloat& a);

rfloat operator+(const rfloat& a, const rfloat& b);
rfloat operator+(double a, const rfloat& b);
rfloat operator+(const rfloat& a, double b);

rfloat operator-(const rfloat& a, const rfloat& b);
rfloat operator-(double a, const rfloat& b);
rfloat operator-(const rfloat& a, double b);

rfloat operator*(const rfloat& a, const rfloat& b);
rfloat operator*(double a, const rfloat& b);
rfloat operator*(const rfloat& a, double b);

rfloat operator/(const rfloat& a, const rfloat& b);
rfloat operator/(double a, const rfloat& b);
rfloat operator/(const rfloat& a, double b);





#endif
CPS_END_NAMESPACE
