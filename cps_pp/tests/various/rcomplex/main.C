#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:47 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/various/rcomplex/main.C,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
//  $Id: main.C,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.9  2002/03/11 22:26:58  anj
//  This should now be the correct, fully merged code from our two versions. Anj
//
//  Revision 1.6.2.1  2002/03/08 16:36:29  anj
//  Checking in the Columbia code branch on tag Columbia4_1_1_test-branch, to be
//  merged with the UKQCD head branch shortly.  Anj
//
//  Revision 1.6  2001/08/17 20:03:39  anj
//  Multiple (extra) changes to make the test suite smaller (16CPUs
//  required, not 64) and faster.  Anj
//
//  Revision 1.5  2001/08/16 10:50:07  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "float".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and float).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.4  2001/07/03 17:00:59  anj
//
//  Multiple minor alterations to change some #include's from referring to
//  files relative to the top-level source directory to referring to files
//  relative to the source-file positions.  This alteration makes the code
//  backwards compatable with the make structure of QCDSP, although this
//  may have to be changed to a more usual form in the future. Anj.
//
//  Revision 1.3  2001/06/21 15:40:15  anj
//  Updated the _TARTAN ifdefs, using PARALLEL instead (where appropriate).Anj
//
//  Revision 1.2  2001/06/19 18:12:32  anj
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
//  Revision 1.2  2001/05/25 06:16:04  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: main.C,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/various/rcomplex/main.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//======================================
// SUI	Oct. 31 1997
//======================================
// test driver for Rcomplex class

CPS_END_NAMESPACE
#include<util/rcomplex.h>
#include <stdio.h>
#include <stdlib.h>
#include<config.h>
CPS_START_NAMESPACE

#ifdef _TARTAN
  #include <matht.h>
#else
  #include <math.h>
#endif

#define BUF_SIZE 1000

IFloat buffer[BUF_SIZE];
IFloat *buffer_p = buffer;

main()
{
  Rcomplex c0(1.98765, 1.23456);
  Rcomplex c1(0.6, 0.610001);
  Rcomplex c2(0.123456, 0.187654);

  IFloat f = 1.123;

  int i = 0;	// buffer write index
  int n = 0;	// buffer read index

  printf(" test + and hence += :\n\n");
  {
    Rcomplex r1;
    Rcomplex r2;
    Rcomplex r3;

    for (int l = 0; l < 100; l++) {
      r1 = r1 + c0;
      r2 = f + r2;
      r3 = r3 + f;
    }
    buffer[i++] = r1.real();
    buffer[i++] = r2.real();
    buffer[i++] = r3.real();
    buffer[i++] = r1.imag();
    buffer[i++] = r2.imag();
    buffer[i++] = r3.imag();
  }
  for (n = 0; n < i ; n++) printf("%e\n", buffer[n]);

  printf(" test operator - and hence -=:\n\n");
  {
    Rcomplex r1;
    Rcomplex r2;
    Rcomplex r3;

    for (int l = 0; l < 100; l++) {
      r1 = r1 - c0;
      r2 = r2 - f;
      r3 = f - r3;
    }
    buffer[i++] = r1.real();
    buffer[i++] = r2.real();
    buffer[i++] = r3.real();
    buffer[i++] = r1.imag();
    buffer[i++] = r2.imag();
    buffer[i++] = r3.imag();
  }
  for (; n < i ; n++) printf("%e\n", buffer[n]);
    
  printf(" test operator *= and * :\n\n");
  {
    Rcomplex r1 = c1;
    Rcomplex r2 = c1;
    Rcomplex r3 = c1;


    for (int l = 0; l < 100; l++) {
      r1 = r1 * c1;
      r2 = r2 * f;
      r3 = f * r3;
    }
    buffer[i++] = r1.real();
    buffer[i++] = r2.real();
    buffer[i++] = r3.real();
    buffer[i++] = r1.imag();
    buffer[i++] = r2.imag();
    buffer[i++] = r3.imag();
  }
  for ( ; n < i ; n++) printf("%e\n", buffer[n]);


  printf(" test operator /= and / :\n\n");
  {
    Rcomplex r1 = c1;
    Rcomplex r2 = c1;
    Rcomplex r3 = c1;
    
    for (int l = 0; l < 100; l++) {
      r1 = r1 / c1;
      r2 = r2 / f;
      r3 = f / r3;
    }
    buffer[i++] = r1.real();
    buffer[i++] = r2.real();
    buffer[i++] = r3.real();
    buffer[i++] = r1.imag();
    buffer[i++] = r2.imag();
    buffer[i++] = r3.imag();
  }
  for ( ; n < i ; n++) printf("%e\n", buffer[n]);

  printf(" test unary operator - (negate), and conjugate:\n\n");
  {
    Rcomplex r1 = c1;
    Rcomplex r2 = -r1;
    Rcomplex r3 = conj(r1);

    buffer[i++] = real(r1);
    buffer[i++] = real(r2);
    buffer[i++] = real(r3);
    buffer[i++] = imag(r1);
    buffer[i++] = imag(r2);
    buffer[i++] = imag(r3);
  }
  for ( ; n < i ; n++) printf("%e\n", buffer[n]);

  printf(" test norm() and abs() in both forms:\n\n");
  {
    Rcomplex r1 = c1;
    Rcomplex r2 = c1;

    buffer[i++] = r1.norm();
    buffer[i++] = norm(r2);
    buffer[i++] = r1.abs();
    buffer[i++] = abs(r2);
  }
  for ( ; n < i ; n++) printf("%e\n", buffer[n]);

  printf("Done!");
}
CPS_END_NAMESPACE
