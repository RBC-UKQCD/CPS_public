#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
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
#include <util/qcdio.h>
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
