#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:46 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/vector/noarch/vector_util.C,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
//  $Id: vector_util.C,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.4  2001/08/16 10:50:40  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.2  2001/06/19 18:13:41  anj
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
//  $RCSfile: vector_util.C,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/vector/noarch/vector_util.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
/*------------------------------------------------------------------*/
/*
   vector_util.c

   This file contains the definitions of some genaral c-style
   functions that perform operations on vectors of general
   length. For these functions there exists optimized assembly 
   code.
*/
/*------------------------------------------------------------------*/

CPS_END_NAMESPACE
#include <string.h>		/* memcpy */
#include<util/vector.h>
CPS_START_NAMESPACE



//---------------------------------------------------------------//
void moveMem(void *b, const void *a, int len)
{
    memcpy(b, a, len);
}
//---------------------------------------------------------------//


//---------------------------------------------------------------//

//
//  Assume that c!=a and c!=b
//
//    c = a*b
//
void mDotMEqual(IFloat* c, const IFloat* a, const IFloat* b)
{
    *c      = *a      * *b      - *(a+1)  * *(b+1)    +
    	      *(a+2)  * *(b+6)  - *(a+3)  * *(b+7)    +
    	      *(a+4)  * *(b+12) - *(a+5)  * *(b+13);
    *(c+1)  = *a      * *(b+1)  + *(a+1)  * *b        +
    	      *(a+2)  * *(b+7)  + *(a+3)  * *(b+6)    +
    	      *(a+4)  * *(b+13) + *(a+5)  * *(b+12);

    *(c+2)  = *a      * *(b+2)  - *(a+1)  * *(b+3)    +
    	      *(a+2)  * *(b+8)  - *(a+3)  * *(b+9)    +
    	      *(a+4)  * *(b+14) - *(a+5)  * *(b+15);
    *(c+3)  = *a      * *(b+3)  + *(a+1)  * *(b+2)    +
    	      *(a+2)  * *(b+9)  + *(a+3)  * *(b+8)    +
    	      *(a+4)  * *(b+15) + *(a+5)  * *(b+14);

    *(c+4)  = *a      * *(b+4)  - *(a+1)  * *(b+5)    +
    	      *(a+2)  * *(b+10) - *(a+3)  * *(b+11)   +
    	      *(a+4)  * *(b+16) - *(a+5)  * *(b+17);
    *(c+5)  = *a      * *(b+5)  + *(a+1)  * *(b+4)    +
    	      *(a+2)  * *(b+11) + *(a+3)  * *(b+10)   +
    	      *(a+4)  * *(b+17) + *(a+5)  * *(b+16);

    *(c+6)  = *(a+6)  * *b      - *(a+7)  * *(b+1)    +
    	      *(a+8)  * *(b+6)  - *(a+9)  * *(b+7)    +
    	      *(a+10) * *(b+12) - *(a+11) * *(b+13);
    *(c+7)  = *(a+6)  * *(b+1)  + *(a+7)  * *b        +
    	      *(a+8)  * *(b+7)  + *(a+9)  * *(b+6)    +
    	      *(a+10) * *(b+13) + *(a+11) * *(b+12);

    *(c+8)  = *(a+6)  * *(b+2)  - *(a+7)  * *(b+3)    +
    	      *(a+8)  * *(b+8)  - *(a+9)  * *(b+9)    +
    	      *(a+10) * *(b+14) - *(a+11) * *(b+15);
    *(c+9)  = *(a+6)  * *(b+3)  + *(a+7)  * *(b+2)    +
    	      *(a+8)  * *(b+9)  + *(a+9)  * *(b+8)    +
    	      *(a+10) * *(b+15) + *(a+11) * *(b+14);

    *(c+10) = *(a+6)  * *(b+4)  - *(a+7)  * *(b+5)    +
    	      *(a+8)  * *(b+10) - *(a+9)  * *(b+11)   +
    	      *(a+10) * *(b+16) - *(a+11) * *(b+17);
    *(c+11) = *(a+6)  * *(b+5)  + *(a+7)  * *(b+4)    +
    	      *(a+8)  * *(b+11) + *(a+9)  * *(b+10)   +
    	      *(a+10) * *(b+17) + *(a+11) * *(b+16);

    *(c+12) = *(a+12) * *b      - *(a+13) * *(b+1)    +
    	      *(a+14) * *(b+6)  - *(a+15) * *(b+7)    +
    	      *(a+16) * *(b+12) - *(a+17) * *(b+13);
    *(c+13) = *(a+12) * *(b+1)  + *(a+13) * *b        +
    	      *(a+14) * *(b+7)  + *(a+15) * *(b+6)    +
    	      *(a+16) * *(b+13) + *(a+17) * *(b+12);

    *(c+14) = *(a+12) * *(b+2)  - *(a+13) * *(b+3)    +
    	      *(a+14) * *(b+8)  - *(a+15) * *(b+9)    +
    	      *(a+16) * *(b+14) - *(a+17) * *(b+15);
    *(c+15) = *(a+12) * *(b+3)  + *(a+13) * *(b+2)    +
    	      *(a+14) * *(b+9)  + *(a+15) * *(b+8)    +
    	      *(a+16) * *(b+15) + *(a+17) * *(b+14);

    *(c+16) = *(a+12) * *(b+4)  - *(a+13) * *(b+5)    +
    	      *(a+14) * *(b+10) - *(a+15) * *(b+11)   +
    	      *(a+16) * *(b+16) - *(a+17) * *(b+17);
    *(c+17) = *(a+12) * *(b+5)  + *(a+13) * *(b+4)    +
    	      *(a+14) * *(b+11) + *(a+15) * *(b+10)   +
    	      *(a+16) * *(b+17) + *(a+17) * *(b+16);
}


//
//  Assume that c!=a and c!=b
//
//    c += a*b
//
void mDotMPlus(IFloat* c, const IFloat* a, const IFloat* b)
{
    *c     += *a      * *b      - *(a+1)  * *(b+1)    +
    	      *(a+2)  * *(b+6)  - *(a+3)  * *(b+7)    +
    	      *(a+4)  * *(b+12) - *(a+5)  * *(b+13);
    *(c+1) += *a      * *(b+1)  + *(a+1)  * *b        +
    	      *(a+2)  * *(b+7)  + *(a+3)  * *(b+6)    +
    	      *(a+4)  * *(b+13) + *(a+5)  * *(b+12);

    *(c+2) += *a      * *(b+2)  - *(a+1)  * *(b+3)    +
    	      *(a+2)  * *(b+8)  - *(a+3)  * *(b+9)    +
    	      *(a+4)  * *(b+14) - *(a+5)  * *(b+15);
    *(c+3) += *a      * *(b+3)  + *(a+1)  * *(b+2)    +
    	      *(a+2)  * *(b+9)  + *(a+3)  * *(b+8)    +
    	      *(a+4)  * *(b+15) + *(a+5)  * *(b+14);

    *(c+4) += *a      * *(b+4)  - *(a+1)  * *(b+5)    +
    	      *(a+2)  * *(b+10) - *(a+3)  * *(b+11)   +
    	      *(a+4)  * *(b+16) - *(a+5)  * *(b+17);
    *(c+5) += *a      * *(b+5)  + *(a+1)  * *(b+4)    +
    	      *(a+2)  * *(b+11) + *(a+3)  * *(b+10)   +
    	      *(a+4)  * *(b+17) + *(a+5)  * *(b+16);

    *(c+6) += *(a+6)  * *b      - *(a+7)  * *(b+1)    +
    	      *(a+8)  * *(b+6)  - *(a+9)  * *(b+7)    +
    	      *(a+10) * *(b+12) - *(a+11) * *(b+13);
    *(c+7) += *(a+6)  * *(b+1)  + *(a+7)  * *b        +
    	      *(a+8)  * *(b+7)  + *(a+9)  * *(b+6)    +
    	      *(a+10) * *(b+13) + *(a+11) * *(b+12);

    *(c+8) += *(a+6)  * *(b+2)  - *(a+7)  * *(b+3)    +
    	      *(a+8)  * *(b+8)  - *(a+9)  * *(b+9)    +
    	      *(a+10) * *(b+14) - *(a+11) * *(b+15);
    *(c+9) += *(a+6)  * *(b+3)  + *(a+7)  * *(b+2)    +
    	      *(a+8)  * *(b+9)  + *(a+9)  * *(b+8)    +
    	      *(a+10) * *(b+15) + *(a+11) * *(b+14);

    *(c+10)+= *(a+6)  * *(b+4)  - *(a+7)  * *(b+5)    +
    	      *(a+8)  * *(b+10) - *(a+9)  * *(b+11)   +
    	      *(a+10) * *(b+16) - *(a+11) * *(b+17);
    *(c+11)+= *(a+6)  * *(b+5)  + *(a+7)  * *(b+4)    +
    	      *(a+8)  * *(b+11) + *(a+9)  * *(b+10)   +
    	      *(a+10) * *(b+17) + *(a+11) * *(b+16);

    *(c+12)+= *(a+12) * *b      - *(a+13) * *(b+1)    +
    	      *(a+14) * *(b+6)  - *(a+15) * *(b+7)    +
    	      *(a+16) * *(b+12) - *(a+17) * *(b+13);
    *(c+13)+= *(a+12) * *(b+1)  + *(a+13) * *b        +
    	      *(a+14) * *(b+7)  + *(a+15) * *(b+6)    +
    	      *(a+16) * *(b+13) + *(a+17) * *(b+12);

    *(c+14)+= *(a+12) * *(b+2)  - *(a+13) * *(b+3)    +
    	      *(a+14) * *(b+8)  - *(a+15) * *(b+9)    +
    	      *(a+16) * *(b+14) - *(a+17) * *(b+15);
    *(c+15)+= *(a+12) * *(b+3)  + *(a+13) * *(b+2)    +
    	      *(a+14) * *(b+9)  + *(a+15) * *(b+8)    +
    	      *(a+16) * *(b+15) + *(a+17) * *(b+14);

    *(c+16)+= *(a+12) * *(b+4)  - *(a+13) * *(b+5)    +
    	      *(a+14) * *(b+10) - *(a+15) * *(b+11)   +
    	      *(a+16) * *(b+16) - *(a+17) * *(b+17);
    *(c+17)+= *(a+12) * *(b+5)  + *(a+13) * *(b+4)    +
    	      *(a+14) * *(b+11) + *(a+15) * *(b+10)   +
    	      *(a+16) * *(b+17) + *(a+17) * *(b+16);
}


//---------------------------------------------------------------//

//
//  y   =  U x
//
void uDotXEqual(IFloat* y, const IFloat* u, const IFloat* x)
{
    *y     =  *u      * *x     - *(u+1)  * *(x+1) + *(u+2)  * *(x+2)
	    - *(u+3)  * *(x+3) + *(u+4)  * *(x+4) - *(u+5)  * *(x+5);
    *(y+1) =  *u      * *(x+1) + *(u+1)  * *x     + *(u+2)  * *(x+3)
	    + *(u+3)  * *(x+2) + *(u+4)  * *(x+5) + *(u+5)  * *(x+4);
    *(y+2) =  *(u+6)  * *x     - *(u+7)  * *(x+1) + *(u+8)  * *(x+2)
	    - *(u+9)  * *(x+3) + *(u+10) * *(x+4) - *(u+11) * *(x+5);
    *(y+3) =  *(u+6)  * *(x+1) + *(u+7)  * *x     + *(u+8)  * *(x+3)
	    + *(u+9)  * *(x+2) + *(u+10) * *(x+5) + *(u+11) * *(x+4);
    *(y+4) =  *(u+12) * *x     - *(u+13) * *(x+1) + *(u+14) * *(x+2)
	    - *(u+15) * *(x+3) + *(u+16) * *(x+4) - *(u+17) * *(x+5);
    *(y+5) =  *(u+12) * *(x+1) + *(u+13) * *x     + *(u+14) * *(x+3)
	    + *(u+15) * *(x+2) + *(u+16) * *(x+5) + *(u+17) * *(x+4);
}


//---------------------------------------------------------------//

IFloat dotProduct(const IFloat *a, const IFloat *b, int len)
{
    IFloat sum = 0.0;
    for(int i = 0; i < len; ++i) {
    	sum += *a++ * *b++;
    }
    return sum;
}

void vecTimesEquFloat(IFloat *a, IFloat b, int len)
{
    for(int i = 0; i < len; ++i) {
    	*a++ *= b;
    }
}

void vecAddEquVec(IFloat *a, const IFloat *b, int len)
{
    for(int i = 0; i < len; ++i) {
    	*a++ += *b++;
    }
}

void vecMinusEquVec(IFloat *a, const IFloat *b, int len)
{
    for(int i = 0; i < len; ++i) {
    	*a++ -= *b++;
    }
}

void fTimesV1PlusV2(IFloat *a, IFloat b, const IFloat *c,
	const IFloat *d, int len)
{
    for(int i = 0; i < len; ++i) {
    	*a++ = b * *c++ + *d++;
    }
}

void fTimesV1MinusV2(IFloat *a, IFloat b, const IFloat *c,
	const IFloat *d, int len)
{
    for(int i = 0; i < len; ++i) {
    	*a++ = b * *c++ - *d++;
    }
}

void oneMinusfTimesMatrix(IFloat *a, IFloat b, const IFloat *c,
	int n)
{
    IFloat *p = a;
    for(int i = 0; i < n; ++i) {
        *p++ = -b * *c++;
    }
    *a += 1.0;    *(a+8) += 1.0;    *(a+16) += 1.0;
}

void vecNegative(IFloat *a, const IFloat *b, int len)
{
    for(int i = 0; i < len; ++i) {
        *a++ = -*b++;
    }
}

void compDotProduct(IFloat *c_r, IFloat *c_i, 
		    const IFloat *a, const IFloat *b, int len)
{
    *c_r = *c_i = 0.0;
    for(int i = 0; i < len; i += 2, a += 2, b += 2) 
    {
      *c_r += *a * *b     + *(a+1) * *(b+1);   // real part
      *c_i += *a * *(b+1) - *(a+1) * *b;       // imag part
    }
}

void cTimesV1PlusV2(IFloat *a, IFloat re, IFloat im, const IFloat *c,
	const IFloat *d, int len)
{
    for(int i = 0; i < len; i += 2, c += 2) 
    {
      *a++ = re * *c     - im * *(c+1) + *d++;   // real part
      *a++ = re * *(c+1) + im * *c     + *d++;   // imag part
    }
}



CPS_END_NAMESPACE
