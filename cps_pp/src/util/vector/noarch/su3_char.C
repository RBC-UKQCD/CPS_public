#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Definitions of functions computing 
  the characters of some SU(3) representations.

  $Id: su3_char.C,v 1.2 2003-07-24 16:53:54 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2003-07-24 16:53:54 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/vector/noarch/su3_char.C,v 1.2 2003-07-24 16:53:54 zs Exp $
//  $Id: su3_char.C,v 1.2 2003-07-24 16:53:54 zs Exp $
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
//  Revision 1.2  2001/06/19 18:13:40  anj
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
//  $RCSfile: su3_char.C,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/vector/noarch/su3_char.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/data_types.h>
CPS_START_NAMESPACE

static IFloat f1o3 = 1.0 / 3.0 ;
static IFloat f2o3 = 2.0 / 3.0 ;

IFloat reChar6(IFloat *p) ;
IFloat imChar6(IFloat *p) ;

IFloat reChar8(IFloat *p) ;

IFloat reChar10(IFloat *p) ;
IFloat imChar10(IFloat *p) ;

IFloat reChar6(IFloat *p)
{
//  return   *(p   ) * *(p   ) - *(p+ 1) * *(p+ 1)
//         + *(p+ 8) * *(p+ 8) - *(p+ 9) * *(p+ 9)
//         + *(p+16) * *(p+16) - *(p+17) * *(p+17)
//         + *(p   ) * *(p+ 8) - *(p+ 1) * *(p+ 9)
//         + *(p   ) * *(p+16) - *(p+ 1) * *(p+17)
//         + *(p+ 8) * *(p+16) - *(p+ 9) * *(p+17)
//         + *(p+ 2) * *(p+ 6) - *(p+ 3) * *(p+ 7)
//         + *(p+ 4) * *(p+12) - *(p+ 5) * *(p+13)
//         + *(p+10) * *(p+14) - *(p+11) * *(p+15) ;

  return   *(p   ) * ( *(p   ) + *(p+ 8) + *(p+16) )
         - *(p+ 1) * ( *(p+ 1) + *(p+ 9) + *(p+17) )
         + *(p+ 8) * ( *(p+ 8) + *(p+16) )
	 - *(p+ 9) * ( *(p+ 9) + *(p+17) )
         + *(p+16) * *(p+16) - *(p+17) * *(p+17)
         + *(p+ 2) * *(p+ 6) - *(p+ 3) * *(p+ 7)
         + *(p+ 4) * *(p+12) - *(p+ 5) * *(p+13)
         + *(p+10) * *(p+14) - *(p+11) * *(p+15) ;
}


IFloat imChar6(IFloat *p)
{
//  return   *(p   ) * *(p+ 1) + *(p   ) * *(p+ 1)
//         + *(p+ 8) * *(p+ 9) + *(p+ 8) * *(p+ 9)
//         + *(p+16) * *(p+17) + *(p+16) * *(p+17)
//         + *(p   ) * *(p+ 9) + *(p+ 1) * *(p+ 8)
//         + *(p   ) * *(p+17) + *(p+ 1) * *(p+16)
//         + *(p+ 8) * *(p+17) + *(p+ 9) * *(p+16)
//         + *(p+ 2) * *(p+ 7) + *(p+ 3) * *(p+ 6)
//         + *(p+ 4) * *(p+13) + *(p+ 5) * *(p+12)
//         + *(p+10) * *(p+15) + *(p+11) * *(p+14) ;

  return   *(p   ) * ( *(p+ 1) + *(p+ 9) + *(p+17) )
         + *(p+ 1) * ( *(p   ) + *(p+ 8) + *(p+16) )
         + *(p+ 8) * ( *(p+ 9) + *(p+17) )
	 + *(p+ 9) * ( *(p+ 8) + *(p+16) )
         + *(p+16) * *(p+17) + *(p+16) * *(p+17)
         + *(p+ 2) * *(p+ 7) + *(p+ 3) * *(p+ 6)
         + *(p+ 4) * *(p+13) + *(p+ 5) * *(p+12)
         + *(p+10) * *(p+15) + *(p+11) * *(p+14) ;
}

IFloat reChar8(IFloat *p) 
{
  return    2.0 * (   *(p   ) * ( *(p+ 8) + *(p+16) )
                    + *(p+ 1) * ( *(p+ 9) + *(p+17) )
                    + *(p+ 8) * *(p+16) + *(p+ 9) * *(p+17) )
	 + f2o3 * (   *(p   ) * *(p   ) + *(p+ 1) * *(p+ 1)
                    + *(p+ 8) * *(p+ 8) + *(p+ 9) * *(p+ 9)
                    + *(p+16) * *(p+16) + *(p+17) * *(p+17) )
         - f1o3 * (   *(p+ 2) * *(p+ 2) + *(p+ 3) * *(p+ 3)
                    + *(p+ 4) * *(p+ 4) + *(p+ 5) * *(p+ 5)
                    + *(p+ 6) * *(p+ 6) + *(p+ 7) * *(p+ 7)
                    + *(p+10) * *(p+10) + *(p+11) * *(p+11)
                    + *(p+12) * *(p+12) + *(p+13) * *(p+13)
                    + *(p+14) * *(p+14) + *(p+15) * *(p+15) ) ;
}

IFloat reChar10(IFloat *p)
{
  return 0.0 ;
}

IFloat imChar10(IFloat *p)
{
  return 0.0 ;
}

CPS_END_NAMESPACE
