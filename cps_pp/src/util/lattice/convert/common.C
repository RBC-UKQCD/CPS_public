#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:46 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/convert/common.C,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
//  $Id: common.C,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.5  2001/08/16 10:50:32  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.2  2001/06/19 18:13:21  anj
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
//  $RCSfile: common.C,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/convert/common.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include<util/data_types.h>
CPS_START_NAMESPACE

#ifdef __cplusplus
extern "C" {
#endif

void negate_link(unsigned link_size, IFloat *link) ;
void site2cram(IFloat *src, IFloat *dst, unsigned site_size) ;
void site2dram(IFloat *src, IFloat *dst, unsigned *link_tbl, unsigned site_size) ;

#ifdef __cplusplus
}
#endif

void negate_link(unsigned link_size, IFloat *link)
{
	unsigned idx ;

	for (idx=0; idx<link_size; idx++) *(link+idx) = - *(link+idx) ;
}


void site2cram(IFloat *src, IFloat *dst, unsigned site_size)
{
	unsigned offset ;

	for (offset=0; offset<site_size; offset++)
		*(dst+offset) = *(src+offset) ;
}


void site2dram(IFloat *src, IFloat *dst, unsigned *link_tbl, unsigned site_size)
{
	unsigned offset ;

	for (offset=0; offset<site_size; offset++)
		*(dst+*(link_tbl+offset)) = *(src+offset) ;
}
CPS_END_NAMESPACE
