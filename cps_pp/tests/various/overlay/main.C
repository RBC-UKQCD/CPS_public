#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:47 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/various/overlay/main.C,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
//  $Id: main.C,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.2  2001/06/19 18:12:31  anj
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
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/various/overlay/main.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
/*
 *  test for alg
 */
CPS_END_NAMESPACE
#include <string.h>  //memmove()
CPS_START_NAMESPACE

extern unsigned int wfm_ram1_start;    // Physical starting address 
extern unsigned int wfm_ram1_end;      // Physical ending address 
extern unsigned int wfm_ram1_dest;     // Virtual starting address


extern "C" {
void test(void);
}

int out_buf[100];

main(int argc,char *argv[])
{

  for(int ii=0; ii<100; ii++)
    out_buf[ii] = ii+1;

  memmove(&wfm_ram1_dest, &wfm_ram1_start, 
	  &wfm_ram1_end - &wfm_ram1_start);

  test();

}






CPS_END_NAMESPACE
