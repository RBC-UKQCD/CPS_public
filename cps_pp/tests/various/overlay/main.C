#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-06-04 21:14:18 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/various/overlay/main.C,v 1.2 2004-06-04 21:14:18 chulwoo Exp $
//  $Id: main.C,v 1.2 2004-06-04 21:14:18 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.2 $
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
