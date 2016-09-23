#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /space/cvs/cps/cps++/tests/various/overlay/main.C,v $
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
