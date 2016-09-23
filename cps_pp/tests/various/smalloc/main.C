#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2008/02/08 18:35:09 $
//  $Header: /space/cvs/cps/cps++/tests/various/smalloc/main.C,v 1.6 2008/02/08 18:35:09 chulwoo Exp $
//  $Id: main.C,v 1.6 2008/02/08 18:35:09 chulwoo Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $RCSfile: main.C,v $
//  $Revision: 1.6 $
//  $Source: /space/cvs/cps/cps++/tests/various/smalloc/main.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
CPS_END_NAMESPACE
#include <util/qcdio.h>
#include<config.h>
#include<util/smalloc.h>
#include<util/verbose.h>
CPS_START_NAMESPACE

#ifdef PARALLEL
CPS_END_NAMESPACE
#include <comms/sysfunc_cps.h>
CPS_START_NAMESPACE
#endif

Verbose VRB;

int
main(int, char **)
{
  
  char * p[10];
  for (int i = 0 ; i < 10; i++) {
    p[i] = (char *)smalloc(1);
    printf("%x\n", int(p[i]) );
  }
  
  sfree(p[1]);
  sfree(p[0]);
  
  printf("%x ....\n", smalloc(5000) );
  sclear();
  for (i = 0 ; i < 15; i++) {
    char * p1 = (char *)smalloc(97);
    printf("%x\n", int(p1) );
  }

}

CPS_END_NAMESPACE
