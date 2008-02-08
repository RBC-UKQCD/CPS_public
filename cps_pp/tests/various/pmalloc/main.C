#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2008-02-08 18:35:09 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/various/pmalloc/main.C,v 1.6 2008-02-08 18:35:09 chulwoo Exp $
//  $Id: main.C,v 1.6 2008-02-08 18:35:09 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: main.C,v $
//  $Revision: 1.6 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/various/pmalloc/main.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
CPS_END_NAMESPACE
#include <util/qcdio.h>
#include<config.h>
#include<util/pmalloc.h>
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
    p[i] = (char *)pmalloc(1);
    printf("%x\n", int(p[i]) );
  }
  
  pfree(p[1]);
  pfree(p[0]);
  
  printf("%x ....\n", pmalloc(5000) );
  pclear();
  for (i = 0 ; i < 15; i++) {
    char * p1 = (char *)pmalloc(97);
    printf("%x\n", int(p1) );
  }

}

CPS_END_NAMESPACE
