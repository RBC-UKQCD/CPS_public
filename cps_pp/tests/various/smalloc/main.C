#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2003-10-31 14:15:34 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/various/smalloc/main.C,v 1.2 2003-10-31 14:15:34 zs Exp $
//  $Id: main.C,v 1.2 2003-10-31 14:15:34 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: main.C,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/various/smalloc/main.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
CPS_END_NAMESPACE
#include <stdio.h>
#include<config.h>
#include<util/smalloc.h>
#include<util/verbose.h>
CPS_START_NAMESPACE

#ifdef PARALLEL
CPS_END_NAMESPACE
#include <comms/sysfunc.h>
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
