#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-08-18 11:58:13 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/various/tools/main.C,v 1.5 2004-08-18 11:58:13 zs Exp $
//  $Id: main.C,v 1.5 2004-08-18 11:58:13 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: main.C,v $
//  $Revision: 1.5 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/various/tools/main.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
CPS_END_NAMESPACE
#include <util/qcdio.h>
#include <stdlib.h>	// exit()
#include<config.h>
CPS_START_NAMESPACE

#ifdef PARALLEL
CPS_END_NAMESPACE
#include <comms/sysfunc.h>
CPS_START_NAMESPACE
#endif

int
main(int, char **)
{
  for(int i=0; i<1000; i++)
    printf("%d\n",0);

}

CPS_END_NAMESPACE
