#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2008/02/08 18:35:09 $
//  $Header: /space/cvs/cps/cps++/tests/various/tools/main.C,v 1.6 2008/02/08 18:35:09 chulwoo Exp $
//  $Id: main.C,v 1.6 2008/02/08 18:35:09 chulwoo Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $RCSfile: main.C,v $
//  $Revision: 1.6 $
//  $Source: /space/cvs/cps/cps++/tests/various/tools/main.C,v $
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
#include <comms/sysfunc_cps.h>
CPS_START_NAMESPACE
#endif

int
main(int, char **)
{
  for(int i=0; i<1000; i++)
    printf("%d\n",0);

}

CPS_END_NAMESPACE
