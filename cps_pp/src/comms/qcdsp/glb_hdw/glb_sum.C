#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-06-04 21:14:04 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp/glb_hdw/glb_sum.C,v 1.3 2004-06-04 21:14:04 chulwoo Exp $
//  $Id: glb_sum.C,v 1.3 2004-06-04 21:14:04 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.3 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp/glb_hdw/glb_sum.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
// glb_sum.C
//
// Wrapper for Ping's optimized global sum. This will
// eventually be part of the qos.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include<comms/glb.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/gjp.h>
#include<global_sum.h>
CPS_START_NAMESPACE

static int is_initialized = 0;

void glb_sum(Float *float_p)
{
  char *cname = " ";
  char *fname = "glb_sum";

  //----------------------------------------------------------------
  // If this is the first call initialize the global sum parameters
  //----------------------------------------------------------------
  if(is_initialized == 0){
    global_sum_init(GJP.GsumFastMode(), GJP.GsumMaxTry());
    VRB.Flow(cname, fname,"Initialized global sum\n");
    is_initialized = 1; 
  }

  //----------------------------------------------------------------
  // Do the global sum over the full machine 
  //----------------------------------------------------------------
  Float sum;
  IFloat sum_f;
  int partitions;
  sum = *float_p;
  sum_f = sum;
  int exp = global_sum(&sum_f);

  //----------------------------------------------------------------
  // Divide the sum by the number of partitions.
  //----------------------------------------------------------------
  sum = sum_f;
  partitions = SizeX() / GJP.Xnodes();
  partitions = partitions * ( SizeY() / GJP.Ynodes() );
  partitions = partitions * ( SizeZ() / GJP.Znodes() );
  partitions = partitions * ( SizeT() / GJP.Tnodes() );
  sum = sum / Float(partitions);
  *float_p = sum;
}


CPS_END_NAMESPACE
