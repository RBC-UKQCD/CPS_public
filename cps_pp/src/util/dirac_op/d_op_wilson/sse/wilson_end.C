#ifdef USE_SSE
#include <config.h>
CPS_START_NAMESPACE
/*! \file
  \brief  Routine used internally in the DiracOpWilson class.

  $Id: wilson_end.C,v 1.2 2011-02-26 00:19:27 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2011-02-26 00:19:27 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/sse/wilson_end.C,v 1.2 2011-02-26 00:19:27 chulwoo Exp $
//  $Id: wilson_end.C,v 1.2 2011-02-26 00:19:27 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/sse/wilson_end.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/wilson.h>
#include <util/smalloc.h>
#include <util/verbose.h>
//for SSEOMP part
#include <util/gjp.h>
CPS_START_NAMESPACE

void wilson_end_free(int dir, Wilson *wilson_p)
{
  int idx;

  idx=dir;

  QMP_free_msghandle(wilson_p->multiple[idx]);
  QMP_free_msgmem(wilson_p->msgmem[idx][0]);
  QMP_free_msgmem(wilson_p->msgmem[idx][1]);
  free(wilson_p->recv_buf[idx]);
  
  idx = dir +4;
  QMP_free_msghandle(wilson_p->multiple[idx]);
  QMP_free_msgmem(wilson_p->msgmem[idx][0]);
  QMP_free_msgmem(wilson_p->msgmem[idx][1]);
  free(wilson_p->recv_buf[idx]);
}
    
void wilson_end( Wilson *wilson_p)
{
  char *cname = " ";
  char *fname = "wilson_end(Wilson*)";
  VRB.Func(cname,fname);

  VRB.Sfree(cname,fname, "af[1]", wilson_p->af[1]);
  sfree(wilson_p->af[1]);

  VRB.Sfree(cname,fname, "af[0]", wilson_p->af[0]);
  sfree(wilson_p->af[0]);


  // for SSEOMP
  if(GJP.Xnodes()!=1){
    int dir = 0;  wilson_end_free(dir, wilson_p);
  }
  if(GJP.Ynodes()!=1){
    int dir = 1;  wilson_end_free(dir, wilson_p);
  }
  if(GJP.Znodes()!=1){
    int dir = 2;  wilson_end_free(dir, wilson_p);
  }
  if(GJP.Tnodes()!=1){
    int dir = 3;  wilson_end_free(dir, wilson_p);
  }

#if 0
#pragma omp parallel for schedule(static,1)  
  for(int i=0; i<num_threads; ++i) free( omp_chi[i] );
#endif


  VRB.Sfree(cname,fname, "ptr", wilson_p->ptr);
  sfree(wilson_p->ptr);




}

CPS_END_NAMESPACE
#endif
