/*! \file
  \brief  Definition of parallel transport definitions for QCDOC.
  
  $Id: pt_vvpd.C,v 1.2 2006/12/14 17:54:28 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2006/12/14 17:54:28 $
//  $Header: /space/cvs/cps/cps++/src/util/parallel_transport/pt_base/qcdoc/pt_vvpd.C,v 1.2 2006/12/14 17:54:28 chulwoo Exp $
//  $Id: pt_vvpd.C,v 1.2 2006/12/14 17:54:28 chulwoo Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $RCSfile: pt_vvpd.C,v $
//  $Revision: 1.2 $
//  $Source: /space/cvs/cps/cps++/src/util/parallel_transport/pt_base/qcdoc/pt_vvpd.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
#include "asq_data_types.h"
#include "pt_int.h"
#include "pt_qcdoc.h"

/*! 
  Computes sum[x] = vect[x] vect[x + hop dir]^dagger
  where the sum is over n_vect vectors and the hop is in a forward direction.
*/
void PT::vvpd(IFloat **vect, int n_vect, const int *dir, 
	     int n_dir, int hop, IFloat **sum){

#define NEW_VVPD
  //#undef NEW_VVPD
  //---------------------------------------------------------------------
  //Adapt old vvpd for usage with new interface

#ifdef NEW_VVPD
  //Convert direction array
  int newdir[n_dir];
  for(int i = 0; i < n_dir; i++)
    newdir[i] = 2*dir[i];
  IFloat ***vect_shift = (IFloat ***)FastAlloc("","","vect_shift",n_vect*sizeof(IFloat**));
  for(int i = 0; i < n_vect; i++)
    {
      vect_shift[i] = (IFloat**)FastAlloc("","","vect_shift[i]",n_dir*sizeof(IFloat *));
      for(int j = 0; j < n_dir; j++)
	vect_shift[i][j] = vect[i];
    }
  vvpd(vect, vect_shift, n_vect, newdir, n_dir, hop, sum, 1);
  for(int i = 0; i < n_vect; i++)
    Free(vect_shift[i]);
  Free(vect_shift);
  //----------------------------------------------------------------------
  //Old vvpd
#else
  char *fname = "pt_vvpd()";
//  VRB.Func("",fname);
  int i, v;
  Float f = 2.0;
  int wire[n_dir];
  for(i=0;i<n_dir;i++) wire[i] = dir[i]; // from (x,y,z,t) to (t,x,y,z)

  #ifdef TESTING
  printf("vvpd called.  Checking input data.\n");
  for(int kk = 0; kk < 1; kk++)
  for(int count = 0; count < 6*vol; count++)
    printf("vect[%d][%d]=%f\n",kk,count,(Float)vect[kk][count]);
  printf("End input data check.\n");
  #endif

  SCUDirArgIR *SCUarg_p[2*n_dir];  
  SCUDirArgIR *SCUarg_p2[2*n_dir];  

  // Only do communication in forward direction
  int comms=0;
  for(i=0;i<n_dir;i++)
  if( !local[wire[i]]) {
    if ( size[wire[i]] <hop)
      fprintf(stderr, 
		"%s::size(%d) in direction %d is smaller than the hop(%d)\n",
		fname,size[wire[i]],wire[i],hop);
    SCUarg_p[2*comms] = SCUarg[hop-1][4*wire[i]];
    SCUarg_p[2*comms+1] = SCUarg[hop-1][4*wire[i]+1];
    SCUarg_p2[2*comms] = SCUarg2[hop-1][4*wire[i]];
    SCUarg_p2[2*comms+1] = SCUarg2[hop-1][4*wire[i]+1];
    comms++;
  }

  for(v=0; v<n_vect; v++){
    SCUDirArgMulti SCUmulti;

    if (v%2==0) {
      comms=0;
      for(i=0;i<n_dir;i++)
        if( !local[wire[i]]){ 
	  SCUarg_p[2*comms+1]->Addr((void *)(vect[v]+VECT_LEN*set_offset(2*wire[i], hop)));
          comms++;
        }

      // Start communication
      if (comms) SCUmulti.Init(SCUarg_p,2*comms);
    } else {
      comms=0;
      for(i=0;i<n_dir;i++)
        if( !local[wire[i]]){ 
	  SCUarg_p2[2*comms+1]->Addr((void *)(vect[v]+VECT_LEN*set_offset(2*wire[i], hop)));
          comms++;
      }
      // Start communication
      if (comms) SCUmulti.Init(SCUarg_p2,2*comms);
    } 

    if (comms) SCUmulti.SlowStartTrans();
    // Finalise communication
    if (comms) SCUmulti.TransComplete();

    // Perform non-local calculation for previous v
    if (v>0)
      if (v==1) {
	for(i=0; i<n_dir; i++) 
	  if(non_local_chi[2*wire[i]] > 0)
	    cross_over_lin(sum[i], &f, vect[v-1],rcv_buf[2*wire[i]], hop*non_local_chi[2*wire[i]],
		src_nl[hop-1][2*wire[i]], dest_nl[hop-1][2*wire[i]]);
      } else if (v%2==1) {
	for(i=0; i<n_dir; i++) 
	  if(non_local_chi[2*wire[i]] > 0)
	    cross_lin(sum[i], &f, vect[v-1],rcv_buf[2*wire[i]], hop*non_local_chi[2*wire[i]],
		src_nl[hop-1][2*wire[i]], dest_nl[hop-1][2*wire[i]]);
      } else {
	for(i=0; i<n_dir; i++) 
	  if(non_local_chi[2*wire[i]] > 0)
	    cross_lin(sum[i], &f,vect[v-1],rcv_buf2[2*wire[i]], hop*non_local_chi[2*wire[i]],
		src_nl[hop-1][2*wire[i]], dest_nl[hop-1][2*wire[i]]);
      }
    
    // Perform local calculation for current v
    if (v==0)
      for(i=0; i<n_dir; i++)
	{
	  #ifdef TESTING
	  printf("wire[%d] = %d\n",i,wire[i]);
	  printf("sum[%d] = %f\n",i,(float)*(sum[i]));
	  printf("vol-hop*non_local_chi[%d]=%d\n",2*wire[i],vol-hop*non_local_chi[2*wire[i]]);
	  for(int count = 0; count < vol-hop*non_local_chi[2*wire[i]]; count++)
	    {
	      printf("src_l[%d][%d][%d] = %ld\n",hop-1,2*wire[i],count,(int)src_l[hop-1][2*wire[i]][count]);
	      printf("dest_l[%d][%d][%d] = %ld\n",hop-1,2*wire[i],count,(int)dest_l[hop-1][2*wire[i]][count]);
	      for(int v_count = 0; v_count < 6; v_count++)
		{
		  printf("vect[%d][%d] = %f\n",v,(6*((int)src_l[hop-1][2*wire[i]][count])+v_count),(Float)(vect[v][6*((int)src_l[hop-1][2*wire[i]][count])+v_count]));
		  printf("vect[%d][%d] = %f\n",v,(6*((int)dest_l[hop-1][2*wire[i]][count])+v_count),(Float)(vect[v][(6*((int)dest_l[hop-1][2*wire[i]][count])+v_count)]));
		}
	    }
	  printf("cross_over_look called\n");
	  #endif
	  if((vol-hop*non_local_chi[2*wire[i]]) > 0)
	    cross_over_look(sum[i], &f, vect[v], vect[v], vol-hop*non_local_chi[2*wire[i]], src_l[hop-1][2*wire[i]], dest_l[hop-1][2*wire[i]]);
        }
    else
      {
      for(i=0; i<n_dir; i++)
	if((vol-hop*non_local_chi[2*wire[i]])> 0)
	  cross_look(sum[i], &f, vect[v], vect[v], vol-hop*non_local_chi[2*wire[i]], src_l[hop-1][2*wire[i]], dest_l[hop-1][2*wire[i]]);
      }

  }

  if (v==1) {
    for(i=0; i<n_dir; i++) 
      if(non_local_chi[2*wire[i]] > 0)
	cross_over_lin(sum[i], &f, vect[v-1],rcv_buf[2*wire[i]], hop*non_local_chi[2*wire[i]],
	    src_nl[hop-1][2*wire[i]], dest_nl[hop-1][2*wire[i]]);
  } else if (v%2==1) {
    for(i=0; i<n_dir; i++) 
      if(non_local_chi[2*wire[i]] > 0)
	cross_lin(sum[i], &f, vect[v-1],rcv_buf[2*wire[i]], hop*non_local_chi[2*wire[i]],
		  src_nl[hop-1][2*wire[i]], dest_nl[hop-1][2*wire[i]]);
  } else {
    for(i=0; i<n_dir; i++)
      if(non_local_chi[2*wire[i]] > 0)
	cross_lin(sum[i], &f,vect[v-1],rcv_buf2[2*wire[i]], hop*non_local_chi[2*wire[i]],
		  src_nl[hop-1][2*wire[i]], dest_nl[hop-1][2*wire[i]]);
  }  
//  ParTrans::PTflops += 90*n_vect*n_dir*vol;
#endif
  //------------------------------------------------------------------------
}

/*! 
  Computes sum[x] = vect2[x] vect[x + hop dir]^dagger
  where the sum is over n_vect vectors and the hop is in forward
  or negative direction.
*/
void PT::vvpd(IFloat **vect2, IFloat ***vect, int n_vect, const int *dir, int n_dir, int hop, IFloat **sum, int overwrite){

  char *fname = "pt_vvpd()";
//  VRB.Func("",fname);
  int i, s, v;
  Float f = 2.0;
  int wire[n_dir];
  for(i=0;i<n_dir;i++) wire[i] = dir[i]; // from (x,y,z,t) to (t,x,y,z)

  #ifdef TESTING
  for(i = 0; i < n_dir; i++)
    {
      printf("wire[%d] = %d\n",i,wire[i]);
      for(int count = 0; count < vol - hop*non_local_chi[wire[i]]; count++)
	{
	  printf("src_l[%d][%d][%ld] = %d\n",hop-1,wire[i],count,src_l[hop-1][wire[i]][count]);
	  printf("dest_l[%d][%d][%ld] = %d\n",hop-1,wire[i],count,dest_l[hop-1][wire[i]][count]);
	}
     for(int count = 0; count < hop*non_local_chi[wire[i]]; count++)
	{
	  printf("src_nl[%d][%d][%ld] = %d\n",hop-1,wire[i],count,src_nl[hop-1][wire[i]][count]);
	  printf("dest_nl[%d][%d][%ld] = %d\n",hop-1,wire[i],count,dest_nl[hop-1][wire[i]][count]);
	}
    }
  #endif

  SCUDirArgIR *SCUarg_p[2*n_dir];  
  SCUDirArgIR *SCUarg_p2[2*n_dir];  

  //Setup communciation
  int comms=0;
  for(i=0;i<n_dir;i++)
  if( !local[wire[i]/2]) {
    if ( size[wire[i]/2] <hop)
      fprintf(stderr, 
		"%s:size(%d) in direction %d is smaller than the hop(%d)\n",
		fname,size[wire[i]],wire[i],hop);
    SCUarg_p[2*comms] = SCUarg[hop-1][2*wire[i]];
    SCUarg_p[2*comms+1] = SCUarg[hop-1][2*wire[i]+1];
    SCUarg_p2[2*comms] = SCUarg2[hop-1][2*wire[i]];
    SCUarg_p2[2*comms+1] = SCUarg2[hop-1][2*wire[i]+1];
    comms++;
  }

  for(v=0; v<n_vect; v++){
    SCUDirArgMulti SCUmulti;

    if (v%2==0) {
      comms=0;
      for(i=0;i<n_dir;i++)
        if( !local[wire[i]/2]){ 
	  SCUarg_p[2*comms+1]->Addr((void *)(vect[v][i]+VECT_LEN*set_offset(wire[i], hop)));
          comms++;
        }

      // Start communication
      if (comms) SCUmulti.Init(SCUarg_p,2*comms);
    } else {
      comms=0;
      for(i=0;i<n_dir;i++)
        if( !local[wire[i]/2]){ 
	  SCUarg_p2[2*comms+1]->Addr((void *)(vect[v][i]+VECT_LEN*set_offset(wire[i], hop)));
          comms++;
      }

      // Start communication
      if (comms) SCUmulti.Init(SCUarg_p2,2*comms);
    }
    if (comms) SCUmulti.SlowStartTrans();
    // Finalise communication
    if (comms) SCUmulti.TransComplete();

    // Perform non-local calculation for previous v
    if (v>0)
      if (v==1 && overwrite==1) {
	for(i=0; i<n_dir; i++)
	  if(non_local_chi[wire[i]]>0)
	    cross_over_lin(sum[i], &f, vect2[v-1],rcv_buf[wire[i]], hop*non_local_chi[wire[i]],
			   src_nl[hop-1][wire[i]], dest_nl[hop-1][wire[i]]);
      } else if (v%2==1) {
	for(i=0; i<n_dir; i++) 
	  if(non_local_chi[wire[i]]>0)
	    cross_lin(sum[i], &f, vect2[v-1],rcv_buf[wire[i]], hop*non_local_chi[wire[i]],
		      src_nl[hop-1][wire[i]], dest_nl[hop-1][wire[i]]);
      } else {
	for(i=0; i<n_dir; i++) 
	  if(non_local_chi[wire[i]]>0)
	    cross_lin(sum[i], &f,vect2[v-1],rcv_buf2[wire[i]], hop*non_local_chi[wire[i]],
		      src_nl[hop-1][wire[i]], dest_nl[hop-1][wire[i]]);
      }
    
    // Perform local calculation for current v
    if (v==0 && overwrite==1)
      {
	for(i=0; i<n_dir; i++)
	  if((vol-hop*non_local_chi[wire[i]])>0)
	    cross_over_look(sum[i], &f, vect2[v], vect[v][i], vol-hop*non_local_chi[wire[i]], src_l[hop-1][wire[i]], dest_l[hop-1][wire[i]]);
      }
    else
      {
	for(i=0; i<n_dir; i++)
	  if((vol-hop*non_local_chi[wire[i]])>0)
	    cross_look(sum[i], &f, vect2[v], vect[v][i], vol-hop*non_local_chi[wire[i]], src_l[hop-1][wire[i]], dest_l[hop-1][wire[i]]);
      }

  }

  if (v==1 && overwrite==1) {
    for(i=0; i<n_dir; i++) 
      if(non_local_chi[wire[i]]>0)
	cross_over_lin(sum[i], &f, vect2[v-1],rcv_buf[wire[i]], hop*non_local_chi[wire[i]],
		       src_nl[hop-1][wire[i]], dest_nl[hop-1][wire[i]]);
  } else if (v%2==1) {
    for(i=0; i<n_dir; i++)
      if(non_local_chi[wire[i]]>0)
	cross_lin(sum[i], &f, vect2[v-1],rcv_buf[wire[i]], hop*non_local_chi[wire[i]],
		  src_nl[hop-1][wire[i]], dest_nl[hop-1][wire[i]]);
  } else {
    for(i=0; i<n_dir; i++) 
      if(non_local_chi[wire[i]]>0)
	cross_lin(sum[i], &f,vect2[v-1],rcv_buf2[wire[i]], hop*non_local_chi[wire[i]],
		  src_nl[hop-1][wire[i]], dest_nl[hop-1][wire[i]]);
  }
  
  //  ParTrans::PTflops += 90*n_vect*n_dir*vol;
}
