/*! \file
  \brief  Definition of parallel transport definitions for QCDOC.
  
  $Id: pt_vec.C,v 1.1 2006/07/05 18:13:49 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2006/07/05 18:13:49 $
//  $Header: /space/cvs/cps/cps++/src/util/parallel_transport/pt_base/qcdoc/pt_vec.C,v 1.1 2006/07/05 18:13:49 chulwoo Exp $
//  $Id: pt_vec.C,v 1.1 2006/07/05 18:13:49 chulwoo Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $RCSfile: pt_vec.C,v $
//  $Revision: 1.1 $
//  $Source: /space/cvs/cps/cps++/src/util/parallel_transport/pt_base/qcdoc/pt_vec.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
#include "asq_data_types.h"
#include "pt_int.h"
#include "pt_qcdoc.h"

//Parallel transport of a vector defined on single parity sites
//
//Parameters
//
//n - number of directions in which to parallel transport
//vout - Transported vector
//vin - Initial vector
//dir - list of directions in which to transport the vectors
//cb - Parity of the sites where the vectors are defined
//gauge - Pointer to block of gauge fields in STAG order

#ifndef SCIDAC 

//Normal parallel transport with normal gauge fields
#undef PROFILE
void PT::vec_cb(int n, IFloat **vout, IFloat **vin, const int *dir, int
 parity)
{
  vec_cb_norm(n,vout,vin,dir,parity,gauge_field_addr);
}

//Normal parallel transport, but with user-specified gauge fields
#undef PROFILE
void PT::vec_cb(int n, IFloat **vout, IFloat **vin, const int *dir, int
parity, IFloat * new_gauge_field)
{
  vec_cb_norm(n,vout,vin,dir,parity,new_gauge_field);
}

//Padded parallel transport with normal gauge fields
#undef PROFILE
void PT::vec_cb(int n, IFloat *vout, IFloat **vin, const int *dir, int
parity, int pad)
{
  vec_cb_pad(n,vout,vin,dir,parity,gauge_field_addr);
}

//Padded parallel transport, but with user-specified gauge fields
#undef PROFILE
void PT::vec_cb(int n, IFloat *vout, IFloat **vin, const int *dir, int
parity, int pad, IFloat * new_gauge_field)
{
  vec_cb_pad(n,vout,vin,dir,parity,new_gauge_field);
}

#define PROFILE
#undef PROFILE
void PT::vec_cb_norm(int n, IFloat **vout, IFloat **vin, const int *dir,int parity, IFloat * gauge)
{
  //List of the different directions
  int wire[n];
  int i;
//  int j,d,s,k;
  //SCUDirArgs for sending and receiving in the n directions
  SCUDirArgIR *SCUarg_p[2*non_local_dirs];
  //SCUDirArgIR *SCUarg_p[2*n];
  SCUDirArgMulti SCUmulti;
  static int call_num = 0;
  int vlen = VECT_LEN;
  int vlen2 = VECT_LEN;
//  printf("gauge=%p\n",gauge);

  call_num++;
  
  //Name our function
//  char *fname="pt_1vec_cb_norm()";
  
  //Set the transfer directions
  //If wire[i] is even, then we have communication in the negative direction
  //If wire[i] is odd, then we have communication in the positive direction
  for(i=0;i<n;i++)
    wire[i]=dir[i];

  Float dtime;


  //If wire[i] is odd, then we have parallel transport in the
  //positive direction.  In this case, the matrix multiplication is
  //done before the field is transferred over to the adjacent node
  //
  //If we have transfer in the negative T direction (wire[i] = 6) then
  //we have to copy the appropriate fields into the send buffer
  if(conjugated)
    for(i=0;i<n;i++)
      {
	if(!local[wire[i]/2])
	  {
	    if(wire[i]%2)
	      {
		//printf("dir = %d, pre-mulitply\n", wire[i]);
#ifdef PROFILE
  dtime  = - dclock();
#endif

  partrans_cmv(non_local_chi_cb[wire[i]]/2,uc_nl_cb_pre[parity][wire[i]/2],gauge,vin[i],snd_buf_cb[wire[i]/2]);

#ifdef PROFILE
  dtime +=dclock();
  print_flops("",fname,66*non_local_chi_cb[wire[i]],dtime);
#endif
	      }
	    else if((wire[i] == 6))
	      {
#ifdef PROFILE
  dtime  = - dclock();
#endif
#if 1
            pt_copy_buffer(non_local_chi_cb[6],(long)vin[i],(long)snd_buf_t_cb,(long)Toffset[parity]);
#else
		for(j = 0; j < non_local_chi_cb[6];j++)
		  for(k = 0; k < VECT_LEN;k++)
		    *(snd_buf_t_cb+j*VECT_LEN+k) = *(vin[i] + *(Toffset[parity]+j)+ k);
		  //moveMem(snd_buf_t_cb + j*VECT_LEN,vin[i] + *(Toffset[parity]+j)*vlen,VECT_LEN*sizeof(IFloat));
#endif
#ifdef PROFILE
  dtime +=dclock();
  print_flops(fname,"pt_copy_buffer()",0,dtime);
#endif
	      }
	  }
      }
  else
    for(i=0;i<n;i++)
      {
	if(!local[wire[i]/2])
	  {
	    if(wire[i]%2)
	      {
#ifdef PROFILE
  dtime  = - dclock();
#endif
 
  partrans_cmv_dag(non_local_chi_cb[wire[i]]/2,uc_nl_cb_pre[parity][wire[i]/2],gauge,vin[i],snd_buf_cb[wire[i]/2]);	 

#ifdef PROFILE
  dtime +=dclock();
  print_flops("",fname,66*non_local_chi_cb[wire[i]],dtime);
#endif
	      }
	    else if(wire[i] == 6)
	      {
#ifdef PROFILE
  dtime  = - dclock();
#endif
#if 1
            pt_copy_buffer(non_local_chi_cb[6],(long)vin[i],(long)snd_buf_t_cb,(long)Toffset[parity]);
#else
		for(j = 0; j < non_local_chi_cb[6];j++)
		  for(k = 0; k < VECT_LEN;k++)
		    *(snd_buf_t_cb+j*VECT_LEN+k) = *(vin[i] + *(Toffset[parity]+j)+ k);
//		  moveMem(snd_buf_t_cb + j*VECT_LEN,vin[i] + *(Toffset[parity]+j),VECT_LEN*sizeof(IFloat));
#endif
#ifdef PROFILE
  dtime +=dclock();
  print_flops(fname,"pt_copy_buffer()",0,dtime);
#endif
	      }
	  }
      }
#ifdef PROFILE
  dtime  = - dclock();
#endif

  #if 1
  int comms = 0;
  for(i=0;i<n;i++)
    {
      if(!local[wire[i]/2])
	{
	  //Calculate the starting address for the data to be sent
	  IFloat *addr = vin[i] + VECT_LEN * offset_cb[wire[i]];
	  //This points to the appropriate SCUDirArg for receiving
	  SCUarg_p[2*comms] = SCUarg_cb[2*wire[i]];
	  //This points to the appropriate SCUDirArg for sending
	  SCUarg_p[2*comms+1] = SCUarg_cb[2*wire[i]+1];
	  
	  //Set the send address
	  if(wire[i]%2)
	    SCUarg_p[2*comms+1]->Addr((void *)snd_buf_cb[wire[i]/2]);
	  else if(wire[i] == 6)
	    SCUarg_p[2*comms+1]->Addr((void *)snd_buf_t_cb);
	  else
	    SCUarg_p[2*comms+1]->Addr((void *)addr);
	  comms++;
	}
    }
  #endif


  if(comms){
    SCUmulti.Init(SCUarg_p,2*comms);
//Begin transmission
    SCUmulti.SlowStartTrans();
  }

  //Do local calculations
  if(conjugated)
    {
      for(i=0;i<n;i++)
	{
#ifdef PROFILE
  dtime  = - dclock();
#endif 
	
  if(wire[i]%2)
    partrans_cmv(local_chi_cb[wire[i]]/2,uc_l_cb[parity][wire[i]],gauge,vin[i],vout[i]);
  else
    partrans_cmv_dag(local_chi_cb[wire[i]]/2,uc_l_cb[parity][wire[i]],gauge,vin[i],vout[i]);
  
#ifdef PROFILE
  dtime +=dclock();
  print_flops("",fname,66*local_chi_cb[wire[i]],dtime);
#endif
	}
    }
  else
    {
    for(i=0;i<n;i++)
      {
#ifdef PROFILE
  dtime  = - dclock();
#endif

  if(!(wire[i]%2))
    partrans_cmv(local_chi_cb[wire[i]]/2,uc_l_cb[parity][wire[i]],gauge,vin[i],vout[i]);
  else
    partrans_cmv_dag(local_chi_cb[wire[i]]/2,uc_l_cb[parity][wire[i]],gauge,vin[i],vout[i]);

#ifdef PROFILE
  dtime +=dclock();
  print_flops("",fname,66*local_chi_cb[wire[i]],dtime);
#endif
      }
    }
  
  //End transmission
  if(comms){ SCUmulti.TransComplete(); }

  //If wire[i] is even, then we have transport in the negative direction.
  //In this case, the vector field is multiplied by the SU(3) link matrix
  //after all communication is complete
  IFloat *fp0, *fp1;
  for(i=0;i<n;i++)
    {
      if(!local[wire[i]/2])
      	{
	  if(!(wire[i]%2))
	    {
#ifdef PROFILE
  dtime  = - dclock();
#endif

	    if(conjugated)
	      partrans_cmv_dag(non_local_chi_cb[wire[i]]/2,uc_nl_cb[parity][wire[i]],gauge,rcv_buf[wire[i]],vout[i]);
	    else
	      partrans_cmv(non_local_chi_cb[wire[i]]/2,uc_nl_cb[parity][wire[i]],gauge,rcv_buf[wire[i]],vout[i]);

#ifdef PROFILE
  dtime +=dclock();
  print_flops("",fname,66*non_local_chi_cb[wire[i]],dtime);
#endif
	    }
	  //Otherwise we have parallel transport in the positive direction.
	  //In this case, the received data has already been pre-multiplied
	  //All we need to do is to put the transported field in the correct place
	  else
	    {
#ifdef PROFILE
  dtime  = - dclock();
#endif
#if 1
              pt_copy(non_local_chi_cb[wire[i]]/2,uc_nl_cb[parity][wire[i]],rcv_buf[wire[i]],vout[i]);
#else
	      //Place the data in the receive buffer into the result vector
	      for(s=0;s<non_local_chi_cb[wire[i]];s++)
		{
		  fp0 = (IFloat *)((long)rcv_buf[wire[i]]+uc_nl_cb[parity][wire[i]][s].src);
		  fp1 = (IFloat *)((long)vout[i]+uc_nl_cb[parity][wire[i]][s].dest);
		  for(d = 0;d<VECT_LEN;d++)
		    *(fp1+d) = *(fp0+d);
		  //moveMem(fp1,fp0,VECT_LEN*sizeof(IFloat));
		}
#endif
#ifdef PROFILE
  dtime +=dclock();
  print_flops(fname,"pt_copy()",0,dtime);
#endif
	    }
	}
    }
//  ParTrans::PTflops +=33*n*vol;
}

#define PROFILE
#undef PROFILE

void PT::vec_cb_pad(int n, IFloat *vout, IFloat **vin, const int *dir,int parity, IFloat * gauge)
{
  //List of the different directions
  int wire[n];
  int i;
  //SCUDirArgs for sending and receiving in the n directions
  SCUDirArgIR *SCUarg_p[2*non_local_dirs];
  //SCUDirArgIR *SCUarg_p[2*n];
  SCUDirArgMulti SCUmulti;
  static int call_num = 0;
//  int vlen = VECT_LEN;
//  int vlen2 = 8;
#ifdef PROFILE
  printf("gauge=%p parity =%d\n",gauge,parity);
  for(i=0;i<n;i++){
    printf("%d: vin=%p vout=%p\n",i,vin[i],vout);
  }
#endif

  call_num++;
  
  //Name our function
//  char *fname="pt_1vec_cb_pad()";
  //VRB.Func("",fname);
  
  //Set the transfer directions
  //If wire[i] is even, then we have communication in the negative direction
  //If wire[i] is odd, then we have communication in the positive direction
  for(i=0;i<n;i++)
    wire[i]=dir[i];

  Float dtime;

  //If wire[i] is odd, then we have parallel transport in the
  //positive direction.  In this case, the matrix multiplication is
  //done before the field is transferred over to the adjacent node
  //
  //If we have transfer in the negative T direction (wire[i] = 6) then
  //we have to copy the appropriate fields into the send buffer
  for(i=0;i<n;i++)
    {
      if(!local[wire[i]/2])
      	{
	  if(wire[i]%2)
	    {
	      //printf("dir = %d, pre-mulitply \n", wire[i]);
#ifdef PROFILE
  dtime  = - dclock();
#endif

  if(conjugated)
    {
      partrans_cmv(non_local_chi_cb[wire[i]]/2,uc_nl_cb_pre[parity][wire[i]/2],gauge,vin[i],snd_buf_cb[wire[i]/2]);
    }
  else
    partrans_cmv_dag(non_local_chi_cb[wire[i]]/2,uc_nl_cb_pre[parity][wire[i]/2],gauge,vin[i],snd_buf_cb[wire[i]/2]);
  
#ifdef PROFILE
  dtime +=dclock();
  print_flops("",fname,66*non_local_chi_cb[wire[i]],dtime);
#endif
	    }
	  else if(wire[i] == 6)
	    {
#ifdef PROFILE
  dtime  = - dclock();
#endif
#if 1
            pt_copy_buffer(non_local_chi_cb[6],(long)vin[i],(long)snd_buf_t_cb,(long)Toffset[parity]);
#else
	      for(int j = 0; j < non_local_chi_cb[6];j++)
		  for(k = 0; k < VECT_LEN;k++)
		    *(snd_buf_t_cb+j*VECT_LEN+k) = *(vin[i] + *(Toffset[parity]+j)+ k);
#endif
#ifdef PROFILE
  dtime +=dclock();
  print_flops(fname,"pt_copy_buffer()",0,dtime);
#endif
	    }
	}
    }

#ifdef PROFILE
  dtime  = - dclock();
#endif
  int comms = 0;
  for(i=0;i<n;i++)
    {
      if(!local[wire[i]/2])
	{
	  //Calculate the starting address for the data to be sent
	  IFloat *addr = vin[i] + VECT_LEN * offset_cb[wire[i]];
	  //This points to the appropriate SCUDirArg for receiving
	  SCUarg_p[2*comms] = SCUarg_cb[2*wire[i]];
	  //This points to the appropriate SCUDirArg for sending
	  SCUarg_p[2*comms+1] = SCUarg_cb[2*wire[i]+1];
	  
	  //Set the send address
	  if(wire[i]%2)
	    SCUarg_p[2*comms+1]->Addr((void *)snd_buf_cb[wire[i]/2]);
	  else if(wire[i] == 6)
	    SCUarg_p[2*comms+1]->Addr((void *)snd_buf_t_cb);
	  else
	    SCUarg_p[2*comms+1]->Addr((void *)addr);
	  comms++;
	}
    }
#ifdef PROFILE
  dtime +=dclock();
  print_flops(fname,"Addr",0,dtime);
#endif

#ifdef PROFILE
  dtime  = - dclock();
#endif
  if(comms){
    SCUmulti.Init(SCUarg_p,2*comms);
//Begin transmission
    SCUmulti.SlowStartTrans();
  }
#ifdef PROFILE
  dtime +=dclock();
  print_flops(fname,"StartTrans()",0,dtime);
#endif


  //Do local calculations
  for(i=0;i<n;i++)
    {
#ifdef PROFILE
  dtime  = - dclock();
#endif

      if((wire[i]%2 && conjugated) || ((wire[i]%2 == 0) && (conjugated == 0)))
	{
	partrans_cmv_pad(local_chi_cb[wire[i]]/2,uc_l_pad_cb[parity][wire[i]],gauge,vin[i],vout);
	}
      else
	partrans_cmv_dag_pad(local_chi_cb[wire[i]]/2,uc_l_pad_cb[parity][wire[i]],gauge,vin[i],vout);

#ifdef PROFILE
  dtime +=dclock();
  print_flops("",fname,66*local_chi_cb[wire[i]],dtime);
#endif
    }
#ifdef PROFILE
  dtime  = - dclock();
#endif
  if(comms){ SCUmulti.TransComplete(); }
#ifdef PROFILE
  dtime +=dclock();
  print_flops(fname,"TransComplete()",0,dtime);
#endif

  //If wire[i] is even, then we have transport in the negative direction.
  //In this case, the vector field is multiplied by the SU(3) link matrix
  //after all communication is complete
//  IFloat *fp0,*fp1;
  for(i=0;i<n;i++)
    {
      if(!local[wire[i]/2])
	{
	  if(!(wire[i]%2))
	    {
#ifdef PROFILE
  dtime  = - dclock();
#endif
 
  if(conjugated)
    partrans_cmv_dag_pad(non_local_chi_cb[wire[i]]/2,uc_nl_pad_cb[parity][wire[i]],gauge,rcv_buf[wire[i]],vout);
  else
    partrans_cmv_pad(non_local_chi_cb[wire[i]]/2,uc_nl_pad_cb[parity][wire[i]],gauge,rcv_buf[wire[i]],vout);

#ifdef PROFILE
  dtime +=dclock();
  print_flops("",fname,66*non_local_chi_cb[wire[i]],dtime);
#endif
	    }
	  //Otherwise we have parallel transport in the positive direction.
	  //In this case, the received data has already been pre-multiplied
	  //All we need to do is to put the transported field in the correct place
	  else
	    {
#ifdef PROFILE
  dtime  = - dclock();
#endif
#if 1
              pt_copy_pad(non_local_chi_cb[wire[i]]/2,uc_nl_pad_cb[parity][wire[i]],rcv_buf[wire[i]],vout);
#else
	      //Place the data in the receive buffer into the result vector
	      for(int s=0;s<non_local_chi_cb[wire[i]];s++)
		{
		  fp0 = (IFloat *)((long)rcv_buf[wire[i]]+uc_nl_pad_cb[parity][wire[i]][s].src);
		  fp1 = (IFloat *)((long)vout +  uc_nl_pad_cb[parity][wire[i]][s].dest);
		  for(int d = 0;d<VECT_LEN;d++)
		    *(fp1+d) = *(fp0+d);
		  //moveMem(fp1,fp0,VECT_LEN*sizeof(IFloat));
		}
#endif
#ifdef PROFILE
  dtime +=dclock();
  print_flops(fname,"pt_copy_pad()",0,dtime);
#endif
	    }
	}
    }
//  ParTrans::PTflops +=33*n*vol;
}

#endif // #ifndef SCIDAC

#undef PROFILE
//Parallel transport of a vector through one hop
void PT::vec(int n, IFloat **vout, IFloat **vin, const int *dir){
  int i;
  static int call_num=0;
  SCUDirArgIR *SCUarg_p[2*n];
  call_num++;
  //for(int s = 0; s < GJP.VolNodeSites(); s++)
  //  {
  //    for(int t = 0; t < 4; t++)
  //	{
  //	  printf("site = %d, direction = %d\n",s,t);
  //	  for(int u = 0; u < 9; u++)
  //	    printf("%e %e\n",*(gauge_field_addr+4*GAUGE_LEN*s + GAUGE_LEN*t + 2*u),*(gauge_field_addr+4*GAUGE_LEN*s + GAUGE_LEN*t + 2*u+1));
  //	}
  //  }

#ifdef PROFILE
  Float dtime  = - dclock();
#endif
  int wire[n];
  SCUDirArgMulti SCUmulti;

  char *fname="pt_1vec";
//  VRB.Func("",fname);
	
  int non_local_dir=0;
  for(i=0;i<n;i++) wire[i] = dir[i]; // from (x,y,z,t) to (t,x,y,z)
//  for(i=0;i<n;i++) printf("wire[%d]=%d\n",i,dir[i]);
  for(i=0;i<n;i++)
  if (!local[wire[i]/2]){
    IFloat * addr = (vin[i]+VECT_LEN*offset[wire[i]]);
    SCUarg_p[2*non_local_dir] = SCUarg[0][2*wire[i]];
    SCUarg_p[2*non_local_dir+1] = SCUarg[0][2*wire[i]+1];
    SCUarg_p[2*non_local_dir+1]->Addr((void *)addr);
    non_local_dir++;
  }
  if(non_local_dir){
    SCUmulti.Init(SCUarg_p,non_local_dir*2);
    SCUmulti.SlowStartTrans();
  }
	
  for(i=0;i<n;i++) 
    partrans_cmv_agg(local_chi[wire[i]],(long)uc_l[wire[i]], (long)vin[i],(long)vout[i]);
	
  if(non_local_dir){ SCUmulti.TransComplete(); }

  for(i=0;i<n;i++) 
    partrans_cmv_agg(non_local_chi[wire[i]],(long)uc_nl[wire[i]], (long)rcv_buf[wire[i]],(long)vout[i]);

#ifdef PROFILE
  dtime +=dclock();
  print_flops("",fname,66*n*vol,dtime);
#endif
  Flops +=66*n*vol;
}
