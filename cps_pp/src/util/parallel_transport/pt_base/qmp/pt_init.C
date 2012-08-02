#ifdef USE_QMP
/*! \file
  \brief  Definition of parallel transport definitions for QCDOC.
  
  $Id: pt_init.C,v 1.6 2012-08-02 21:20:01 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2012-08-02 21:20:01 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_base/qmp/pt_init.C,v 1.6 2012-08-02 21:20:01 chulwoo Exp $
//  $Id: pt_init.C,v 1.6 2012-08-02 21:20:01 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: pt_init.C,v $
//  $Revision: 1.6 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_base/qmp/pt_init.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
#include <strings.h>
#include <string.h>
#include "asq_data_types.h"
#include "pt_int.h"

int PT::size[NDIM];
int PT::vol;
int PT::evenodd;

//Initialization of Parallel Transport class
void PT::init(PTArg *pt_arg)
{
  char *cname = "";
  char *fname = "pt_init()";
//printf("%s\n",fname);
  int i, j, x[NDIM],nei[NDIM];
  int local_count[2*NDIM];
  int non_local_count[2*NDIM];
  int vlen = VECT_LEN*sizeof(IFloat); //size of incoming vector
  int vlen2 =VECT_LEN_OUT*sizeof(IFloat); //size of outgoing vector (maybe different to optimize for QCDOC PEC)

  //---------------------------------------------------------------------------
  int local_count_cb[2][2*NDIM];
  int non_local_count_cb[2][2*NDIM];
  //---------------------------------------------------------------------------

  size[0] = pt_arg->size[0];
  size[1] = pt_arg->size[1];
  size[2] = pt_arg->size[2];
  size[3] = pt_arg->size[3];
  local[0] = pt_arg->local[0];
  local[1] = pt_arg->local[1];
  local[2] = pt_arg->local[2];
  local[3] = pt_arg->local[3];
  non_local_dirs = 2*(4-local[0]-local[1]-local[2]-local[3]);
#ifdef UNIFORM_SEED_NO_COMMS
  if( non_local_dirs>0){
    fprintf(stderr,"PT::non_local_dirs=%d\n",non_local_dirs);
    exit(-33);
  }
#endif
  //printf("Local directions = %d %d %d %d\n", local[0],local[1],local[2],local[3]);

  gauge_field_addr = pt_arg->gauge_field_addr;
  g_str_ord = pt_arg->g_str_ord;
  g_conj = pt_arg->g_conj;
  v_str_ord = pt_arg->v_str_ord;
  v_str_ord_cb = pt_arg->v_str_ord_cb;
  evenodd = pt_arg->evenodd;
  prec = pt_arg->prec;

//  printf("g_str_ord=%d v_str_ord=%d v_str_ord_cb=%d g_conj=%d\n", g_str_ord, v_str_ord, v_str_ord_cb,g_conj);

  switch(g_str_ord){
    case PT_XYZT:
      LexGauge = lex_g_xyzt;
      LexGauge2 = lex_g_txyz_cb;
      break;
    case PT_XYZT_CB_O:
      LexGauge = lex_g_xyzt_cb_o;
      LexGauge2 = lex_g_txyz_cb;
      break;
    case PT_XYZT_CB_E:
      LexGauge = lex_g_xyzt_cb_e;
      LexGauge2 = lex_g_txyz_cb;
      break;
    default:
      fprintf(stderr,"PT::init got invalid g_str_ord\n");
      break;
  }

  switch(v_str_ord){
    case PT_XYZT:
      LexVector = lex_xyzt;
      break;
    case PT_XYZT_CB_O:
      LexVector = lex_xyzt_cb_o;
      break;
    default:
      fprintf(stderr,"PT::init got invalid v_str_ord\n");
      break;
  }

  switch(v_str_ord_cb){
    case PT_TXYZ:
      LexVector_cb = lex_txyz_cb;
      break;
    default:
      fprintf(stderr,"PT::init got invalid v_str_ord_cb\n");
      break;
  }
  
  if (g_conj) { Copy = &dag_cpy; DagCopy = &cpy; conjugated = PT_DAG_YES;}
  else        { Copy = &cpy; DagCopy = &dag_cpy; conjugated = PT_DAG_NO;}


  //For the fastest changing index, data must be sent in many short messages
  //For the slowest changing index, the boundaries of the hypersurface are 
  //stored together in large blocks, so a few long messages can be sent.

  blklen[0] = blklen[1]= vlen;
  for(i=1;i<NDIM;i++) {blklen[2*i+1] = blklen[2*i] = blklen[2*i-1]*size[i-1]; }

  numblk[2*NDIM-1]=numblk[2*NDIM-2]=1;
  for(i=NDIM-2;i>=0;i--) {numblk[i*2+1] = numblk[2*i] = numblk[2*i+2]*size[i+1]; }

  //The stride length is longer when the blocks are large
  for(i=0;i<NDIM*2;i++)  stride[i] = blklen[i]* (size[i/2]-1);

  //Calculate the local volume
  vol = 1;
  for(i=0; i<NDIM;i++) vol *= size[i];

  //Calculate the number of local and non-local parallel transports
  //are needed in each direction
  for(i=0; i<NDIM;i++) {
    if (local[i])
      non_local_chi[2*i+1] = non_local_chi[2*i] = 0;
    else 
      non_local_chi[2*i+1] = non_local_chi[2*i] = vol/size[i];
    local_chi[2*i+1] = local_chi[2*i] = vol - non_local_chi[2*i];
//printf("local_chi[%d]=%d non_local_chi[%d]=%d\n", i*2,local_chi[i*2],i*2,non_local_chi[i*2]);
//printf("local_chi[%d]=%d non_local_chi[%d]=%d\n", i*2+1,local_chi[i*2+1],i*2+1,non_local_chi[i*2+1]);
  }

  //---------------------------------------------------------------------------
  //Calculation of block length, number of blocks, and stride for 
  //checkerboarded storage

  //Block length for checkerboarded scheme is similar to canonical scheme
  //However, there must be a few modifications.
  //
  //For the T (3rd) direction, a block strided communciation will not work
  //As a result, pointers to the appropriate fields will be aggregated
  //into a send buffer before being sent as a single block.
  //
  //For transfers in all other directions, a block-strided move is allowed.
  //In these cases, the fastest changing index (X) requires many short messages
  //while the slowest changing index (Z) can be transfered in one block
  
  blklen_cb[6] = blklen_cb[7] = vol*vlen/(2*size[3]);
  numblk_cb[6] = numblk_cb[7] = 1;

  blklen_cb[4] = blklen_cb[5] = vol*vlen/(2*size[2]);
  numblk_cb[4] = numblk_cb[5] = 1;

  blklen_cb[2] = blklen_cb[3] = vol*vlen/(2*size[1]*size[2]);
  numblk_cb[2] = blklen_cb[3] = size[2];

  blklen_cb[0] = blklen_cb[1] = vol*vlen/(2*size[0]*size[1]*size[2]);
  numblk_cb[0] = numblk_cb[1] = size[1]*size[2];

  //The stride is also similar
  stride_cb[0] = stride_cb[1] = 0;
  for(i = 0; i < NDIM*2; i++)
      stride_cb[i] = blklen_cb[i] * (size[i/2]-1);
  stride_cb[6] = stride_cb[7] = 0;

  //Calculate the number of local and non-local parallel transports
  for(i = 0; i < NDIM; i++)
    {
      if(!local[i])
	non_local_chi_cb[2*i+1] = non_local_chi_cb[2*i] = vol/(2*size[i]);
      else
	non_local_chi_cb[2*i+1] = non_local_chi_cb[2*i] = 0;
      local_chi_cb[2*i+1] = local_chi_cb[2*i] = vol/2 - non_local_chi_cb[2*i];
  //    printf("non_local_chi_cb[2*%d] = %d; local_chi_cb[2*%d] = %d\n",i,non_local_chi_cb[2*i],i,local_chi_cb[2*i]);
    }

  //---------------------------------------------------------------------------

  for(i=0; i<2*NDIM;i++){
    local_count[i]=non_local_count[i]=0;
    //-------------------------------------------------------------------------
    for(int parity = 0; parity < 2; parity++)
      local_count_cb[parity][i] = non_local_count_cb[parity][i] = 0;
    //------------------------------------------------------------------------


      uc_l[i] = (gauge_agg *)Alloc(cname,fname,"uc_l[i]",
				     sizeof(gauge_agg)*(1+local_chi[i]));
      uc_nl[i] = (gauge_agg *)Alloc(cname,fname,"uc_nl[i]",
				      sizeof(gauge_agg)*(1+non_local_chi[i]),0,non_local_chi[i]);

    //-------------------------------------------------------------------------
    //Allocate memory for gauge_agg_cb
    for(int parity = 0; parity < 2; parity++)
    {
      uc_l_cb[parity][i] = (gauge_agg_cb *)FastAlloc(cname,fname,"uc_l_cb",sizeof(gauge_agg_cb)*(1+local_chi_cb[i]),local_chi_cb[i]);
      uc_nl_cb[parity][i] = (gauge_agg_cb *)FastAlloc(cname,fname,"uc_nl_cb",sizeof(gauge_agg_cb)*(1+non_local_chi_cb[i]),non_local_chi_cb[i]);

      uc_l_pad_cb[parity][i] = (gauge_agg_cb *)(unsigned long)FastAlloc(cname,fname,"uc_l_pad_cb",sizeof(gauge_agg_cb)*(1+local_chi_cb[i]),local_chi_cb[i]);
      uc_nl_pad_cb[parity][i] = (gauge_agg_cb *)(unsigned long)FastAlloc(cname,fname,"uc_nl_pad_cb",sizeof(gauge_agg_cb)*(1+non_local_chi_cb[i]),non_local_chi_cb[i]);
//      printf("uc_pad_cb = %p %p\n",uc_l_pad_cb[parity][i],uc_nl_pad_cb[parity][i]);
    }

    //-------------------------------------------------------------------------

    // This buffer is actually overkill, but ensures will work if
    // shift_field is called with hop>1
    if(non_local_chi[i]>0){
      rcv_buf[i] = (IFloat *)FastAlloc(3*MAX_HOP*non_local_chi[i]*vlen);
    if(!rcv_buf[i])PointerErr("",fname,"rcv_buf[i]");

    //Used buffer used in vvpd
      rcv_buf2[i] = (IFloat *)FastAlloc(MAX_HOP*non_local_chi[i]*vlen);
      if(!rcv_buf2[i])PointerErr("",fname,"rcv_buf2[i]");
    } else{
      rcv_buf[i] = rcv_buf2[i] = NULL;
    }
  }

  //---------------------------------------------------------------------------
  //Allocate memory for send buffer

  for(i=0; i<NDIM;i++)
      snd_buf_cb[i] = (IFloat *)FastAlloc(cname,fname,"snd_buf_cb[i]",3*non_local_chi_cb[2*i+1]*vlen,non_local_chi_cb[2*i+1]);
    snd_buf_t_cb = (IFloat *)FastAlloc("",fname,"snd_buf_t_cb",3*non_local_chi_cb[6]*vlen,non_local_chi_cb[6]);

  for(i = 0; i < 2;i++)
    Toffset[i] = (int *)FastAlloc(cname,fname,"Toffset[parity]",non_local_chi_cb[6]*sizeof(int),non_local_chi_cb[6]);

  //Allocate memory for the gauge_agg_cb used for matrix pre-multiplication

  for(i = 0; i< NDIM;i++)
    for(int parity = 0; parity<2;parity++)
      uc_nl_cb_pre[parity][i] = (gauge_agg_cb *)FastAlloc(cname,fname,"uc_nl_cb_pre[parity][i]",sizeof(gauge_agg_cb)*(1+non_local_chi_cb[2*i+1]),non_local_chi_cb[2*i+1]);

  int parity = 0;
  //---------------------------------------------------------------------------
 
  //Calculate source and destination indices for gauge aggregates (see set_hop_pointer())
  //Only for one hop
   for(x[3]=0,nei[3]=0;x[3]<size[3];x[3]++,nei[3]++)
    for(x[2]=0,nei[2]=0;x[2]<size[2];x[2]++,nei[2]++)
      for(x[1]=0,nei[1]=0;x[1]<size[1];x[1]++,nei[1]++)
	for(x[0]=0,nei[0]=0;x[0]<size[0];x[0]++,nei[0]++){
	  for(i=0;i<NDIM;i++){
	    
//	    printf("%d %d %d %d %d\n",x[0],x[1],x[2],x[3],i);
	    // positive direction
	    //This is for transport of a vector in the negative direction
	    //An even index for uc_nl, uc_l, uc_nl_cb, uc_l_cb corresponds
	    //to parallel transport in the negative direction
	    
	    if( (x[i] == 0) && (!local[i])){
	      nei[i] = size[i]-1;
	      (uc_nl[2*i]+non_local_count[2*i])->src = non_local_count[2*i]*vlen;
	      (uc_nl[2*i]+non_local_count[2*i])->dest = LexVector(nei)*vlen2;
	      non_local_count[i*2]++;
	      if (non_local_count[i*2]>non_local_chi[i*2])
		fprintf(stderr,"%s:non_local_count[%d](%d)>non_local_chi[%d](%d)\n",
			fname,2*i,non_local_count[2*i],2*i,non_local_chi[2*i]);
	    } 
	    else {
	      nei[i] = (x[i]-1+size[i])%size[i];
	      if(local_count[2*i]<0) fprintf(stderr,"%s:local_count[%d]=%d]n",
			fname,2*i,local_count[2*i]);
	      (uc_l[2*i]+local_count[2*i])->src = LexVector(x)*vlen;
	      (uc_l[2*i]+local_count[2*i])->dest = LexVector(nei)*vlen2;
	      local_count[i*2]++;
	      if (local_count[i*2]>local_chi[i*2])
		fprintf(stderr,"%s:local_count[%d](%d)>local_chi[%d](%d)\n",
			fname,2*i,local_count[2*i],2*i,local_chi[2*i]);
	    }
	    // negative direction
	    //This is parallel transport in the positive direction
	    //An odd index for uc_l, uc_nl, uc_l_cb,uc_nl_cb corresponds to
	    //transport in the positive direction
	    if((x[i] == (size[i] -1))  && (!local[i])){
	      nei[i] = 0;
	      (uc_nl[2*i+1]+non_local_count[2*i+1])->src = non_local_count[2*i+1]*vlen;
	      (uc_nl[2*i+1]+non_local_count[2*i+1])->dest = LexVector(nei)*vlen2;
	      non_local_count[i*2+1]++;
	      if (non_local_count[i*2+1]>non_local_chi[i*2+1])
		fprintf(stderr,"%s:non_local_count[%d](%d)>non_local_chi[%d](%d)\n",
			fname,2*i+1,non_local_count[2*i+1],2*i+1,non_local_chi[2*i+1]);
	    } else {
	      nei[i] = (x[i]+1)%size[i];
	      if(local_count[2*i+1]<0) fprintf(stderr,"%s:local_count[%d]=%d]n",
					       fname,2*i+local_count[2*i+1]);
	      (uc_l[2*i+1]+local_count[2*i+1])->src = LexVector(x)*vlen;
	      (uc_l[2*i+1]+local_count[2*i+1])->dest = LexVector(nei)*vlen2;
	      local_count[i*2+1]++;
	      if (local_count[i*2+1]>local_chi[i*2+1])
		fprintf(stderr,"%s:local_count[%d](%d)>local_chi[%d](%d)\n",
			fname,2*i+1,local_count[2*i+1],2*i+1,local_chi[2*i+1]);
	    }
	    nei[i] = x[i];
	  }
	}

   //--------------------------------------------------------------------------
  //Calculate source and destination indices for gauge aggregates
  //Only for one hop
   for(x[2]=0,nei[2]=0;x[2]<size[2];x[2]++,nei[2]++)
     for(x[1]=0,nei[1]=0;x[1]<size[1];x[1]++,nei[1]++)
       for(x[0]=0,nei[0]=0;x[0]<size[0];x[0]++,nei[0]++)
	 for(x[3]=0,nei[3]=0;x[3]<size[3];x[3]++,nei[3]++)
	   {

	    parity = (x[0]+x[1]+x[2]+x[3])%2;
	    //Calculate offsets for transfers in the negative T direction
	    if((x[3] == 0) && !local[3])
	      *(Toffset[parity] +non_local_count_cb[parity][6]) = LexVector_cb(x)*VECT_LEN;

	    for(i=0;i<NDIM;i++){

	    // positive direction
	    //This is for transport of a vector in the negative direction
	    //An even index for uc_nl, uc_l, uc_nl_cb, uc_l_cb corresponds
	    //to parallel transport in the negative direction

	      if((x[i] == 0) && !local[i])
	      {
	      nei[i] = size[i]-1;

	      //The src and dest indexes index the Vector, and do not include 
	      //information for the size of the vector, nor
	      //the size of the IFloat object.  This is to allow runt-time 
	      //adjustment of these parameters
	      //
	      //src - Source index in the receive buffer coming from a 
	      //      positive adjacent node
	      //dest - Destination index for the transported vector field
	      //dest2 - Destination index for the padded vector field
	      //gauge_index - Index of the SU(3) gauge link needed to 
	      //              transport the field
	      //dagger - determines if the gauge link needs to be conjugated

	      (uc_nl_cb[parity][2*i]+non_local_count_cb[parity][2*i])->src = non_local_count_cb[parity][2*i]*6*sizeof(IFloat);
	      (uc_nl_cb[parity][2*i]+non_local_count_cb[parity][2*i])->dest = LexVector_cb(nei)*6*sizeof(IFloat);
	      (uc_nl_cb[parity][2*i]+non_local_count_cb[parity][2*i])->gauge = LexGauge2(nei,i)*GAUGE_LEN*sizeof(IFloat);

	      (uc_nl_pad_cb[parity][2*i]+non_local_count_cb[parity][2*i])->src = non_local_count_cb[parity][2*i]*6*sizeof(IFloat);
	      (uc_nl_pad_cb[parity][2*i]+non_local_count_cb[parity][2*i])->dest = (LexVector_cb(nei)*8+2*i)*8*sizeof(IFloat);
	      (uc_nl_pad_cb[parity][2*i]+non_local_count_cb[parity][2*i])->gauge = LexGauge2(nei,i)*GAUGE_LEN*sizeof(IFloat);


	      non_local_count_cb[parity][i*2]++;
	      if(non_local_count_cb[parity][i*2]>non_local_chi_cb[i*2])
		fprintf(stderr,
			"%s:non_local_count_cb[%d][%d](%d)>non_local_chi_cb[%d](%d)\n",
			fname,parity,2*i,non_local_count_cb[parity][2*i],2*i,non_local_chi[2*i]);
	    } 
	    else 
	      {
	      nei[i] = (x[i]+size[i]-1)%size[i];

	      (uc_l_cb[parity][2*i]+local_count_cb[parity][2*i])->src = LexVector_cb(x)*6*sizeof(IFloat);
	      (uc_l_cb[parity][2*i]+local_count_cb[parity][2*i])->dest = LexVector_cb(nei)*6*sizeof(IFloat);
	      (uc_l_cb[parity][2*i]+local_count_cb[parity][2*i])->gauge = LexGauge2(nei,i)*GAUGE_LEN*sizeof(IFloat);

	      (uc_l_pad_cb[parity][2*i]+local_count_cb[parity][2*i])->src = LexVector_cb(x)*6*sizeof(IFloat);
	      (uc_l_pad_cb[parity][2*i]+local_count_cb[parity][2*i])->dest = (LexVector_cb(nei)*8+2*i)*8*sizeof(IFloat);
	      (uc_l_pad_cb[parity][2*i]+local_count_cb[parity][2*i])->gauge = LexGauge2(nei,i)*GAUGE_LEN*sizeof(IFloat);

	      local_count_cb[parity][i*2]++;
	      if(local_count_cb[parity][i*2]>local_chi_cb[i*2])
		fprintf(stderr,"%s:local_count_cb[%d][%d](%d)>local_chi_cb[%d](%d)\n",fname,parity,2*i,local_count_cb[parity][2*i],2*i,local_chi[2*i]);
	      }

	    // negative direction
	    //This is parallel transport in the positive direction
	    //An odd index for uc_l, uc_nl, uc_l_cb,uc_nl_cb corresponds to
	    //transport in the positive direction

	      if((x[i] == (size[i] -1)) && !local[i])
	      {
	      nei[i] = 0;

	      //src - Source index in the receive buffer
	      //dest - Destination index for the transported vector field
	      //In only this case, the field is transported pre-multiplied by
	      //the gauge link.  As a result, gauge_index and dagger are not 
	      //strictly necessary.
	      //
	      //However, we do need to specify another gauge aggregate that 
	      //will contain the information necessary
	      //for the pre-multiplication of the SU(3) link matrix
	      
	      (uc_nl_cb[parity][2*i+1]+non_local_count_cb[parity][2*i+1])->src = non_local_count_cb[parity][2*i+1]*6*sizeof(IFloat);
	      (uc_nl_cb[parity][2*i+1]+non_local_count_cb[parity][2*i+1])->dest = LexVector_cb(nei)*6*sizeof(IFloat);
	      (uc_nl_cb[parity][2*i+1]+non_local_count_cb[parity][2*i+1])->gauge = LexGauge2(x,i)*GAUGE_LEN*sizeof(IFloat);

	      
	      (uc_nl_pad_cb[parity][2*i+1]+non_local_count_cb[parity][2*i+1])->src = non_local_count_cb[parity][2*i+1]*6*sizeof(IFloat);
	      (uc_nl_pad_cb[parity][2*i+1]+non_local_count_cb[parity][2*i+1])->dest = (LexVector_cb(nei)*8+2*i+1)*8*sizeof(IFloat);
	      (uc_nl_pad_cb[parity][2*i+1]+non_local_count_cb[parity][2*i+1])->gauge = LexGauge2(x,i)*GAUGE_LEN*sizeof(IFloat);

	      (uc_nl_cb_pre[parity][i]+non_local_count_cb[parity][2*i+1])->src = LexVector_cb(x)*6*sizeof(IFloat);
	      (uc_nl_cb_pre[parity][i]+non_local_count_cb[parity][2*i+1])->dest = non_local_count_cb[parity][2*i+1]*6*sizeof(IFloat);
	      (uc_nl_cb_pre[parity][i]+non_local_count_cb[parity][2*i+1])->gauge = LexGauge2(x,i)*GAUGE_LEN*sizeof(IFloat);

	      non_local_count_cb[parity][i*2+1]++;
	      if(non_local_count_cb[parity][i*2+1]>non_local_chi_cb[i*2+1])
		fprintf(stderr,"%s:non_local_count_cb[%d][%d](%d)>non_local_chi_cb[%d](%d)\n",fname,parity,2*i+1,non_local_count_cb[parity][2*i+1],2*i+1,non_local_chi[2*i+1]);
	      } 
	    else 
	      {
	      nei[i] = (x[i]+1)%size[i];

	      (uc_l_cb[parity][2*i+1]+local_count_cb[parity][2*i+1])->src = LexVector_cb(x)*6*sizeof(IFloat);
	      (uc_l_cb[parity][2*i+1]+local_count_cb[parity][2*i+1])->dest = LexVector_cb(nei)*6*sizeof(IFloat);
	      (uc_l_cb[parity][2*i+1]+local_count_cb[parity][2*i+1])->gauge = LexGauge2(x,i)*GAUGE_LEN*sizeof(IFloat);
	      (uc_l_pad_cb[parity][2*i+1]+local_count_cb[parity][2*i+1])->src = LexVector_cb(x)*6*sizeof(IFloat);
	      (uc_l_pad_cb[parity][2*i+1]+local_count_cb[parity][2*i+1])->dest = (LexVector_cb(nei)*8+2*i+1)*8*sizeof(IFloat);
	      (uc_l_pad_cb[parity][2*i+1]+local_count_cb[parity][2*i+1])->gauge = LexGauge2(x,i)*GAUGE_LEN*sizeof(IFloat);

	      local_count_cb[parity][i*2+1]++;
	      if(local_count_cb[parity][i*2+1]>local_chi_cb[i*2+1])
		fprintf(stderr,"%s:local_count_cb[%d][%d](%d)>local_chi_cb[%d](%d)\n",fname,parity,2*i+1,local_count_cb[parity][2*i+1],2*i+1,local_chi[2*i+1]);
	    }
	    nei[i] = x[i];
	  }
	}
   //--------------------------------------------------------------------------

  //Sets bits in uc_l and uc_nl to zero
  for(i=0;i<NDIM*2;i++){
    gauge_agg *tmp = uc_l[i]+local_chi[i];
    memcpy(tmp,tmp-1,sizeof(gauge_agg));
	if(non_local_chi[i]){
      tmp = uc_nl[i]+non_local_chi[i];
      memcpy(tmp,tmp-1,sizeof(gauge_agg));
	}
  }

  //Calculate offsets?
  //For even array index (transfer in the negative  direction) the
  //offset is 0
  //For odd indexes, the offsets are:
  //offset[1] = size[0]-1
  //offset[3] = size[0]*(size[1]-1)
  //offset[5] = size[0]*size[1]*(size[2]-1)
  //offset[7] = size[0]*size[1]*size[2]*(size[3]-1)
  //These offsets correspond to the starting index for data transfer
  //in the positive direction
  int temp=1;
  for(i=0;i<NDIM;i++){
    offset[2*i]  = 0;
    offset[2*i+1] = temp*(size[i]-1);
    temp *= size[i];
  }

  //-------------------------------------------------------------------
  temp = 1;
  /*
  for(i = 0;i<NDIM;i++)
    {
      offset_cb[2*i] = 0;
      offset_cb[2*i+1] = stride_cb[2*i+1]/sizeof(IFloat);
    }
  */
  //-------------------------------------------------------------------

  // Allocate memory for hop pointers
  for (j=0; j<MAX_HOP; j++) {
  
    for(i=0; i<2*NDIM; i++){
      
      //Calculate the number of local and non-local sites needed 
      //j+1 is the length of the hop
      //i indicates the communication direction
      int nl_size = (j+1)*non_local_chi[i] ;
      int l_size = vol - nl_size ;

      if (l_size>0){
        hp_l[j][i] = (hop_pointer*) FastAlloc(cname,fname,"hp_l[j][i]",(1+l_size)*sizeof(hop_pointer));
        src_l[j][i] = (unsigned long*)FastAlloc(cname,fname,"src_l[j][i]",(1+l_size)*sizeof(unsigned long));
        dest_l[j][i] = (unsigned long*)FastAlloc(cname,fname,"dest_l[j][i]",(1+l_size)*sizeof(unsigned long));
      } else {
        hp_l[j][i] = NULL;
        src_l[j][i] = NULL;
        dest_l[j][i] = NULL;
      }
      if (nl_size>0){
  	hp_nl[j][i] = (hop_pointer*) FastAlloc(cname,fname,"hp_nl[j][i]",(1+nl_size)*sizeof(hop_pointer));
  	src_nl[j][i] = (unsigned long*)FastAlloc(cname,fname,"src_nl[j][i]",(1+nl_size)*sizeof(unsigned long));
  	dest_nl[j][i] = (unsigned long*)FastAlloc(cname,fname,"dest_nl[j][i]",(1+nl_size)*sizeof(unsigned long));
      } else {
        hp_nl[j][i] = NULL;
        src_nl[j][i] = NULL;
        dest_nl[j][i] = NULL;
      }
    }
  }
  
  set_hop_pointer();

  //Calculate the indices for the source and destination
  for (j=0; j<MAX_HOP; j++) {
    for(i=0; i<2*NDIM; i++){
      int nl_size = (j+1)*non_local_chi[i];
      int l_size = vol - nl_size;
      if (l_size>0){
        for (int s=0; s<l_size; s++) {
          src_l[j][i][s] = hp_l[j][i][s].src/(VECT_LEN*sizeof(IFloat));
          dest_l[j][i][s] = hp_l[j][i][s].dest/(VECT_LEN_OUT*sizeof(IFloat));
        }
	src_l[j][i][l_size] = src_l[j][i][l_size-1];
	dest_l[j][i][l_size] = dest_l[j][i][l_size-1] ;
      }
      if (nl_size>0){
        for (int s=0; s<nl_size; s++) {
          src_nl[j][i][s] = hp_nl[j][i][s].src/(VECT_LEN*sizeof(IFloat));
          dest_nl[j][i][s] = hp_nl[j][i][s].dest/(VECT_LEN_OUT*sizeof(IFloat));
        }
	src_nl[j][i][nl_size] = src_nl[j][i][nl_size-1];
	dest_nl[j][i][nl_size] = dest_nl[j][i][nl_size-1] ;
      }
    }
  }
//printf("pt_init() done\n");
	
}

//Free memory associated with the parallel transport of the fermions
void PT::delete_buf(){
  char *fname = "pt_delete()";
	
  for(int i = 0; i < 2*NDIM; i++){
    Free(uc_l[i]);
    Free(uc_nl[i]);
    //--------------------------------------------------------------------

    for(int parity = 0; parity < 2; parity++)
      {
	Free(uc_l_cb[parity][i]);
	Free(uc_nl_cb[parity][i]);
	Free( (void *) ((unsigned long)uc_l_pad_cb[parity][i]) );
	Free( (void *) ((unsigned long)uc_nl_pad_cb[parity][i]) );
      }

    //-------------------------------------------------------------------
//    sfree(rcv_buf[i],cname,fname,"rcv_buf[i]");
//    sfree(rcv_buf2[i],cname,fname,"rcv_buf2[i]");
    if(non_local_chi[i] > 0)
      {
	Free(rcv_buf[i]);
	Free(rcv_buf2[i]);
      }
  }

  //-----------------------------------------------------------------------

  for(int i = 0; i < NDIM; i++)
    {
//      sfree(snd_buf_cb[i],cname,fname,"snd_buf_cb[i]");
      if(non_local_chi_cb[2*i+1]>0)
	Free(snd_buf_cb[i]);
      for(int parity = 0; parity < 2; parity++)
	Free(uc_nl_cb_pre[parity][i]);
    }
  if(non_local_chi_cb[6] > 0)
    {
      Free(snd_buf_t_cb);
      
      for(int i = 0; i < 2; i++)
	Free(Toffset[i]);    
    }

  //-----------------------------------------------------------------------

  for (int hop=0; hop<MAX_HOP; hop++) {
    for(int i = 0; i < 2*NDIM; i++){
      int nl_size = (hop+1)*non_local_chi[i];
      int l_size = vol - nl_size ;
      if (l_size>0){
        Free(hp_l[hop][i]);
        Free(src_l[hop][i]);
        Free(dest_l[hop][i]);
      }
      if (nl_size>0){
	  Free(hp_nl[hop][i]);
	  Free(src_nl[hop][i]);
	  Free(dest_nl[hop][i]);
      }
    }
  }
}

//CPS_END_NAMESPACE
#endif
