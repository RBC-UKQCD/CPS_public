/*! \file
  \brief  Definition of ParTransAsqtad class methods for QCDOC.
  
  $Id: pt.C,v 1.12 2004-08-09 07:47:25 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-08-09 07:47:25 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_asqtad/qcdoc/pt.C,v 1.12 2004-08-09 07:47:25 chulwoo Exp $
//  $Id: pt.C,v 1.12 2004-08-09 07:47:25 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: pt.C,v $
//  $Revision: 1.12 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_asqtad/qcdoc/pt.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
#if 0
#include <config.h>
#include <util/gjp.h>
#include <util/pt.h>
#include <util/time.h>
#include <sysfunc.h>
#include <comms/scu.h>
#include <stdio.h>
#include <qalloc.h>

#undef CPP
CPS_START_NAMESPACE

void dirac_cmv_jcw_agg_cpp( int sites, long chi, long u,long in, long out);
extern "C"{
  void cmm_agg(gauge_agg *chi, Matrix *phi,Matrix *result, int counter);
  void cmm_agg_cpp( int sites, long chi, long u,long in, long out);
  void cmv_agg_cpp( int sites, long u,long in, long out);
  void pt_asqtad_agg( int sites, long chi, long u,long in, long out);
  void copy_buffer(int n, long src, long dest, long ptable);
}

static const int NDIM=4;
static const int MAX_HOP=3;
static int size[NDIM];

enum {GAUGE_LEN=18,VECT_LEN=6, VECT_LEN2=6};
int vol = 1;

static gauge_agg *uc_l[2*SCUMachDim];
static gauge_agg *uc_nl[2*SCUMachDim];

static hop_pointer *hp_l[MAX_HOP][2*SCUMachDim];
static hop_pointer *hp_nl[MAX_HOP][2*SCUMachDim];

//static gauge_agg *uc_l_txyz[2][2*SCUMachDim];
//static gauge_agg *uc_nl_txyz[2][2*SCUMachDim];
//static int *Toffset[2][2];
static int blklen[2*SCUMachDim];
static int numblk[2*SCUMachDim];
static int stride[2*SCUMachDim];
static int local_chi[2*SCUMachDim];
static int non_local_chi[2*SCUMachDim];
static int offset[2*SCUMachDim];

static SCUDirArgIR *SCUarg[MAX_HOP][4*SCUMachDim];
static SCUDirArgIR *SCUarg_mat[MAX_HOP][4*SCUMachDim];

//static SCUDirArgIR *SCUarg_txyz[4*SCUMachDim];

static void Copy (IFloat *dest, IFloat *src){
  for(int i=0;i<18;i++)
    dest[i]=src[i];
}

static void DagCopy (IFloat *dest, IFloat *src){
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++){
      dest[2*(3*i+j)]=src[2*(3*j+i)];
      dest[2*(3*i+j)+1]=-src[2*(3*j+i)+1];
    }
}

static int LexVector(int *x){
  int result = x[0] + size[0]*(x[1]+size[1]*(x[2]+size[2]*x[3]));
  if (result<0 or result>=vol)
    ERR.General("","LexVector","index out of bounds %d %d %d %d\n",
		x[0],x[1],x[2],x[3]);
  return result;
}

static int LexVector_txyz(int *x){
  return  (x[3] + size[3]*(x[0]+size[0]*(x[1]+size[1]*x[2])))/2 ;
}

static int LexSurface(int *x){
  return  (x[0]+size[0]*(x[1]+size[1]*x[2]))/2 ;
}

static int LexGauge(int *x, int mu){
  int temp =  (x[0] + size[0]*(x[1]+size[1]*(x[2]+size[2]*x[3])));
  return (temp*NDIM + mu);
}

// Calculate the required offset given the direction and hop
int ParTransAsqtad::Offset(int dir, int hop) {

  // if positive direction then start at 0
  if (dir%2 == 0) return 0;

  int temp=1;
  int offset;
  for(int i=0;i<dir/2+1;i++){
    offset = temp*(size[i]-hop);
    temp *= size[i];
  }
  return offset;

}

void ParTransAsqtad::pt_init(const void *gauge_u)
{
  char *fname = "pt_init()";
  VRB.Func(cname,fname);
  //	IFloat *u_fl = (IFloat *)gauge_u;
  int i, j, x[NDIM],nei[NDIM];
  int local_count[2*NDIM];
  int non_local_count[2*NDIM];
  //	int local_count_txyz[2][2*NDIM];
  //	int non_local_count_txyz[2][2*NDIM];
  //	int Tcount[2][2];
  int vlen = VECT_LEN*sizeof(IFloat); //size of incoming vector
  int vlen2 =VECT_LEN2*sizeof(IFloat); //size of outgoing vector (maybe different to optimize for QCDOC PEC)

  size[0] = GJP.XnodeSites();
  size[1] = GJP.YnodeSites();
  size[2] = GJP.ZnodeSites();
  size[3] = GJP.TnodeSites();
  gauge_field_addr = (IFloat *) gauge_u;

  //	local_chi = ( int *)smalloc(sizeof(int)*2*NDIM);
  //	non_local_chi = ( int *)smalloc(sizeof(int)*2*NDIM);
  //	blklen = ( int *)smalloc(sizeof(int)*2*NDIM);
  //	numblk = ( int *)smalloc(sizeof(int)*2*NDIM);
  //	stride = ( int *)smalloc(sizeof(int)*2*NDIM);
  //	printf("stride=%p\n",stride);
  //	offset = ( int *)smalloc(sizeof(int)*2*NDIM);
  //	printf("offset=%p\n",offset);

  blklen[0] = blklen[1]= vlen;
  for(i=1;i<NDIM;i++) {blklen[2*i+1] = blklen[2*i] = blklen[2*i-1]*size[i-1]; }
  numblk[2*NDIM-1]=numblk[2*NDIM-2]=1;
  for(i=NDIM-2;i>=0;i--) {numblk[i*2+1] = numblk[2*i] = numblk[2*i+2]*size[i+1]; }
  for(i=0;i<NDIM*2;i++)  stride[i] = blklen[i]* (size[i/2]-1);
  vol = 1;
  for(i=0; i<NDIM;i++) vol *= size[i];
  //	printf("vol=%d\n",vol);
  for(i=0; i<NDIM;i++){
    non_local_chi[2*i+1] = non_local_chi[2*i] = vol/size[i];
    local_chi[2*i+1] = local_chi[2*i] = vol - non_local_chi[2*i];
    //		printf("local_chi[%d]=%d non_local_chi[%d]=%d\n",2*i,local_chi[2*i],2*i,non_local_chi[2*i]);
  }
  for(i=0; i<2*NDIM;i++){
    local_count[i]=non_local_count[i]=0;
    if(vol> 4096){
      uc_l[i]=uc_nl[i]=NULL;
    } else {	
      uc_l[i] = (gauge_agg *)qalloc(QFAST,sizeof(gauge_agg)*(1+local_chi[i]));
      uc_nl[i] = (gauge_agg *)qalloc(QFAST,sizeof(gauge_agg)*(1+non_local_chi[i]));
    }
    if(uc_l[i]==NULL) {
      uc_l[i] = (gauge_agg *)smalloc(sizeof(gauge_agg)*(1+local_chi[i]));
      printf("%s:%s: uc_l[%d] = %p\n",cname,fname,i,uc_l[i]);
    }
    if(uc_nl[i]==NULL) {
      uc_nl[i] = (gauge_agg *)smalloc(sizeof(gauge_agg)*(1+non_local_chi[i]));
      printf("%s:%s: uc_nl[%d] = %p\n",cname,fname,i,uc_nl[i]);
    }
    if(uc_l[i]==NULL)ERR.Pointer(cname,fname,"uc_l[i]");
    if(uc_nl[i]==NULL)ERR.Pointer(cname,fname,"uc_nl[i]");
#if 0
    for(int j = 0;j<2;j++){
      uc_l_txyz[j][i] = (gauge_agg *)smalloc(sizeof(gauge_agg)*local_chi[i]/2);
      uc_nl_txyz[j][i] = (gauge_agg *)smalloc(sizeof(gauge_agg)*non_local_chi[i]/2);
      local_count_txyz[j][i]=non_local_count_txyz[j][i]=0;
    }
#endif
#if 0
    rcv_buf[i] = (IFloat *)qalloc(QFAST,non_local_chi[i]*vlen*3);
    tmp_buf[i] = (IFloat *)qalloc(QFAST,vol*vlen2*3);
    if(rcv_buf[i]==NULL)
      rcv_buf[i] = (IFloat *)qalloc(QCOMMS,non_local_chi[i]*vlen*3);
    if(tmp_buf[i]==NULL)
      tmp_buf[i] = (IFloat *)qalloc(QCOMMS,vol*vlen2*3);
#endif
    rcv_buf[i] = (IFloat *)smalloc(MAX_HOP*3*non_local_chi[i]*vlen);
    //		tmp_buf[i] = (IFloat *)smalloc(vol*vlen2*3);
    if(rcv_buf[i]==NULL)ERR.Pointer(cname,fname,"rcv_buf[i]");
    //		if(tmp_buf[i]==NULL)ERR.Pointer(cname,fname,"tmp_buf[i]");
    //		printf("rcv_buf[%d]=%p\n",i,rcv_buf[i]);
    //		printf("tmp_buf[%d]=%p\n",i,tmp_buf[i]);
  }

  // Allocate memory for hop pointers
  for (j=0; j<MAX_HOP; j++) {
  
    for(i=0; i<2*NDIM; i++){
      
      int nl_size = (j+1)*non_local_chi[i];
      int l_size = vol - nl_size;
      
      hp_l[j][i] = (hop_pointer*) smalloc(l_size*sizeof(hop_pointer));
      VRB.Smalloc(cname,fname,"hp_l[j][i]",hp_l[j][i],l_size*sizeof(hop_pointer));
      if (hp_l[j][i] == NULL) ERR.Pointer(cname,fname,"hp_l[j][i]");

      hp_nl[j][i] = (hop_pointer*) smalloc(nl_size*sizeof(hop_pointer));
      VRB.Smalloc(cname,fname,"hp_nl[j][i]",hp_nl[j][i],nl_size*sizeof(hop_pointer));
      if (hp_nl[j][i] == NULL) ERR.Pointer(cname,fname,"hp_nl[j][i]");
    }
  }

  setHopPointer();

#if 0
  Toffset[0][0] = (int *)smalloc(sizeof(int)*size[0]*size[1]*size[2]/2);
  Toffset[1][0] = (int *)smalloc(sizeof(int)*size[0]*size[1]*size[2]/2);
  Toffset[0][1] = (int *)smalloc(sizeof(int)*size[0]*size[1]*size[2]/2);
  Toffset[1][1] = (int *)smalloc(sizeof(int)*size[0]*size[1]*size[2]/2);
#endif


  for(x[3]=0,nei[3]=0;x[3]<size[3];x[3]++,nei[3]++)
    for(x[2]=0,nei[2]=0;x[2]<size[2];x[2]++,nei[2]++)
      for(x[1]=0,nei[1]=0;x[1]<size[1];x[1]++,nei[1]++)
	for(x[0]=0,nei[0]=0;x[0]<size[0];x[0]++,nei[0]++){
	  for(i=0;i<NDIM;i++){
	    // positive direction
	    if(x[i] == 0){
	      nei[i] = size[i]-1;  
	      (uc_nl[2*i]+non_local_count[2*i])->src = non_local_count[2*i]*vlen;
	      (uc_nl[2*i]+non_local_count[2*i])->dest = LexVector(nei)*vlen2;
	      non_local_count[i*2]++;
	      if (non_local_count[i*2]>non_local_chi[i*2])
		ERR.General(cname,fname,"non_local_count[%d](%d)>non_local_chi[%d](%d)\n",2*i,non_local_count[2*i],2*i,non_local_chi[2*i]);
	    } else {
	      nei[i] = x[i]-1;
	      if(local_count[2*i]<0) ERR.General(cname,fname,"local_count[%d]=%d]n",2*i,local_count[2*i]);
	      (uc_l[2*i]+local_count[2*i])->src = LexVector(x)*vlen;
	      (uc_l[2*i]+local_count[2*i])->dest = LexVector(nei)*vlen2;
	      local_count[i*2]++;
	      if (local_count[i*2]>local_chi[i*2])
		ERR.General(cname,fname,"local_count[%d](%d)>local_chi[%d](%d)\n",2*i,local_count[2*i],2*i,local_chi[2*i]);
	    }
	    // negative direction
	    if(x[i] == (size[i] -1)){
	      nei[i] = 0;
	      (uc_nl[2*i+1]+non_local_count[2*i+1])->src = non_local_count[2*i+1]*vlen;
	      (uc_nl[2*i+1]+non_local_count[2*i+1])->dest = LexVector(nei)*vlen2;
	      non_local_count[i*2+1]++;
	      if (non_local_count[i*2+1]>non_local_chi[i*2+1])
		ERR.General(cname,fname,"non_local_count[%d](%d)>non_local_chi[%d](%d)\n",2*i+1,non_local_count[2*i+1],2*i+1,non_local_chi[2*i+1]);
	    } else {
	      nei[i] = x[i]+1;
	      if(local_count[2*i+1]<0) ERR.General(cname,fname,"local_count[%d]=%d]n",2*i+local_count[2*i+1]);
	      (uc_l[2*i+1]+local_count[2*i+1])->src = LexVector(x)*vlen;
	      (uc_l[2*i+1]+local_count[2*i+1])->dest = LexVector(nei)*vlen2;
	      local_count[i*2+1]++;
	      if (local_count[i*2+1]>local_chi[i*2+1])
		ERR.General(cname,fname,"local_count[%d](%d)>local_chi[%d](%d)\n",2*i+1,local_count[2*i+1],2*i+1,local_chi[2*i+1]);
	    }
	    nei[i] = x[i];
	  }
	}
  for(i=0;i<NDIM;i++){
    memset(uc_l[2*i]+local_count[2*i],0,sizeof(gauge_agg));
    memset(uc_l[2*i+1]+local_count[2*i+1],0,sizeof(gauge_agg));
    memset(uc_nl[2*i]+non_local_count[2*i],0,sizeof(gauge_agg));
    memset(uc_nl[2*i+1]+non_local_count[2*i+1],0,sizeof(gauge_agg));
  }

  int temp=1;
  for(i=0;i<NDIM;i++){
    offset[2*i] = 0;
    offset[2*i+1] = temp*(size[i]-1);
    temp *= size[i];
    //printf("offset[%d]=%d offset[%d]=%d\n",2*i,offset[2*i],2*i+1,offset[2*i+1]);
  }
#if 0
  for(x[2]=0,nei[2]=0;x[2]<size[2];x[2]++,nei[2]++)
    for(x[1]=0,nei[1]=0;x[1]<size[1];x[1]++,nei[1]++)
      for(x[0]=0,nei[0]=0;x[0]<size[0];x[0]++,nei[0]++)
	for(x[3]=0,nei[3]=0;x[3]<size[3];x[3]++,nei[3]++){
	  int odd = (x[0]+x[1]+x[2]+x[3])%2;
	  for(i=0;i<NDIM;i++){
	    // positive direction
	    if(x[i] == 0){
	      nei[i] = size[i]-1;  
	      (uc_nl_txyz[odd][2*i]+non_local_count_txyz[odd][2*i])->src = non_local_count_txyz[odd][2*i]*vlen;
	      (uc_nl_txyz[odd][2*i]+non_local_count_txyz[odd][2*i])->dest = LexVector_txyz(nei)*vlen2;
	      non_local_count_txyz[odd][i*2]++;
	      if( i == 3){ //t
		Toffset[odd][0][Tcount[odd][0]] = LexVector(x)*vlen;
		Tcount[odd][0]++;
	      }
	    } else {
	      nei[i] = x[i]-1;
	      (uc_l_txyz[odd][2*i]+local_count_txyz[odd][2*i])->src = LexVector_txyz(x)*vlen;
	      (uc_l_txyz[odd][2*i]+local_count_txyz[odd][2*i])->dest = LexVector_txyz(nei)*vlen2;
	      local_count_txyz[odd][i*2]++;
	    }
	    // negative direction
	    if(x[i] == (size[i] -1)){
	      nei[i] = 0;
	      (uc_nl_txyz[odd][2*i+1]+non_local_count_txyz[odd][2*i+1])->src = non_local_count_txyz[odd][2*i+1]*vlen;
	      (uc_nl_txyz[odd][2*i+1]+non_local_count_txyz[odd][2*i+1])->dest = LexVector_txyz(nei)*vlen2;
	      non_local_count_txyz[odd][i*2+1]++;
	      if( i == 3){ //t
		Toffset[odd][1][Tcount[odd][1]] = LexVector(x)*vlen;
		Tcount[odd][1]++;
	      }
	    } else {
	      nei[i] = x[i]+1;
	      (uc_l_txyz[odd][2*i+1]+local_count_txyz[odd][2*i+1])->src = LexVector_txyz(x)*vlen;
	      (uc_l_txyz[odd][2*i+1]+local_count_txyz[odd][2*i+1])->dest = LexVector_txyz(nei)*vlen2;
	      local_count_txyz[odd][i*2+1]++;
	    }
	    nei[i] = x[i];
	  }
	}
#endif
}

void ParTransAsqtad::setHopPointer() {

  char *fname = "setHopPointer()";

  int vlen = VECT_LEN*sizeof(IFloat);
  int vlen2 =VECT_LEN2*sizeof(IFloat);

  int size[NDIM], x[NDIM], nei[NDIM];
  
  int hp_local_count[MAX_HOP][2*NDIM];
  int hp_non_local_count[MAX_HOP][2*NDIM];
  int hop, i;

  size[0] = GJP.XnodeSites();
  size[1] = GJP.YnodeSites();
  size[2] = GJP.ZnodeSites();
  size[3] = GJP.TnodeSites();

  for (hop=0; hop<MAX_HOP; hop++) {
    for (i=0; i<2*NDIM; i++) {
      hp_non_local_count[hop][i] = 0;
      hp_local_count[hop][i] = 0;
    }
  }
  
  for (hop = 1; hop <= MAX_HOP; hop++) {
    hop_pointer **h_l = hp_l[hop-1];
    hop_pointer **h_nl = hp_nl[hop-1];

    int *local_count = hp_local_count[hop-1];
    int *non_local_count = hp_non_local_count[hop-1];

    for (i=0; i<NDIM; i++) {

      int non_local_check = hop*non_local_chi[i*2];
      int local_check = vol - non_local_check;;

      for(x[3]=0,nei[3]=0;x[3]<size[3];x[3]++,nei[3]++)
	for(x[2]=0,nei[2]=0;x[2]<size[2];x[2]++,nei[2]++)
	  for(x[1]=0,nei[1]=0;x[1]<size[1];x[1]++,nei[1]++)
	    for(x[0]=0,nei[0]=0;x[0]<size[0];x[0]++,nei[0]++){

	      // positive direction
	      if(x[i] < hop){
		nei[i] = size[i]-hop+x[i];  
		(h_nl[2*i]+non_local_count[2*i])->src = non_local_count[2*i]*vlen;
		(h_nl[2*i]+non_local_count[2*i])->dest = LexVector(nei)*vlen2;
		non_local_count[i*2]++;
		if (non_local_count[i*2]>non_local_check)
		  ERR.General(cname,fname,"non_local_count[%d](%d)>non_local_check[%d](%d)\n",
			      2*i,non_local_count[2*i],2*i,non_local_check);
	      } else {
		nei[i] = x[i]-hop;
		(h_l[2*i]+local_count[2*i])->src = LexVector(x)*vlen;
		(h_l[2*i]+local_count[2*i])->dest = LexVector(nei)*vlen2;
		local_count[i*2]++;
		if (local_count[i*2]>local_check)
		  ERR.General(cname,fname,"local_count[%d](%d)>local_check[%d](%d)\n",
			      2*i,local_count[2*i],2*i,local_check);
	      }
	      
	      // negative direction
	      if(x[i] >= (size[i]-hop)){
		nei[i] = x[i]+hop-size[i];
		(h_nl[2*i+1]+non_local_count[2*i+1])->src = non_local_count[2*i+1]*vlen;
		(h_nl[2*i+1]+non_local_count[2*i+1])->dest = LexVector(nei)*vlen2;
		non_local_count[i*2+1]++;
		if (non_local_count[i*2]>non_local_check)
		  ERR.General(cname,fname,"non_local_count[%d](%d)>non_local_check[%d](%d)\n",
			      2*i,non_local_count[2*i],2*i,non_local_check);
	      } else {
		nei[i] = x[i]+hop;
		(h_l[2*i+1]+local_count[2*i+1])->src = LexVector(x)*vlen;
		(h_l[2*i+1]+local_count[2*i+1])->dest = LexVector(nei)*vlen2;
		local_count[i*2+1]++;
		if (local_count[i*2]>local_check)
		  ERR.General(cname,fname,"local_count[%d](%d)>local_check[%d](%d)\n",
			      2*i,local_count[2*i],2*i,local_check);
	      }
	      // Need to reset the neighbour pointer
	      nei[i] = x[i];
	    }
    }
  }

}

void ParTransAsqtad::pt_delete(){
  char *fname = "pt_delete()";
  VRB.Func(cname,fname);
	
  for(int i = 0; i < 2*NDIM; i++){
    qfree(uc_l[i]);
    qfree(uc_nl[i]);
#if 0
    for(int j = 0; j < 2; j++){
      sfree(uc_l_txyz[j][i]);
      sfree(uc_nl_txyz[j][i]);
    }
#endif
    sfree(rcv_buf[i]);
    //		qfree(tmp_buf[i]);
  }
#if 0
  for(int i = 0; i < 2; i++)
    for(int j = 0; j < 2; j++)
      sfree(Toffset[j][i]);
#endif
  for (int hop=0; hop<MAX_HOP; hop++) {
    for(int i = 0; i < 2*NDIM; i++){
      sfree(hp_l[hop][i]);
      sfree(hp_nl[hop][i]);
    }
  }
}

void ParTransAsqtad::pt_delete_g(){
  char *fname = "pt_delete_g()";
  VRB.Func(cname,fname);
  for (int hop=0; hop<MAX_HOP; hop++) {
    for(int i = 0; i < 4*NDIM; i++){
      delete SCUarg[hop][i];
      delete SCUarg_mat[hop][i];
      //		delete SCUarg_txyz[i];
    }
  }
}

static IFloat rcv[GAUGE_LEN];
void ParTransAsqtad::pt_init_g(void){
  int x[NDIM], nei[NDIM];
  int local_count[2*NDIM];
  int non_local_count[2*NDIM];
  int i;

  char *fname = "pt_init_g()";
  VRB.Func(cname,fname);
  for(i=0;i<18;i++) rcv[i]=0.;
  for(i=0; i<2*NDIM;i++){
    local_count[i]=non_local_count[i]=0;
  }
  IFloat *u = gauge_field_addr;

  SCUDir rcv_dir[]={SCU_XP, SCU_XM, SCU_YP, SCU_YM, SCU_ZP, SCU_ZM,SCU_TP,SCU_TM};
  SCUDir snd_dir[]={SCU_XM, SCU_XP, SCU_YM, SCU_YP, SCU_ZM, SCU_ZP,SCU_TM,SCU_TP};
  for(x[3]=0,nei[3]=0;x[3]<size[3];x[3]++,nei[3]++)
    for(x[2]=0,nei[2]=0;x[2]<size[2];x[2]++,nei[2]++)
      for(x[1]=0,nei[1]=0;x[1]<size[1];x[1]++,nei[1]++)
	for(x[0]=0,nei[0]=0;x[0]<size[0];x[0]++,nei[0]++){
	  for(i=0;i<NDIM;i++){
	    // positive direction
	    if(x[i] == 0){
	      nei[i] = size[i]-1;  
	      DagCopy((uc_nl[2*i]+non_local_count[2*i])->mat, u+LexGauge(nei,i)*GAUGE_LEN);
	      non_local_count[i*2]++;
	    } else {
	      nei[i] = x[i]-1;
	      DagCopy((uc_l[2*i]+local_count[2*i])->mat, u+LexGauge(nei,i)*GAUGE_LEN);
	      local_count[i*2]++;
	    }
	    // negative direction
	    if(x[i] == (size[i] -1)){
	      nei[i] = 0;
#if 1
	      getMinusData( rcv, u+LexGauge(x,i)*GAUGE_LEN, GAUGE_LEN, i);
	      Copy((uc_nl[2*i+1]+non_local_count[2*i+1])->mat, rcv);
#else
	      Copy((uc_nl[2*i+1]+local_count[2*i+1])->mat, u+LexGauge(x,i)*GAUGE_LEN);
#endif
	      non_local_count[i*2+1]++;
	    } else {
	      nei[i] = x[i]+1;
	      Copy((uc_l[2*i+1]+local_count[2*i+1])->mat, u+LexGauge(x,i)*GAUGE_LEN);
	      local_count[i*2+1]++;
	    }
	    nei[i] = x[i];
	  }
	}
  for(i=0;i<2*NDIM;i++) {
    for (int hop=1; hop<=MAX_HOP; hop++) {
      SCUarg[hop-1][i*2] = new SCUDirArgIR;
      SCUarg[hop-1][i*2]->Init((void *)rcv_buf[i],rcv_dir[i],SCU_REC,
			hop*non_local_chi[i]*VECT_LEN*sizeof(IFloat),1,0,IR_9);
      SCUarg[hop-1][i*2+1] = new SCUDirArgIR;
      SCUarg[hop-1][i*2+1]->Init((void *)rcv_buf[i],snd_dir[i],SCU_SEND,hop*blklen[i],
				 numblk[i],stride[i]+(1-hop)*blklen[i],IR_9);
    
      SCUarg_mat[hop-1][i*2] = new SCUDirArgIR;
      SCUarg_mat[hop-1][i*2]->Init((void *)rcv_buf[i],rcv_dir[i],SCU_REC,
				   3*hop*non_local_chi[i]*VECT_LEN*sizeof(IFloat),1,0,IR_9);
      SCUarg_mat[hop-1][i*2+1] = new SCUDirArgIR;
      SCUarg_mat[hop-1][i*2+1]->Init((void *)rcv_buf[i],snd_dir[i],SCU_SEND,
				     3*hop*blklen[i],numblk[i],
				     3*(stride[i]+(1-hop)*blklen[i]),IR_9);
    }
  }
  //	printf("pt_init_g() done \n");
}

#undef PROFILE
void ParTransAsqtad::run(int n, Matrix **mout, Matrix **min, const int *dir){
    
  int wire[n];
  int i;
  SCUDirArgIR *SCUarg_p[2*n];
  SCUDirArgMulti SCUmulti;
  static int call_num = 0;

  call_num++;
  //	printf("run(i,M**,M**,i):call_num=%d\n",call_num);
  char *fname="run(i,M**,M**,i)";
  VRB.Func(cname,fname);
	
  //	for(i=0;i<n;i++) wire[i] = (dir[i]+2)%(2*NDIM); // from (x,y,z,t) to (t,x,y,z)
  for(i=0;i<n;i++) wire[i] = dir[i]; 
#ifdef PROFILE
  Float dtime  = - dclock();
  //	int nflops=0;
#endif
  for(i=0;i<n;i++) {
    Matrix * addr = (min[i]+offset[wire[i]]);
    //		printf("min[%d]=%p wire[%d]=%d addr=%p\n",i,min[i],i,wire[i],addr);
    SCUarg_p[2*i] = SCUarg_mat[0][2*wire[i]];
    SCUarg_p[2*i+1] = SCUarg_mat[0][2*wire[i]+1];
    SCUarg_p[2*i+1]->Addr((void *)addr);
  }
  SCUmulti.Init(SCUarg_p,n*2);
  SCUmulti.SlowStartTrans();
			
  for(i=0;i<n;i++){
    //		nflops +=local_chi[wire[i]]*198;
#ifdef CPP
    cmm_agg_cpp(local_chi[wire[i]],0, (long)uc_l[wire[i]], (long)min[i],(long)mout[i]);
#else
    cmm_agg(uc_l[wire[i]],min[i],mout[i],local_chi[wire[i]]/2);
#endif
  }
  //	dtime +=dclock();
  //	print_flops(nflops,dtime);
  SCUmulti.TransComplete();
  //	dtime = -dclock();
  //	nflops = 0;
  for(i=0;i<n;i++){ 
    //		nflops +=non_local_chi[wire[i]]*198;
#ifdef CPP
    cmm_agg_cpp(non_local_chi[wire[i]],0, (long)uc_nl[wire[i]], (long)rcv_buf[wire[i]],(long)mout[i]);
#else
    cmm_agg(uc_nl[wire[i]],(Matrix *)rcv_buf[wire[i]],mout[i],non_local_chi[wire[i]]/2);
#endif
  }
#ifdef PROFILE
  dtime +=dclock();
  //	print_flops(nflops,dtime);
  print_flops(cname,fname,198*vol*n,dtime);
#endif
  PTflops +=198*n*vol;
#if 0
  for(i=0;i<n;i++){
    printf("mout[%d]= ",i);
    IFloat * f_tmp = (IFloat *)(mout[i]);
    for(j=0;j<18*vol;j++){
      printf("%e ",*f_tmp++);
      if (j%6==5) printf("\n");
    }
  }
  if(call_num==20) exit(6);
#endif
}

#undef PROFILE
void ParTransAsqtad::run(int n, Vector **vout, Vector **vin, const int *dir){
  int i;
  static int call_num=0;
  SCUDirArgIR *SCUarg_p[2*n];
  call_num++;
#ifdef PROFILE
  Float dtime  = - dclock();
#endif
  int wire[n];
  SCUDirArgMulti SCUmulti;

  //char *fname="run(i,V**,V**,i)";
  //	VRB.Func(cname,fname);
	
  for(i=0;i<n;i++) wire[i] = dir[i]; // from (x,y,z,t) to (t,x,y,z)
  for(i=0;i<n;i++) {
    Vector * addr = (vin[i]+offset[wire[i]]);
    SCUarg_p[2*i] = SCUarg[0][2*wire[i]];
    SCUarg_p[2*i+1] = SCUarg[0][2*wire[i]+1];
    SCUarg_p[2*i+1]->Addr((void *)addr);
  }
  SCUmulti.Init(SCUarg_p,n*2);
  SCUmulti.SlowStartTrans();
	
#if 1
  for(i=0;i<n;i++) pt_asqtad_agg(local_chi[wire[i]],0, (long)uc_l[wire[i]], (long)vin[i],(long)vout[i]);
#else
  for(i=0;i<n;i++) cmv_agg_cpp(local_chi[wire[i]],(long)uc_l[wire[i]], (long)vin[i],(long)vout[i]);
#endif
  SCUmulti.TransComplete();
	
#if 1
  for(i=0;i<n;i++) pt_asqtad_agg(non_local_chi[wire[i]],0, (long)uc_nl[wire[i]], (long)rcv_buf[wire[i]],(long)vout[i]);
#else
  for(i=0;i<n;i++) cmv_agg_cpp(non_local_chi[wire[i]],(long)uc_nl[wire[i]], (long)rcv_buf[wire[i]],(long)vout[i]);
#endif
#ifdef PROFILE
  dtime +=dclock();
  print_flops(cname,fname,66*n*vol,dtime);
#endif
  PTflops +=66*n*vol;
}

/*! 
  Computes sum[x] = vect[x] vect[x + hop dir]^dagger
  where the sum is over n_vect vectors and the hop is in a forward direction.
*/
void ParTransAsqtad::vvpd(Vector **vect, int n_vect, const int *dir, 
			  int n_dir, int hop, Matrix **sum){
  int i, s, v;
  int wire[n_dir];
  for(i=0;i<n_dir;i++) wire[i] = dir[i]; // from (x,y,z,t) to (t,x,y,z)

  SCUDirArgIR *SCUarg_p[2*n_dir];  

  // Only do communication in forward direction
  for(i=0;i<n_dir;i++) {
    SCUarg_p[2*i] = SCUarg[hop-1][4*wire[i]];
    SCUarg_p[2*i+1] = SCUarg[hop-1][4*wire[i]+1];
  }

  Vector *v1, *v2;
  Matrix *s1;
  Matrix m;

  for (i=0; i<n_dir; i++)
    for (s=0; s<vol; s++) sum[i][s].ZeroMatrix();

  for(v=0; v<n_vect; v++){
    SCUDirArgMulti SCUmulti;

    for(i=0;i<n_dir;i++)
      SCUarg_p[2*i+1]->Addr((void *)(vect[v]+Offset(2*wire[i], hop)));

    // Start communication
    SCUmulti.Init(SCUarg_p,2*n_dir);
    SCUmulti.SlowStartTrans();

    // Perform local calculation
    for(i=0; i<n_dir; i++){ 
      for(s=0; s<vol-hop*non_local_chi[2*wire[i]]; s++) {
	v1 = (Vector*)((int)vect[v] + hp_l[hop-1][2*wire[i]][s].dest);
	v2 = (Vector*)((int)vect[v] + hp_l[hop-1][2*wire[i]][s].src);
	s1 = (Matrix*)((int)sum[i] + 3*hp_l[hop-1][2*wire[i]][s].dest);
	m.Cross2(*v1, *v2);
	*s1 += m;
      }
    }
    
    // Finalise communication
    SCUmulti.TransComplete();

    // Perform non-local calculation
    for(i=0; i<n_dir; i++){ 
      for(s=0; s<hop*non_local_chi[2*wire[i]]; s++) {
	v1 = (Vector*)((int)vect[v] + hp_nl[hop-1][2*wire[i]][s].dest);
	v2 = (Vector*)((int)rcv_buf[2*wire[i]] + hp_nl[hop-1][2*wire[i]][s].src);
	s1 = (Matrix*)((int)sum[i] + 3*hp_nl[hop-1][2*wire[i]][s].dest);
	m.Cross2(*v1, *v2);
	*s1 += m;
      }
    }  

  }
  PTflops += 90*n_vect*n_dir*vol;

}

//! u[x] = v[x+dir] for n_dir forward or backward directions dir.
void ParTransAsqtad::shift_field(Matrix **v, const int *dir, int n_dir,
				 int hop, Matrix **u){

  int i, s;
  int wire[n_dir];
  for (i=0; i<n_dir;i++) wire[i] = dir[i];
  SCUDirArgMulti SCUmulti;
  SCUDirArgIR *SCUarg_p[2*n_dir];

  for (i=0; i<n_dir; i++) {
    SCUarg_p[2*i] = SCUarg_mat[hop-1][2*wire[i]];
    SCUarg_p[2*i+1] = SCUarg_mat[hop-1][2*wire[i]+1];
    SCUarg_p[2*i+1]->Addr((void *)(v[i]+Offset(wire[i], hop)));
  }

  SCUmulti.Init(SCUarg_p,2*n_dir);
  SCUmulti.SlowStartTrans();

  Matrix *m1, *m2;

  for (i=0; i<n_dir; i++)
    for (s=0; s<vol-hop*non_local_chi[wire[i]]; s++) {
      m1 = (Matrix*)((int)u[i] + 3*hp_l[hop-1][wire[i]][s].dest);
      m2 = (Matrix*)((int)v[i] + 3*hp_l[hop-1][wire[i]][s].src);
      *m1 = *m2;
    }

  SCUmulti.TransComplete();

  for (i=0; i<n_dir; i++)
    for (s=0; s<hop*non_local_chi[wire[i]]; s++) {
      m1 = (Matrix*)((int)u[i] + 3*hp_nl[hop-1][wire[i]][s].dest);
      m2 = (Matrix*)((int)rcv_buf[wire[i]] + 3*hp_nl[hop-1][wire[i]][s].src);
      *m1 = *m2;
    }
}

//! u[-/+nu](x) = U_[-/+nu](x) 
void ParTransAsqtad::shift_link(Matrix **u, const int *dir, int n_dir){

  for (int i=0; i<n_dir; i++) {
    for (int s=0; s<local_chi[dir[i]]; s++) {
      Copy((IFloat*)((int)u[i] + 3*hp_l[0][dir[i]][s].dest),
	   (IFloat*)(uc_l[dir[i]][s].mat)); 
    }

    for (int s=0; s<non_local_chi[dir[i]]; s++) {
      Copy((IFloat*)((int)u[i] + 3*hp_nl[0][dir[i]][s].dest),
	   (IFloat*)(uc_nl[dir[i]][s].mat));
    }
  }

}


CPS_END_NAMESPACE
#endif
