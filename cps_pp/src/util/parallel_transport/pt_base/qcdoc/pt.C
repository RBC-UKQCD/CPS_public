/*! \file
  \brief  Definition of parallel transport definitions for QCDOC.
  
  $Id: pt.C,v 1.5 2004-08-18 11:58:07 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-08-18 11:58:07 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_base/qcdoc/pt.C,v 1.5 2004-08-18 11:58:07 zs Exp $
//  $Id: pt.C,v 1.5 2004-08-18 11:58:07 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: pt.C,v $
//  $Revision: 1.5 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_base/qcdoc/pt.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
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
void cmm_agg(gauge_agg *chi, IFloat *phi,IFloat *result, int counter);
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
static int blklen[2*SCUMachDim];
static int numblk[2*SCUMachDim];
static int stride[2*SCUMachDim];
static int local_chi[2*SCUMachDim];
static int non_local_chi[2*SCUMachDim];
static int offset[2*SCUMachDim];
static   IFloat * rcv_buf[2*6];
static   IFloat *gauge_field_addr;

static SCUDirArgIR *SCUarg[MAX_HOP][4*SCUMachDim];
static SCUDirArgIR *SCUarg_mat[MAX_HOP][4*SCUMachDim];


static void (*Copy) (IFloat *dest, IFloat *src);
static void (*DagCopy) (IFloat *dest, IFloat *src);
static int (*LexVector)(int *x);
static int (*LexGauge) (int *x,int mu);

static void cpy (IFloat *dest, IFloat *src){
	for(int i=0;i<18;i++)
		dest[i]=src[i];
}

static void dag_cpy (IFloat *dest, IFloat *src){
	for(int i=0;i<3;i++)
	for(int j=0;j<3;j++){
		dest[2*(3*i+j)]=src[2*(3*j+i)];
		dest[2*(3*i+j)+1]=-src[2*(3*j+i)+1];
	}
}

static int lex_xyzt(int *x){
	static char *fname = "lex_xyzt()";
//	printf("%s: %d %d %d %d\n",fname,x[0],x[1],x[2],x[3]);
	int result = x[0] + size[0]*(x[1]+size[1]*(x[2]+size[2]*x[3]));
	return result;
}

static int lex_xyzt_cb_o(int *x){
	int result = x[0] + size[0]*(x[1]+size[1]*(x[2]+size[2]*x[3]));
  if ( (x[0]+x[1]+x[2]+x[3])%2 == 0) result = result/2+vol/2;
  else result = result/2;
  return result;
}

static int lex_txyz(int *x){
	return  (x[3] + size[3]*(x[0]+size[0]*(x[1]+size[1]*x[2])))/2 ;
}

static int LexSurface(int *x){
	return  (x[0]+size[0]*(x[1]+size[1]*x[2]))/2 ;
}

static int lex_g_xyzt(int *x, int mu){
	int temp =  lex_xyzt(x);
	return (temp*NDIM + mu);
}

static int lex_g_xyzt_cb_o(int *x, int mu){
	int temp =  lex_xyzt_cb_o(x);
	return (temp*NDIM + mu);
}

// Calculate the required offset given the direction and hop
int pt_offset(int dir, int hop) {

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

void pt_set_hop_pointer() {

  char *fname = "set_hop_pointer()";

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
//			printf("%d %d %d %d %d hop=%d\n",x[0],x[1],x[2],x[3],i,hop);

	      // positive direction
	      if(x[i] < hop){
		nei[i] = size[i]-hop+x[i];  
		(h_nl[2*i]+non_local_count[2*i])->src = non_local_count[2*i]*vlen;
		(h_nl[2*i]+non_local_count[2*i])->dest = LexVector(nei)*vlen2;
		non_local_count[i*2]++;
		if (non_local_count[i*2]>non_local_check)
		  ERR.General("",fname,"non_local_count[%d](%d)>non_local_check[%d](%d)\n",
			      2*i,non_local_count[2*i],2*i,non_local_check);
	      } else {
		nei[i] = x[i]-hop;
		(h_l[2*i]+local_count[2*i])->src = LexVector(x)*vlen;
		(h_l[2*i]+local_count[2*i])->dest = LexVector(nei)*vlen2;
		local_count[i*2]++;
		if (local_count[i*2]>local_check)
		  ERR.General("",fname,"local_count[%d](%d)>local_check[%d](%d)\n",
			      2*i,local_count[2*i],2*i,local_check);
	      }
	      
	      // negative direction
	      if(x[i] >= (size[i]-hop)){
		nei[i] = x[i]+hop-size[i];
		(h_nl[2*i+1]+non_local_count[2*i+1])->src = non_local_count[2*i+1]*vlen;
		(h_nl[2*i+1]+non_local_count[2*i+1])->dest = LexVector(nei)*vlen2;
		non_local_count[i*2+1]++;
		if (non_local_count[i*2]>non_local_check)
		  ERR.General("",fname,"non_local_count[%d](%d)>non_local_check[%d](%d)\n",
			      2*i,non_local_count[2*i],2*i,non_local_check);
	      } else {
		nei[i] = x[i]+hop;
		(h_l[2*i+1]+local_count[2*i+1])->src = LexVector(x)*vlen;
		(h_l[2*i+1]+local_count[2*i+1])->dest = LexVector(nei)*vlen2;
		local_count[i*2+1]++;
		if (local_count[i*2]>local_check)
		  ERR.General("",fname,"local_count[%d](%d)>local_check[%d](%d)\n",
			      2*i,local_count[2*i],2*i,local_check);
	      }
	      // Need to reset the neighbour pointer
	      nei[i] = x[i];
	    }
    }
  }

}

void pt_init(Lattice &lat)
{
	char *cname = "";
	char *fname = "pt_init()";
	VRB.Func("",fname);
	int i, j, x[NDIM],nei[NDIM];
	int local_count[2*NDIM];
	int non_local_count[2*NDIM];
	int vlen = VECT_LEN*sizeof(IFloat); //size of incoming vector
	int vlen2 =VECT_LEN2*sizeof(IFloat); //size of outgoing vector (maybe different to optimize for QCDOC PEC)

	size[0] = GJP.XnodeSites();
	size[1] = GJP.YnodeSites();
	size[2] = GJP.ZnodeSites();
	size[3] = GJP.TnodeSites();
	gauge_field_addr = (IFloat *) lat.GaugeField();
	StrOrdType str_ord = lat.StrOrd();
	if (str_ord == CANONICAL){
		Copy = cpy; DagCopy = dag_cpy;
		LexVector = lex_xyzt;
		LexGauge = lex_g_xyzt;
	} else if (str_ord == WILSON){
		Copy = cpy; DagCopy = dag_cpy;
		LexVector = lex_xyzt_cb_o;
		LexGauge = lex_g_xyzt_cb_o;
	} else if (str_ord == STAG){
		Copy = dag_cpy; DagCopy = cpy;
		LexVector = lex_xyzt;
		LexGauge = lex_g_xyzt;
	} else
    ERR.General(cname,fname,"storage ordering not implemented");

	blklen[0] = blklen[1]= vlen;
	for(i=1;i<NDIM;i++) {blklen[2*i+1] = blklen[2*i] = blklen[2*i-1]*size[i-1]; }
	numblk[2*NDIM-1]=numblk[2*NDIM-2]=1;
	for(i=NDIM-2;i>=0;i--) {numblk[i*2+1] = numblk[2*i] = numblk[2*i+2]*size[i+1]; }
  for(i=0;i<NDIM*2;i++)  stride[i] = blklen[i]* (size[i/2]-1);
	vol = 1;
	for(i=0; i<NDIM;i++) vol *= size[i];
	for(i=0; i<NDIM;i++){
		non_local_chi[2*i+1] = non_local_chi[2*i] = vol/size[i];
		local_chi[2*i+1] = local_chi[2*i] = vol - non_local_chi[2*i];
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
		}
		if(uc_nl[i]==NULL) {
			uc_nl[i] = (gauge_agg *)smalloc(sizeof(gauge_agg)*(1+non_local_chi[i]));
		}
		if(uc_l[i]==NULL)ERR.Pointer("",fname,"uc_l[i]");
		if(uc_nl[i]==NULL)ERR.Pointer("",fname,"uc_nl[i]");
		rcv_buf[i] = (IFloat *)smalloc(non_local_chi[i]*vlen*3);
		if(rcv_buf[i]==NULL)ERR.Pointer("",fname,"rcv_buf[i]");
  }

	for(x[3]=0,nei[3]=0;x[3]<size[3];x[3]++,nei[3]++)
	for(x[2]=0,nei[2]=0;x[2]<size[2];x[2]++,nei[2]++)
	for(x[1]=0,nei[1]=0;x[1]<size[1];x[1]++,nei[1]++)
	for(x[0]=0,nei[0]=0;x[0]<size[0];x[0]++,nei[0]++){
		for(i=0;i<NDIM;i++){
//			printf("%d %d %d %d %d\n",x[0],x[1],x[2],x[3],i);
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
  
  pt_set_hop_pointer();
	
}

void pt_delete(){
	char *fname = "pt_delete()";
	VRB.Func("",fname);
	
	for(int i = 0; i < 2*NDIM; i++){
		qfree(uc_l[i]);
		qfree(uc_nl[i]);
		qfree(rcv_buf[i]);
	}
  for (int hop=0; hop<MAX_HOP; hop++) {
    for(int i = 0; i < 2*NDIM; i++){
      sfree(hp_l[hop][i]);
      sfree(hp_nl[hop][i]);
    }
  }
}

void pt_delete_g(){
	char *fname = "pt_delete_g()";
	VRB.Func("",fname);
    for(int hop = 0; hop < MAX_HOP; hop++)
	for(int i = 0; i < 4*NDIM; i++){
		delete SCUarg[hop][i];
		delete SCUarg_mat[hop][i];
	}
}

void pt_init_g(void){
	int x[NDIM], nei[NDIM];
	int local_count[2*NDIM];
	int non_local_count[2*NDIM];
	int i;

	char *fname = "pt_init_g()";
	VRB.Func("",fname);
	for(i=0; i<2*NDIM;i++){
		local_count[i]=non_local_count[i]=0;
	}
	IFloat *u = gauge_field_addr;

	SCUDir rcv_dir[]={SCU_XP, SCU_XM, SCU_YP, SCU_YM, SCU_ZP, SCU_ZM,SCU_TP,SCU_TM};
	SCUDir snd_dir[]={SCU_XM, SCU_XP, SCU_YM, SCU_YP, SCU_ZM, SCU_ZP,SCU_TM,SCU_TP};

	IFloat *rcv_mat = (IFloat *)qalloc(QFAST|QNONCACHE,18*sizeof(IFloat));
    sys_cacheflush(0);
	for(i=0;i<NDIM;i++){
		SCUDirArgIR snd(u,snd_dir[i*2+1],SCU_SEND,sizeof(Matrix));
		SCUDirArgIR rcv(rcv_mat,rcv_dir[i*2+1],SCU_REC,sizeof(Matrix));
		for(x[3]=0,nei[3]=0;x[3]<size[3];x[3]++,nei[3]++)
		for(x[2]=0,nei[2]=0;x[2]<size[2];x[2]++,nei[2]++)
		for(x[1]=0,nei[1]=0;x[1]<size[1];x[1]++,nei[1]++)
		for(x[0]=0,nei[0]=0;x[0]<size[0];x[0]++,nei[0]++){
// positive direction
			if(x[i] == 0){
				nei[i] = size[i]-1;  
				Copy((uc_nl[2*i]+non_local_count[2*i])->mat, u+LexGauge(nei,i)*GAUGE_LEN);
				non_local_count[i*2]++;
			} else {
				nei[i] = x[i]-1;
				Copy((uc_l[2*i]+local_count[2*i])->mat, u+LexGauge(nei,i)*GAUGE_LEN);
				local_count[i*2]++;
			}
// negative direction
			if(x[i] == (size[i] -1)){
				nei[i] = 0;
#if 0
				getMinusData( rcv_mat, u+LexGauge(x,i)*GAUGE_LEN, GAUGE_LEN, i);
#else
				snd.Addr(u+LexGauge(x,i)*GAUGE_LEN);
				snd.StartTrans();rcv.StartTrans();
				snd.TransComplete();rcv.TransComplete();
#endif
				DagCopy((uc_nl[2*i+1]+non_local_count[2*i+1])->mat, rcv_mat);
				non_local_count[i*2+1]++;
			} else {
				nei[i] = x[i]+1;
				DagCopy((uc_l[2*i+1]+local_count[2*i+1])->mat, u+LexGauge(x,i)*GAUGE_LEN);
				local_count[i*2+1]++;
			}
			nei[i] = x[i];
		} // x[]
	} // i
#if 0
	for(i=0;i<2*NDIM;i++) {
		SCUarg[i*2] = new SCUDirArgIR;
		SCUarg[i*2]->Init((void *)rcv_buf[i],rcv_dir[i],SCU_REC,non_local_chi[i]*VECT_LEN*sizeof(IFloat),1,0,IR_9);
		SCUarg[i*2+1] = new SCUDirArgIR;
//
//  inputs a dummy but valid address to pass syscall test and changed later, CJ
//
		SCUarg[i*2+1]->Init((void *)rcv_buf[i],snd_dir[i],SCU_SEND,blklen[i],numblk[i],stride[i],IR_9);
		SCUarg_mat[i*2] = new SCUDirArgIR;
		SCUarg_mat[i*2]->Init((void *)rcv_buf[i],rcv_dir[i],SCU_REC,non_local_chi[i]*VECT_LEN*sizeof(IFloat)*3,1,0,IR_9);
		SCUarg_mat[i*2+1] = new SCUDirArgIR;
		SCUarg_mat[i*2+1]->Init((void *)rcv_buf[i],snd_dir[i],SCU_SEND,blklen[i]*3,numblk[i],stride[i]*3,IR_9);
	}
#else
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
#endif
    qfree(rcv_mat);
}

#undef PROFILE
void pt_mat(int n, IFloat **mout, IFloat **min, const int *dir){
    
	int wire[n];
	int i;
	SCUDirArgIR *SCUarg_p[2*n];
	SCUDirArgMulti SCUmulti;
	static int call_num = 0;

	call_num++;
	char *fname="pt_mat()";
	VRB.Func("",fname);
	
	for(i=0;i<n;i++) wire[i] = dir[i]; 
#ifdef PROFILE
	Float dtime  = - dclock();
#endif
	for(i=0;i<n;i++) {
		IFloat * addr = (min[i]+GAUGE_LEN*offset[wire[i]]);
		SCUarg_p[2*i] = SCUarg_mat[0][2*wire[i]];
		SCUarg_p[2*i+1] = SCUarg_mat[0][2*wire[i]+1];
		SCUarg_p[2*i+1]->Addr((void *)addr);
	}
	SCUmulti.Init(SCUarg_p,n*2);
	SCUmulti.SlowStartTrans();
			
	for(i=0;i<n;i++){
#ifdef CPP
	 cmm_agg_cpp(local_chi[wire[i]],0, (long)uc_l[wire[i]], (long)min[i],(long)mout[i]);
#else
	cmm_agg(uc_l[wire[i]],min[i],mout[i],local_chi[wire[i]]/2);
#endif
	}
	SCUmulti.TransComplete();
	for(i=0;i<n;i++){ 
#ifdef CPP
		cmm_agg_cpp(non_local_chi[wire[i]],0, (long)uc_nl[wire[i]], (long)rcv_buf[wire[i]],(long)mout[i]);
#else
	cmm_agg(uc_nl[wire[i]],rcv_buf[wire[i]],mout[i],non_local_chi[wire[i]]/2);
#endif
	}
#ifdef PROFILE
	dtime +=dclock();
	print_flops("",fname,198*vol*n,dtime);
#endif
	ParTrans::PTflops +=198*n*vol;
}

#undef PROFILE
void pt_1vec(int n, IFloat **vout, IFloat **vin, const int *dir){
	int i;
 	static int call_num=0;
	SCUDirArgIR *SCUarg_p[2*n];
 	call_num++;
#ifdef PROFILE
	Float dtime  = - dclock();
#endif
	int wire[n];
	SCUDirArgMulti SCUmulti;

	char *fname="pt_1vec";
	VRB.Func("",fname);
	
	for(i=0;i<n;i++) wire[i] = dir[i]; // from (x,y,z,t) to (t,x,y,z)
	for(i=0;i<n;i++) {
		IFloat * addr = (vin[i]+VECT_LEN*offset[wire[i]]);
		SCUarg_p[2*i] = SCUarg[0][2*wire[i]];
		SCUarg_p[2*i+1] = SCUarg[0][2*wire[i]+1];
		SCUarg_p[2*i+1]->Addr((void *)addr);
	}
	SCUmulti.Init(SCUarg_p,n*2);
	SCUmulti.SlowStartTrans();
	
#ifndef CPP
	for(i=0;i<n;i++) pt_asqtad_agg(local_chi[wire[i]],0, (long)uc_l[wire[i]], (long)vin[i],(long)vout[i]);
#else
	for(i=0;i<n;i++) cmv_agg_cpp(local_chi[wire[i]],(long)uc_l[wire[i]], (long)vin[i],(long)vout[i]);
#endif
	SCUmulti.TransComplete();
	
#ifndef CPP
	for(i=0;i<n;i++) pt_asqtad_agg(non_local_chi[wire[i]],0, (long)uc_nl[wire[i]], (long)rcv_buf[wire[i]],(long)vout[i]);
#else
	for(i=0;i<n;i++) cmv_agg_cpp(non_local_chi[wire[i]],(long)uc_nl[wire[i]], (long)rcv_buf[wire[i]],(long)vout[i]);
#endif
#ifdef PROFILE
	dtime +=dclock();
	print_flops("",fname,66*n*vol,dtime);
#endif
	ParTrans::PTflops +=66*n*vol;
}

/*! 
  Computes sum[x] = vect[x] vect[x + hop dir]^dagger
  where the sum is over n_vect vectors and the hop is in a forward direction.
*/
void pt_vvpd(Vector **vect, int n_vect, const int *dir, 
			  int n_dir, int hop, Matrix **sum){
char *fname = "pt_vvpd()";
VRB.Func("",fname);
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
      SCUarg_p[2*i+1]->Addr((void *)(vect[v]+pt_offset(2*wire[i], hop)));

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
  ParTrans::PTflops += 90*n_vect*n_dir*vol;

}

//! u[x] = v[x+dir] for n_dir forward or backward directions dir.
void pt_shift_field(Matrix **v, const int *dir, int n_dir,
				 int hop, Matrix **u){

  int i, s;
  int wire[n_dir];
  for (i=0; i<n_dir;i++) wire[i] = dir[i];
  SCUDirArgMulti SCUmulti;
  SCUDirArgIR *SCUarg_p[2*n_dir];

  for (i=0; i<n_dir; i++) {
    SCUarg_p[2*i] = SCUarg_mat[hop-1][2*wire[i]];
    SCUarg_p[2*i+1] = SCUarg_mat[hop-1][2*wire[i]+1];
    SCUarg_p[2*i+1]->Addr((void *)(v[i]+pt_offset(wire[i], hop)));
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
void pt_shift_link(Matrix **u, const int *dir, int n_dir){

char *fname = "pt_shift_link()";
VRB.Func("",fname);
  for (int i=0; i<n_dir; i++) {
    for (int s=0; s<local_chi[dir[i]]; s++) {
      cpy((IFloat*)((int)u[i] + 3*hp_l[0][dir[i]][s].dest),
	   (IFloat*)(uc_l[dir[i]][s].mat)); 
    }

    for (int s=0; s<non_local_chi[dir[i]]; s++) {
      cpy((IFloat*)((int)u[i] + 3*hp_nl[0][dir[i]][s].dest),
	   (IFloat*)(uc_nl[dir[i]][s].mat));
    }
  }

}

CPS_END_NAMESPACE
