#include <config.h>
#include <util/gjp.h>
#include <util/pt.h>
#include <sysfunc.h>
#include <comms/scu.h>
#include <stdio.h>
CPS_START_NAMESPACE
void dirac_cmv_jcw_agg_cpp( int sites, long chi, long u,long in, long out);
void dirac_cmm_jcw_agg_cpp( int sites, long chi, long u,long in, long out);
extern "C" void dirac_cmv_jcw_agg( int sites, long chi, long u,long in, long out);
static const int NDIM=4;
static int size[NDIM];
enum {GAUGE_LEN=18,VECT_LEN=6, VECT_LEN2=8};
int vol = 1;

gauge_agg *uc_l[2*SCUMachDim];
gauge_agg *uc_nl[2*SCUMachDim];
int *local_chi;
int *non_local_chi;
int *blklen;
int *numblk;
int *stride;
int *offset;

static SCUDirArgIR *SCUarg[4*SCUMachDim];
static SCUDirArgIR *SCUarg_mat[4*SCUMachDim];

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
	return  (x[0] + size[0]*(x[1]+size[1]*(x[2]+size[2]*x[3])));
}

static int LexGauge(int *x, int mu){
#if 0
	int temp =  (x[1] + size[1]*(x[2]+size[2]*(x[3]+size[3]*x[0])));
	return (temp*NDIM + (mu+3)%NDIM);
#else
	int temp =  (x[0] + size[0]*(x[1]+size[1]*(x[2]+size[2]*x[3])));
	return (temp*NDIM + mu);
#endif
}

void ParTransAsqtad::pt_init(const void *gauge_u)
{
	char *fname = "pt_init()";
	VRB.Func(cname,fname);
	IFloat *u_fl = (IFloat *)gauge_u;
//	printf("pt_init():gauge_u=%p\t u_fl[0]=%e\n",gauge_u,u_fl[0]);
	int i, x[NDIM],nei[NDIM];
	int local_count[2*NDIM];
	int non_local_count[2*NDIM];
	int vlen = VECT_LEN*sizeof(IFloat); //size of incoming vector
	int vlen2 =VECT_LEN2*sizeof(IFloat); //size of outgoing vector (maybe different to optimize for QCDOC PEC)

	size[0] = GJP.XnodeSites();
	size[1] = GJP.YnodeSites();
	size[2] = GJP.ZnodeSites();
	size[3] = GJP.TnodeSites();
	gauge_field_addr = (IFloat *) gauge_u;

	local_chi = ( int *)smalloc(sizeof(int)*2*NDIM);
	non_local_chi = ( int *)smalloc(sizeof(int)*2*NDIM);
	blklen = ( int *)smalloc(sizeof(int)*2*NDIM);
	numblk = ( int *)smalloc(sizeof(int)*2*NDIM);
	stride = ( int *)smalloc(sizeof(int)*2*NDIM);
	offset = ( int *)smalloc(sizeof(int)*2*NDIM);

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
		uc_l[i] = (gauge_agg *)smalloc(sizeof(gauge_agg)*local_chi[i]);
		uc_nl[i] = (gauge_agg *)smalloc(sizeof(gauge_agg)*non_local_chi[i]);
		rcv_buf[i] = (IFloat *)smalloc(non_local_chi[i]*vlen*3);
		tmp_buf[i] = (IFloat *)smalloc(vol*vlen2*3);
  }

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
//				printf("%d %d %d %d non_local_count[%d]=%d LexVector(x)=%d LexVector(nei)=%d\n", x[0],x[1],x[2],x[3],i*2,non_local_count[2*i],LexVector(x),LexVector(nei));
				non_local_count[i*2]++;
			} else {
				nei[i] = x[i]-1;
				(uc_l[2*i]+local_count[2*i])->src = LexVector(x)*vlen;
				(uc_l[2*i]+local_count[2*i])->dest = LexVector(nei)*vlen2;
//				printf("%d %d %d %d local_count[%d]=%d LexVector(x)=%d LexVector(nei)=%d\n", x[0],x[1],x[2],x[3],i*2,local_count[2*i],LexVector(x),LexVector(nei));
				local_count[i*2]++;
			}
// negative direction
			if(x[i] == (size[i] -1)){
				nei[i] = 0;
				(uc_nl[2*i+1]+non_local_count[2*i+1])->src = non_local_count[2*i+1]*vlen;
				(uc_nl[2*i+1]+non_local_count[2*i+1])->dest = LexVector(nei)*vlen2;
				non_local_count[i*2+1]++;
			} else {
				nei[i] = x[i]+1;
				(uc_l[2*i+1]+local_count[2*i+1])->src = LexVector(x)*vlen;
				(uc_l[2*i+1]+local_count[2*i+1])->dest = LexVector(nei)*vlen2;
//				printf("%d %d %d %d local_count[%d]=%d LexVector(x)=%d LexVector(nei)=%d\n", x[0],x[1],x[2],x[3],i*2+1,local_count[2*i+1],LexVector(x),LexVector(nei));
				local_count[i*2+1]++;
			}
			nei[i] = x[i];
		}
	}
	int temp=1;
	for(i=0;i<NDIM;i++){
		offset[2*i] = 0;
		offset[2*i+1] = temp*(size[i]-1);
		temp *= size[i];
	}
//	printf("  blklen numblk stride local_chi non_local_chi offset\n");
  for(i=0;i<NDIM*2;i++){
//	printf("%d: %d %d %d %d %d %d\n",i,blklen[i],numblk[i],stride[i],local_chi[i],non_local_chi[i], offset[i]);
  }
}

void ParTransAsqtad::pt_delete(){
	char *fname = "pt_delete()";
	VRB.Func(cname,fname);
	
	sfree (local_chi);
	sfree (non_local_chi);
	sfree (blklen);
	sfree (numblk);
	sfree (stride);
	sfree (offset);
	for(int i = 0; i < 2*NDIM; i++){
		sfree(uc_l[i]);
		sfree(uc_nl[i]);
		sfree(rcv_buf[i]);
		sfree(tmp_buf[i]);
	}
}

void ParTransAsqtad::pt_delete_g(){
	char *fname = "pt_delete_g()";
	VRB.Func(cname,fname);
	for(int i = 0; i < 4*NDIM; i++){
		delete SCUarg[i];
		delete SCUarg_mat[i];
	}
}

static IFloat rcv[GAUGE_LEN];
void ParTransAsqtad::pt_init_g(void){
	int x[NDIM], nei[NDIM];
	int local_count[2*NDIM];
	int non_local_count[2*NDIM];
	int i;
//	FILE *fp = fopen("/host/chulwoo/pt_init_g.out","w");

	char *fname = "pt_init_g()";
	VRB.Func(cname,fname);
	for(i=0;i<18;i++) rcv[i]=0.;
	for(i=0; i<2*NDIM;i++){
		local_count[i]=non_local_count[i]=0;
	}
	IFloat *u = gauge_field_addr;
//	fprintf(fp,"rcv=%p\n",rcv);fflush(fp);
//	fprintf(fp,"u=%p\n",u);fflush(fp);

	SCUDir rcv_dir[]={SCU_XP, SCU_XM, SCU_YP, SCU_YM, SCU_ZP, SCU_ZM,SCU_TP,SCU_TM};
	SCUDir snd_dir[]={SCU_XM, SCU_XP, SCU_YM, SCU_YP, SCU_ZM, SCU_ZP,SCU_TM,SCU_TP};
	for(x[3]=0,nei[3]=0;x[3]<size[3];x[3]++,nei[3]++)
	for(x[2]=0,nei[2]=0;x[2]<size[2];x[2]++,nei[2]++)
	for(x[1]=0,nei[1]=0;x[1]<size[1];x[1]++,nei[1]++)
	for(x[0]=0,nei[0]=0;x[0]<size[0];x[0]++,nei[0]++){
//		fprintf(fp,"x=%d %d %d %d ",x[0],x[1],x[2],x[3]);fflush(fp);
//		fprintf(fp,"nei=%d %d %d %d ",nei[0],nei[1],nei[2],nei[3]);fflush(fp);
		for(i=0;i<NDIM;i++){
//		fprintf(fp,"dir=%d\n",i);fflush(fp);
// positive direction
			if(x[i] == 0){
				nei[i] = size[i]-1;  
				DagCopy((uc_nl[2*i]+non_local_count[2*i])->mat, u+LexGauge(nei,i)*GAUGE_LEN);
//fprintf(fp,"uc_nl[%d][%d]->mat=%p\n",2*i,non_local_count[2*i],(uc_nl[2*i]+non_local_count[2*i])->mat); fflush(fp);
//fprintf(fp,"u[%d]=%e\n",LexGauge(nei,i),*(u+LexGauge(nei,i)*GAUGE_LEN)); fflush(fp);
//				fprintf(fp,"non_local_count[%d]=%d\n",2*i,non_local_count[2*i]);fflush(fp);
				non_local_count[i*2]++;
			} else {
				nei[i] = x[i]-1;
				DagCopy((uc_l[2*i]+local_count[2*i])->mat, u+LexGauge(nei,i)*GAUGE_LEN);
//				fprintf(fp,"local_count[%d]=%d\n",2*i,local_count[2*i]);fflush(fp);
				local_count[i*2]++;
			}
// negative direction
			if(x[i] == (size[i] -1)){
				nei[i] = 0;
#if 1
//				fprintf(fp,"rcv[0]=%e, u[%d][0]=%e\n",rcv[0],LexGauge(x,i),u[LexGauge(x,i)*GAUGE_LEN]);
				getMinusData( rcv, u+LexGauge(x,i)*GAUGE_LEN, GAUGE_LEN*sizeof(IFloat), i);
//				fprintf(fp,"rcv[0]=%e, u[%d][0]=%e\n",rcv[0],LexGauge(x,i),u[LexGauge(x,i)*GAUGE_LEN]);
				Copy((uc_nl[2*i+1]+non_local_count[2*i+1])->mat, rcv);
#else
				Copy((uc_nl[2*i+1]+local_count[2*i+1])->mat, u+LexGauge(x,i)*GAUGE_LEN);
#endif
//				fprintf(fp,"non_local_count[%d]=%d\n",2*i+1,non_local_count[2*i+1]);fflush(fp);
				non_local_count[i*2+1]++;
			} else {
				nei[i] = x[i]+1;
				Copy((uc_l[2*i+1]+local_count[2*i+1])->mat, u+LexGauge(x,i)*GAUGE_LEN);
//				fprintf(fp,"local_count[%d]=%d\n",2*i+1,local_count[2*i+1]);fflush(fp);
				local_count[i*2+1]++;
			}
			nei[i] = x[i];
		}
	}
	for(i=0;i<2*NDIM;i++){
		VRB.Flow(cname,fname,"local_count[%d]=%d\n",i,local_count[i]);
		VRB.Flow(cname,fname,"non_local_count[%d]=%d\n",i,non_local_count[i]);
	}
	for(i=0;i<2*NDIM;i++) {
		SCUarg[i*2] = new SCUDirArgIR;
		SCUarg[i*2]->Init((void *)rcv_buf[i],rcv_dir[i],SCU_REC,non_local_chi[i]*VECT_LEN*sizeof(IFloat),1,0,IR_9);
		SCUarg[i*2+1] = new SCUDirArgIR;
		SCUarg[i*2+1]->Init((void *)0,snd_dir[i],SCU_SEND,blklen[i],numblk[i],stride[i],IR_9);
		SCUarg_mat[i*2] = new SCUDirArgIR;
		SCUarg_mat[i*2]->Init((void *)rcv_buf[i],rcv_dir[i],SCU_REC,non_local_chi[i]*VECT_LEN*sizeof(IFloat)*3,1,0,IR_9);
		SCUarg_mat[i*2+1] = new SCUDirArgIR;
		SCUarg_mat[i*2+1]->Init((void *)0,snd_dir[i],SCU_SEND,blklen[i]*3,numblk[i],stride[i]*3,IR_9);
	}
//	fclose(fp);
}

#if 0
void ParTransAsqtad::run(int n, Vector **vout, Vector **vin, const SCUDir *scudir){
	int dir[n];
	for(int i=0;i<n;i++){
#if  0
		switch(scudir[i]){
			case SCU_XP: dir[i]=0;break;
			case SCU_XM: dir[i]=1;break;
			case SCU_YP: dir[i]=2;break;
			case SCU_YM: dir[i]=3;break;
			case SCU_ZP: dir[i]=4;break;
			case SCU_ZM: dir[i]=5;break;
			case SCU_TP: dir[i]=6;break;
			case SCU_TM: dir[i]=7;break;
			default: printf("wrong scudir"); exit(3);
		}
#else
		dir[i] = (int)scudir[i];
#endif
		printf("dir[%d]=%d\n",i,dir[i]);
	}
	run(n,vout,vin,dir);
}
#endif

void ParTransAsqtad::run(int n, Matrix **mout, Matrix **min, const int *dir){
	int wire[n];
	int i,j;
	SCUDirArgIR *SCUarg_p[2*n];
	SCUDirArgMulti SCUmulti;

	char *fname="run(i,M**,M**,i)";
	VRB.Func(cname,fname);
	
//	for(i=0;i<n;i++) wire[i] = (dir[i]+2)%(2*NDIM); // from (x,y,z,t) to (t,x,y,z)
	for(i=0;i<n;i++) wire[i] = dir[i]; 
	for(i=0;i<n;i++) {
		printf("wire[%d]=%d\n",i,wire[i]);
		Matrix * addr = (min[i]+offset[wire[i]]);
		printf("min[%d]=%p addr=%p\n",i,min[i],addr);
		SCUarg_p[2*i] = SCUarg_mat[2*wire[i]];
		SCUarg_p[2*i+1] = SCUarg_mat[2*wire[i]+1];
		SCUarg_p[2*i+1]->Addr((void *)addr);
	}
	SCUmulti.Init(SCUarg_p,n*2);
	SCUmulti.SlowStartTrans();
	
	for(i=0;i<n;i++) dirac_cmm_jcw_agg_cpp(local_chi[wire[i]],0, (long)uc_l[wire[i]], (long)min[i],(long)tmp_buf[i]);
	SCUmulti.TransComplete();
	
	for(i=0;i<n;i++) dirac_cmm_jcw_agg_cpp(non_local_chi[wire[i]],0, (long)uc_nl[wire[i]], (long)rcv_buf[wire[i]],(long)tmp_buf[i]);
	for(i=0;i<n;i++){
		for(j=0;j<vol;j++){
			IFloat * temp = (IFloat *)(mout[i]+j);
			for(int k = 0;k<VECT_LEN*3;k += 3){
			temp[0] = *(tmp_buf[i]+(j*VECT_LEN2*3)+0);
			temp[1] = *(tmp_buf[i]+(j*VECT_LEN2*3)+1);
			temp[2] = *(tmp_buf[i]+(j*VECT_LEN2*3)+2);
      }
		}
	}
}

void ParTransAsqtad::run(int n, Vector **vout, Vector **vin, const int *dir){
	int wire[n];
	int i,j;
	SCUDirArgIR *SCUarg_p[2*n];
	SCUDirArgMulti SCUmulti;

	char *fname="run(i,V**,V**,i)";
	VRB.Func(cname,fname);
	
	for(i=0;i<n;i++) wire[i] = dir[i]; // from (x,y,z,t) to (t,x,y,z)
	for(i=0;i<n;i++) {
//		printf("wire[%d]=%d\n",i,wire[i]);
		Vector * addr = (vin[i]+offset[wire[i]]);
//		printf("vin[%d]=%p addr=%p\n",i,vin[i],addr);
		SCUarg_p[2*i] = SCUarg[2*wire[i]];
		SCUarg_p[2*i+1] = SCUarg[2*wire[i]+1];
		SCUarg_p[2*i+1]->Addr((void *)addr);
	}
	SCUmulti.Init(SCUarg_p,n*2);
	SCUmulti.SlowStartTrans();
	
	for(i=0;i<n;i++) dirac_cmv_jcw_agg(local_chi[wire[i]],0, (long)uc_l[wire[i]], (long)vin[i],(long)tmp_buf[i]);
	SCUmulti.TransComplete();
	
	for(i=0;i<n;i++) dirac_cmv_jcw_agg(non_local_chi[wire[i]],0, (long)uc_nl[wire[i]], (long)rcv_buf[wire[i]],(long)tmp_buf[i]);
	for(i=0;i<n;i++){
		for(j=0;j<vol;j++){
			IFloat * temp = (IFloat *)(vout[i]+j);
			temp[0] = *(tmp_buf[i]+(j*VECT_LEN2)+0);
			temp[1] = *(tmp_buf[i]+(j*VECT_LEN2)+1);
			temp[2] = *(tmp_buf[i]+(j*VECT_LEN2)+2);
			temp[3] = *(tmp_buf[i]+(j*VECT_LEN2)+3);
			temp[4] = *(tmp_buf[i]+(j*VECT_LEN2)+4);
			temp[5] = *(tmp_buf[i]+(j*VECT_LEN2)+5);
		}
	}
}

CPS_END_NAMESPACE
