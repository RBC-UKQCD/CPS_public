/*
  $Id: main.C,v 1.2 2005-04-05 06:44:52 chulwoo Exp $
*/

#include<config.h>
#include <util/qcdio.h>
#include <util/asqtad_int.h>
#include <util/time.h>
#include <math.h>
#ifdef HAVE_STRINGS_H
#include<strings.h>
#endif

static const char *f_asqtad_test_filename = "f_asqtad_test";
static const char *psi_filename = "psi";
static const char *input_filename = "f_asqtad_inv.in";

static int f_offset (int *size, int *coor, int odd){
  int index = coor[3];
  index = index*size[2]+ coor[2];
  index = index*size[1]+ coor[1];
  index =  index*size[0]+ coor[0];
  if ((coor[0]+coor[1]+coor[2]+coor[3]+odd)%2==0 )
  index = index/2;
  else
  index = index/2 + (size[0]*size[1]*size[2]*size[3])/2;
  
  return index;
}

static int g_offset (int *size, int *coor){
  int index = coor[0];
  index = index*size[3]+ coor[3];
  index = index*size[2]+ coor[2];
  index =  index*size[1]+ coor[1];
  return index;
}

int main(int argc,char *argv[]){

    Start();
    printf("Sizes = %d %d %d %d %d %d\n",SizeX(),SizeY(),SizeZ(),SizeT(),SizeS(),SizeW());
    printf("Coors = %d %d %d %d %d %d\n",CoorX(),CoorY(),CoorZ(),CoorT(),CoorS(),CoorW());

    FILE *fp;
    double dtime;

    //----------------------------------------------------------------
    // Initializes all Global Job Parameters
    //----------------------------------------------------------------
    AsqDArg asq_arg;
    int nx,ny,nz,nt;
    int pos[4];

    if (argc < 9){
        printf("f_asqtad_test::main(): usage: %s nx ny nz nt px py pz pt\n",argv[0]);
	exit(4);
    }
    sscanf(argv[1],"%d",&nx);
    sscanf(argv[2],"%d",&ny);
    sscanf(argv[3],"%d",&nz);
    sscanf(argv[4],"%d",&nt);
    for(int i = 0;i<4;i++)sscanf(argv[5+i],"%d",pos+i);

    int size[4];
    int sites[4];
    asq_arg.NP[0] = SizeT(); asq_arg.size[0] =nt/SizeT();
    asq_arg.NP[1] = SizeX(); asq_arg.size[1] =nx/SizeX();
    asq_arg.NP[2] = SizeY(); asq_arg.size[2] =ny/SizeY();
    asq_arg.NP[3] = SizeZ(); asq_arg.size[3] =nz/SizeZ();

    asq_arg.coor[0] = CoorT();
    asq_arg.coor[1] = CoorX();
    asq_arg.coor[2] = CoorY();
    asq_arg.coor[3] = CoorZ();

#if 0
    do_arg.x_bc = BND_CND_PRD;
    do_arg.y_bc = BND_CND_PRD;
    do_arg.z_bc = BND_CND_PRD;
    do_arg.t_bc = BND_CND_APRD;
#endif

    asq_arg.c1 = (1.0/8.0)+(6.0/16.0)+(1.0/8.0);
    asq_arg.c2 = -1.0/24.0;
    asq_arg.c3 = (1.0/8.0)*0.5;
    asq_arg.c5 = ( 1.0/8.0)*0.25*0.5;
    asq_arg.c7 = (1.0/8.0)*0.125*(1.0/6.0);
    asq_arg.c6 = -1.0/16;
	asq_arg.Fat=NULL;
	asq_arg.Naik=NULL;
	asq_arg.NaikM=NULL;

    InvArg inv_arg;

    inv_arg.mass = 1.;
    inv_arg.stop_rsd = 1e-12;
    inv_arg.evenodd = 0;
    inv_arg.niter = 1000;

    printf("total sites = %d %d %d %d\n",nx,ny,nz,nt);
    
    fp = Fopen(ADD_ID,"asqd_test.out","w");

    AsqD dirac;
    dirac.init(&asq_arg);

    int vec_len = dirac.Vol()*6*sizeof(Float);
    int gauge_len = dirac.Vol()*72*sizeof(Float);
    Float *result = (Float *)qalloc(QCOMMS, vec_len);
    bzero((char *)result,vec_len);
    Float *X_out = (Float *)qalloc(QCOMMS, vec_len);
    bzero((char *)X_out,vec_len);
    Float *X_out2 = (Float *)qalloc(QCOMMS, vec_len);
    Float *X_in = (Float *)qalloc(QCOMMS, vec_len);
    bzero((char *)X_in,vec_len);

    int s[4];

//    Matrix *gf = lat.GaugeField();
//    IFloat *gf_p = (IFloat *)lat.GaugeField();
    Float *fatlink =  (Float *)qalloc(QCOMMS, gauge_len);
    Float *longlink = (Float *)qalloc(QCOMMS, gauge_len);
    bzero((char *)fatlink,gauge_len);
    bzero((char *)longlink,gauge_len);

    int offset[4], odd=0;
    for(int i=0;i<4;i++){
    offset[i] = asq_arg.size[i]*asq_arg.coor[i];
    odd += offset[i];
    }
    odd = odd%2;
 



    for(s[3]=0; s[3]<asq_arg.size[3]; s[3]++)
	for(s[2]=0; s[2]<asq_arg.size[2]; s[2]++)
	for(s[1]=0; s[1]<asq_arg.size[1]; s[1]++)
	for(s[0]=0; s[0]<asq_arg.size[0]; s[0]++) {
		int n = f_offset(asq_arg.size,s,odd);
		int g_coor[4];
 		for(int i = 0;i<4;i++) g_coor[i] = s[i]+offset[i];
		
	    IFloat crd = 1.0*g_coor[0]
			+0.1*g_coor[1]
			+0.01*g_coor[2]
			+0.001*g_coor[3];
#if 1
	    if((offset[0]+s[0])==pos[3] &&
	    (offset[1]+s[1])==pos[0] &&
	    (offset[2]+s[2])==pos[1] &&
	    (offset[3]+s[3])==pos[2] )
        crd=1.0; else crd = 0.0;
#endif				
	    for(int v=0; v<6; v++){ 
			if (v==0)
			    *(X_in+6*n+v) = crd;
//			    *(X_in+6*n+v) = 0.;
			else
			    *(X_in+6*n+v) = 0.;
	    }

	    crd = 1.0 +1e-1*g_coor[0]
			+1e-2*g_coor[1]
			+1e-3*g_coor[2]
			+1e-4*g_coor[3];
	    for(int dir=0; dir<4; dir++){ 
		    n = g_offset(asq_arg.size,s);
			n += ((dir+1)%4)*dirac.Vol();
		    for(int v=0; v<18; v++){ 
				if (v%8==0){
#if 0
		  			*(fatlink+18*n+v) = crd+1e-5*(dir);
			    	*(longlink+18*n+v) = crd+1e-5*(dir+4);
#else
		  			*(fatlink+18*n+v) = 1.;
			    	*(longlink+18*n+v) = 1.;
#endif
				} else{
			   		*(fatlink+18*n+v) = 0.;
			    	*(longlink+18*n+v) = 0.;
				}
			}
		}
	}

    double maxdiff =0.;
    Float *out;
//    AsqD dirac(lat,X_out->Vec(),X_in->Vec(),&cg_arg,CNV_FRM_NO);
    Float *frm_tmp = (Float *)qalloc(QCOMMS, vec_len);
    dirac.init_g(frm_tmp,fatlink,longlink, NULL);

    for(int k = 0; k< 1; k++){
	printf("k=%d ",k);
	    out = result;
	bzero((char *)out, dirac.Vol()*6*sizeof(Float));
    Float mass_sq = inv_arg.mass*inv_arg.mass*4;
#if 1
        double true_res;
	int iter = dirac.InvCg(&inv_arg,result,X_in,&true_res);
	dirac.MdagM(&mass_sq,X_out,result);
	printf("iter=%d\n",iter);
#else
	dirac.Dslash(result,X_in);
//	dirac.MdagM(&mass_sq,result,X_in);
#endif

#if 0
	    bzero((char *)X_out2, dirac.Vol()*6*sizeof(IFloat));
	    dirac.Dslash(X_out2->Vec(),out->Vec(offset),CHKB_ODD,DAG_NO);
	    dirac.Dslash(X_out2->Vec(offset),out->Vec(),CHKB_EVEN,DAG_NO);
	    lat.Fconvert(X_out2->Vec(),CANONICAL,STAG);
	//X_out2->FTimesV1PlusV2(2*cg_arg.mass,out,X_out2);
#endif
    
	Float dummy;
	Float dt = 2;

    
	for(s[3]=0; s[3]<asq_arg.size[3]; s[3]++) 
	for(s[2]=0; s[2]<asq_arg.size[2]; s[2]++)
	for(s[1]=0; s[1]<asq_arg.size[1]; s[1]++)
	for(s[0]=0; s[0]<asq_arg.size[0]; s[0]++) {
	    int n = f_offset(asq_arg.size,s,odd);
		Float tmp1,tmp2,tmp3,tmp4;
	    for(int i=0; i<3; i++){
		    for(int mu = 0;mu<4;mu++)
		    Fprintf(fp," %d",offset[mu]+s[mu]);
		    Fprintf(fp," %d",i);
		    Fprintf(fp,"  (%0.7e %0.7e) (%0.3e %0.3e)",
		    tmp1=*(out+n*6+i*2), tmp2=*(out+n*6+i*2+1),
		    tmp3=*(X_in+n*6+i*2), tmp4=*(X_in+n*6+i*2+1) );
		    Fprintf(fp,"\n");
        if ((tmp1*tmp1+tmp2*tmp2+tmp3*tmp3+tmp4*tmp4)>1e-10){
		    for(int mu = 0;mu<4;mu++)
		    Fprintf(stderr," %d",offset[mu]+s[mu]);
		    Fprintf(stderr," %d",i);
		    Fprintf(stderr,"  (%0.7e %0.7e) (%0.3e %0.3e)",
		    tmp1=*(out+n*6+i*2), tmp2=*(out+n*6+i*2+1),
		    tmp3=*(X_in+n*6+i*2), tmp4=*(X_in+n*6+i*2+1) );
		    Fprintf(stderr," (%0.2e %0.2e)",
		    *(X_out+n*6+i*2),
		    *(X_out+n*6+i*2+1) );
		    Fprintf(stderr,"\n");
		}
#if 0
		    double diff = (*(X_out2->Field(n,0,i*2)))-*(X_in->Field(n,0,i*2));
		    if (fabs(diff)>maxdiff){
				 maxdiff = fabs(diff);
			}
		    diff = (*(X_out2->Field(n,0,i*2+1)))-*(X_in->Field(n,0,i*2+1));
		    if (fabs(diff)>maxdiff){
				maxdiff = fabs(diff);
			}
#endif
	    }
	}
    }
    Fclose(fp);
//    printf("Max diff between X_in and M*X_out = %0.2e\n", maxdiff);
    
    qfree( X_in);
    qfree( result);
    qfree( X_out);
    qfree( X_out2);
//    End();
    return 0; 
}
