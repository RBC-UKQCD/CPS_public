#include<config.h>
#include <stdio.h>
#include<util/lattice.h>
#include<util/gjp.h>
#include<util/verbose.h>
#include<util/dirac_op.h>
#include<util/error.h>
#include<alg/alg_hmd.h>
#include<alg/do_arg.h>
#include<alg/common_arg.h>
#include<alg/cg_arg.h>
#include<alg/hmd_arg.h>
#include<alg/ghb_arg.h>

#include "Gauge_edram_random.h"
#include "Source_edram.h"
CPS_START_NAMESPACE
GlobalJobParameter GJP;
LatRanGen LRG;
Verbose VRB;
Error ERR;
CPS_END_NAMESPACE

USING_NAMESPACE_CPS

const char *f_asqtad_test_filename = CWDPREFIX("f_asqtad_test");
const char *psi_filename = CWDPREFIX("psi");

const int nx=4;
const int ny=4;
const int nz=4;
const int nt=4;


int main(int argc,char *argv[]){

#if TARGET == QCDOC
    DefaultSetup();
    printf("Sizes = %d %d %d %d %d %d\n",SizeX(),SizeY(),SizeZ(),SizeT(),SizeS(),SizeW());
    printf("Coors = %d %d %d %d %d %d\n",CoorX(),CoorY(),CoorZ(),CoorT(),CoorS(),CoorW());
#endif
    printf("start\n");
    FILE *fp, *fp2;

    //----------------------------------------------------------------
    // Initializes all Global Job Parameters
    //----------------------------------------------------------------
    DoArg do_arg;

#if TARGET == QCDOC
    do_arg.x_node_sites = nx/SizeX();
    do_arg.y_node_sites = ny/SizeY();
    do_arg.z_node_sites = nz/SizeZ();
    do_arg.t_node_sites = nt/SizeT();
    do_arg.s_node_sites = 0;
    do_arg.x_nodes = SizeX();
    do_arg.y_nodes = SizeY();
    do_arg.z_nodes = SizeZ();
    do_arg.t_nodes = SizeT();
    do_arg.s_nodes = 1;
#else
    do_arg.x_node_sites = nx;
    do_arg.y_node_sites = ny;
    do_arg.z_node_sites = nz;
    do_arg.t_node_sites = nt;
    do_arg.s_node_sites = 0;
    do_arg.x_nodes = 1;
    do_arg.y_nodes = 1;
    do_arg.z_nodes = 1;
    do_arg.t_nodes = 1;
    do_arg.s_nodes = 1;
#endif 
    do_arg.x_bc = BND_CND_PRD;
    do_arg.y_bc = BND_CND_PRD;
    do_arg.z_bc = BND_CND_PRD;
    do_arg.t_bc = BND_CND_PRD;
#if 1
#if 1
    do_arg.start_conf_kind = START_CONF_ORD;
#else
    do_arg.start_conf_kind = START_CONF_DISORD;
#endif
#else
    do_arg.start_conf_kind = START_CONF_LOAD;
    do_arg.start_conf_load_addr = (Matrix *)Gauge_edram;
#endif
    do_arg.start_seed_kind = START_SEED_FIXED;
    do_arg.colors = 3;
    do_arg.beta = 5.5;
    do_arg.dwf_height = 0.9;
    do_arg.clover_coeff = 2.0171;
//    do_arg.verbose_level = -120403;
    do_arg.verbose_level = -12;

    do_arg.asqtad_KS = (1.0/8.0)+(6.0/16.0)+(1.0/8.0);
    do_arg.asqtad_naik = -1.0/24.0;
//    do_arg.asqtad_naik = 0.;
    do_arg.asqtad_lepage = -1.0/16;
    do_arg.asqtad_3staple = (-1.0/8.0)*0.5;
    do_arg.asqtad_5staple = ( 1.0/8.0)*0.25*0.5;
    do_arg.asqtad_7staple = (-1.0/8.0)*0.125*(1.0/6.0);
    
    CgArg cg_arg;

    cg_arg.mass = 1.0;
    cg_arg.stop_rsd = 1e-12;
    cg_arg.max_num_iter = 1000;
    GJP.Initialize(do_arg);

    VRB.Level(GJP.VerboseLevel());

#if TARGET == QCDOC
    char filename [200];
    sprintf(filename,"%s%0.2d%0.2d%0.2d%0.2d%0.2d%0.2d.out",f_asqtad_test_filename,CoorX(),CoorY(),CoorZ(),CoorT(),CoorS(),CoorW());
   fp = fopen(filename,"w");
    sprintf(filename,"%s%0.2d%0.2d%0.2d%0.2d%0.2d%0.2d.h",psi_filename,CoorX(),CoorY(),CoorZ(),CoorT(),CoorS(),CoorW());
    fp2 = fopen(filename,"w");
#else
    fp = fopen("f_asqtad_test.out","w");
    fp2 = fopen("psi.out","w");
#endif
    fprintf(fp2,"unsigned long psi[]={\n");

    GwilsonFasqtad lat;

	Vector *result = 
	(Vector*)smalloc(GJP.VolNodeSites()*lat.FsiteSize()*sizeof(IFloat));

    Vector *X_out =
	(Vector*)smalloc(GJP.VolNodeSites()*lat.FsiteSize()*sizeof(IFloat));
    if(!X_out) ERR.Pointer("","","X_out");
    Vector *X_out2 =
	(Vector*)smalloc(GJP.VolNodeSites()*lat.FsiteSize()*sizeof(IFloat));
    if(!X_out2) ERR.Pointer("","","X_out2");
    Vector *X_out_odd = &(X_out[GJP.VolNodeSites()/2]);

    int s[4];
#if 0
#if 0
    Vector *X_in = (Vector *) Source_edram;
#else
    Vector *X_in =
	(Vector*)smalloc(GJP.VolNodeSites()*lat.FsiteSize()*sizeof(IFloat));
    if(!X_in) ERR.Pointer("","","X_in");
	lat.RandGaussVector(X_in,1.0);
#endif
#else
    Vector *X_in =
	(Vector*)smalloc(GJP.VolNodeSites()*lat.FsiteSize()*sizeof(IFloat));
    if(!X_in) ERR.Pointer("","","X_in");
    Vector *X_in_odd = &(X_in[GJP.VolNodeSites()/2]);

    Matrix *gf = lat.GaugeField();
    IFloat *gf_p = (IFloat *)lat.GaugeField();

    for(s[3]=0; s[3]<GJP.NodeSites(3); s[3]++)
	for(s[2]=0; s[2]<GJP.NodeSites(2); s[2]++)
	    for(s[1]=0; s[1]<GJP.NodeSites(1); s[1]++)
		for(s[0]=0; s[0]<GJP.NodeSites(0); s[0]++) {

		    int n = lat.FsiteOffset(s);
			IFloat *temp_p = (IFloat *)(gf+4*n+3);

		    IFloat crd = 1.0*s[0]+0.1*s[1]+0.01*s[2]+0.001*s[3];
#if TARGET==QCDOC
		  if(CoorX()==0 && CoorY()==0 && CoorZ()==0 && CoorT()==0 &&n==0) crd=1.0; else crd = 0.0;
#else
	if(n==0) crd = 1.0; else crd = 0.0;
#endif
					
		    for(int v=0; v<6; v+=2){ 
			if (v==0)
			*((IFloat*)&X_in[n]+v) = crd;
			else
			*((IFloat*)&X_in[n]+v) = 0;
			*((IFloat*)&X_in[n]+v+1) = 0.0;
		    }
		}
#endif

    Vector *out;
    DiracOpAsqtad dirac(lat,X_out,X_in,&cg_arg,CNV_FRM_NO);

	for(int k = 0; k< 2; k++){
		fprintf(fp, "k=%d\n",k);
		printf("k=%d ",k);
		if (k ==0)
			out = result;
		else
			out = X_out;
		bzero((char *)out, GJP.VolNodeSites()*lat.FsiteSize()*sizeof(IFloat));
		lat.Fconvert(out,STAG,CANONICAL);
		lat.Fconvert(X_in,STAG,CANONICAL);
   		int iter = dirac.MatInv(out,X_in);
		printf("iter=%d\n",iter);
		if (k == 0){
		int offset = GJP.VolNodeSites()*lat.FsiteSize()/ (2*6);
			bzero((char *)X_out2, GJP.VolNodeSites()*lat.FsiteSize()*sizeof(IFloat));
			dirac.Dslash(X_out2,out+offset,CHKB_ODD,DAG_NO);
			dirac.Dslash(X_out2+offset,out,CHKB_EVEN,DAG_NO);
			lat.Fconvert(X_out2,CANONICAL,STAG);
		}
		lat.Fconvert(out,CANONICAL,STAG);
		lat.Fconvert(X_in,CANONICAL,STAG);
		X_out2->FTimesV1PlusV2(2*cg_arg.mass,out,X_out2,GJP.VolNodeSites
()*lat.FsiteSize());
    
    
    Float dummy;
    Float dt = 2;

if (k==0)    fprintf(fp," x y z t\n");
    
    for(s[3]=0; s[3]<GJP.NodeSites(3); s[3]++) 
	for(s[2]=0; s[2]<GJP.NodeSites(2); s[2]++)
	    for(s[1]=0; s[1]<GJP.NodeSites(1); s[1]++)
		for(s[0]=0; s[0]<GJP.NodeSites(0); s[0]++) {

		    int n = lat.FsiteOffset(s);
			unsigned long  *pt = (unsigned long *)&X_out[n];
			unsigned long  *pt2 = (unsigned long *)&result[n];
			for(int i=0; i<3; i++){
#if TARGET == QCDOC
		    if ( k==0 )
				fprintf(fp," %d %d %d %d %d ", CoorX()*GJP.NodeSites(0)+s[0], CoorY()*GJP.NodeSites(1)+s[1], CoorZ()*GJP.NodeSites(2)+s[2], CoorT()*GJP.NodeSites(3)+s[3], i);
#else
		    if ( k==0 )
				fprintf(fp," %d %d %d %d %d ", s[0], s[1], s[2], s[3], i);
#endif
		    if ( k==0 )
				fprintf(fp," (%0.7e %0.7e) (%0.7e %0.7e) (%0.2e %0.2e)\n",
				*((IFloat*)&result[n]+i*2), *((IFloat*)&result[n]+i*2+1),
				*((IFloat*)&X_in[n]+i*2), *((IFloat*)&X_in[n]+i*2+1),
#if 0
				*((IFloat*)&X_out2[n]+i*2), *((IFloat*)&X_out2[n]+i*2+1));
#else
				*((IFloat*)&X_out2[n]+i*2)-*((IFloat*)&X_in[n]+i*2), *((IFloat*)&X_out2[n]+i*2+1)-*((IFloat*)&X_in[n]+i* 2+1));
#endif
				for (int j = 0; j<4;j++)
				if( k !=0 )
				if(*(pt2+j) != *(pt +j)){
					fprintf(fp, "ERROR! X_out[%d][%d](0x%x) != psi[%d][%d](0x%x)\n",n,i*4+j,*(pt+j),n,i*4+j,*(pt2+j));
					printf("ERROR! X_out[%d][%d](0x%x) != psi[%d][%d](0x%x)\n",n,i*4+j,*(pt+j),n,i*4+j,*(pt2+j));
				}
				pt2 +=4;
		    if ( k==0 )
			    fprintf(fp2,"0x%0.8x, 0x%0.8x, 0x%0.8x, 0x%0.8x,\n", *(pt2),*(pt2+1),*(pt2+2),*(pt2+3));
				pt +=4;
			}
		}
}
    fclose(fp);
		fprintf(fp2,"};\n");
    fclose(fp2);
    
    return 0; 
}
