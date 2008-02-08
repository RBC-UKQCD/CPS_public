/*
  $Id: main.C,v 1.19 2008-02-08 18:35:08 chulwoo Exp $
*/

#include<config.h>
#include <util/qcdio.h>
#include <math.h>
#include<util/lattice.h>
#include<util/gjp.h>
#include<util/verbose.h>
#include<util/lat_data.h>
#include<util/dirac_op.h>
#include<util/error.h>
#include<util/time_cps.h>
#include<comms/scu.h>
#include<comms/sysfunc_cps.h>
#include<alg/alg_hmd.h>
#include<alg/do_arg.h>
#ifdef HAVE_STRINGS_H
#include<strings.h>
#endif
#if TARGET == QCDOC
#include <qalloc.h>
extern "C" void _mcleanup(void);
#endif

#define CG


USING_NAMESPACE_CPS

static const char *f_asqtad_test_filename = "f_asqtad_test";
static const char *psi_filename = "psi";
static const char *input_filename = "f_asqtad_inv.in";


int main(int argc,char *argv[]){

    Start();
    printf("Sizes = %d %d %d %d %d %d\n",SizeX(),SizeY(),SizeZ(),SizeT(),SizeS(),SizeW());
    printf("Coors = %d %d %d %d %d %d\n",CoorX(),CoorY(),CoorZ(),CoorT(),CoorS(),CoorW());

    FILE *fp;
    double dtime;

    //----------------------------------------------------------------
    // Initializes all Global Job Parameters
    //----------------------------------------------------------------
    DoArg do_arg;
    int nx,ny,nz,nt;

    if (argc < 5){
        ERR.General("f_asqtad_test","main()","usage: %s nx ny nz nt\n",argv[0]);
	exit(4);
    }
    sscanf(argv[1],"%d",&nx);
    sscanf(argv[2],"%d",&ny);
    sscanf(argv[3],"%d",&nz);
    sscanf(argv[4],"%d",&nt);

    int size[4];
    int sites[4];
    size[0] = SizeX(); sites[0] =nx;
    size[1] = SizeY(); sites[1] =ny;
    size[2] = SizeZ(); sites[2] =nz;
    size[3] = SizeT(); sites[3] =nt;
    for (int i = 0;i<4;i++)
      if( (size[i]*2) > sites[i] ) size[i] = sites[i]/2;
    for (int i = 0;i<4;i++)
      if( sites[i] % size[i] !=0) size[i] = 1;
    do_arg.x_node_sites = nx/size[0];
    do_arg.y_node_sites = ny/size[1];
    do_arg.z_node_sites = nz/size[2];
    do_arg.t_node_sites = nt/size[3];
    do_arg.x_nodes = size[0];
    do_arg.y_nodes = size[1];
    do_arg.z_nodes = size[2];
    do_arg.t_nodes = size[3];

    do_arg.x_bc = BND_CND_PRD;
    do_arg.y_bc = BND_CND_PRD;
    do_arg.z_bc = BND_CND_PRD;
    do_arg.t_bc = BND_CND_APRD;

    do_arg.start_conf_kind = START_CONF_DISORD;

    do_arg.start_seed_kind = START_SEED_FIXED;
    do_arg.beta = 5.5;
    do_arg.dwf_height = 0.9;
    do_arg.clover_coeff = 2.0171;

    do_arg.asqtad_KS = (1.0/8.0)+(6.0/16.0)+(1.0/8.0);
    do_arg.asqtad_naik = -1.0/24.0;
    do_arg.asqtad_3staple = (1.0/8.0)*0.5;
    do_arg.asqtad_5staple = ( 1.0/8.0)*0.25*0.5;
    do_arg.asqtad_7staple = (1.0/8.0)*0.125*(1.0/6.0);
    do_arg.asqtad_lepage = -1.0/16;

#if 0
    do_arg.asqtad_naik = 1e-16;
    do_arg.asqtad_KS = 1.;
    do_arg.asqtad_3staple = 1e-16;
    do_arg.asqtad_5staple = 1e-16;
    do_arg.asqtad_7staple = 1e-16;
    do_arg.asqtad_lepage = 1e-16;
#endif
    
    CgArg cg_arg;

    cg_arg.mass = 0.1;
    cg_arg.stop_rsd = 1e-12;
    cg_arg.max_num_iter = 500;

#if TARGET==cpsMPI
    MPISCU::set_pe_grid(do_arg.x_nodes, do_arg.y_nodes, do_arg.z_nodes, do_arg.t_nodes);    
    using MPISCU::fprintf;
    using MPISCU::printf;
#endif

    printf("total sites = %d %d %d %d\n",nx,ny,nz,nt);
    
    GJP.Initialize(do_arg);

    VRB.Level(0);
//    VRB.ActivateLevel(VERBOSE_SMALLOC_LEVEL);
//    VRB.ActivateLevel(VERBOSE_FLOW_LEVEL);
//    VRB.ActivateLevel(VERBOSE_FUNC_LEVEL);
    VRB.ActivateLevel(VERBOSE_RNGSEED_LEVEL);

    fp = Fopen(ADD_ID,"f_asqtad_test.out","w");

    GwilsonFasqtad lat;

    LatVector *result = new LatVector(1); // 1 = number of spinors per site
    LatVector *X_out  = new LatVector(1);
    LatVector *X_out2  = new LatVector(1);
    LatVector *X_in  = new LatVector(1);

    int s[4];
#if 1
    lat.RandGaussVector(X_in->Vec(),1.0);
#else

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
			    *(X_in->Field(n,0,v)) = crd;
			else
			    *(X_in->Field(n,0,v)) = 0.;
		        *(X_in->Field(n,0,v+1)) = 0.;
		    }
		}
#endif

    double maxdiff =0.;
    LatVector *out;
    DiracOpAsqtad dirac(lat,X_out->Vec(),X_in->Vec(),&cg_arg,CNV_FRM_NO);

    for(int k = 0; k< 1; k++){
	printf("k=%d ",k);
	if (k ==0)
	    out = result;
	else
	    out = X_out;
	bzero((char *)out->Field(), out->Size()*sizeof(IFloat));
	lat.Fconvert(out->Vec(),STAG,CANONICAL);
	lat.Fconvert(X_in->Vec(),STAG,CANONICAL);
	int offset = GJP.VolNodeSites()/2;
#ifdef CG
	int vol = GJP.VolNodeSites();
	dtime = -dclock();
	int iter = dirac.MatInv(out->Vec(),X_in->Vec());
	dtime +=dclock();
	print_flops(1182*iter*vol,dtime);
	printf("iter=%d\n",iter);
#else
	dirac.Dslash(out->Vec(),X_in->Vec(offset),CHKB_ODD,DAG_NO);
	dirac.Dslash(out->Vec(offset),X_in->Vec(),CHKB_EVEN,DAG_NO);
#endif

	if (k == 0){
	    bzero((char *)X_out2->Vec(), X_out2->Size()*sizeof(IFloat));
	    dirac.Dslash(X_out2->Vec(),out->Vec(offset),CHKB_ODD,DAG_NO);
	    dirac.Dslash(X_out2->Vec(offset),out->Vec(),CHKB_EVEN,DAG_NO);
	    lat.Fconvert(X_out2->Vec(),CANONICAL,STAG);
	}
	lat.Fconvert(out->Vec(),CANONICAL,STAG);
	lat.Fconvert(X_in->Vec(),CANONICAL,STAG);
	X_out2->FTimesV1PlusV2(2*cg_arg.mass,out,X_out2);
    
	Float dummy;
	Float dt = 2;

    
	for(s[3]=0; s[3]<GJP.NodeSites(3); s[3]++) 
	for(s[2]=0; s[2]<GJP.NodeSites(2); s[2]++)
	for(s[1]=0; s[1]<GJP.NodeSites(1); s[1]++)
	for(s[0]=0; s[0]<GJP.NodeSites(0); s[0]++) {
	    int n = lat.FsiteOffset(s);
	    for(int i=0; i<3; i++){
		if ( k==0 ){
		    for(int mu = 0;mu<4;mu++)
		    Fprintf(fp," %d",GJP.NodeCoor(mu)*GJP.NodeSites(mu)+s[mu]);
		    Fprintf(fp," %d",i);
		    Fprintf(fp,"  (%0.7e %0.7e) (%0.7e %0.7e)",
		    *(out->Field(n,0,i*2)), *(out->Field(n,0,i*2+1)),
		    *(X_in->Field(n,0,i*2)), *(X_in->Field(n,0,i*2+1)) );
#if 1
		    Fprintf(fp,"\n");
#else
		    Fprintf(fp," (%0.2e %0.2e)\n",
		    (*(X_out2->Field(n,0,i*2))),
		    (*(X_out2->Field(n,0,i*2+1))) );
//		    (*(X_out2->Field(n,0,i*2)))-(*(X_in->Field(n,0,i*2))), 
//		    (*(X_out2->Field(n,0,i*2+1)))-*(X_in->Field(n,0,i* 2+1)) );
#endif
		    double diff = (*(X_out2->Field(n,0,i*2)))-*(X_in->Field(n,0,i*2));
		    if (fabs(diff)>maxdiff){
				 maxdiff = fabs(diff);
			}
		    diff = (*(X_out2->Field(n,0,i*2+1)))-*(X_in->Field(n,0,i*2+1));
		    if (fabs(diff)>maxdiff){
				maxdiff = fabs(diff);
			}
		}
	    }
	}
    }
    Fclose(fp);
    printf("Max diff between X_in and M*X_out = %0.2e\n", maxdiff);
    
    delete X_in;
    delete result;
    delete X_out;
    delete X_out2;
    End();
    return 0; 
}
