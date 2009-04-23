#include<config.h>
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2009-04-23 03:33:25 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/f_p4_inv/main.C,v 1.5 2009-04-23 03:33:25 chulwoo Exp $
//  $Id: main.C,v 1.5 2009-04-23 03:33:25 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: main.C,v $
//  $Revision: 1.5 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/f_p4_inv/main.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------


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
#include<alg/alg_rnd_gauge.h>
#include<alg/do_arg.h>
#ifdef HAVE_STRINGS_H
#include<strings.h>
#endif
#if TARGET == QCDOC
#include <qalloc.h>
extern "C" void _mcleanup(void);
#endif



USING_NAMESPACE_CPS

int main(int argc,char *argv[]){

    Start(&argc, &argv);
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
    }
    sscanf(argv[1],"%d",&nx);
    sscanf(argv[2],"%d",&ny);
    sscanf(argv[3],"%d",&nz);
    sscanf(argv[4],"%d",&nt);
    printf("size(int)=%d size(long)=%d\n",sizeof(int),sizeof(long));

    do_arg.x_node_sites = nx/SizeX();
    do_arg.y_node_sites = ny/SizeY();
    do_arg.z_node_sites = nz/SizeZ();
    do_arg.t_node_sites = nt/SizeT();
    do_arg.x_nodes = SizeX();
    do_arg.y_nodes = SizeY();
    do_arg.z_nodes = SizeZ();
    do_arg.t_nodes = SizeT();
    do_arg.x_bc = BND_CND_PRD;
    do_arg.y_bc = BND_CND_PRD;
    do_arg.z_bc = BND_CND_PRD;
    do_arg.t_bc = BND_CND_APRD;

    #if TARGET==QCDOC
    do_arg.start_conf_alloc_flag = QFAST;
    #endif

    do_arg.start_conf_kind = START_CONF_DISORD;
//    do_arg.start_conf_kind = START_CONF_ORD;
    do_arg.start_seed_kind = START_SEED_FIXED;
    do_arg.beta = 5.5;
    do_arg.dwf_height = 0.9;
    do_arg.clover_coeff = 2.0171;

//    do_arg.p4_KS = 1.;
    do_arg.p4_KS = 3./4.;
//    do_arg.p4_KS = 1e-15;
//    do_arg.p4_knight = 1.;
    do_arg.p4_knight = 1.0/48.0;
//    do_arg.p4_knight = 1e-100;
    do_arg.p4_3staple = (-1.0/8.0)*0.5;
//    do_arg.p4_3staple = 1e-100;//(-1.0/8.0)*0.5;
    do_arg.p4_5staple = 1e-100;//(1.0/8.0)*.25*.5;
    do_arg.p4_7staple = 1e-100;//(-1.0/8.0)*.125*(1.0/6.0);
    do_arg.p4_lepage = 1e-100;//-1.0/16.0;
    
    CgArg cg_arg;

    cg_arg.mass = 0.1 ;
    cg_arg.stop_rsd = 1e-8;
    cg_arg.max_num_iter = 5000;

#if TARGET==cpsMPI
    MPISCU::set_pe_grid(do_arg.x_nodes, do_arg.y_nodes, do_arg.z_nodes, do_arg.t_nodes);    
    using MPISCU::fprintf;
    using MPISCU::printf;
#endif

    printf("total sites = %d %d %d %d\n",nx,ny,nz,nt);
    
    GJP.Initialize(do_arg);

    VRB.Level(0);
    VRB.ActivateLevel(VERBOSE_RNGSEED_LEVEL);
    VRB.ActivateLevel(VERBOSE_RESULT_LEVEL);
//    VRB.ActivateLevel(VERBOSE_FUNC_LEVEL);
//    VRB.ActivateLevel(VERBOSE_FLOW_LEVEL);
    VRB.ActivateLevel(VERBOSE_INPUT_LEVEL);
    VRB.ActivateLevel(VERBOSE_DEBUG_LEVEL);
//    VRB.DeactivateAll();

    FILE *fp2 = Fopen(ADD_ID,"f_p4_inv.out","w");
    fp = fp2;
//     fp = stdout;

{
    GnoneFp4 lat;

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
    IFloat crd;

    for(s[3]=0; s[3]<GJP.NodeSites(3); s[3]++)
	for(s[2]=0; s[2]<GJP.NodeSites(2); s[2]++)
	    for(s[1]=0; s[1]<GJP.NodeSites(1); s[1]++)
		for(s[0]=0; s[0]<GJP.NodeSites(0); s[0]++) {

		    long n = lat.FsiteOffset(s);
      //             printf("lat.GsiteOffset(s)=%d\n",lat.GsiteOffset(s));
		   Matrix *g = gf+lat.GsiteOffset(s);
                    double index = s[0];
		    for(int v=1; v<4; v++) index = index*0.1+s[v];
		    for(int v=0; v<4; v++){ 
//                       *(g+v) = (Float)(index*0.1+v);
			Float *temp_p = (Float *)(g+v);
//			printf("index=%d g+v(%p)=%e\n",index*10+v,g+v,*(temp_p));
//                       *(g+v) = s[0]+v+1;
                    }

		    crd=0.;
#if PARALLEL
		    if(CoorX()==0 && CoorY()==0 && CoorZ()==0 && CoorT()==0)
#endif
		    if (s[0]==3 &&s[1]==3 &&s[2]==3&&s[3]==3)  crd = 1.0; 
					
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
{
    DiracOpP4 dirac(lat,X_out->Vec(),X_in->Vec(),&cg_arg,CNV_FRM_NO);

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
#if 1
	int vol = GJP.VolNodeSites();
	dtime = -dclock();
	int iter = dirac.MatInv(out->Vec(),X_in->Vec());
	dtime +=dclock();
	print_flops(DiracOp::CGflops,dtime);
	printf("iter=%d\n",iter);
#else
	dirac.Dslash(out->Vec(),X_in->Vec(offset),CHKB_ODD,DAG_NO);
	dirac.Dslash(out->Vec(offset),X_in->Vec(),CHKB_EVEN,DAG_NO);
#endif

#if 1
	if (k == 0){
	    bzero((char *)X_out2->Vec(), X_out2->Size()*sizeof(IFloat));
	    dirac.Dslash(X_out2->Vec(),out->Vec(offset),CHKB_ODD,DAG_NO);
	    dirac.Dslash(X_out2->Vec(offset),out->Vec(),CHKB_EVEN,DAG_NO);
	    lat.Fconvert(X_out2->Vec(),CANONICAL,STAG);
	}
#endif
	lat.Fconvert(out->Vec(),CANONICAL,STAG);
	lat.Fconvert(X_in->Vec(),CANONICAL,STAG);
	X_out2->FTimesV1PlusV2(2*cg_arg.mass,out,X_out2);


	for(s[3]=0; s[3]<GJP.TnodeSites(); s[3]++) 
	for(s[2]=0; s[2]<GJP.ZnodeSites(); s[2]++)
	for(s[1]=0; s[1]<GJP.YnodeSites(); s[1]++)
	for(s[0]=0; s[0]<GJP.XnodeSites(); s[0]++) {
	    int n = lat.FsiteOffset(s);
	    //invert gauge transformation on the output field
	    for(int i=0; i<1; i++){
		if ( k==0 ){
		  if(fabs(*(out->Field(n,0,2*i)))> 1e-12)
		  {
		    for(int mu = 0;mu<4;mu++)
		    Fprintf(fp," %d",GJP.NodeCoor(mu)*GJP.NodeSites(mu)+s[mu]);
		    Fprintf(fp," %d",i);
		    /*Fprintf(fp,"  (%0.7e %0.7e) (%0.7e %0.7e)",
		    *(out->Field(n,0,i*2)), *(out->Field(n,0,i*2+1)),
		    *(X_in->Field(n,0,i*2)), *(X_in->Field(n,0,i*2+1)) );*/
		    Fprintf(fp,"  (%0.7e %0.7e)",
		     *(out->Field(n,0,i*2)), *(out->Field(n,0,i*2+1)));
#if 0
		    Fprintf(fp,"\n");
#else
		    Fprintf(fp," (%0.2e %0.2e)\n",
		    *(X_out2->Field(n,0,i*2))-*(X_in->Field(n,0,i*2)), 
		    *(X_out2->Field(n,0,i*2+1))-*(X_in->Field(n,0,i* 2+1)));
#endif
		    }
		    double diff = *(X_out2->Field(n,0,i*2))-*(X_in->Field(n,0,i*2));
		    if (fabs(diff)>maxdiff) maxdiff = fabs(diff);
		    diff = *(X_out2->Field(n,0,i*2+1))-*(X_in->Field(n,0,i*2+1));
		    if (fabs(diff)>maxdiff) maxdiff = fabs(diff);
		}
	    }
	}
    }
}
    Fclose(fp2);
    printf("Max diff between X_in and M*X_out = %0.2e\n", maxdiff);
    
    delete X_in;
    delete result;
    delete X_out;
    delete X_out2;
}
    End();
    return 0; 
}
