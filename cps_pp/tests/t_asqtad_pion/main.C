///  \brief
///  Test the asqtad inverter by computing the pion correlator.
///

/*!----------------------------------------------------------------------
  $Id: main.C,v 1.8 2004-09-17 17:56:21 chulwoo Exp $
  Test Asqtad dirac operator code.	

   This is a simple regression test of the Asqtad inverter.
   The pseudo-goldstone pion is computed and compared
   against a known result.

----------------------------------------------------------------------*/



#include<config.h>

#include <util/qcdio.h>
#include <math.h>

#include<util/lattice.h>
#include<util/gjp.h>
#include<util/verbose.h>
#include<util/dirac_op.h>
#include<util/error.h>
#include<alg/alg_hmd.h>
#include<alg/do_arg.h>
#include <util/qcdio.h>
#include<util/testing_framework.h>



CPS_START_NAMESPACE
GlobalJobParameter GJP;
LatRanGen LRG;
Verbose VRB;
Error ERR;
CPS_END_NAMESPACE

USING_NAMESPACE_CPS



int main(int argc,char *argv[])
{
    
    const int nx=4;
    const int ny=4;
    const int nz=4;
    const int nt=4;

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
#else
    do_arg.x_node_sites = nx;
    do_arg.y_node_sites = ny;
    do_arg.z_node_sites = nz;
    do_arg.t_node_sites = nt;
    do_arg.s_node_sites = 0;
    do_arg.x_nodes = 1;
    do_arg.y_nodes = 1;
    do_arg.z_nodes = 1;
#endif 
    do_arg.x_bc = BND_CND_PRD;
    do_arg.y_bc = BND_CND_PRD;
    do_arg.z_bc = BND_CND_PRD;
    do_arg.t_bc = BND_CND_APRD;
    //    do_arg.t_bc = BND_CND_PRD;
    do_arg.start_conf_kind = START_CONF_DISORD;
    do_arg.start_seed_kind = START_SEED_FIXED;
    do_arg.beta = 5.5;
    do_arg.u0 = 1.0  ;

    Float u0 = do_arg.u0 ;
    do_arg.asqtad_KS = (1.0/8.0)+(6.0/16.0)+(1.0/8.0);
    do_arg.asqtad_naik = -1.0/(24.0*pow(u0,2.0))  ;
    do_arg.asqtad_lepage = -1.0/(16) * (1.0/pow(u0,4.0)) ;
    do_arg.asqtad_3staple = (-1.0/8.0)*0.5 * (1.0/pow(u0,2.0)) ;
    do_arg.asqtad_5staple = ( 1.0/8.0)*0.25*0.5 * (1.0/pow(u0,4.0)) ;
    do_arg.asqtad_7staple = (-1.0/8.0)*0.125*(1.0/6.0)* (1.0/(pow(u0,6.0))) ;
    
#if TARGET==cpsMPI
    MPISCU::set_pe_grid(do_arg.x_nodes, do_arg.y_nodes, do_arg.z_nodes, do_arg.t_nodes);
    using MPISCU::printf;
    using MPISCU::fprintf;
#endif

    printf("Computation of staggered pion correlator\n") ; 
    printf("Physical lattice volume [xyzt] = %d %d %d %d\n",nx,ny,nz,nt); 
#if TARGET == QCDOC
    DefaultSetup();
    printf("Sizes = %d %d %d %d %d %d\n",SizeX(),SizeY(),SizeZ(),SizeT(),SizeS(),SizeW());
    printf("Coors = %d %d %d %d %d %d\n",CoorX(),CoorY(),CoorZ(),CoorT(),CoorS(),CoorW());
#endif
    printf("Boundary conditions: [x,y,z,t] = [%d,%d,%d,%d] \n",
	   do_arg.x_bc,do_arg.y_bc,do_arg.z_bc,do_arg.t_bc) ; 
    printf("Tadpole factor = %f\n",do_arg.u0) ;
    printf("HOT gauge configuration created\n") ;

    CgArg cg_arg;

    cg_arg.mass = 0.03;
    printf("Staggered mass parameter = %f\n",cg_arg.mass) ;

    cg_arg.stop_rsd = 1e-12;
    cg_arg.max_num_iter = 1000;

    GJP.Initialize(do_arg);



    // known output
    Float pion_corr_good[] = { 7.563369 , 6.546032 , 5.474610 , 5.255596 } ;
    int good_size = sizeof(pion_corr_good) / sizeof(Float) ; 

    // compute the time sliced pion correlator
    int time_size=GJP.Tnodes()*GJP.TnodeSites();

    if( time_size != good_size )
	ERR.General("", "main", "time_size = %d != good_size = %d\n",
		    time_size,good_size);

    
    Float* pion_corr = (Float*)smalloc(time_size*sizeof(IFloat));

    GwilsonFasqtad lat;
    // char gauge_name[50] = "../gauge/hot";
    // qload_gauge(gauge_name, lat,sizeof(Float),1,1);

    //    char gauge_name[50] = "../gauge/Q5600U003500";
    //qload_gauge(gauge_name, lat,sizeof(double),1,1);

    staggered_local_pion(lat,cg_arg.mass,pion_corr,time_size) ; 


    Float tol = 0.00001 ; 
    compare_array_relative(pion_corr, pion_corr_good,tol,
			   time_size) ;


    sfree(pion_corr); 
    
    return 0; 
}
