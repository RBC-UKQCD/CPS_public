#include<config.h>
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-09-21 20:16:50 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/PlaqTest/main.C,v 1.10 2004-09-21 20:16:50 chulwoo Exp $
//  $Id: main.C,v 1.10 2004-09-21 20:16:50 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: main.C,v $
//  $Revision: 1.10 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/PlaqTest/main.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

// Simple plaquette measurement.

#include <util/qcdio.h>
#include<util/lattice.h>
#include<util/gjp.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/random.h>
#include<alg/alg_plaq.h>
#include <alg/do_arg.h>




USING_NAMESPACE_CPS

int main(int argc,char *argv[]) {


    DoArg do_arg;

    do_arg.x_node_sites = 4;
    do_arg.y_node_sites = 4;
    do_arg.z_node_sites = 4;
    do_arg.t_node_sites = 4;
#ifdef PARALLEL
    do_arg.x_nodes = 1;
    do_arg.y_nodes = 1;
    do_arg.z_nodes = 2;
    do_arg.t_nodes = 2;
#else
    do_arg.x_nodes = 1;
    do_arg.y_nodes = 1;
    do_arg.z_nodes = 1;
    do_arg.t_nodes = 1;
#endif 
    do_arg.x_bc = BND_CND_PRD;
    do_arg.y_bc = BND_CND_PRD;
    do_arg.z_bc = BND_CND_PRD;
    do_arg.t_bc = BND_CND_APRD;
    do_arg.start_conf_kind = START_CONF_DISORD;
    do_arg.start_seed_kind = START_SEED_FIXED_UNIFORM;

    
#if TARGET==cpsMPI
    MPISCU::set_pe_grid(do_arg.x_nodes, do_arg.y_nodes,
			do_arg.z_nodes, do_arg.t_nodes);
    MPI_Init(&argc, &argv);
    using MPISCU::printf;
#endif
    
    GJP.Initialize(do_arg);

    VRB.Level(VERBOSE_FLOW_LEVEL);
    VRB.DeactivateLevel(VERBOSE_SMALLOC_LEVEL);    

    GJP.Initialize(do_arg);

    CommonArg common;
    NoArg none;

    common.set_filename("plaqdump");

    GwilsonFnone lat;
    AlgPlaq plaquette(lat, &common, &none);
    plaquette.run();
    printf(" plaquette = %f\n", lat.SumReTrPlaqNode()/(GJP.VolNodeSites()*3.0*6.0));
    return 0;

}





  


