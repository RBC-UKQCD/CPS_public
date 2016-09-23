#include<config.h>
<<<<<<< HEAD
=======
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004/09/21 20:16:50 $
//  $Header: /space/cvs/cps/cps++/tests/PlaqTest/main.C,v 1.10 2004/09/21 20:16:50 chulwoo Exp $
//  $Id: main.C,v 1.10 2004/09/21 20:16:50 chulwoo Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $RCSfile: main.C,v $
//  $Revision: 1.10 $
//  $Source: /space/cvs/cps/cps++/tests/PlaqTest/main.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
>>>>>>> 8d384c144ad6ce719f894a1b46153b32abf2e66a

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

    Start(&argc,&argv);
    DoArg do_arg;

    do_arg.x_node_sites = 4;
    do_arg.y_node_sites = 4;
    do_arg.z_node_sites = 4;
    do_arg.t_node_sites = 4;
    do_arg.x_nodes = SizeX();
    do_arg.y_nodes = SizeY();
    do_arg.z_nodes = SizeZ();
    do_arg.t_nodes = SizeT();
    do_arg.x_bc = BND_CND_PRD;
    do_arg.y_bc = BND_CND_PRD;
    do_arg.z_bc = BND_CND_PRD;
    do_arg.t_bc = BND_CND_APRD;
    do_arg.start_conf_kind = START_CONF_DISORD;
    do_arg.start_seed_kind = START_SEED_FIXED_UNIFORM;

    
    
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
    End();
    return 0;

}





  


