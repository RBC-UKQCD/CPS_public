#include<config.h>

//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-05-10 15:26:55 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/PlaqTest/main.C,v 1.3 2004-05-10 15:26:55 zs Exp $
//  $Id: main.C,v 1.3 2004-05-10 15:26:55 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: main.C,v $
//  $Revision: 1.3 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/PlaqTest/main.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

#include<util/lattice.h>
#include<util/gjp.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/random.h>
#include<alg/alg_plaq.h>
#include <alg/do_arg.h>


CPS_START_NAMESPACE
GlobalJobParameter GJP;
LatRanGen LRG;
Verbose VRB;
Error ERR;
CPS_END_NAMESPACE


USING_NAMESPACE_CPS

int main(int argc,char *argv[])
{

    //----------------------------------------------------------------
    // Initializes all Global Job Parameters
    //----------------------------------------------------------------
    DoArg do_arg;

    do_arg.x_node_sites = 4;
    do_arg.y_node_sites = 4;
    do_arg.z_node_sites = 4;
    do_arg.t_node_sites = 4;
    do_arg.s_node_sites = 1;

#ifdef PARALLEL
    do_arg.x_nodes = 2; 
    do_arg.y_nodes = 2; 
    do_arg.z_nodes = 2; 
    do_arg.t_nodes = 2; 
    do_arg.s_nodes = 1;
#else
    do_arg.x_nodes = 1;
    do_arg.y_nodes = 1;
    do_arg.z_nodes = 1;
    do_arg.t_nodes = 1;
    do_arg.s_nodes = 1;
#endif

    do_arg.x_bc = BND_CND_PRD;
    do_arg.y_bc = BND_CND_PRD;
    do_arg.z_bc = BND_CND_PRD;
    do_arg.t_bc = BND_CND_APRD;
    do_arg.start_conf_kind = START_CONF_DISORD;
    do_arg.start_seed_kind = START_SEED_FIXED;

    VRB.DeactivateAll();
    VRB.Level(VERBOSE_FLOW_LEVEL);
    VRB.DeactivateLevel(VERBOSE_SMALLOC_LEVEL);    

    GJP.Initialize(do_arg);

    CommonArg common;
    NoArg none;

    common.set_filename("plaqdump");

    GwilsonFnone lat;
    AlgPlaq plaquette(lat, &common, &none);
    plaquette.run();

    return 0;

}








