#include<config.h>
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-02-09 14:30:07 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/asqtad_update/main.C,v 1.3 2004-02-09 14:30:07 zs Exp $
//  $Id: main.C,v 1.3 2004-02-09 14:30:07 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: main.C,v $
//  $Revision: 1.3 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/asqtad_update/main.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------


#include <stdio.h>
#include<util/lattice.h>
#include<util/gjp.h>
#include<util/verbose.h>
#include<util/error.h>
#include<alg/alg_hmd.h>
#include<alg/do_arg.h>
#include<alg/common_arg.h>
#include<alg/cg_arg.h>
#include<alg/hmd_arg.h>
#include<alg/ghb_arg.h>

CPS_START_NAMESPACE
GlobalJobParameter GJP;
LatRanGen LRG;
Verbose VRB;
Error ERR;
CPS_END_NAMESPACE

USING_NAMESPACE_CPS


int main(int argc,char *argv[]){

    //----------------------------------------------------------------
    // Initializes all Global Job Parameters
    //----------------------------------------------------------------
    DoArg do_arg;

    do_arg.x_node_sites = 2;
    do_arg.y_node_sites = 2;
    do_arg.z_node_sites = 2;
    do_arg.t_node_sites = 2;
    do_arg.s_node_sites = 0;
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
    do_arg.t_bc = BND_CND_PRD;
    do_arg.start_conf_kind = START_CONF_ORD;
    do_arg.start_seed_kind = START_SEED_FIXED;
    do_arg.colors = 3;
    do_arg.beta = 5.5;
    do_arg.dwf_height = 0.9;
    do_arg.clover_coeff = 2.0171;
    do_arg.verbose_level = DEFAULT_VERBOSE_LEVEL;

    
    do_arg.asqtad_KS = (1.0/8.0)+(6.0/16.0)+(1.0/8.0);
    do_arg.asqtad_naik = -1.0/24.0;
    do_arg.asqtad_lepage = -1.0/16;
    do_arg.asqtad_3staple = (-1.0/8.0)*0.5;
    do_arg.asqtad_5staple = ( 1.0/8.0)*0.25*0.5;
    do_arg.asqtad_7staple = (-1.0/8.0)*0.125*(1.0/6.0);

//     do_arg.asqtad_KS = 0.0;
//     do_arg.asqtad_naik = 0.0;
//     do_arg.asqtad_3staple = 0.0;
//     do_arg.asqtad_lepage = 0.0; 
//     do_arg.asqtad_5staple = 0.0;
//    do_arg.asqtad_7staple = 0.0;
    
    GJP.Initialize(do_arg);

    VRB.Level(GJP.VerboseLevel());


    GnoneFasqtad lat;

    Matrix *mom =
	(Matrix*)smalloc(GJP.VolNodeSites()*lat.GsiteSize()*sizeof(IFloat));
    if(!mom) ERR.Pointer("","","mom");

    Vector *X =
	(Vector*)smalloc(GJP.VolNodeSites()*lat.FsiteSize()*sizeof(IFloat));
    if(!X) ERR.Pointer("","","X");

    Matrix *gf = lat.GaugeField();

    int s[4];
    for(s[3]=0; s[3]<GJP.NodeSites(3); s[3]++)
	for(s[2]=0; s[2]<GJP.NodeSites(2); s[2]++)
	    for(s[1]=0; s[1]<GJP.NodeSites(1); s[1]++)
		for(s[0]=0; s[0]<GJP.NodeSites(0); s[0]++) {

		    int n = lat.FsiteOffset(s);

		    IFloat crd = 1.0*s[0]+0.1*s[1]+0.01*s[2]+0.001*s[3];
					
		    for(int v=0; v<6; v++) *((IFloat*)&X[n]+v) = crd;

		    for(int d=0; d<4; d++){
			Complex icrd(0, crd+0.0001*d);
			*(mom+4*n+d) = *(gf+4*n+d) = icrd;
			(*(mom+4*n+d))(2,2) = -2.0*icrd;
		    }
		    if(s[3]==GJP.NodeSites(3)-1)  // MILC bound. conds.
			*(gf+4*n+3) *= -1.0;
		}


    Float dummy;
    Float dt = 2;
    
    lat.EvolveMomFforce(mom, X, dummy, dt);


    FILE *fp;
    if( (fp = fopen("update.dat", "a")) == NULL ) 
	ERR.FileA(" ","main", "info.dat");
    

    fprintf(fp, " x y z t\n");
    
    for(s[3]=0; s[3]<GJP.NodeSites(3); s[3]++) 
	for(s[2]=0; s[2]<GJP.NodeSites(2); s[2]++)
	    for(s[1]=0; s[1]<GJP.NodeSites(1); s[1]++)
		for(s[0]=0; s[0]<GJP.NodeSites(0); s[0]++) {

		    fprintf(fp, "\n %d %d %d %d", s[0], s[1], s[2], s[3]);
		    int n = lat.FsiteOffset(s);

		    for(int d=0; d<4; d++){
			fprintf(fp, "\n         %d\n\n", d);
			for(int i=0; i<3; i++){
			    for(int j=0; j<3; j++)
				fprintf(fp, " (%+7e %+7e)",
				       (*(mom+4*n+d))(i,j).real(),
				       (*(mom+4*n+d))(i,j).imag());
			    fprintf(fp, "\n");
			}
		    }
		}
    
    
    return 0; 
}






  


