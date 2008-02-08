/*
  $Id: main.C,v 1.14 2008-02-08 18:35:08 chulwoo Exp $
*/

#include<config.h>
#include <util/qcdio.h>
#include<util/lattice.h>
#include<util/gjp.h>
#include<util/verbose.h>
#include<util/error.h>
#include<alg/do_arg.h>
#include <comms/sysfunc_cps.h>


USING_NAMESPACE_CPS


int main(int argc,char *argv[]){

    //----------------------------------------------------------------
    // Initializes all Global Job Parameters
    //----------------------------------------------------------------
    DoArg do_arg;

    do_arg.x_node_sites = 4;
    do_arg.y_node_sites = 4;
    do_arg.z_node_sites = 4;
    do_arg.t_node_sites = 4;

#if TARGET==cpsMPI
    using MPISCU::printf;
    using MPISCU::fprintf;
    MPISCU::set_pe_grid(2,2,1,1);    
#endif

    do_arg.x_nodes = SizeX();
    do_arg.y_nodes = SizeY();    
    do_arg.z_nodes = SizeZ();
    do_arg.t_nodes = SizeT();

    do_arg.x_bc = BND_CND_PRD;
    do_arg.y_bc = BND_CND_PRD;
    do_arg.z_bc = BND_CND_PRD;
    do_arg.t_bc = BND_CND_PRD;

    do_arg.start_conf_kind = START_CONF_ORD;
    do_arg.start_seed_kind = START_SEED_FIXED;
    
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
//     do_arg.asqtad_7staple = 0.0;
    

    GJP.Initialize(do_arg);

    GnoneFasqtad lat;

    Matrix *mom =
	(Matrix*)smalloc(GJP.VolNodeSites()*lat.GsiteSize()*sizeof(IFloat));
    if(!mom) ERR.Pointer("","","mom");

    const int degree = 2; 

    Vector *X = (Vector*)smalloc(degree*GJP.VolNodeSites()*lat.FsiteSize()*sizeof(IFloat));
    
    Matrix *gf = lat.GaugeField();

    int s[4];
    for(s[3]=0; s[3]<GJP.NodeSites(3); s[3]++)
	for(s[2]=0; s[2]<GJP.NodeSites(2); s[2]++)
	    for(s[1]=0; s[1]<GJP.NodeSites(1); s[1]++)
		for(s[0]=0; s[0]<GJP.NodeSites(0); s[0]++) {

		    int oe = lat.FsiteOffsetChkb_all(s);
		    int lex = lat.GsiteOffset(s);

		    IFloat crd = 1.0*s[0]+0.1*s[1]+0.01*s[2]+0.001*s[3];

		    for(int d=0; d<degree; d++)
			for(int v=0; v<6; v++)
			    *((IFloat*)(X+d*GJP.VolNodeSites()+oe)+v) = crd;

		    for(int d=0; d<4; d++){
			Complex icrd(0, crd+0.0001*d);
			*(mom+lex+d) = *(gf+lex+d) = icrd;
			(*(mom+lex+d))(2,2) = -2.0*icrd;
		    }
		    if(s[3]==GJP.NodeSites(3)-1)  // MILC bound. conds.
			*(gf+lex+3) *= -1.0;
		}

    Float dummy;
    Float dt = 2;


    Float alpha[degree];
    for(int d=0;  d<degree; d++) alpha[d] = 1.0;

     for(int d=0; d<degree; d++)
  	lat.EvolveMomFforce(mom, X+d*GJP.VolNodeSites(), dummy, dt);

//    lat.RHMC_EvolveMomFforce(mom, X, degree, alpha, dummy, dt, X);
		

    FILE *fp;
    if( (fp = Fopen("update.dat", "w")) == NULL ) 
	ERR.FileA(" ","main", "update.dat");
    

    Fprintf(fp, " x y z t\n");
    
    for(s[3]=0; s[3]<GJP.NodeSites(3); s[3]++) 
	for(s[2]=0; s[2]<GJP.NodeSites(2); s[2]++)
	    for(s[1]=0; s[1]<GJP.NodeSites(1); s[1]++)
		for(s[0]=0; s[0]<GJP.NodeSites(0); s[0]++) {
		    int n = lat.GsiteOffset(s);		    

 		    Fprintf(fp, "\nsite (%d %d %d %d) = %d",
 			    s[0], s[1], s[2], s[3], n/4);

		    for(int d=0; d<4; d++){
			Fprintf(fp, "\n         %d\n\n", d);
			for(int i=0; i<3; i++){
			    for(int j=0; j<3; j++)
				Fprintf(fp, " (%+7e %+7e)",
					(*(mom+n+d))(i,j).real(),
					(*(mom+n+d))(i,j).imag());
			    Fprintf(fp, "\n");
			}
			
		    }
		}
    
    
    return 0; 
}






  


