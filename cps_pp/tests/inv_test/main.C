/*
  $Id: main.C,v 1.2 2006-02-21 21:14:13 chulwoo Exp $
*/

#include<stdlib.h>
#include<stdio.h>
#include<unistd.h>
#include<config.h>
#include <util/qcdio.h>
#include <math.h>
#include<util/lattice.h>
#include<util/gjp.h>
#include<util/verbose.h>
#include<util/dirac_op.h>
#include<util/error.h>
#include<util/time.h>
#include<comms/scu.h>
#include<alg/alg_hmd.h>
#include<alg/do_arg.h>
#if TARGET == QCDOC
#include <qalloc.h>
extern "C"{
void _mcleanup(void);
}
#endif




USING_NAMESPACE_CPS



int main(int argc,char *argv[]){

  char *fname = "main()";
  char *cname = "blah";

#if TARGET == QCDOC
    DefaultSetup();
    printf("Sizes = %d %d %d %d %d %d\n",SizeX(),SizeY(),SizeZ(),SizeT(),SizeS(),SizeW());
    printf("Coors = %d %d %d %d %d %d\n",CoorX(),CoorY(),CoorZ(),CoorT(),CoorS(),CoorW());
#endif
    FILE *fp;
    double dtime;

    //----------------------------------------------------------------
    // Initializes all Global Job Parameters
    //----------------------------------------------------------------
    DoArg do_arg;
    int nx,ny,nz,nt;
    CgArg cg_arg;

    if ( argc!=5 ) { 
      printf("Args: doarg-file cgarg_file initial-directory\n");
      exit(-1);
    }
    
    chdir(argv[4]);
    
    if ( !do_arg.Decode(argv[1],"do_arg") ) { 
      do_arg.Encode("bum_arg","bum_arg");
      printf("Bum do_arg\n"); 
      exit(-1);
    }
    
    //Layout the lattice on the machine (without regard to even-odd)
    do_arg.x_nodes = SizeX(); 
    do_arg.y_nodes = SizeY();
    do_arg.z_nodes = SizeZ();
    do_arg.t_nodes = SizeT();
    do_arg.s_nodes = SizeS();
    
    do_arg.x_node_sites = do_arg.x_sites/do_arg.x_nodes; 
    do_arg.y_node_sites = do_arg.y_sites/do_arg.y_nodes;
    do_arg.z_node_sites = do_arg.z_sites/do_arg.z_nodes;
    do_arg.t_node_sites = do_arg.t_sites/do_arg.t_nodes;
    do_arg.s_node_sites = do_arg.s_sites/do_arg.s_nodes;
    
    if (do_arg.x_sites!=do_arg.x_node_sites*do_arg.x_nodes) 
      {printf("Lattice does not fit\n");exit(-1);}
    if (do_arg.y_sites!=do_arg.y_node_sites*do_arg.y_nodes) 
      {printf("Lattice does not fit\n");exit(-1);}
    if (do_arg.z_sites!=do_arg.z_node_sites*do_arg.z_nodes) 
      {printf("Lattice does not fit\n");exit(-1);}
    if (do_arg.t_sites!=do_arg.t_node_sites*do_arg.t_nodes) 
      {printf("Lattice does not fit\n");exit(-1);}
    if (do_arg.s_sites!=do_arg.s_node_sites*do_arg.s_nodes) 
      {printf("Lattice does not fit\n");exit(-1);}
    
    if ( !cg_arg.Decode(argv[2],"cg_arg")){printf("Bum cg_arg\n"); exit(-1);}
    
    GJP.Initialize(do_arg);

    VRB.Level(0);
    VRB.ActivateLevel(VERBOSE_RESULT_LEVEL);

    VRB.ActivateLevel(VERBOSE_FLOW_LEVEL);

    GwilsonFclover lat;

    int f_size = GJP.VolNodeSites() * lat.FsiteSize();
    Vector *X_out = 
      (Vector*)smalloc(f_size*sizeof(IFloat),"x_out",fname,cname);
    Vector *X_in =(Vector*)smalloc(f_size*sizeof(IFloat),"X_in",fname,cname);

    lat.RandGaussVector(X_in,1.0);

    DiracOpClover dirac(lat,X_out,X_in,&cg_arg,CNV_FRM_NO);
    
    int iter = dirac.MatInv(X_out,X_in);
    
    sfree(X_in);
    sfree(X_out);
    return 0; 
}
