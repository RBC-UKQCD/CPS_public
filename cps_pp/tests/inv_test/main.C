/*
  $Id: main.C,v 1.4 2012-12-05 16:39:19 chulwoo Exp $
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
#include<util/time_cps.h>
#include<comms/scu.h>
#include<alg/alg_hmd.h>
#include<alg/do_arg.h>
#include <alg/qpropw.h>

//#include<alg/common_arg.h>
//#include<alg/no_arg.h>
//#include<alg/alg_plaq.h>

#ifdef USE_QUDA
//#include <alg/quda_arg.h>
#endif

#if TARGET == QCDOC
#include <qalloc.h>
extern "C"{
void _mcleanup(void);
}
#endif

USING_NAMESPACE_CPS


int main(int argc,char *argv[]){

  char *fname = "main()";
  char *cname = "none";

//  double dtime;
#ifdef USE_QUDA
//  QudaArg QudaParam;
  if ( !QudaParam.Decode(argv[4],"QudaParam") ) { printf("Bum quda_arg\n"); exit(-1);}  
  //Start(&argc, &argv, 1);
  //Start(&argc, &argv, 0);
#else
  //Start(&argc,&argv);
#endif

//  QudaParam = quda_arg;

  //----------------------------------------------------------------
  // Initializes all Global Job Parameters
  //----------------------------------------------------------------
  DoArg do_arg;
  int nx,ny,nz,nt;
  CgArg cg_arg;
  //CommandLine::is(argc,argv);

  //if ( argc!=4 ) { 
  //  printf("Args: argc = %d doarg-file cgarg_file initial-directory\n",argc);
  //  exit(-1);
  //}

  chdir(argv[1]);
  if ( !do_arg.Decode(argv[2],"do_arg") ) { printf("Bum do_arg\n"); exit(-1);}
  if ( !cg_arg.Decode(argv[3],"cg_arg") ) { printf("Bum cg_arg\n"); exit(-1);}
  
  Start(&argc, &argv);

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

  //quda_arg.gauge_prec = CUDA_DOUBLE_PRECISION;
  //quda_arg.spinor_prec = CUDA_DOUBLE_PRECISION;
  //quda_arg.reconstruct = CUDA_RECONSTRUCT_NO;
  //quda_arg.gauge_prec_sloppy = CUDA_DOUBLE_PRECISION;
  //quda_arg.spinor_prec_sloppy = CUDA_DOUBLE_PRECISION;
  //quda_arg.reconstruct_sloppy = CUDA_RECONSTRUCT_NO;
  //quda_arg.max_num_iter = 100;
  //quda_arg.tolerance = 1e-8;
  //quda_arg.reliable_delta = 0.1;
  //quda_arg.max_restart = 100;
  //quda_arg.device = 0;
  //printf("Used GPU # : %d\t tolerance : %e \n",QudaParam.device,QudaParam.tolerance);

  GJP.Initialize(do_arg);

  VRB.Level(0);
  VRB.ActivateLevel(VERBOSE_RESULT_LEVEL);

  VRB.ActivateLevel(VERBOSE_FLOW_LEVEL);

  GwilsonFwilson lat;
  
  int x[4] = {1, 0, 0, 0};
  int mu = 0;
  Float *Mat = (Float*)lat.GetLink(x,mu);

  printf("(%d,%d,%d,%d) pt, %d direc. Link\n",x[0],x[1],x[2],x[3],mu);
  for(int m_i = 0; m_i < 3 ; m_i++)
  {
    for(int m_j = 0; m_j <3; m_j++)
    {
      printf("\t%e\t%e",Mat[2*(3*m_i+m_j)],Mat[2*(3*m_i+m_j)+1]);
    }
    printf("\n");
  }

  int f_size = GJP.VolNodeSites() * lat.FsiteSize();
  printf("Lattice site size = %d, Fermion size = %d\n",GJP.VolNodeSites(),lat.FsiteSize());
  
//  Vector *X_out = (Vector*)smalloc(f_size*sizeof(IFloat),"x_out",fname,cname);
//  Vector *X_in = (Vector*)smalloc(f_size*sizeof(IFloat),"X_in",fname,cname);
//  
//  X_out->VecZero(f_size_5d);
//  X_in->VecZero(f_size_5d);

    FermionVectorTp src;
    FermionVectorTp sol;
    src.ZeroSource();
    sol.ZeroSource();
//   src.SetVolSource(color,spin);
    src.SetVolSource(0,0);
//    sol.SetVolSource(0,0);
//    src.SetPointSource(0, mu, x[0], x[1], x[2], x[3]);
//    sol.SetPointSource(0, mu, x[0], x[1], x[2], x[3]);
//    src.SetPointSource(0, 0, 1, 0, 0, 0);
//    src.SetPointSource(0, 0, 0, 1, 0, 0);
//    src.SetPointSource(0, 0, 0, 2, 0, 0);
//    src.SetPointSource(0, 0, 0, 3, 0, 0);
//    src.SetPointSource(0, 0, 0, 4, 0, 0);
  //lat.RandGaussVector(X_in,1.0);
  
  //------------------windy-----------------------
//  int i_x,i_y,i_z,i_t,j,offset;
//  Float *Tmp_in = (Float*)X_in;
//  Float *Tmp_out = (Float*)X_out;;
//  for(i_t = 0; i_t < SizeT(); i_t++)
//  for(i_x = 0; i_x < SizeX(); i_x++)
//  for(i_y = 0; i_y < SizeY(); i_y++)
//  for(i_z = 0; i_z < SizeZ(); i_z++)
//  {
//    for(j = 0; j < 24; j++)
//    {
//      offset = i_x + GJP.XnodeSites()*i_y + GJP.XnodeSites()*GJP.YnodeSites()*i_z +
//               GJP.XnodeSites()*GJP.YnodeSites()*GJP.ZnodeSites()*i_t;
//      if((i_t+i_x+i_y+i_z)%2 == 0)
//      {
//        Tmp_in[offset*6+j] = 1.0;
//        Tmp_out[offset*6+j] = 1.0;
//      }
//      else
//      {
//        Tmp_in[offset*6+j] = 0.0;
//        Tmp_out[offset*6+j] = 0.0;
//      }
//
//    }
//  }
//  
//  for(i_z = 0; i_z < f_size; i_z++)
//  {
//    Tmp_out[i_z] = 1.0;
//  }
//
  //------------------windy-----------------------

  //DiracOpWilson dirac(lat,(Vector*)X_out,(Vector*)X_in,&cg_arg,CNV_FRM_NO);
  DiracOpWilson dirac(lat,(Vector*)sol.data(),(Vector*)src.data(),&cg_arg,CNV_FRM_YES);

  //int iter = dirac.MatInv(X_out,X_in);
  int iter = dirac.MatInv((Vector*)sol.data(),(Vector*)src.data());

//  sfree(X_in);
//  sfree(X_out);
  End();
  return 0; 
}
