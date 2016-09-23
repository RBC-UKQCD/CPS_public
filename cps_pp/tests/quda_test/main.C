/*
  $Id: main.C,v 1.3 2008/02/08 18:35:08 chulwoo Exp $
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
#include <util/eigen_container.h>
#if TARGET == QCDOC
#include <qalloc.h>
extern "C"{
void _mcleanup(void);
}
#endif

USING_NAMESPACE_CPS

using namespace std;

// needed to declare globally
std::vector<EigenCache*> cps::EigenCacheList(0);

//Search contents that match to arguments, return 0 if not found
EigenCache* cps::EigenCacheListSearch( char* fname_root_bc, int neig )
{
  EigenCache* ecache=0;

  for(int i=0; i< EigenCacheList.size(); ++i  )
    if( EigenCacheList[i]-> is_cached( fname_root_bc, neig ) )
      ecache = EigenCacheList[i];

  return ecache;
}

// Cleanup list EigenCache, it also destroies contents pointed by the elements.
void cps::EigenCacheListCleanup( )
{
  for(size_t i=0;i< EigenCacheList. size(); ++i){
    EigenCacheList[i]-> dealloc();
  }
  EigenCacheList.clear() ;
}

int main(int argc,char *argv[]){

  char *fname = "main()";
  char *cname = "blah";

#ifdef USE_QUDA
//  QudaArg QudaParam;
  if ( !QudaParam.Decode(argv[5],"QudaParam") ) { printf("Bum quda_arg\n"); exit(-1);}  
  printf("device %d\n", QudaParam.device);
  //Start(&argc, &argv, 1);
  //Start(&argc, &argv, 0);
#else
  //Start(&argc,&argv);
#endif


    //----------------------------------------------------------------
    // Initializes all Global Job Parameters
    //----------------------------------------------------------------
    DoArg do_arg;
    DoArgExt do_ext;
    int nx,ny,nz,nt;
    CgArg cg_arg;

    //if ( argc!=5 ) { 
    //  printf("Args: doarg-file cgarg_file initial-directory\n");
    //  exit(-1);
    //}
    
    chdir(argv[1]);
    if ( !do_arg.Decode(argv[2],"do_arg") ) { printf("Bum do_arg\n"); exit(-1);}
    if ( !do_ext.Decode(argv[3],"do_ext") ) { printf("Bum do_ext\n"); exit(-1);}
    if ( !cg_arg.Decode(argv[4],"cg_arg") ) { printf("Bum cg_arg\n"); exit(-1);}
    
    Start(&argc,&argv);
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
    
    
    
    GJP.Initialize(do_arg);
    GJP.InitializeExt(do_ext);

#if 0
    VRB.Level(0);
    VRB.ActivateLevel(VERBOSE_RESULT_LEVEL);

    VRB.ActivateLevel(VERBOSE_FLOW_LEVEL);
#endif

    GwilsonFdwf lat;
    int x[4] = {0, 0, 0, 0};
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
    int ls_glb = GJP.SnodeSites()*GJP.Snodes();

    int f_size = GJP.VolNodeSites() * lat.FsiteSize()/GJP.SnodeSites();
    int f_size_5d = GJP.SnodeSites() * f_size;
    
    Vector *X_out = (Vector*)smalloc(f_size_5d*sizeof(IFloat),"x_out",fname,cname);
    Vector *X_in = (Vector*)smalloc(f_size_5d*sizeof(IFloat),"X_in",fname,cname);
    Vector *X_tmp = (Vector*)smalloc(f_size_5d*sizeof(IFloat),"X_in",fname,cname);

    X_out->VecZero(f_size_5d);
    X_in->VecZero(f_size_5d);
    X_tmp->VecZero(f_size_5d);
  
    FermionVectorTp src;
    FermionVectorTp sol;
    src.ZeroSource();
    sol.ZeroSource();
    //----------------------------------
    // src.SetVolSource(color,spin)
    //----------------------------------
    src.SetVolSource(0,0);
//    sol.SetVolSource(0,0);
//    src.SetPointSource(0, mu, x[0], x[1], x[2], x[3]);
//    sol.SetPointSource(0, mu, x[0], x[1], x[2], x[3]);
//    src.SetPointSource(0, 0, 1, 0, 0, 0);
//    src.SetPointSource(0, 0, 0, 1, 0, 0);
//    src.SetPointSource(0, 0, 0, 2, 0, 0);
//    src.SetPointSource(0, 0, 0, 3, 0, 0);
//    src.SetPointSource(0, 0, 0, 4, 0, 0);
//    lat.RandGaussVector((Vector*)sol.data(),1.0);
//    lat.RandGaussVector(X_out,1.0);
    
    Vector *src_4d = (Vector*)src.data();
    Vector *sol_4d = (Vector*)sol.data();
//    lat.RandGaussVector(src_4d,1.0);
    //lat.RandGaussVector(sol_4d,1.0);

  // Convert: 4d --> 5d
    //lat.Ffour2five(X_in, src_4d, 0, ls_glb-1);
    //lat.Ffour2five(X_out, sol_4d, 0, ls_glb-1);
    //lat.Ffour2five(X_in, src_4d, 0, ls_glb-1);
    lat.Ffour2five(X_in, src_4d, 0, ls_glb-1);
    lat.Ffour2five(X_out, sol_4d, ls_glb-1, 0);
/*
    srand(1823);

    Float *tmp_v1 = (Float*)X_in;
    Float *tmp_v2 = (Float*)X_out;
    int k,j;
    for(k = 0; k < f_size; k++)
    {
      //printf("site[%d]:\t",k);
      for(j = 0; j <GJP.SnodeSites(); j++ )
      {
        tmp_v1[k + j*f_size] = (double)rand()/(double)RAND_MAX;
        tmp_v2[k + j*f_size] = 0.0;
        //printf("%e\t",tmp_v1[k + j*f_size]);
      }
      //printf("\n");
    }
*/
    //Set Mobius Argument
    MdwfArg mobius_dwf;
    mobius_dwf.M5 = GJP.DwfHeight();
    mobius_dwf.use_single_precision = 0;
    mobius_dwf.use_mdwf_for_dwf = 0;
    const int dwf_ls = GJP.SnodeSites();
    mobius_dwf.b5.b5_len = dwf_ls;
    mobius_dwf.c5.c5_len = dwf_ls;
    mobius_dwf.rsd_vec.rsd_vec_len = dwf_ls;
    mobius_dwf.b5.b5_val = (Float *)smalloc(cname, fname, "b5_val", sizeof(Float)*dwf_ls);
    mobius_dwf.c5.c5_val = (Float *)smalloc(cname, fname, "c5_val", sizeof(Float)*dwf_ls);
    mobius_dwf.rsd_vec.rsd_vec_val = (Float *)smalloc(cname, fname, "rsd_vec", sizeof(Float)*dwf_ls);
    for(int i=0;i<dwf_ls; ++i){
      mobius_dwf.b5.b5_val[i] = GJP.Mobius_b();
      mobius_dwf.c5.c5_val[i] = GJP.Mobius_c();
      mobius_dwf.rsd_vec.rsd_vec_val[i] = 0.0;
    }
    mobius_dwf.cg_arg = cg_arg;
    
    GJP.SetMdwfArg(&mobius_dwf);

    double kappa_b = 1.0/(GJP.Mobius_b()*(4.0-GJP.DwfHeight())+1.0);

    Float *input = (Float*) X_in;
//    for(int k = 0 ; k < f_size_5d/2 ; k++)
//    {
//      printf("In[%d] = %e + j %e\n",k,input[2*k],input[2*k+1]);
//    }

///    lat.Convert( DWF_4D_EOPREC_EE, X_out, X_in); 
//    lat.Convert(CANONICAL, X_out, X_in);

    {
      DiracOpMobius dirac(lat, X_out, X_in, &cg_arg, CNV_FRM_YES);

      printf("MDWF operator is allocated: b_5:%e\t c_5:%e\n",GJP.Mobius_b(),GJP.Mobius_c());
      int iter = dirac.MatInv(X_out,X_in, &cg_arg.stop_rsd , PRESERVE_YES);
    }
#if 0
    dirac.Mat(X_out, X_in);
    dirac.Mat(X_tmp, X_out);

    lat.Convert(CANONICAL, X_out, X_in);
    Float *output = (Float*) X_out;
    Float *recal = (Float*) X_tmp;
//    for(int k = 0 ; k < f_size_5d/2 ; k++)
//    {
//      printf("Out[%d] = %e + j %e\n",k,output[2*k]*kappa_b,output[2*k+1]*kappa_b);
//    }
    for(int k = 0 ; k < f_size_5d/2 ; k++)
    {
      if((fabs(input[2*k] - recal[2*k])>1e-8) || (fabs(input[2*k+1] - recal[2*k+1])>1e-8))
      {
        printf("X_in[%d] = %e + j %e\n",k,input[2*k],input[2*k+1]);
        printf("recon[%d] = %e + j %e\n",k,recal[2*k],recal[2*k+1]);
      }
    }
//    lat.Convert(CANONICAL, X_out, X_in);
#endif

    if(X_in != NULL)
    {
      sfree(X_in); 
      X_in = NULL;
    }
    if(X_out != NULL)
    {
      sfree(X_out); 
      X_out = NULL;
    }
    if(X_tmp != NULL)
    {
      sfree(X_tmp);
      X_tmp = NULL;
    }
    End();
    return 0; 
}
