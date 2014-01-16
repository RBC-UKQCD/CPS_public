#include<config.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------
/*!\file
  \brief Methods of the AlgLanczos class.
  
  $Id: alg_lanczos.C,v 1.2 2013-04-05 17:46:30 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2013-04-05 17:46:30 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_lanczos/alg_lanczos.C,v 1.2 2013-04-05 17:46:30 chulwoo Exp $
//  $Id: alg_lanczos.C,v 1.2 2013-04-05 17:46:30 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: alg_lanczos.C,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_lanczos/alg_lanczos.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/qcdio.h>
#include <util/enum.h>
#include <math.h>
#include <alg/alg_lanczos.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/smalloc.h>
#include <util/time_cps.h>
#include <util/vector.h>
#include <util/verbose.h>
#include <util/error.h>
//#include <util/time_cps.h>
//#include <util/qcdio.h>
//#include <util/qio_writeGenericFields.h>
//#include <util/qio_readGenericFields.h>

#include <util/eigen_container.h>


#include <alg/alg_eig.h>



CPS_START_NAMESPACE

#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

extern void gamma_5(Float *v_out, Float *v_in, int num_sites);

//------------------------------------------------------------------
// Constructor 
/*!
  \param latt The lattice on which to compute the condensate.
  \param c_arg The common argument structure for all algorithms.
  \param arg The algorithm parameters.
 */
//------------------------------------------------------------------
AlgLanczos::AlgLanczos(Lattice& latt, 
		       CommonArg *c_arg,
		       LanczosArg *arg,
		       EigenCache *acache) : 
  Alg(latt, c_arg) 
{
  cname = "AlgLanczos";
  char *fname = "AlgLanczos(L&,CommonArg*,LanczosArg*)";
  VRB.Func(cname,fname);

  // Initialize the argument pointer
  //----------------------------------------------------------------
  if(arg == 0)
    ERR.Pointer(cname,fname, "arg");
  alg_lanczos_arg = arg;

  Ncb = NumChkb(alg_lanczos_arg->RitzMat_lanczos);

  ecache=acache;

  VRB.FuncEnd(cname,fname);
}



//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
AlgLanczos::~AlgLanczos() {
  char *fname = "~AlgLanczos()";
  VRB.Func(cname,fname);

}



void AlgLanczos::run(int init_flag, int ncompress, char* comp_file ){

#if TARGET==cpsMPI
  using MPISCU::fprintf;
#endif
  
  Float time = -dclock();
  int iter=0;
  LanczosArg *lanczos_arg;
  char *fname = "run(int init_flag, int ncompress, char* comp_file)";
  VRB.Func(cname,fname);
  
  // Set the Lattice pointer lanczos_arg
  //----------------------------------------------------------------
  Lattice& lat = AlgLattice(  );
  lanczos_arg = alg_lanczos_arg;

  int nk= lanczos_arg-> nk_lanczos_vectors;
  int np= lanczos_arg-> np_lanczos_vectors;
  int m = nk+np;
  int f_size = GJP.VolNodeSites() * lat.FsiteSize() * Ncb / 2;

#if 0
  if( (ecache=EigenCacheListSearch(cachename))==0 ){
 
    lambda = (Float *) smalloc(cname,fname, "eval", m * sizeof(Float));
    eigenv = (Vector **) smalloc(cname,fname, "evec", m * sizeof(Vector*));
    // We have to use the sequencial area of memory to be consistent with the lanczos routine in the dirac operator.
    int f_totsize=m* f_size; // m*GJP.VolNodeSites() * lattice.FsiteSize() / 2 ;
    Float* tmpalloc= (Float*) smalloc(cname,fname,"evec[...]",
  				    f_totsize* sizeof(Float));
    for(int i=0;i<m;i++)
      eigenv[i] = (Vector*)( tmpalloc + i* (f_totsize/m));
  }else{
    printf("AlgLanczos::run using cache n=%d name=%s\n",m,filename);
    lambda = ecache->eval_address();
    eigenv = ecache->evec_address();
  }
#endif

  lambda = ecache->eval_address();
  eigenv = ecache->evec_address();

  if(!(  (lanczos_arg -> RitzMat_lanczos == MATPCDAG_MATPC 
	  || lanczos_arg -> RitzMat_lanczos == MATPCDAG_MATPC_SHIFT )
	 &&
	 (lanczos_arg -> RitzMat_convcheck == MATPCDAG_MATPC
	  || lanczos_arg -> RitzMat_convcheck == MATPC_HERM )
	 ))
    ERR.NotImplemented(cname,fname, "This RitzMat_{lanczos, convcheck} combination may or maynot work. We didn't check it yet. Likely there needs a (minor) modification/check/debugs. Good luck\n");

  //==============================================
  // print out mass here in case we want to
  // do any output from within the FeigSolv call
  //==============================================
  
  if (lanczos_arg->results!=0){
    FILE* filep;
    filep=Fopen(lanczos_arg->results,"a");
    Fprintf(filep,"mass = %g\n",(Float)lanczos_arg->mass);
    Fclose(filep); // close file
  }
  
  VRB.Result(cname,fname, "mass = %g\n", (Float)lanczos_arg->mass);
  
  //initialization is out of this routine
  if( init_flag==0 ){
    // set the initial vector 
    int nodes = GJP.Nodes(0)*GJP.Nodes(1)*GJP.Nodes(2)*GJP.Nodes(3)*GJP.Nodes(4);
    for(int n=0;n<f_size;n+=2){
      *((Float*)(eigenv[0])+n) = sqrt(2.0/(Float)(nodes*f_size));
      *((Float*)(eigenv[0])+1+n) = 0.0;
    } 
  }else
    // ----------- experiment ------------------------
    // try the linear conbination of eigenvector as an initial vector
    {
      ERR.NotImplemented(cname,fname,"This  init_flag is not supported yet\n");
      Lattice& lattice= lat;
      const int n_fields =  GJP.SnodeSites();  //   *nk ; 
      const int f_size_per_site = lattice.FsiteSize() / GJP.SnodeSites()  / (lattice.FchkbEvl()+1);
      
#if 1
      int neig= 20;
      char* comp_file = "eigsave/eig4dee.mass0.01.traj2200.bc0001";
#else
      int neig= 100;
      char* comp_file = "/scratch1/izubuchi/cps_eigen/DW_b2.13_16x32_ms0.032-mu0.01/eig4dee.mass0.01.traj2200.bc0001";
#endif
      
      
      EigenContainer eigcon( lattice, comp_file, neig, f_size_per_site, n_fields );
      double mass=0.01; 

      int step_eig= (neig - 1)/(ncompress-1);
      
      int f_size = f_size_per_site * n_fields *GJP.VolNodeSites();
      Vector* sol = eigenv[0];

      Float* eig1 = (Float*) smalloc(cname,fname,"eig1",
				     sizeof(Float)*ncompress );

      eigcon. compress(sol, mass, step_eig, eig1);


      if(init_flag>1) 
	eigcon. decompress(eigenv, mass, step_eig, eig1);

	// Try to improve using Ritz
      if(init_flag>2) 
	{
	  EigArg eig_arg;
	  //eig_arg.Encode("eig_arg.tmp","eig_arg");
	  if ( !eig_arg.Decode("ritz_arg.vml", "eig_arg") ) 
	    { 
	      printf("Decoding of ritz_arg failed\n"); exit(-1);
	    }
	  
	  eig_arg. N_eig = 1; //ncompress;
	  eig_arg. N_eigacc = -1;
	  eig_arg.fname = "ritz";
	  
	  CommonArg carg; carg.set_filename("ritz.dat");
	  AlgEig  eig(lattice, &carg, &eig_arg);
	  eig.run(0, eigenv);
	  
	  for(int i=0; i<ncompress;++i)
	    eigcon. nev_check( eigenv[i],  mass );	
	  
	}
    }
  
  // ----------- experiment ends ------------------------

  // Solve for eigenvectors and eigenvalues.
  if(Ncb==2)
    iter = lat.FeigSolv(eigenv, lambda, lanczos_arg, CNV_FRM_YES);
  else if(Ncb==1)
    iter = lat.FeigSolv(eigenv, lambda, lanczos_arg, CNV_FRM_NO);
 
  // Now Let's save them

  // FIXME : pass the argument contains the ensemble informations
  char* ensemble_id = "n/a";
  char* ensemble_label = "n/a";
  int seqNum = -777; // untill GJP has the sequencial number ...
  char* field_type_label = "n/a";

  const int n_fields =  GJP.SnodeSites();  //   *nk ; 
  const int f_size_per_site = lat.FsiteSize() / GJP.SnodeSites()  * Ncb / 2;

  // write to disk if desired, confirm save in cache
  char filename[1024];
  snprintf(filename,1024, "%s.bc%d%d%d%d", alg_lanczos_arg->file, GJP.Bc(0),GJP.Bc(1),GJP.Bc(2),GJP.Bc(3));
  EigenContainer eigcon( lat, filename, nk, f_size_per_site, n_fields, ecache);
  eigcon. save_eval( lambda );
  ecache->eval_cached=1;

  for(int iev=0; iev < nk; iev++)
    ecache->index[iev]=iev;

  int save_stride = GJP.SaveStride();
  if(lanczos_arg->save){
    for(int iev=0; iev < nk; iev+= save_stride){
      // save in "nev" format
      eigcon.nev_save( iev, eigenv[iev], 
		       field_type_label, ensemble_id, ensemble_label, seqNum );
    }
  }

#if 0
  // Free memory if not cached
  //----------------------------------------------------------------
  if( ecache==0 ){
    sfree(cname,fname, "eigenv[0]", eigenv[0]);
    sfree(cname,fname, "eigenv", eigenv);
    sfree(cname,fname,"lambda", lambda);
  }
#endif

  VRB.Result(cname,fname,"Lanczos iterations %d\n",iter);
}


CPS_END_NAMESPACE
