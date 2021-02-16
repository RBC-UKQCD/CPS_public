#include<config.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------
/*!\file
  \brief Methods of the AlgLanczos class.
  
  $Id: alg_lanczos.C,v 1.25 2009/03/23 19:13:32 chulwoo Exp $
*/

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
//#include <iomanip>
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
  int nt= lanczos_arg-> nt_lanczos_vectors;
  int np= lanczos_arg-> np_lanczos_vectors;
  if (nt>nk)
  ERR.General(cname,fname,"nk(%d) cannot be smaller than nt(%d)\n",nk,nt);
  int m = nk+np;
  size_t f_size = GJP.VolNodeSites() * lat.FsiteSize() * Ncb / 2;

  lambda = ecache->eval_address();
  // eigenvectors are stored in cache. Pass along addresses.
  eigenv = (Vector**)smalloc(m*sizeof(Vector*));
  for(int i=0;i<m;++i){
	 eigenv[i] = ecache->vec_ptr(i);
	VRB.Result (cname, fname, "eigenv[%d]=%p\n", i, eigenv[i]);
  }


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
      *((float*)(eigenv[0])+n) = (float)(std::sqrt(2.0)/std::sqrt((float)nodes)/std::sqrt((float)f_size));
      *((float*)(eigenv[0])+1+n) = 0.0;
    }
#if 0
    if(!UniqueID()){printf("Nodes: %d %d %d %d %d %d\n",GJP.Nodes(0),GJP.Nodes(2),GJP.Nodes(2),GJP.Nodes(3),GJP.Nodes(4),nodes);};
    if(!UniqueID()){printf("f_zise: %d\n",f_size);};
    for(int n=0;n<f_size;n+=2){
	printf("evec0 %d %d %e %e\n",
	UniqueID(),n,
	*((float*)(eigenv[0])+n),*((float*)(eigenv[0])+n+1));
    } 
#endif
  }else
    // ----------- experiment ------------------------
    // try the linear conbination of eigenvector as an initial vector
    {
      ERR.NotImplemented(cname,fname,"This  init_flag is not supported yet\n");
#if 0
      Lattice& lattice= lat;
      const int n_fields =  GJP.SnodeSites();  //   *nk ; 
      const size_t f_size_per_site = lattice.FsiteSize() / GJP.SnodeSites()  / (lattice.FchkbEvl()+1);
      
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
      
      size_t f_size = f_size_per_site * n_fields *GJP.VolNodeSites();
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
#endif
    }
  
  // ----------- experiment ends ------------------------

  // Solve for eigenvectors and eigenvalues.
  if(Ncb==2)
    iter = lat.FeigSolv(eigenv, lambda, lanczos_arg, CNV_FRM_YES);
  else if(Ncb==1)
    iter = lat.FeigSolv(eigenv, lambda, lanczos_arg, CNV_FRM_NO);
  for(int i=0;i<m;++i){
	 ecache->set_ptr(i,eigenv[i]);
	VRB.Result (cname, fname, "eigenv[%d]=%p ecache %p\n", i, eigenv[i],ecache->vec_ptr(i));
  }
 
#if 0
//EigenContainer deprecated
  // Now Let's save them

  // FIXME : pass the argument contains the ensemble informations
  char* ensemble_id = "n/a";
  char* ensemble_label = "n/a";
  int seqNum = -777; // untill GJP has the sequencial number ...
  char* field_type_label = "n/a";

  const int n_fields =  GJP.SnodeSites();  //   *nk ; 
  const size_t f_size_per_site = lat.FsiteSize() / GJP.SnodeSites()  * Ncb / 2 / 2; // 2nd two is for float!

  // write to disk if desired, confirm save in cache
  char filename[1024];
//  snprintf(filename,1024, "%s.bc%d%d%d%d", alg_lanczos_arg->file, GJP.Bc(0),GJP.Bc(1),GJP.Bc(2),GJP.Bc(3));
  snprintf(filename,1024, "%s", alg_lanczos_arg->file);
  EigenContainer eigcon( lat, filename, nk, f_size_per_site, n_fields, ecache);
//  if(lanczos_arg->save) {
// I'm not sure why the eigenvectors are saved only when 'save' flag is on?
  eigcon. save_eval( lambda );
  ecache->eval_cached=1;
//  }
  for(int iev=0; iev < nt; iev++)
    ecache->index[iev]=iev;

  int save_stride = GJP.SaveStride();
  if(lanczos_arg->save){
    for(int iev=0; iev < nt; iev+= save_stride){
      // save in "nev" format
      eigcon.nev_save( iev, eigenv[iev], 
		       field_type_label, ensemble_id, ensemble_label, seqNum );
    }
  }
#endif

#if 0
  // Free memory if not cached
  //----------------------------------------------------------------
  if( ecache==0 ){
    sfree(cname,fname, "eigenv[0]", eigenv[0]);
    sfree(cname,fname, "eigenv", eigenv);
    sfree(cname,fname,"lambda", lambda);
  }
#endif

  // free temp pointer to evecs. Still stored in cache.
  sfree(cname,fname, "eigenv", eigenv);

  VRB.Result(cname,fname,"Lanczos iterations %d\n",iter);
}


CPS_END_NAMESPACE
