#include<config.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------
/*!\file
  \brief Methods of the AlgEig class.
  
  $Id: alg_eig.C,v 1.27 2013-04-05 17:46:30 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2013-04-05 17:46:30 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_eig/alg_eig.C,v 1.27 2013-04-05 17:46:30 chulwoo Exp $
//  $Id: alg_eig.C,v 1.27 2013-04-05 17:46:30 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: alg_eig.C,v $
//  $Revision: 1.27 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_eig/alg_eig.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/qcdio.h>
#include <util/enum.h>
#include <math.h>
#include <alg/alg_eig.h>
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


CPS_START_NAMESPACE

using namespace std;

extern void gamma_5(Float *v_out, Float *v_in, int num_sites);

int NumChkb( RitzMatType ritz_mat)
{
  // Determine the number of checkerboards
  int Ncb=0;
  char *fname = "NumChkb(RitzType)";
  switch(ritz_mat)
  {
  case MAT_HERM:
  case MATDAG_MAT:
  case NEG_MATDAG_MAT:
  case MATDAG_MAT_NORM:
  case NEG_MATDAG_MAT_NORM:
    Ncb = 2;
    break;
    
  case MATPC_HERM:
  case MATPCDAG_MATPC:
  case NEG_MATPCDAG_MATPC:
    Ncb = 1;
    break;

  default:
    ERR.General("",fname,"RitzMatOper %d not implemented",
		ritz_mat);
  }
  return Ncb;
}

//------------------------------------------------------------------
// Constructor 
/*!
  \param latt The lattice on which to compute the condensate.
  \param c_arg The common argument structure for all algorithms.
  \param arg The algorithm parameters.
 */
//------------------------------------------------------------------
AlgEig::AlgEig(Lattice& latt, 
	       CommonArg *c_arg,
	       EigArg *arg) : 
	       Alg(latt, c_arg) 
{
  cname = "AlgEig";
  char *fname = "AlgEig(L&,CommonArg*,EigArg*)";
  VRB.Func(cname,fname);

  // Initialize the argument pointer
  //----------------------------------------------------------------
  if(arg == 0)
    ERR.Pointer(cname,fname, "arg");
  alg_eig_arg = arg;

  Ncb = NumChkb(alg_eig_arg->RitzMatOper);

  // Set the node size of the full (non-checkerboarded) fermion field
  // NOTE: at this point we must know on what lattice size the operator 
  // will act.
  //----------------------------------------------------------------
  int f_size = GJP.VolNodeSites() * latt.FsiteSize() * Ncb / 2;
  VRB.Flow(cname,fname,"f_size=%d\n",0);

  //  int f_size = GJP.VolNodeSites() * Ncb / 2;
  //  exit(1);

  
  //VRB.Flow(cname,fname,"Doing Ritz");
  int N_eig = alg_eig_arg->N_eig;
  
  // Allocate memory for the eigenvectors and eigenvalues
  //----------------------------------------------------------------
  eigenv = (Vector **) smalloc (cname,fname, "eigenv", N_eig * sizeof(Vector *));
  
  for(int n = 0; n < N_eig; ++n)
    {
      eigenv[n] = (Vector *) smalloc(cname,fname, "eigenv[n]", (f_size)* sizeof(Float));
    }
  
  lambda = (Float *) smalloc(cname,fname, "lambda", 2*N_eig * sizeof(Float));
  
  chirality = (Float *) smalloc(cname,fname,"chirality", N_eig * sizeof(Float));
  
  valid_eig = (int *) smalloc(cname,fname,"valid_eig",N_eig * sizeof(int));
  
  // Print out input parameters
  //----------------------------------------------------------------
  VRB.Input(cname,fname,
	    "N_eig = %d\n",int(N_eig));
  VRB.Input(cname,fname,
	    "MaxCG = %d\n",alg_eig_arg->MaxCG);
  VRB.Input(cname,fname,
	    "Mass_init = %g\n",IFloat(alg_eig_arg->Mass_init));
  VRB.Input(cname,fname,
	    "Mass_final = %g\n",IFloat(alg_eig_arg->Mass_final));
  VRB.Input(cname,fname,
	    "Mass_step = %g\n",IFloat(alg_eig_arg->Mass_step));
  
  // Calculate n_masses if necessary
  VRB.Flow(cname,fname,"alg_eig_arg->pattern_kind=%d\n",alg_eig_arg->pattern_kind);
  switch( alg_eig_arg->pattern_kind ) {
  case ARRAY: 
    n_masses = alg_eig_arg->Mass.Mass_len;
    break;
  case LOG:
    n_masses = (int) ((log(alg_eig_arg->Mass_final - alg_eig_arg->Mass_init)
		       / log(alg_eig_arg->Mass_step)) + 1.000001);
    break;
  case LIN:   
    n_masses = (int) (fabs((alg_eig_arg->Mass_final 
			    - alg_eig_arg->Mass_init)/alg_eig_arg->Mass_step) + 1.000001); 
    break;
  default: 
    ERR.General(cname, fname,
		"alg_eig_arg->pattern_kind = %d is unrecognized\n", 
		alg_eig_arg->pattern_kind);
    break;
  }  
  
  VRB.FuncEnd(cname,fname);
}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
AlgEig::~AlgEig() {
  char *fname = "~AlgEig()";
  VRB.Func(cname,fname);

  // Free memory
  //----------------------------------------------------------------
  sfree(cname,fname, "valid_eig", valid_eig);
    
  sfree(cname,fname, "chirality", chirality);
    
  sfree(cname,fname, "lambda", lambda);
    
  for(int n = alg_eig_arg->N_eig - 1; n >= 0; --n)
    {
      sfree(cname,fname, "eigenv[n] ",eigenv[n]);
    }
    
  sfree(cname,fname, "eigenv", eigenv);
  //???
}

//!< Overloaded for backwards compatibility
void AlgEig::run()
{
  run((Float**)0);
}

//------------------------------------------------------------------
//! Performs the computation.
/*!
  \post The results are written to the file specified in the common_arg
  structure.
*/
//------------------------------------------------------------------
void AlgEig::run(Float **evalues, Vector **in_eigv)
{
#if TARGET==cpsMPI
    using MPISCU::fprintf;
#endif

  Float time = -dclock();
  int iter=0;
  EigArg *eig_arg;
  char *fname = "run()";
  VRB.Func(cname,fname);

  // Set the Lattice pointer eig_arg
  //----------------------------------------------------------------
  Lattice& lat = AlgLattice();
  eig_arg = alg_eig_arg;
  Float **hsum;
  const int N_eig = eig_arg->N_eig;
  const int f_size = GJP.VolNodeSites() * lat.FsiteSize() * Ncb / 2;
  int hsum_len = 0;
  int n;

  if (eig_arg->print_hsum)
  {
    switch(eig_arg->hsum_dir)
    {
    case 0:
      hsum_len = GJP.Xnodes()*GJP.XnodeSites();
      break;
    case 1:
      hsum_len = GJP.Ynodes()*GJP.YnodeSites();
      break;
    case 2:
      hsum_len = GJP.Znodes()*GJP.ZnodeSites();
      break;
    case 3:
      hsum_len = GJP.Tnodes()*GJP.TnodeSites();
      break;
    case 4:
      if (lat.Fclass() == F_CLASS_DWF) 
        hsum_len = GJP.Snodes()*GJP.SnodeSites();
      else
        ERR.General(cname,fname,"Invalid direction\n");
      break;
     default:
      ERR.General(cname,fname,"Invalid direction\n");
    }

    hsum = (Float **) smalloc(cname,fname,"hsum",N_eig * sizeof(Float*)); // surely Float* ?
  
    for(n = 0; n < N_eig; ++n)
    {
      hsum[n] = (Float *) smalloc(cname,fname,"hsum[n]",hsum_len * sizeof(Float));
    }
  }
  else
  {
    hsum = (Float **) 0;
  }

  // allocate memory to store the eigenvectors
  Vector** eig_store=0;
  if(eig_arg->ncorr) 
  {
    eig_store = (Vector**)smalloc(cname,fname, "eig_store",N_eig * sizeof(Vector*));
    for(n=0;n<N_eig;++n)
    {
      eig_store[n] = (Vector*) smalloc(cname,fname, "eig_store",f_size * sizeof(Float));
    }
  }


  // Initialize eigenvectors to gaussian
  // and compute eigenvectors
  //----------------------------------------------------------------
  //for(int n = 0; n < eig_arg->N_eig; ++n)
  //  lat.RandGaussVector(eigenv[n], 0.5, Ncb);
  
  int sign_dm=1;

  // Initialize the cg_arg mass, with the first mass we
  // want to compute for:
  switch( eig_arg->pattern_kind ) {
  case ARRAY: 
    eig_arg->mass = eig_arg->Mass.Mass_val[0]; 
    break;
  case LIN:   
    // Loop over mass values
    sign_dm = (eig_arg->Mass_step < 0.0) ? -1 : 1;
    eig_arg->mass = eig_arg->Mass_init; 
    if (sign_dm*eig_arg->Mass_init > sign_dm*eig_arg->Mass_final)
      ERR.General(cname,fname,"initial and final mass not valid\n");
    break;
  case LOG:   
    eig_arg->mass = eig_arg->Mass_init; 
    break;
  default: 
    ERR.General(cname, fname,
		"eig_arg->pattern_kind = %d is unrecognized\n", 
		eig_arg->pattern_kind);
    break;
  }

  // Loop over masses
  for(int m=0; m<n_masses; m++){
    
    //for(Float mass = eig_arg->Mass_init; 
    //sign_dm*mass <= sign_dm*eig_arg->Mass_final;
    //mass += eig_arg->Mass_step){

    //eig_arg->mass = mass;
    // count the number of masses 
    int count(0);


    // store eigenvectors from previous mass
    if(eig_arg->ncorr && ( count > 0 ) ){
      for ( n=0; n<N_eig; n++ )
      {
        eig_store[n]->CopyVec(eigenv[n],f_size); 
      }
    }


    // TIZB use the input eigenvector as starting vectors
    //if(in_eigv){
     // for(n = 0; n<N_eig; ++n)
//	eigenv[n]->CopyVec(in_eigv[n],f_size); 
 //   } else {
      // random guess every time; do *not* put in the old solution
      // TIZB why ? Ask Chulwoo !
      for(n = 0; n<N_eig; ++n)
	{
	  lat.RandGaussVector(eigenv[n], 0.5, Ncb);
	}
  //  }



/*
    // DEBUG, dumping start vector from QCDOC to be read in QCDSP
    cout << "Dump DEBUGGING info...startvector.dat" << endl;
    int f_size = GJP.VolNodeSites() * lat.FsiteSize() * Ncb / 2;
    ofstream fout ( "startvector.dat");
    fout << "f_size=" << f_size << endl;
    for(n=0;n<N_eig;++n) {
      fout << "n=" << n << endl;
      for(int i=0;i<f_size;i++) {
	fout << setprecision(15) << ((Float*)eigenv[n])[i] << endl;
      }
      fout << endl;
    }
    fout.close();
    //    exit(0);
    //DEBUG end
*/
   
    //==============================================
    // print out mass here in case we want to
    // do any output from within the FeigSolv call
    //==============================================
			      
    if (eig_arg->fname!=0)
    {
      FILE* filep;
      filep=Fopen(eig_arg->fname,"a");
      Fprintf(filep,"mass = %g\n",(Float)eig_arg->mass);
      Fclose(filep); // close file
    }
         
    VRB.Result(cname,fname, "mass = %g\n", (Float)eig_arg->mass);

    // Solve for eigenvectors and eigenvalues.
    // Use eigenv as initial guess. Lambda is not used initially.

    if(Ncb==2)
      iter = lat.FeigSolv(eigenv, lambda, chirality, valid_eig, 
			hsum, eig_arg, CNV_FRM_YES);
    else if(Ncb==1)
      iter = lat.FeigSolv(eigenv, lambda, chirality, valid_eig, 
      			hsum, eig_arg, CNV_FRM_NO);

    //------------------------------------------------------------
    // Solve for eigenvectors and eigenvalues.
    // Lambda is not used initially.
    // This call will return a negative value for the iteration number
    // if either the solver maxes out on one of the limits.
#if 0
    iter = lat.FeigSolv(eigenv,lambda,chirality, valid_eig,
    			hsum, eig_arg, CNV_FRM_YES);
#endif
   
    //!< Copy over eigenvalues to return them
    //if (evalues != 0) {
    //  for (int eig=0; eig<eig_arg->N_eig; eig++) {
//	evalues[eig][m] = lambda[eig];
 //     }
  //  }

    if ( iter < 0 )
    {
      FILE* filep;
      filep=Fopen(eig_arg->fname,"a");
      if ( iter == -1 )
      {
        Fprintf(filep, "maxed out in ritz\n");
      }
      else
      {
        Fprintf(filep, "maxed out for KS steps\n");
      }
      Fclose(filep);
    }
    else
    {   
      //----------------------------------------
      //  spatial correlations of eigen vectors
      //---------------------------------------
      
      if ( eig_arg->fname != 0x0 )
	{
	  
	  FILE* filep;
	  filep=Fopen(eig_arg->fname,"a");
              
	  int i_eig,j_eig;
              
	  // GM5Correlation
              
	  //tmp vector v1
	  Vector* v1 = (Vector *)smalloc(cname, fname, "v1",f_size*sizeof(Float));
              
	  for(i_eig=0;i_eig<N_eig;i_eig++){
	    for(j_eig=i_eig;j_eig<N_eig;j_eig++){
	      gamma_5((Float*)v1, 
		      (Float*)eigenv[j_eig], 
		      f_size/24);
	      Complex cr = eigenv[i_eig]->CompDotProductGlbSum(v1,f_size);
	      Fprintf(filep,"GM5CORR: %d %d %g %g\n",
		      i_eig, j_eig, (Float)cr.real(), (Float)cr.imag());
	    }
	  }
	  sfree(cname,fname, "v1", v1);
              
	  // Correlation with previous eigen vector
	  if(count >0 && eig_arg->ncorr ){
	    for(i_eig=0;i_eig<N_eig;i_eig++){
	      for(j_eig=0;j_eig<N_eig;j_eig++){
		Complex cr = eig_store[i_eig]
		  ->CompDotProductGlbSum(eigenv[j_eig],f_size);
		Fprintf(filep,"NeibCORR: %d %d %g %g\n",
			i_eig, j_eig, (Float)cr.real(), (Float)cr.imag());
	      }
	    }
	  }


	  //------------------------------------------
	  // Print out number of iterations  and hsum
	  //------------------------------------------
	  
	  int i;
	  Fprintf(filep, "  iter = %d\n", iter);
	  if (eig_arg->print_hsum)
	    {
	      for(n = 0; n < eig_arg->N_eig; ++n)
		{
		  for(i = 0; i < hsum_len; ++i)
		    {
		      Fprintf(filep, "  hsum[%d][%d] = %g\n",n,i,
			      (Float)hsum[n][i]);
		    }
		}
	    }
	  Fclose(filep);
	} // output
    } // solver worked

    // If there is another mass loop iteration ahead, we should
    // set the eig_arg->mass to it's next desired value
    if( m < n_masses - 1 ) {
      switch( eig_arg->pattern_kind ) {
      case ARRAY: 
	eig_arg->mass = eig_arg->Mass.Mass_val[m+1]; 
	break;
      case LIN:   
	eig_arg->mass += eig_arg->Mass_step; 
	break;
      case LOG:   
	eig_arg->mass *= eig_arg->Mass_step; 
	break;
      }
    }
    count++;
  } // mass loop

  //==============================
  // deallocate eigenvalue store
  //==============================
  
  if ( eig_arg->ncorr )
    {
      for(n = 0; n < N_eig; ++n)
        {
          sfree(cname,fname,"eig_store[n]",eig_store[n]);
        }
      sfree(cname,fname,"eig_store",eig_store);
    }
  time +=dclock();
  print_flops(cname,fname,0,time);

  // TIZB
  //if(in_eigv){
   // for(n = 0; n<N_eig; ++n)
    //  in_eigv[n]->CopyVec(eigenv[n],f_size); 
  //}

}



CPS_END_NAMESPACE
