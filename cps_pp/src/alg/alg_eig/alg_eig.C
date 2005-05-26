#include<config.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------
/*!\file
  \brief Methods of the AlgEig class.
  
  $Id: alg_eig.C,v 1.12 2005-05-26 14:32:26 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2005-05-26 14:32:26 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_eig/alg_eig.C,v 1.12 2005-05-26 14:32:26 chulwoo Exp $
//  $Id: alg_eig.C,v 1.12 2005-05-26 14:32:26 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: alg_eig.C,v $
//  $Revision: 1.12 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_eig/alg_eig.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/qcdio.h>
#include <math.h>
#include <alg/alg_eig.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/smalloc.h>
#include <util/vector.h>
#include <util/verbose.h>
#include <util/error.h>
CPS_START_NAMESPACE


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

  // Determine the number of checkerboards
  switch(alg_eig_arg->RitzMatOper)
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
    ERR.General(cname,fname,"RitzMatOper %d not implemented",
		alg_eig_arg->RitzMatOper);
  }


  // Set the node size of the full (non-checkerboarded) fermion field
  // NOTE: at this point we must know on what lattice size the operator 
  // will act.
  //----------------------------------------------------------------
  int f_size = GJP.VolNodeSites() * latt.FsiteSize() * Ncb / 2;
  int N_eig = alg_eig_arg->N_eig;

  // Allocate memory for the eigenvectors and eigenvalues
  //----------------------------------------------------------------
  eigenv = (Vector **) smalloc (cname,fname, "eigenv", N_eig * sizeof(Vector *));
  
  for(int n = 0; n < N_eig; ++n)
  {
    eigenv[n] = (Vector *) smalloc(cname,fname, "eigenv[n]", f_size * sizeof(Float));
  }

  lambda = (Float *) smalloc(cname,fname, "lambda", N_eig * sizeof(Float));

  chirality = (Float *) smalloc(N_eig * sizeof(Float));
  if(chirality == 0)
    ERR.Pointer(cname,fname, "chirality");
  VRB.Smalloc(cname,fname, "chirality", chirality, N_eig * sizeof(Float));

  valid_eig = (int *) smalloc(N_eig * sizeof(int));
  if(valid_eig == 0)
    ERR.Pointer(cname,fname, "valid_eig");
  VRB.Smalloc(cname,fname, "valid_eig", valid_eig, N_eig * sizeof(int));

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
  switch( alg_eig_arg->pattern_kind ) {
  case ARRAY: 
  case LOG:   
    // Should be specified already
    break;
  case LIN:   
    alg_eig_arg->n_masses = (int) (fabs((alg_eig_arg->Mass_final 
      - alg_eig_arg->Mass_init)/alg_eig_arg->Mass_step) + 1.000001); 
    break;
  default: 
    ERR.General(cname, fname,
		"alg_eig_arg->pattern_kind = %d is unrecognized\n", 
		alg_eig_arg->pattern_kind);
    break;
  }  

}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
AlgEig::~AlgEig() {
  char *fname = "~AlgEig()";
  VRB.Func(cname,fname);

  // Free memory
  //----------------------------------------------------------------
  VRB.Sfree(cname,fname, "valid_eig", valid_eig);
  sfree(valid_eig);

  VRB.Sfree(cname,fname, "chirality", chirality);
  sfree(chirality);

  VRB.Sfree(cname,fname, "lambda", lambda);
  sfree(lambda);

  for(int n = alg_eig_arg->N_eig - 1; n >= 0; --n)
  {
    VRB.Sfree(cname,fname, "eigenv[n] ",eigenv[n]);
    sfree(eigenv[n]);
  }

  VRB.Sfree(cname,fname, "eigenv", eigenv);
  sfree(eigenv);

  //???
}


//------------------------------------------------------------------
//! Performs the computation.
/*!
  \post The results are written to the file specified in the common_arg
  structure.
*/
//------------------------------------------------------------------
void AlgEig::run()
{
#if TARGET==cpsMPI
    using MPISCU::fprintf;
#endif

  int iter=0;
  EigArg *eig_arg;
  char *fname = "run()";
  VRB.Func(cname,fname);

  // Set the Lattice pointer eig_arg
  //----------------------------------------------------------------
  Lattice& lat = AlgLattice();
  eig_arg = alg_eig_arg;
  Float **hsum;
  int N_eig = eig_arg->N_eig;
//  int f_size = GJP.VolNodeSites() * lat.FsiteSize() * Ncb / 2;
  int hsum_len = 0;

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

    hsum = (Float **) smalloc(N_eig * sizeof(Float*)); // surely Float* ?
    if(hsum == 0)
      ERR.Pointer(cname,fname, "hsum");
    VRB.Smalloc(cname,fname, "hsum", hsum, N_eig * sizeof(Float*));
  
    for(int n = 0; n < N_eig; ++n)
    {
      hsum[n] = (Float *) smalloc(hsum_len * sizeof(Float));
      if(hsum[n] == 0)
	ERR.Pointer(cname,fname, "hsum[n]");
      VRB.Smalloc(cname,fname, "hsum[n]", hsum[n], hsum_len*sizeof(Float));
    }
  }
  else
  {
    hsum = (Float **) 0;
  }

  // Initialize eigenvectors to gaussian
  // and compute eigenvectors
  //----------------------------------------------------------------
  for(int n = 0; n < eig_arg->N_eig; ++n)
    lat.RandGaussVector(eigenv[n], 0.5, Ncb);
  
  // Loop over mass values
  int sign_dm = (eig_arg->Mass_step < 0.0) ? -1 : 1;

  // Initialize the cg_arg mass, with the first mass we
  // want to compute for:
  switch( eig_arg->pattern_kind ) {
  case ARRAY: 
    eig_arg->mass = eig_arg->Mass[0]; 
    break;
  case LIN:   
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
  for(int m=0; m<eig_arg->n_masses; m++){
    
    //for(Float mass = eig_arg->Mass_init; 
    //sign_dm*mass <= sign_dm*eig_arg->Mass_final;
    //mass += eig_arg->Mass_step){

    //eig_arg->mass = mass;

    // Solve for eigenvectors and eigenvalues.
    // Use eigenv as initial guess. Lambda is not used initially.
    if(Ncb==2)
      iter = lat.FeigSolv(eigenv, lambda, chirality, valid_eig, 
			hsum, eig_arg, CNV_FRM_YES);
    else if(Ncb==1)
      iter = lat.FeigSolv(eigenv, lambda, chirality, valid_eig, 
      			hsum, eig_arg, CNV_FRM_NO);
    // Print out number of iterations and eigs
    //----------------------------------------------------------------
    if(common_arg->results != 0)
    {
      FILE *fp;
      if( (fp = Fopen((char *)common_arg->results, "a")) == NULL ) {
	ERR.FileA(cname,fname, (char *)common_arg->results);
      }

      Fprintf(fp, "mass = %g\n", (IFloat)(eig_arg->mass));
      Fprintf(fp, "  iter = %d\n", iter);
      for(int n = 0; n < eig_arg->N_eig; ++n)
      {
	Fprintf(fp, "  lambda[%d] = %g  chirality = %g  valid = %d\n", 
		n, (IFloat)lambda[n], (IFloat)chirality[n], valid_eig[n]);
      }
      
      if (eig_arg->print_hsum)
      {
	for(int n = 0; n < eig_arg->N_eig; ++n)
	  for(int i = 0; i < hsum_len; ++i)
	    Fprintf(fp, "  hsum[%d][%d] = %g\n",n,i,(IFloat)hsum[n][i]);
      }
      
      Fclose(fp);
    }

    // If there is another mass loop iteration ahead, we should
    // set the eig_arg->mass to it's next desired value
    if( m < eig_arg->n_masses - 1 ) {
      switch( eig_arg->pattern_kind ) {
      case ARRAY: 
	eig_arg->mass = eig_arg->Mass[m+1]; 
	break;
      case LIN:   
	eig_arg->mass += eig_arg->Mass_step; 
	break;
      case LOG:   
	eig_arg->mass *= eig_arg->Mass_step; 
	break;
      }
    }
  }

}

CPS_END_NAMESPACE
