#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief Methods of the AlgPbp class.
  
  $Id: alg_pbp.C,v 1.14 2012-07-06 20:22:08 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2012-07-06 20:22:08 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_pbp/alg_pbp.C,v 1.14 2012-07-06 20:22:08 chulwoo Exp $
//  $Id: alg_pbp.C,v 1.14 2012-07-06 20:22:08 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: alg_pbp.C,v $
//  $Revision: 1.14 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_pbp/alg_pbp.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// alg_pbp.C
//
// AlgPbp is derived from Alg and is relevant to the 
// stochastic measurement of PsiBar Psi using the
// Conjugate Gradient algorithm. The type of fermion is
// determined by the argument to the constructor.
//
// PsiBarPsi is normalized so that for large values of the
// PbpArg.mass  PsiBarPsi =  1 / mass for any fermion type.
// This normalization results to the following small mass
// behavior for a trivial background gauge field with periodic
// boundary conditions:
// Staggered = 16 / ( Volume * mass )
// Wilson    =  1 / ( Volume * mass )
// Dwf       =  1 / ( Volume * mass )
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/qcdio.h>
#include <alg/alg_pbp.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/smalloc.h>
#include <util/vector.h>
#include <util/verbose.h>
#include <util/error.h>
CPS_START_NAMESPACE

#define POINT
#undef POINT
#define Z2
#undef Z2

//------------------------------------------------------------------
/*!
  \param latt The lattice on which to compute the condensate.
  \param c_arg The common argument structure for all algorithms.
  \param arg The algorithm parameters.
 */
//------------------------------------------------------------------
AlgPbp::AlgPbp(Lattice& latt, 
	       CommonArg *c_arg,
	       PbpArg *arg) : 
	       Alg(latt, c_arg) 
{
  cname = "AlgPbp";
  char *fname = "AlgPbp(L&,CommonArg*,PbpArg*)";
  VRB.Func(cname,fname);

  // Initialize the argument pointer
  //----------------------------------------------------------------
  if(arg == 0)
    ERR.Pointer(cname,fname, "arg");
  alg_pbp_arg = arg;


  // Set the node size of the full (non-checkerboarded) fermion field
  //----------------------------------------------------------------
  f_size = GJP.VolNodeSites() * latt.FsiteSize();


  // Allocate memory for the source.
  //----------------------------------------------------------------
  src = (Vector *) smalloc(f_size * sizeof(Float));
  if(src == 0)
    ERR.Pointer(cname,fname, "src");
  VRB.Smalloc(cname,fname, "src", src, f_size * sizeof(Float));


  // Allocate memory for the solution
  //----------------------------------------------------------------
  sol = (Vector *) smalloc(f_size * sizeof(Float));
  if(sol == 0)
    ERR.Pointer(cname,fname, "sol");
  VRB.Smalloc(cname,fname, "sol", sol, f_size * sizeof(Float));


}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
AlgPbp::~AlgPbp() {
  char *fname = "~AlgPbp()";
  VRB.Func(cname,fname);

  // Free memory
  //----------------------------------------------------------------
  VRB.Sfree(cname,fname, "sol", sol);
  sfree(sol);
  VRB.Sfree(cname,fname, "src", src);
  sfree(src);
}


//------------------------------------------------------------------
//! Performs the computation.
/*!
  \post The results are written to the file specified in the common_arg
  structure.
*/
//------------------------------------------------------------------
void AlgPbp::run()
{
#if TARGET==cpsMPI
    using MPISCU::fprintf;
#endif
  int iter;
  int ls;
  int ls_glb;
  Float pbp= 0., pbd0p, pbdip;
  Float pbg5p= 0.;
  Float pbp_norm;
  Float true_res;
  PbpArg *pbp_arg;
  CgArg cg_arg_struct;
  CgArg *cg_arg = &cg_arg_struct;
  char *fname = "run()";
  VRB.Func(cname,fname);

/////  printf("HERE HERE... \n");


  // Set the Lattice pointer pbp_arg and cg_arg
  //----------------------------------------------------------------
  Lattice& lat = AlgLattice();
  pbp_arg = alg_pbp_arg;
  cg_arg->max_num_iter = pbp_arg->max_num_iter;
  cg_arg->stop_rsd = pbp_arg->stop_rsd;
  
  // Make a Float pointer to sol
  //----------------------------------------------------------------
  Float *f_sol = (Float *) sol;

  //----------------------------------------------------------------
  // Initialize source and solution (initial guess)
  // and calculate PsiBarPsi for:
  //----------------------------------------------------------------


  //----------------------------------------------------------------
  // Domain Wall fermions
  //----------------------------------------------------------------
  if(lat.Fclass() == F_CLASS_DWF || lat.Fclass() == F_CLASS_BFM){
    ls = GJP.SnodeSites();
    ls_glb = GJP.Snodes() * GJP.SnodeSites();

    // Allocate memory for the 4-dimensional source, solution
    Vector *src_4d = (Vector *) smalloc(f_size * sizeof(Float) / ls);
    if(src_4d == 0)
      ERR.Pointer(cname,fname, "src_4d");
    VRB.Smalloc(cname,fname, "src_4d", src_4d, f_size * sizeof(Float) / ls);
    Vector *sol_4d = (Vector *) smalloc(f_size * sizeof(Float) / ls);
    if(sol_4d == 0)
      ERR.Pointer(cname,fname, "sol_4d");
    VRB.Smalloc(cname,fname, "sol_4d", sol_4d, f_size * sizeof(Float) / ls);

    // allocate space for pbp, pb_g5_p array
    Float * pbp_all = (Float *) smalloc(ls_glb * sizeof(Float));
    if(pbp_all == 0)
      ERR.Pointer(cname,fname, "pbp_all");
    VRB.Smalloc(cname,fname, "pbp_all", pbp_all, ls_glb * sizeof(Float));

    Float * pbg5p_all = (Float *) smalloc(ls_glb * sizeof(Float));
    if(pbg5p_all == 0)
      ERR.Pointer(cname,fname, "pbg5p_all");
    VRB.Smalloc(cname,fname, "pbg5p_all", pbg5p_all, ls_glb * sizeof(Float));


    // initialize 4-dimensional source
    lat.RandGaussVector(src_4d, 0.5, FOUR_D);

    // set the 5-dimensional source
    lat.Ffour2five(src, src_4d, pbp_arg->src_u_s, pbp_arg->src_l_s);

    // Initialize the cg_arg mass, with the first mass we
    // want to compute for:
    switch( pbp_arg->pattern_kind ) {
    case ARRAY: 
      cg_arg->mass = pbp_arg->mass[0]; 
      break;
    case LIN:   
      cg_arg->mass = pbp_arg->mass_start; 
      break;
    case LOG:   
      cg_arg->mass = pbp_arg->mass_start; 
      break;
    default: 
      ERR.General(cname, fname,
		  "pbp_arg->pattern_kind = %d is unrecognized\n", 
		  pbp_arg->pattern_kind);
      break;
    }

    // Loop over masses
    for(int m=0; m<pbp_arg->n_masses; m++){

      // initialize 5-dimensional solution (initial guess) to 1
      for(int i=0; i< f_size/2; i++){
	f_sol[2*i] = 1.0;     // real part
	f_sol[2*i+1] = 0.0;   // imaginary part
      }

      // do inversion
      iter = lat.FmatInv(sol, src, cg_arg, &true_res, CNV_FRM_YES);

      // Calculate the pbp normalization factor 
      pbp_norm = GJP.VolSites() * lat.Colors() * lat.SpinComponents();
      if(lat.Fclass() == F_CLASS_DWF) {
          pbp_norm *= 4.0 + GJP.DwfA5Inv() - GJP.DwfHeight();
      }
      // pbp_norm = (4.0 + GJP.DwfA5Inv() - GJP.DwfHeight())
      //   * GJP.VolSites() 
      //   * ( lat.FsiteSize() / (2 * GJP.SnodeSites()) );  

      VRB.Result(cname, fname, "pbp_norm = %17.10e\n", pbp_norm);

      if (pbp_arg->snk_loop) {
	// Loop over sink - source separation
        for (int i = 0; i < ls_glb; i++) {
          // set the 4-dimensional solution
          lat.Ffive2four(sol_4d, sol, i, ls_glb - i - 1);
	  
          // Calculate pbp
          pbp_all[i] = sol_4d->ReDotProductGlbSum4D(src_4d, f_size/ls)
	             / pbp_norm;

          // Calculate pbg5p = Tr[ PsiBar * Gamma5 * Psi]
	  lat.Gamma5(sol_4d, sol_4d, GJP.VolNodeSites());
          pbg5p_all[i] = sol_4d->ReDotProductGlbSum4D(src_4d, f_size/ls)
	               / pbp_norm;
        }
      } 
      else {
        // set the 4-dimensional solution
        lat.Ffive2four(sol_4d, sol, pbp_arg->snk_u_s, pbp_arg->snk_l_s);

        // Calculate pbp
        pbp = sol_4d->ReDotProductGlbSum4D(src_4d, f_size/ls) / pbp_norm;

	// Calculate pbg5p = Tr[ PsiBar * Gamma5 * Psi]
	lat.Gamma5(sol_4d, sol_4d, GJP.VolNodeSites());
	pbg5p = sol_4d->ReDotProductGlbSum4D(src_4d, f_size/ls) / pbp_norm;
      }

      // Open file for results
      if(common_arg->results != 0){
	FILE *fp;
	if( (fp = Fopen((char *)common_arg->results, "a")) == NULL ) {
	  ERR.FileA(cname,fname, (char *)common_arg->results);
        }
        if (pbp_arg->snk_loop) {
          // Print out mass, s-slice, pbp ,pbg5p, number of iterations
          // and true residual
          for (int i = 0; i < ls_glb; i++) {
	    Fprintf(fp,"%0.16e %d %0.16e %0.16e %d %0.16e\n", 
		    IFloat(cg_arg->mass), 
		    i,
		    IFloat(pbp_all[i]), 
		    IFloat(pbg5p_all[i]), 
		    iter,
		    IFloat(true_res));
          }
        } 
	else {
          // Print out mass, pbp, pbg5p, number of iterations and true residual
          Fprintf(fp, "%0.16e %0.16e %0.16e %d %0.16e\n", 
		  IFloat(cg_arg->mass),
                  IFloat(pbp), 
                  IFloat(pbg5p), 
		  iter, 
		  IFloat(true_res));
        }
	Fclose(fp);
      }

      // If there is another mass loop iteration ahead, we should
      // set the cg_arg->mass to it's next desired value
      if( m < pbp_arg->n_masses - 1 ) {
        switch( pbp_arg->pattern_kind ) {
	case ARRAY: 
	  cg_arg->mass = pbp_arg->mass[m+1]; 
	  break;
	case LIN:   
	  cg_arg->mass += pbp_arg->mass_step; 
	  break;
	case LOG:   
	  cg_arg->mass *= pbp_arg->mass_step; 
	  break;
        }
      }
    }

    // free the 4-dimensional source, solution, pbp, and pbg5p
    VRB.Sfree(cname,fname, "pbg5p_all", pbg5p_all);
    sfree(pbg5p_all);
    VRB.Sfree(cname,fname, "pbp_all", pbp_all);
    sfree(pbp_all);
    VRB.Sfree(cname,fname, "sol_4d", sol_4d);
    sfree(sol_4d);
    VRB.Sfree(cname,fname, "src_4d", src_4d);
    sfree(src_4d);

  }

  //----------------------------------------------------------------
  // Wilson or Clover fermions
  //----------------------------------------------------------------
  else if (   (lat.Fclass() == F_CLASS_WILSON) 
	   || (lat.Fclass() == F_CLASS_CLOVER)   ) {  

    // Allocate memory for gamma_5 * solution
    Vector *sol_g5 = (Vector *) smalloc(f_size * sizeof(Float));
    if(sol_g5 == 0)
      ERR.Pointer(cname,fname, "sol_g5");
    VRB.Smalloc(cname,fname, "sol_g5", sol_g5, f_size * sizeof(Float));

    // set source
    lat.RandGaussVector(src, 0.5);

    // Initialize the cg_arg mass, with the first mass we
    // want to compute for:
    switch( pbp_arg->pattern_kind ) {
    case ARRAY: 
      cg_arg->mass = pbp_arg->mass[0]; 
      break;
    case LIN:   
      cg_arg->mass = pbp_arg->mass_start; 
      break;
    case LOG:   
      cg_arg->mass = pbp_arg->mass_start; 
      break;
    default: 
      ERR.General(cname, fname,
		  "pbp_arg->kind = %d is unrecognized\n", 
		  pbp_arg->pattern_kind);
      break;
    }

    // Loop over masses
    for(int m=0; m<pbp_arg->n_masses; m++){
      
      // initialize solution (initial guess) to 1
      for(int i=0; i< f_size/2; i++){
	f_sol[2*i] = 1.0;     // real part
	f_sol[2*i+1] = 0.0;   // imaginary part
      }

      // do inversion
      iter = lat.FmatInv(sol, src, cg_arg, &true_res, CNV_FRM_YES);

      // Calculate the pbp normalization factor
      pbp_norm = ( (4.0 + cg_arg->mass) )
	         * (GJP.VolSites() 
	         * ( lat.FsiteSize()) / 2 );

      // calculate pbp
      pbp = sol->ReDotProductGlbSum4D(src, f_size) / pbp_norm;

      // Calculate pbg5p = Tr[ PsiBar * Gamma5 * Psi]
      lat.Gamma5(sol_g5, sol, GJP.VolNodeSites());
      pbg5p = sol_g5->ReDotProductGlbSum4D(src, f_size) / pbp_norm;

      // Print out mass, pbp, pbg5p number of iterations and true residual
      if(common_arg->results != 0){
	FILE *fp;
	if( (fp = Fopen((char *)common_arg->results, "a")) == NULL ) {
	  ERR.FileA(cname,fname, (char *)common_arg->results);
	}
	Fprintf(fp, "%0.16e %0.16e %0.16e %d %0.16e\n", 
		IFloat(cg_arg->mass), 
		IFloat(pbp), 
		IFloat(pbg5p), 
		iter, 
		IFloat(true_res));
	Fclose(fp);
      }

      // If there is another mass loop iteration ahead, we should
      // set the cg_arg->mass to it's next desired value
      if( m < pbp_arg->n_masses - 1 ) {
        switch( pbp_arg->pattern_kind ) {
	case ARRAY: 
	  cg_arg->mass = pbp_arg->mass[m+1]; 
	  break;
	case LIN:   
	  cg_arg->mass += pbp_arg->mass_step; 
	  break;
	case LOG:   
	  cg_arg->mass *= pbp_arg->mass_step; 
	  break;
        }
      }

    }

    // free the gamma_5 * solution
    VRB.Sfree(cname,fname, "sol_g5", sol_g5);
    sfree(sol_g5);
  }

  //----------------------------------------------------------------
  // Staggered fermions
  //----------------------------------------------------------------
  else if ( (lat.Fclass() == F_CLASS_STAG)
	   || (lat.Fclass() == F_CLASS_P4)
	   || (lat.Fclass() == F_CLASS_ASQTAD)   ) {  

    Vector*  alt_sol = (Vector *) smalloc(f_size * sizeof(Float));
    if(alt_sol == 0)
      ERR.Pointer(cname,fname, "alt_sol");
    VRB.Smalloc(cname,fname, "alt_sol", alt_sol, f_size * sizeof(Float));

    // set source
#ifdef POINT
    bzero((char *)src, f_size*sizeof(Float));
    if(CoorX()==0 && CoorY()==0 && CoorZ()==0 && CoorT()==0)
      {
	*((IFloat *) &src[0]) = 1.0;
	printf("Setting point source at origin\n");
      }
    for(int zz = 0; zz < GJP.VolNodeSites(); zz++)
      for(int yy = 0; yy < 6; yy++)
	if(*(((IFloat *)&src[zz])+yy) > 1e-15)
	  printf("zz = %d, yy = %d, *(src+zz) = %0.16e\n", zz, yy, *(((IFloat *)&src[zz])+yy));
#else
    lat.RandGaussVector(src, 0.5);
#endif

#ifdef Z2
    IFloat * tmp1;
    for(int zz = 0; zz < GJP.VolNodeSites(); zz++)
      for(int yy = 0; yy < 6; yy+=2)
	{
	  tmp1 = (IFloat *)(src + zz);
	  if(*(tmp1+yy) > 0)
	    *(tmp1+yy) = 1.0;
	  else
	    *(tmp1+yy) = -1.0;
	  *(tmp1 +yy+ 1) = 0.0;
	}
    //printf("Setting Z_2 source.\n");
    //int check_site = 29;
    //printf("Site %d =\n",check_site);
    //tmp1 = (IFloat *)(src + check_site);
    //for(int yy = 0; yy < 6; yy++)
    //  printf("%0.16e  ",*(tmp1+yy));
    //printf("\n");
#endif

    // Initialize the cg_arg mass, with the first mass we
    // want to compute for:
    switch( pbp_arg->pattern_kind ) {
    case ARRAY: 
      cg_arg->mass = pbp_arg->mass[0]; 
      break;
    case LIN:   
      cg_arg->mass = pbp_arg->mass_start; 
      break;
    case LOG:   
      cg_arg->mass = pbp_arg->mass_start; 
      break;
    default: 
      ERR.General(cname, fname,
		  "pbp_arg->kind = %d is unrecognized\n", 
		  pbp_arg->pattern_kind);
      break;
    }

    // Loop over masses
    for(int m=0; m<pbp_arg->n_masses; m++){
      
      // initialize solution (initial guess) to 1
      for(int i=0; i< f_size/2; i++){
	f_sol[2*i] = 1.0;     // real part
	f_sol[2*i+1] = 0.0;   // imaginary part
      }

      // do inversion
      iter = lat.FmatInv(sol, src, cg_arg, &true_res, CNV_FRM_YES);

      // Modified for anisotropic lattices
      // Calculate the pbp normalization factor
#ifdef POINT
      pbp_norm = 1.0;
#else
      pbp_norm = GJP.VolSites() * ( lat.FsiteSize() / 2 );
#endif

      Float norm = src->ReDotProductGlbSum4D(src, f_size) / pbp_norm;

      lat.Fdslash(alt_sol, sol, cg_arg, CNV_FRM_YES, 1);
      pbd0p = alt_sol->ReDotProductGlbSum4D(src, f_size) 
	/ (pbp_norm * GJP.XiVXi());

      pbp_norm *= GJP.XiV()/GJP.XiBare(); 
      
      // calculate pbp and pbdip
      pbp = sol->ReDotProductGlbSum4D(src, f_size) / pbp_norm * 2;
      pbdip = (norm- (cg_arg->mass)*pbp 
	       - GJP.XiVXi()*pbd0p)*GJP.XiBare()/GJP.XiV();
        
      // Print out mass, pbp number of iterations and true residual
      if(common_arg->results != 0){
	FILE *fp;
	if( (fp = Fopen((char *)common_arg->results, "a")) == NULL ) {
	  ERR.FileA(cname,fname, (char *)common_arg->results);
	}
	Fprintf(fp, "%0.16e %0.16e %0.16e %0.16e %d %0.16e\n", 
		IFloat(cg_arg->mass), 
		IFloat(pbp), IFloat(pbd0p), IFloat(pbdip), 
		iter, 
		IFloat(true_res));
	Fclose(fp);
      }

      // If there is another mass loop iteration ahead, we should
      // set the cg_arg->mass to it's next desired value
      if( m < pbp_arg->n_masses - 1 ) {
        switch( pbp_arg->pattern_kind ) {
	case ARRAY: 
	  cg_arg->mass = pbp_arg->mass[m+1]; 
	  break;
	case LIN:   
	  cg_arg->mass += pbp_arg->mass_step; 
	  break;
	case LOG:   
	  cg_arg->mass *= pbp_arg->mass_step; 
	  break;
        }
      }

    }
    VRB.Sfree(cname,fname, "alt_sol", alt_sol);
    sfree(alt_sol);
    
  }

  //----------------------------------------------------------------
  // Unknown fermion type
  //----------------------------------------------------------------
  else {
    ERR.General(cname,fname,"Unknown class type %d\n",int(lat.Fclass()));
  }

}




//------------------------------------------------------------------
//! Run the algorithm using point sources and domain wall fermions.
/*!
  The algorithm runs roughly as follows:

-#  Create a vector R zero everywhere except at the site specified where
  only the spin \a s and colour \a c component is 1.
-#  Solve M psi = R for psi where M is the fermion matrix. The initial guess
  for psi is a vector with every complex component equal to 1.
-#  Compute the normalised real part of the dot product <R,psi>.
-#  Compute the normalised real part of the dot product <R, gamma_5 psi>.
-# Repeat from step 1 for all other possible values of \a s and \a c.
-# Sum the dot products over all spin and colour values.

\param x the x coordinate of the point source.
\param y the y coordinate of the point source.
\param z the z coordinate of the point source.
\param t the t coordinate of the point source.
\post The results are written to the file specified in the common_arg
structure,
 */
//------------------------------------------------------------------
void AlgPbp::runPointSource(int x, int y, int z, int t)
{
#if TARGET==cpsMPI
    using MPISCU::fprintf;
#endif
  int i;
  int iter;
  int ls;
  int ls_glb;
  Float pbp, pbptmp;
  Float pbg5p, pbg5ptmp;
  Float pbp_norm;
  Float true_res;
  PbpArg *pbp_arg;
  CgArg cg_arg_struct;
  CgArg *cg_arg = &cg_arg_struct;
  char *fname = "runPointSource()";
  VRB.Func(cname,fname);

  // Set the Lattice pointer pbp_arg and cg_arg
  //----------------------------------------------------------------
  Lattice& lat = AlgLattice();
  pbp_arg = alg_pbp_arg;
  cg_arg->max_num_iter = pbp_arg->max_num_iter;
  cg_arg->stop_rsd = pbp_arg->stop_rsd;
  
  // Make a Float pointer to sol
  //----------------------------------------------------------------
  Float *f_sol = (Float *) sol;

  //----------------------------------------------------------------
  // Initialize source and solution (initial guess)
  // and calculate PsiBarPsi for:
  //----------------------------------------------------------------


  //----------------------------------------------------------------
  // Domain Wall fermions
  //----------------------------------------------------------------
  if(lat.Fclass() == F_CLASS_DWF || lat.Fclass() == F_CLASS_BFM){
    ls = GJP.SnodeSites();
    ls_glb = GJP.Snodes() * GJP.SnodeSites();

    // Allocate memory for the 4-dimensional source, solution
    Vector *src_4d = (Vector *) smalloc(f_size * sizeof(Float) / ls);
    if(src_4d == 0)
      ERR.Pointer(cname,fname, "src_4d");
    VRB.Smalloc(cname,fname, "src_4d", src_4d, f_size * sizeof(Float) / ls);
    Vector *sol_4d = (Vector *) smalloc(f_size * sizeof(Float) / ls);
    if(sol_4d == 0)
      ERR.Pointer(cname,fname, "sol_4d");
    VRB.Smalloc(cname,fname, "sol_4d", sol_4d, f_size * sizeof(Float) / ls);

    // allocate space for pbp, pb_g5_p array
    Float * pbp_all = (Float *) smalloc(ls_glb * sizeof(Float));
    if(pbp_all == 0)
      ERR.Pointer(cname,fname, "pbp_all");
    VRB.Smalloc(cname,fname, "pbp_all", pbp_all, ls_glb * sizeof(Float));

    Float * pbg5p_all = (Float *) smalloc(ls_glb * sizeof(Float));
    if(pbg5p_all == 0)
      ERR.Pointer(cname,fname, "pbg5p_all");
    VRB.Smalloc(cname,fname, "pbg5p_all", pbg5p_all, ls_glb * sizeof(Float));

    // allocate space for pbp, pb_g5_p temporary array
    Float * pbp_tmp = (Float *) smalloc(ls_glb * sizeof(Float));
    if(pbp_tmp == 0)
      ERR.Pointer(cname,fname, "pbp_tmp");
    VRB.Smalloc(cname,fname, "pbp_tmp", pbp_tmp, ls_glb * sizeof(Float));

    Float * pbg5p_tmp = (Float *) smalloc(ls_glb * sizeof(Float));
    if(pbg5p_all == 0)
      ERR.Pointer(cname,fname, "pbg5p_tmp");
    VRB.Smalloc(cname,fname, "pbg5p_tmp", pbg5p_tmp, ls_glb * sizeof(Float));

    // Initialize the cg_arg mass, with the first mass we
    // want to compute for:
    switch( pbp_arg->pattern_kind ) {
    case ARRAY: 
      cg_arg->mass = pbp_arg->mass[0]; 
      break;
    case LIN:   
      cg_arg->mass = pbp_arg->mass_start; 
      break;
    case LOG:   
      cg_arg->mass = pbp_arg->mass_start; 
      break;
    default: 
      ERR.General(cname, fname,
		  "pbp_arg->pattern_kind = %d is unrecognized\n", 
		  pbp_arg->pattern_kind);
      break;
    }
    
    // Loop over masses
    for(int m=0; m<pbp_arg->n_masses; m++){
	  
      pbptmp = (Float) 0.0;
      pbg5ptmp = (Float) 0.0;

      for (i = 0; i < ls_glb; i++) {
	pbp_tmp[i] = (Float) 0.0;
	pbg5p_tmp[i] = (Float) 0.0;
      }
      
      // initialize 4-dimensional source 
      for (int color = 0; color < GJP.Colors(); color++)
	  for (int spin = 0; spin < 4; spin++) {
	    
	    // trap for wrong arguments
	    if (x < 0 || x >= GJP.Xnodes() * GJP.XnodeSites() ||
		y < 0 || y >= GJP.Ynodes() * GJP.YnodeSites() ||
		z < 0 || z >= GJP.Znodes() * GJP.ZnodeSites() ||
		t < 0 || t >= GJP.Tnodes() * GJP.TnodeSites())
	      ERR.General(cname, fname, 
		 "Coordonate arguments out of range: x=%d, y=%d, z=%d, t=%d\n",
			  x, y, z, t);
	
	    if (color < 0 || color >= GJP.Colors())
	      ERR.General(cname, fname,
			  "Color index out of range: color = %d\n", color);
	  
	    if (spin < 0 || spin > 3)
	      ERR.General(cname, fname,
			  "Spin index out of range: spin = %d\n", spin);
	  
	    // zero the vector
	    int fv_size = GJP.VolNodeSites() * GJP.Colors() * 8;
	    for ( i = 0; i < fv_size; i++)
	      *((IFloat *)src_4d + i) = 0;
	    
	    // set point source
	    int procCoorX = x / GJP.XnodeSites();
	    int procCoorY = y / GJP.YnodeSites();
	    int procCoorZ = z / GJP.ZnodeSites();
	    int procCoorT = t / GJP.TnodeSites();
	    int localX = x % GJP.XnodeSites();
	    int localY = y % GJP.YnodeSites();
	    int localZ = z % GJP.ZnodeSites();
	    int localT = t % GJP.TnodeSites();
	    
	    int coor_x = GJP.XnodeCoor();
	    int coor_y = GJP.YnodeCoor();
	    int coor_z = GJP.ZnodeCoor();
	    int coor_t = GJP.TnodeCoor();
	    
	    if (coor_x == procCoorX &&
		coor_y == procCoorY &&
		coor_z == procCoorZ &&
		coor_t == procCoorT)
	      *((IFloat *)src_4d + 2 * (color + GJP.Colors() * (spin + 4 * (
		  localX + GJP.XnodeSites() * (
		  localY + GJP.YnodeSites() * (
		  localZ + GJP.ZnodeSites() * localT)))))) = 1.0;
   
	    // set the 5-dimensional source
	    lat.Ffour2five(src, src_4d, pbp_arg->src_u_s, pbp_arg->src_l_s);
	    
	    // initialize 5-dimensional solution (initial guess) to 1
	    for(i=0; i< f_size/2; i++){
	      f_sol[2*i] = 1.0;     // real part
	      f_sol[2*i+1] = 0.0;   // imaginary part
	    }
	    
	    // do inversion
	    iter = lat.FmatInv(sol, src, cg_arg, &true_res, CNV_FRM_YES);
	    
	    // Calculate the pbp normalization factor
	    pbp_norm = (5.0 - GJP.DwfHeight());
	  
	    if (pbp_arg->snk_loop) {
	      // Loop over sink - source separation
	      for ( i = 0; i < ls_glb; i++) {
		// set the 4-dimensional solution
		lat.Ffive2four(sol_4d, sol, i, ls_glb - i - 1);
		
		// Calculate pbp
		pbp_all[i] = sol_4d->ReDotProductGlbSum4D(src_4d, f_size/ls)
		  / pbp_norm;
		
		pbp_tmp[i] += pbp_all[i];
		
		// Calculate pbg5p = Tr[ PsiBar * Gamma5 * Psi]
		lat.Gamma5(sol_4d, sol_4d, GJP.VolNodeSites());
		pbg5p_all[i] = sol_4d->ReDotProductGlbSum4D(src_4d, f_size/ls)
		  / pbp_norm;
		
		pbg5p_tmp[i] += pbg5p_all[i];
		
	      }
	    } 
	    else {
	      // set the 4-dimensional solution
	      lat.Ffive2four(sol_4d, sol, pbp_arg->snk_u_s, pbp_arg->snk_l_s);
	      
	      // Calculate pbp
	      pbp = sol_4d->ReDotProductGlbSum4D(src_4d, f_size/ls) / pbp_norm;
	      pbptmp += pbp;
	      
	      // Calculate pbg5p = Tr[ PsiBar * Gamma5 * Psi]
	      lat.Gamma5(sol_4d, sol_4d, GJP.VolNodeSites());
	      pbg5p = sol_4d->ReDotProductGlbSum4D(src_4d, f_size/ls) / pbp_norm;
	      pbg5ptmp += pbg5p;
	      
	    }
	  }
      // Open file for results
      if(common_arg->results != 0){
	FILE *fp;
	if( (fp = Fopen((char *)common_arg->results, "a")) == NULL ) {
	  ERR.FileA(cname,fname, (char *)common_arg->results);
	}
	if (pbp_arg->snk_loop) {
	  // Print out mass, s-slice, pbp ,pbg5p, number of iterations
	  // and true residual
	  for ( i = 0; i < ls_glb; i++) {
	    Fprintf(fp,"%0.16e %d %0.16e %0.16e %d %0.16e\n", 
		    IFloat(cg_arg->mass), 
		    i,
		    IFloat(pbp_tmp[i]/4/GJP.Colors()), 
		    IFloat(pbg5p_tmp[i]/4/GJP.Colors()), 
		    iter,
		    IFloat(true_res));
	  }
	} 
	else {
	  // Print out mass, pbp, pbg5p, number of iterations and true residual
	  Fprintf(fp, "%0.16e %0.16e %0.16e %d %0.16e\n", 
		  IFloat(cg_arg->mass),
		  IFloat(pbptmp/4/GJP.Colors()), 
		  IFloat(pbg5ptmp/4/GJP.Colors()), 
		  iter, 
		  IFloat(true_res));
	}
	Fclose(fp);
      }
      
      // If there is another mass loop iteration ahead, we should
      // set the cg_arg->mass to it's next desired value
      if( m < pbp_arg->n_masses - 1 ) {
	switch( pbp_arg->pattern_kind ) {
	case ARRAY: 
	  cg_arg->mass = pbp_arg->mass[m+1]; 
	  break;
	case LIN:   
	  cg_arg->mass += pbp_arg->mass_step; 
	  break;
	case LOG:   
	  cg_arg->mass *= pbp_arg->mass_step; 
	  break;
	}
      }
    }
    
    // free the 4-dimensional source, solution, pbp, and pbg5p
    VRB.Sfree(cname,fname, "pbg5p_tmp", pbg5p_tmp);
    sfree(pbg5p_tmp);
    VRB.Sfree(cname,fname, "pbg5p_all", pbg5p_all);
    sfree(pbg5p_all);
    VRB.Sfree(cname,fname, "pbp_tmp", pbp_tmp);
    sfree(pbp_tmp);
    VRB.Sfree(cname,fname, "pbp_all", pbp_all);
    sfree(pbp_all);
    VRB.Sfree(cname,fname, "sol_4d", sol_4d);
    sfree(sol_4d);
    VRB.Sfree(cname,fname, "src_4d", src_4d);
    sfree(src_4d);

  }

  //----------------------------------------------------------------
  // Other fermion type
  //----------------------------------------------------------------
  else {
    ERR.General(cname,fname,"Not implemented for class type %d\n",int(lat.Fclass()));
  }
  
}

CPS_END_NAMESPACE
