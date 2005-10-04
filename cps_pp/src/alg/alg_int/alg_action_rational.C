#include<config.h>
#include<math.h>
CPS_START_NAMESPACE 
//------------------------------------------------------------------
//
// alg_action_rational.C
//
// AlgActionRational is a bilinear action where the matrix is represented
// by a rational approximation.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include<alg/alg_hmd.h>
#include<util/lattice.h>
#include<util/vector.h>
#include<util/gjp.h>
#include<util/smalloc.h>
#include<util/verbose.h>
#include<util/error.h>
#include<alg/alg_int.h>
#include<alg/alg_remez.h>
CPS_START_NAMESPACE

//!< Dummy contructor - does nothing
AlgActionRational::AlgActionRational()
  : AlgActionBilinear()
{

}

//!< Need to add restart method or similar so do not need to
//reconstruct actions - add calls to restart methods to hmc class

AlgActionRational::AlgActionRational(AlgMomentum &mom,
				     ActionRationalArg &r_arg)
				     
  : AlgActionBilinear(mom, r_arg.bi_arg)
{

  cname = "AlgActionRational";
  char *fname = "AlgActionRational()";

  rat_arg = &r_arg;

  //!< First check n_masses bilinear = n_masses rational
  if (rat_arg->rationals.rationals_len != 
      rat_arg->bi_arg.bilinears.bilinears_len)
    ERR.General(cname, fname,
		"Inconsistency between RationalArg and BilinearArg n_masses");

  // Rational force term not implemented dim < 4
  if( fermion == F_CLASS_ASQTAD && 
      (GJP.XnodeSites()==2 || GJP.YnodeSites()==2 ||
       GJP.ZnodeSites()==2 || GJP.TnodeSites()==2 ) )
    ERR.General(cname,fname,
		"Force not implemented for Fasqtad with any dimension < 4\n");

  //!< Allocate memory for the fermion CG arguments.
  if(n_masses > 0){
    remez_arg_md = (RemezArg*) smalloc(cname, fname, "remez_arg_md",
				       n_masses * sizeof(RemezArg));
    
    remez_arg_mc = (RemezArg*) smalloc(cname, fname, "remez_arg_mc",
				       n_masses * sizeof(RemezArg));
    
    for(int i=0; i<n_masses; i++) {
      remez_arg_md[i].degree = 
	rat_arg->rationals.rationals_val[i].md_approx.md_approx_len;
      remez_arg_md[i].field_type = 
	rat_arg->rationals.rationals_val[i].field_type;
      remez_arg_md[i].lambda_low = 
	rat_arg->rationals.rationals_val[i].lambda_low;
      remez_arg_md[i].lambda_high = 
	rat_arg->rationals.rationals_val[i].lambda_high;
      remez_arg_md[i].power_num = 
	rat_arg->rationals.rationals_val[i].power_num;
      remez_arg_md[i].power_den = 
	rat_arg->rationals.rationals_val[i].power_den;
      remez_arg_md[i].precision = rat_arg->precision;
      remez_arg_md[i].valid_approx = 0;
      
      remez_arg_mc[i].degree = 
	rat_arg->rationals.rationals_val[i].mc_approx.mc_approx_len;
      remez_arg_mc[i].field_type = 
	rat_arg->rationals.rationals_val[i].field_type;
      remez_arg_mc[i].lambda_low = 
	rat_arg->rationals.rationals_val[i].lambda_low;
      remez_arg_mc[i].lambda_high = 
	rat_arg->rationals.rationals_val[i].lambda_high;
      remez_arg_mc[i].power_num = 
	rat_arg->rationals.rationals_val[i].power_num;
      remez_arg_mc[i].power_den = 
	2*rat_arg->rationals.rationals_val[i].power_den;
      remez_arg_mc[i].precision = rat_arg->precision;
      remez_arg_mc[i].valid_approx = 0;
    }
    
    //!< construct approximation if necessary
    generateApprox();

    frm_cg_arg_md = (CgArg ***) smalloc (cname,fname, "frm_cg_arg_md", 
					 n_masses * sizeof(CgArg**));

    frm_cg_arg_mc = (CgArg ***) smalloc (cname,fname, "frm_cg_arg_mc", 
					 n_masses * sizeof(CgArg**));

    for(int i=0; i<n_masses; i++) {
      frm_cg_arg_md[i] = (CgArg **)smalloc(cname,fname, "frm_cg_arg_md[i]",
					   remez_arg_md[i].degree*
					   sizeof(CgArg*));
      frm_cg_arg_mc[i] = (CgArg **)smalloc(cname,fname, "frm_cg_arg_mc[i]",
					   remez_arg_mc[i].degree*
					   sizeof(CgArg*));
      for(int j=0; j<remez_arg_md[i].degree; j++)
	frm_cg_arg_md[i][j] = (CgArg *) 
	  smalloc(cname,fname, "frm_cg_arg_md[i][j]",sizeof(CgArg));

      for(int j=0; j<remez_arg_mc[i].degree; j++)
	frm_cg_arg_mc[i][j] = (CgArg *) 
	  smalloc(cname,fname, "frm_cg_arg_mc[i][j]",sizeof(CgArg));
    }

    //!< Initialize the fermion CG arguments
    for(int i=0; i<n_masses; i++){
      for (int j=0; j<remez_arg_md[i].degree; j++) {
	frm_cg_arg_md[i][j]->mass = mass[i];
	frm_cg_arg_md[i][j]->max_num_iter = max_num_iter[i];
	frm_cg_arg_md[i][j]->stop_rsd = 
	  rat_arg->rationals.rationals_val[i].md_approx.md_approx_val[j].stop_rsd;
      }
      for (int j=0; j<remez_arg_mc[i].degree; j++) {
	frm_cg_arg_mc[i][j]->mass = mass[i];
	frm_cg_arg_mc[i][j]->max_num_iter = max_num_iter[i];
	frm_cg_arg_mc[i][j]->stop_rsd = 
	  rat_arg->rationals.rationals_val[i].mc_approx.mc_approx_val[j].stop_rsd;      }
    }
  
    max_size = 0;
    total_size = 0;
    for (int i=0; i<n_masses; i++) {
      total_size += remez_arg_md[i].degree;
      if (max_size < remez_arg_md[i].degree) max_size = remez_arg_md[i].degree;
    }
    
    //!< Allocate memory for the frmn solution fermion fields.
    if (fermion == F_CLASS_DWF) {
      
      //!< For dwf we need the solution vector contiguous in memory
      frmn = (Vector**) smalloc(max_size*sizeof(Vector*), cname, fname, "frmn");
      
      frmn[0] = (Vector*) smalloc(f_size*max_size*sizeof(Float), 
				  cname, fname, "frmn[0]");
      
      for (int i=1; i<max_size; i++) frmn[i] = frmn[0] + i*f_vec_count;
      
      frmn_d = 0;
    } else {
      //!< For asqtad we need them checkerboarded with dslash applied
      frmn = (Vector**)smalloc(total_size*sizeof(Vector*), cname, fname, "frmn");
      frmn_d = (Vector**)
	smalloc(total_size*sizeof(Vector*), cname, fname, "frmn_d");
      
      for (int i=0; i<total_size; i++) {
	frmn[i] = (Vector*) smalloc(2*f_size*sizeof(Float));
	frmn_d[i] = frmn[i] + f_vec_count;
      }
      
    }

    all_res = (Float *)smalloc(total_size * sizeof(Float),
			       cname,fname, "all_res");
    
    frmn_tmp = (Vector**)smalloc(total_size*sizeof(Vector*),
				 cname,fname,"frmn_tmp");


    //!< Allocate fractionSplit
    fractionSplit = (int**)smalloc(2*sizeof(int*),cname,fname,"fractionSplit");
    fractionSplit[0] = 
      (int*)smalloc(n_masses*sizeof(int),cname,fname,"fractionSplit[0]");
    fractionSplit[1] = 
      (int*)smalloc(n_masses*sizeof(int),cname,fname,"fractionSplit[1]");

    //!< This is just a dummy parameter when called through AlgActionRational
    for (int i=0; i<n_masses; i++) {
      fractionSplit[0][i] = 0;
      fractionSplit[1][i] = remez_arg_md[i].degree;
    }

  }

  init();

}

void AlgActionRational::init() {

  AlgActionBilinear::init();
  evolved = 1;
  heatbathEval = 0;
  energyEval = 0;
    
}

AlgActionRational::~AlgActionRational() {

  char *fname = "~AlgActionRational()" ;
  VRB.Func(cname,fname);

  //!< Free memory for timescale split partial fraction
  if (n_masses > 0) {
    //!< Free dummy fractionSplit parameters
    sfree(fractionSplit[0], cname, fname, "fractionSplit[0]");
    sfree(fractionSplit[1], cname, fname, "fractionSplit[1]");
    sfree(fractionSplit, cname, fname, "fractionSplit");

    //!< Free memory for the residue coefficients
    sfree(frmn_tmp, cname,fname, "frmn_tmp");
    //!< Free memory for the residue coefficients
    sfree(all_res, cname,fname, "all_res");

    //!< Free memory for the frmn (pseudo fermion) solution fields.
    if (fermion == F_CLASS_DWF) {
      sfree(frmn[0], cname, fname, "frmn[0]");
    } else {
      sfree(frmn_d, cname,fname, "frmn_d");    
      for (int i=0; i<total_size; i++) sfree(frmn[i], cname, fname, "frmn[i]");
    }
    sfree(frmn, cname, fname, "frmn");

    //!< Free memory for the fermion CG arguments
    for(int i=0; i<n_masses; i++){
      for (int j=0; j<remez_arg_md[i].degree; j++)
	sfree(frm_cg_arg_md[i][j], cname,fname, "frm_cg_arg_md[i][j]");

      for (int j=0; j<remez_arg_mc[i].degree; j++)
	sfree(frm_cg_arg_mc[i][j], cname,fname, "frm_cg_arg_mc[i][j]");

      sfree(frm_cg_arg_md[i], cname, fname, "frm_cg_arg_md[i]");
      sfree(frm_cg_arg_mc[i], cname, fname, "frm_cg_arg_mc[i]");
    }

    sfree(frm_cg_arg_md, cname,fname, "frm_cg_arg_md");
    sfree(frm_cg_arg_mc, cname,fname, "frm_cg_arg_mc");

    //!< Must not free these until CG args are freed.
    sfree(remez_arg_md,cname,fname,"remez_arg_md");
    sfree(remez_arg_mc,cname,fname,"remez_arg_mc");

  }

}

//!< Heat Bath for the pseudo-fermions (phi)
void AlgActionRational::heatbath() {

  char *fname = "heatbath()";
  Float trueMass;

  //!< Only evaluate heatbath if necessary
  if (!heatbathEval) {

    //!< Create an appropriate lattice
    Lattice &lat = LatticeFactory::Create(fermion, G_CLASS_NONE);  
    h_init = 0.0;
    
    for(int i=0; i<n_masses; i++){

      //!< Potentially can merge all these three functions
      //!, Certainly can for 2 and 3
      lat.RandGaussVector(frmn[0], 0.5, Ncb);
      h_init += lat.FhamiltonNode(frmn[0],frmn[0]);
      phi[i] -> 
	VecEqualsVecTimesEquFloat(frmn[0], remez_arg_mc[i].norm_inv, f_size);
      
      massRenormalise(&(frm_cg_arg_mc[i][0]->mass), &trueMass, 
		      remez_arg_mc[i].degree, remez_arg_mc[i].pole_inv, 
		      RENORM_FORWARDS);
      
      cg_iter = lat.FmatEvlMInv(phi+i, frmn[0], remez_arg_mc[i].pole_inv, 
				remez_arg_mc[i].degree, 0, 
				frm_cg_arg_mc[i], CNV_FRM_NO, SINGLE,
				remez_arg_mc[i].residue_inv);
      
      massRenormalise(&(frm_cg_arg_mc[i][0]->mass), &trueMass, 
		      remez_arg_mc[i].degree, remez_arg_mc[i].pole_inv, 
		      RENORM_BACKWARDS);
      
      updateCgStats(frm_cg_arg_mc[i][0]);
    }

    LatticeFactory::Destroy();

    evolved = 0;
    heatbathEval = 1;
    energyEval = 0;

  }

}

// Calculate rhmc fermion contribution to the Hamiltonian
Float AlgActionRational::energy() {

  char *fname="energy()";

  if (energyEval) {
    return 0.0;
  } else if (!evolved) {
    energyEval = 1;
    return h_init;
  } else {
    int shift = 0;
    Float trueMass;
    Float h = 0.0;

    //!< Create an appropriate lattice
    Lattice &lat = LatticeFactory::Create(fermion, G_CLASS_NONE);  

    for (int i=0; i<n_masses; i++) {
      massRenormalise(&(frm_cg_arg_mc[i][0]->mass), &trueMass, 
		      remez_arg_mc[i].degree, remez_arg_mc[i].pole, 
		      RENORM_FORWARDS);
      
      frmn[0] -> VecEqualsVecTimesEquFloat(phi[i],remez_arg_mc[i].norm,f_size);
      
      cg_iter = lat.FmatEvlMInv(frmn, phi[i], remez_arg_mc[i].pole, 
				remez_arg_mc[i].degree, 0, 
				frm_cg_arg_mc[i], CNV_FRM_NO, SINGLE, 
				remez_arg_mc[i].residue);
      
      massRenormalise(&(frm_cg_arg_mc[i][0]->mass), &trueMass, 
		      remez_arg_mc[i].degree, remez_arg_mc[i].pole, 
		      RENORM_BACKWARDS);
      
      updateCgStats(frm_cg_arg_mc[i][0]);

      // shift this evaluation into minvcg?
      h += lat.FhamiltonNode(frmn[0], frmn[0]);
    }

    LatticeFactory::Destroy();

    energyEval = 1;

    return h;
  }

  
}

void AlgActionRational::evolve(Float dt, int nsteps) {
  evolve(dt, nsteps, fractionSplit);
}

//!< run method evolves the integrator
void AlgActionRational::evolve(Float dt, int nsteps, int **fractionSplit) 
{

  char *fname = "evolve(Float,int)";

  if (n_masses > 0) {
    Float trueMass;
    
    //!< Variables required for ASQTAD partial fraction splitting
    int total_split_degree = 0;
    for (int i=0; i<n_masses; i++)
      total_split_degree += fractionSplit[1][i] - fractionSplit[0][i];

    //!< Create an appropriate lattice
    Lattice &lat = LatticeFactory::Create(fermion, G_CLASS_NONE);  
    
    int shift = 0;
    int split_shift = 0;
    
    for(int steps = 0; steps<nsteps; steps++) {
      for(int i=0; i<n_masses; i++){
	
	int deg = fractionSplit[1][i] - fractionSplit[0][i];
	int isz = fractionSplit[0][i];
	
	if (deg > 0) {

	  massRenormalise(&(frm_cg_arg_md[i][isz]->mass), &trueMass, 
			  deg, remez_arg_md[i].pole+isz,
			  RENORM_FORWARDS);
	  
	  cg_iter = lat.FmatEvlMInv(frmn+shift+isz, phi[i], 
				    remez_arg_md[i].pole+isz, deg, isz, 
				    frm_cg_arg_md[i]+isz, CNV_FRM_NO, 
				    frmn_d+shift+isz);
	  
	  massRenormalise(&(frm_cg_arg_md[i][isz]->mass), &trueMass, 
			  deg, remez_arg_md[i].pole+isz, RENORM_BACKWARDS);
	  
	  updateCgStats(frm_cg_arg_md[i][isz]);
	  
	  if (fermion != F_CLASS_ASQTAD) {
	    lat.RHMC_EvolveMomFforce(mom, frmn+shift+isz, deg, isz,
				     remez_arg_md[i].residue+isz, mass[i],
				     dt, frmn_d+shift+isz);
	  } else {
	    //!< Do appropriate pointer arithmetic for asqtad
	    for (int j=0; j<deg; j++) {
	      all_res[split_shift+j] = remez_arg_md[i].residue[isz+j];
	      frmn_tmp[split_shift+j] = frmn[shift+isz+j];
	    }
	  }
	}

	if (fermion != F_CLASS_DWF) shift += remez_arg_md[i].degree;
	split_shift += deg;
      }
      
      //!< Only for the case of asqtad fermions do we perform this optimisation
      
      if (fermion == F_CLASS_ASQTAD && total_split_degree > 0)
	lat.RHMC_EvolveMomFforce(mom, frmn_tmp, total_split_degree, 0,
				 all_res, 0.0, dt, frmn_d);
      
      evolved = 1;
      heatbathEval = 0;
      energyEval = 0;
      md_steps++;
    }
    
    LatticeFactory::Destroy();

  }

}

/*!< Renormalise the smallest shift into the mass parameters - reduces
  linear algebra in the multi-mass solver.  This optimisation only
  works for staggered type fermions.*/
void AlgActionRational::massRenormalise(Float *mass, Float *trueMass, 
				 int degree, Float *shift, 
				 MassRenormaliseDir direction) 
{

  //!< Can only renormalise mass for staggered or asqtad cases
  if (fermion == F_CLASS_ASQTAD || fermion == F_CLASS_STAG) {
    if (direction == RENORM_FORWARDS) {
      *trueMass = *mass;
      *mass = sqrt((*trueMass)*(*trueMass) + shift[0]/4.0);
      Float zeroPole = shift[0];
      for (int j=0; j<degree; j++) shift[j] -= zeroPole;
    } else if (direction == RENORM_BACKWARDS) {
      Float zeroPole = 4.0*((*mass)*(*mass) - (*trueMass)*(*trueMass));
      for (int j=0; j<degree; j++) shift[j] += zeroPole;
      *mass = *trueMass;
    }
  }

}

//!< Generate the optimal rational approximation for rational representation
// Need to add boson/fermion optimisation - requires extra space for FIRat.
void AlgActionRational::generateApprox()
{

  //!< Construct approximations
  for (int i=0; i<n_masses; i++) {
    for (int j=0; j<i; j++) {
      //!< Avoid unnecessary reconstruction
      if (mass[j] == mass[i]) {
	remez_arg_md[i].valid_approx = 
	  compareApprox(remez_arg_md[i], remez_arg_md[j]);

	remez_arg_mc[i].valid_approx = 
	  compareApprox(remez_arg_mc[i], remez_arg_mc[j]);
      }
    }

    //!< First generate the md appropximation
    if (!remez_arg_md[i].valid_approx) {
      AlgRemez remez(remez_arg_md[i]);
      remez.generateApprox();
      remez_arg_md[i].valid_approx = 1;
    }

    //!< Now generate the mc approximation
    if (!remez_arg_mc[i].valid_approx) {
      AlgRemez remez(remez_arg_mc[i]);
      remez.generateApprox();
      remez_arg_mc[i].valid_approx = 1;
    }
    
  }
  
}

int AlgActionRational::compareApprox(RemezArg &arg1, RemezArg &arg2) {

  //!< no need to recalculate approximation if same mass
  if (arg1.field_type == arg2.field_type &&
      arg1.power_num == arg2.power_num &&
      arg1.power_den == arg2.power_den &&
      arg1.degree == arg2.degree) {
    
    arg1.norm = arg2.norm;
    arg1.norm_inv = arg2.norm_inv;
    
    for (int k=0; k<arg1.degree; k++) {
      arg1.residue[k] = arg2.residue[k];
      arg1.pole[k] = arg2.pole[k];
      arg1.residue_inv[k] = arg2.residue_inv[k];
      arg1.pole_inv[k] = arg2.pole_inv[k];
    }
    return 1;
  } else {
    return 0;
  }
  
}

CPS_END_NAMESPACE
