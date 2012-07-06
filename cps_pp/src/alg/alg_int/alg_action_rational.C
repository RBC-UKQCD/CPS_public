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
  cname = "AlgActionRational";
}

//!< Constructor called by AlgActionRationalQuotient
AlgActionRational::AlgActionRational(AlgMomentum &mom, 
				     ActionBilinearArg &bi_arg)
  : AlgActionBilinear(mom,bi_arg)
{
  cname = "AlgActionRational";  
}

//!< Need to add restart method or similar so do not need to
//reconstruct actions - add calls to restart methods to hmc class

AlgActionRational::AlgActionRational(AlgMomentum &mom,
				     ActionRationalArg &r_arg, int traj_num)
				     
  : AlgActionBilinear(mom, r_arg.bi_arg)
{

  cname = "AlgActionRational";
  char *fname = "AlgActionRational()";
  VRB.Func(cname,fname);
    
  int_type = INT_RATIONAL;
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
    
    //!< construct approximation if necessary
    generateApprox(mass, &remez_arg_md, &remez_arg_mc, 
		   rat_arg->rationals.rationals_val);

    //!< Allocate memory for the fermion CG arguments.
    generateCgArg(mass, &frm_cg_arg_fg, &frm_cg_arg_md, &frm_cg_arg_mc, "frm_cg_arg",
		  rat_arg->rationals.rationals_val);
  
    max_size = 0;
    total_size = 0;
    for (int i=0; i<n_masses; i++) {
      total_size += remez_arg_md[i].degree;
      if (max_size < remez_arg_md[i].degree) max_size = remez_arg_md[i].degree;
    }
    
    //!< Allocate memory for the frmn solution fermion fields.
    if (fermion == F_CLASS_DWF) {
      
      //!< For dwf we need the solution vector contiguous in memory
      frmn = (Vector**) smalloc(max_size*sizeof(Vector*), "frmn", fname, cname);
      
      frmn[0] = (Vector*) smalloc(f_size*max_size*sizeof(Float), 
				  "frmn[0]", fname, cname);
      
      for (int i=1; i<max_size; i++) frmn[i] = frmn[0] + i*f_vec_count;
      
      frmn_d = 0;
    } else {
      //!< For asqtad we need them checkerboarded with dslash applied
      frmn = (Vector**)smalloc(total_size*sizeof(Vector*), "frmn", fname, cname);
      frmn_d = (Vector**)
	smalloc(total_size*sizeof(Vector*), "frmn_d", fname, cname);
      
      for (int i=0; i<total_size; i++) {
	frmn[i] = (Vector*) smalloc(2*f_size*sizeof(Float), "frmn[i]", fname, cname);
	frmn_d[i] = frmn[i] + f_vec_count;
      }
      
    }

    all_res = (Float *)smalloc(total_size*sizeof(Float),"all_res",fname,cname);
    frmn_tmp = (Vector**)smalloc(total_size*sizeof(Vector*),"frmn_tmp",fname,cname);


    //!< Allocate fractionSplit
    fractionSplit = (int**)smalloc(2*sizeof(int*),"fractionSplit",fname,cname);
    fractionSplit[0] = 
      (int*)smalloc(n_masses*sizeof(int),"fractionSplit[0]",fname,cname);
    fractionSplit[1] = 
      (int*)smalloc(n_masses*sizeof(int),"fractionSplit[1]",fname,cname);

    //!< This is just a dummy parameter when called through AlgActionRational
    for (int i=0; i<n_masses; i++) {
      fractionSplit[0][i] = 0;
      fractionSplit[1][i] = remez_arg_md[i].degree;
    }

    //!< Allocate splitCheck
    splitCheck = 
      (int**)smalloc(n_masses*sizeof(int*), "splitCheck", fname, cname);
    
    for (int i=0; i<n_masses; i++) {
      splitCheck[i] = (int*)smalloc(remez_arg_md[i].degree*sizeof(int), 
				    "splitCheck[i]",fname,cname);
      for (int j=0; j<remez_arg_md[i].degree; j++) splitCheck[i][j] = 0;
    }

  }

  init(traj_num);

  if (rat_arg->eigen.eigen_measure == EIGEN_MEASURE_YES) 
    generateEigArg(rat_arg->eigen);
}

void AlgActionRational::init(int traj_num) {

  char *fname = "init(int)" ;
  VRB.Func(cname,fname);
  AlgActionBilinear::init();
  evolved = 1;
  heatbathEval = 0;
  energyEval = 0;
  traj = traj_num-1;  
}

AlgActionRational::~AlgActionRational() {

  char *fname = "~AlgActionRational()" ;
  VRB.Func(cname,fname);

  //!< Free memory for timescale split partial fraction
  if (n_masses > 0  && int_type == INT_RATIONAL) {
    //!< Free splitCheck parameters
    for (int i=0; i<n_masses; i++)
      sfree(splitCheck[i], "splitCheck[i]", fname, cname);
    sfree(splitCheck, "splitCheck", fname, cname);

    //!< Free dummy fractionSplit parameters
    sfree(fractionSplit[0], "fractionSplit[0]", fname, cname);
    sfree(fractionSplit[1], "fractionSplit[1]", fname, cname);
    sfree(fractionSplit, "fractionSplit", fname, cname);

    //!< Free memory for the residue coefficients
    sfree(frmn_tmp, "frmn_tmp",fname,cname);
    //!< Free memory for the residue coefficients
    sfree(all_res, "all_res", fname, cname);

    //!< Free memory for the frmn (pseudo fermion) solution fields.
    if (fermion == F_CLASS_DWF) {
      sfree(frmn[0], "frmn[0]", fname, cname);
    } else {
      sfree(frmn_d, "frmn_d", fname, cname);
      for (int i=0; i<total_size; i++) sfree(frmn[i], "frmn[i]", fname, cname);
    }
    sfree(frmn, "frmn", fname, cname);

    //!< Free memory for the fermion CG arguments
    destroyCgArg(frm_cg_arg_fg, frm_cg_arg_md, frm_cg_arg_mc, "frm_cg_arg",
		 remez_arg_md, remez_arg_mc);

    //!< Must not free these until CG args are freed.
    destroyApprox(remez_arg_md, remez_arg_mc);

    if (rat_arg->eigen.eigen_measure == EIGEN_MEASURE_YES) destroyEigArg();
  }

}

//!< Heat Bath for the pseudo-fermions (phi)
void AlgActionRational::heatbath() {

  char *fname = "heatbath()";

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
      
      cg_iter = lat.FmatEvlMInv(phi+i, frmn[0], remez_arg_mc[i].pole_inv, 
				remez_arg_mc[i].degree, 0, 
				frm_cg_arg_mc[i], CNV_FRM_NO, SINGLE,
				remez_arg_mc[i].residue_inv);
      
      updateCgStats(frm_cg_arg_mc[i][0]);
    }

    LatticeFactory::Destroy();

    evolved = 0;
    heatbathEval = 1;
    energyEval = 0;
  }
  traj++;

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
    Float h = 0.0;

    //!< Before energy is measured, do we want to check bounds?
    if (rat_arg->eigen.eigen_measure == EIGEN_MEASURE_YES) 
      checkApprox(mass, remez_arg_mc, rat_arg->eigen);

    //!< Create an appropriate lattice
    Lattice &lat = LatticeFactory::Create(fermion, G_CLASS_NONE);  

    for (int i=0; i<n_masses; i++) {
      frmn[0] -> VecEqualsVecTimesEquFloat(phi[i],remez_arg_mc[i].norm,f_size);
      
      cg_iter = lat.FmatEvlMInv(frmn, phi[i], remez_arg_mc[i].pole, 
				remez_arg_mc[i].degree, 0, 
				frm_cg_arg_mc[i], CNV_FRM_NO, SINGLE, 
				remez_arg_mc[i].residue);
      
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

void AlgActionRational::prepare_fg(Matrix * force, Float dt_ratio)
{
  char * fname = "prepare_fg(M*,F)";
  Lattice &lat = LatticeFactory::Create(fermion, G_CLASS_NONE);  

  //!< Variables required for ASQTAD partial fraction splitting
  int total_split_degree = 0;
  for (int i=0; i<n_masses; i++)
    total_split_degree += fractionSplit[1][i] - fractionSplit[0][i];

  int shift = 0;
  int split_shift = 0;
    
  for(int i=0; i<n_masses; i++){
	
    int deg = fractionSplit[1][i] - fractionSplit[0][i];
    int isz = fractionSplit[0][i];
	
    if (deg > 0) {
      cg_iter = lat.FmatEvlMInv(frmn+shift+isz, phi[i], 
                                remez_arg_md[i].pole+isz, deg, isz, 
                                frm_cg_arg_fg[i]+isz, CNV_FRM_NO, 
                                frmn_d+shift+isz);
	  
      updateCgStats(frm_cg_arg_fg[i][isz]);

      if (force_measure == FORCE_MEASURE_YES ||
          (fermion != F_CLASS_ASQTAD && fermion != F_CLASS_P4) ) {
        Fdt = lat.RHMC_EvolveMomFforce(force, frmn+shift+isz, deg, isz,
                                       remez_arg_md[i].residue+isz, 
                                       mass[i], dt_ratio, frmn_d+shift+isz,
                                       force_measure);

        if (force_measure == FORCE_MEASURE_YES) {
          char label[200];
          sprintf(label, "%s total, mass = %e:", force_label, mass[i]);
          Fdt.print(dt_ratio, label);
        }
      } else {
        //!< Do appropriate pointer arithmetic for asqtad/p4
        for (int j=0; j<deg; j++) {
          all_res[split_shift+j] = remez_arg_md[i].residue[isz+j];
          frmn_tmp[split_shift+j] = frmn[shift+isz+j];
        }
      }
    }

    if (fermion != F_CLASS_DWF) shift += remez_arg_md[i].degree;
    split_shift += deg;
  }

  //!< Only for the case of asqtad/p4 fermions do we perform this optimisation
  if ( (fermion == F_CLASS_ASQTAD || fermion == F_CLASS_P4) && 
       total_split_degree > 0 && force_measure == FORCE_MEASURE_NO ) {
    Fdt = lat.RHMC_EvolveMomFforce(force, frmn_tmp, total_split_degree, 0, all_res, 
                                   0.0, dt_ratio, frmn_d, force_measure);
    if (force_measure == FORCE_MEASURE_YES) {
      char label[200];
      sprintf(label, "%s total:", force_label);
      Fdt.print(dt_ratio, label);
    }
  }
  LatticeFactory::Destroy();
}

//!< run method evolves the integrator
void AlgActionRational::evolve(Float dt, int nsteps, int **fractionSplit) 
{
  char *fname = "evolve(Float, int, int **)";
  if (n_masses <= 0) return;
  Float trueMass;
    
  //!< Variables required for ASQTAD partial fraction splitting
  int total_split_degree = 0;
  for (int i=0; i<n_masses; i++)
    total_split_degree += fractionSplit[1][i] - fractionSplit[0][i];

  //!< Create an appropriate lattice
  Lattice &lat = LatticeFactory::Create(fermion, G_CLASS_NONE);  
    
  for(int steps = 0; steps<nsteps; steps++) {
    int shift = 0;
    int split_shift = 0;
    
    for(int i=0; i<n_masses; i++){
	
      int deg = fractionSplit[1][i] - fractionSplit[0][i];
      int isz = fractionSplit[0][i];
	
      if (deg > 0) {

        cg_iter = lat.FmatEvlMInv(frmn+shift+isz, phi[i], 
                                  remez_arg_md[i].pole+isz, deg, isz, 
                                  frm_cg_arg_md[i]+isz, CNV_FRM_NO, 
                                  frmn_d+shift+isz);
	  
        updateCgStats(frm_cg_arg_md[i][isz]);

        if (force_measure == FORCE_MEASURE_YES ||
            (fermion != F_CLASS_ASQTAD && fermion != F_CLASS_P4) ) {
          Fdt = lat.RHMC_EvolveMomFforce(mom, frmn+shift+isz, deg, isz,
                                         remez_arg_md[i].residue+isz, 
                                         mass[i], dt, frmn_d+shift+isz,
                                         force_measure);

          if (force_measure == FORCE_MEASURE_YES) {
            char label[200];
            sprintf(label, "%s total, mass = %e:", force_label, mass[i]);
            Fdt.print(dt, label);
          }
        } else {
          //!< Do appropriate pointer arithmetic for asqtad/p4
          for (int j=0; j<deg; j++) {
            all_res[split_shift+j] = remez_arg_md[i].residue[isz+j];
            frmn_tmp[split_shift+j] = frmn[shift+isz+j];
          }
        }
      }

      if (fermion != F_CLASS_DWF) shift += remez_arg_md[i].degree;
      split_shift += deg;
    }
      
    //!< Only for the case of asqtad/p4 fermions do we perform this optimisation
    if ( (fermion == F_CLASS_ASQTAD || fermion == F_CLASS_P4) && 
         total_split_degree > 0 && force_measure == FORCE_MEASURE_NO ) {
      Fdt = lat.RHMC_EvolveMomFforce(mom, frmn_tmp, total_split_degree, 0, all_res, 
                                     0.0, dt, frmn_d, force_measure);
      if (force_measure == FORCE_MEASURE_YES) {
        char label[200];
        sprintf(label, "%s total:", force_label);
        Fdt.print(dt, label);
      }
    }
      
    evolved = 1;
    heatbathEval = 0;
    energyEval = 0;
    md_steps++;
  }
    
  LatticeFactory::Destroy();
}

//!< Generate the optimal rational approximation for rational representation
// Need to add boson/fermion optimisation - requires extra space for FIRat.
void AlgActionRational::generateApprox(Float *mass, 
				       RemezArg **remez_arg_md, 
				       RemezArg **remez_arg_mc, 
				       RationalDescr *rat)
{
  char *fname = "generateApprox(F*,RA*,RA*,RD*)";

  *remez_arg_md = new RemezArg[n_masses];
  *remez_arg_mc = new RemezArg[n_masses];
  
  for(int i=0; i<n_masses; i++) {
    (*remez_arg_md)[i].approx_type = rat[i].md_approx.approx_type;
    switch (rat[i].md_approx.bounds_type) {
    case RATIONAL_BOUNDS_MANUAL:
      (*remez_arg_md)[i].lambda_low = rat[i].md_approx.lambda_low;
      (*remez_arg_md)[i].lambda_high = rat[i].md_approx.lambda_high;
      break;
    case RATIONAL_BOUNDS_AUTOMATIC:
      if (fermion != F_CLASS_STAG && fermion != F_CLASS_ASQTAD && fermion != F_CLASS_P4)
	ERR.General(cname,fname,"RationalApproxType %d not implemented for FclassType %d\n",
		    RATIONAL_BOUNDS_AUTOMATIC, fermion);
      (*remez_arg_md)[i].lambda_low = 4.0*mass[i]*mass[i];
      if (fermion == F_CLASS_STAG)
	(*remez_arg_md)[i].lambda_high = 16+4.*mass[i]*mass[i];
      else if (fermion == F_CLASS_ASQTAD)
	(*remez_arg_md)[i].lambda_high = 196./9.+4.*mass[i]*mass[i];
      else if (fermion == F_CLASS_P4)
	(*remez_arg_md)[i].lambda_high = 50./9.+4.*mass[i]*mass[i];
      break;
    default:
      ERR.General(cname,fname,"RationalBoundsType %d not implemented\n", rat[i].md_approx.bounds_type);
    }

    (*remez_arg_md)[i].degree = rat[i].md_approx.stop_rsd.stop_rsd_len;
    (*remez_arg_md)[i].field_type = rat[i].field_type;
    (*remez_arg_md)[i].power_num = rat[i].power_num;
    if (rat[i].field_type == BOSON && int_type == INT_RATIONAL_QUOTIENT)
      (*remez_arg_md)[i].power_den = 2*rat[i].power_den;
    else
      (*remez_arg_md)[i].power_den = rat[i].power_den;
    (*remez_arg_md)[i].precision = rat[i].precision;
    (*remez_arg_md)[i].valid_approx = 0;
    (*remez_arg_md)[i].delta_m = 4.0*rat[i].stag_bsn_mass*rat[i].stag_bsn_mass - 
      (*remez_arg_md)[i].lambda_low;
    
    (*remez_arg_mc)[i].approx_type = rat[i].mc_approx.approx_type;
    switch (rat[i].mc_approx.bounds_type) {
    case RATIONAL_BOUNDS_MANUAL:
      (*remez_arg_mc)[i].lambda_low = rat[i].mc_approx.lambda_low;
      (*remez_arg_mc)[i].lambda_high = rat[i].mc_approx.lambda_high;
      break;
    case RATIONAL_BOUNDS_AUTOMATIC:
      if (fermion != F_CLASS_STAG && fermion != F_CLASS_ASQTAD && fermion != F_CLASS_P4)
	ERR.General(cname,fname,"RationalApproxType %d not implemented for FclassType %d\n",
		    RATIONAL_BOUNDS_AUTOMATIC, fermion);
      (*remez_arg_mc)[i].lambda_low = 4.0*mass[i]*mass[i];
      if (fermion == F_CLASS_STAG)
	(*remez_arg_mc)[i].lambda_high = 16+4.*mass[i]*mass[i];
      else if (fermion == F_CLASS_ASQTAD)
	(*remez_arg_mc)[i].lambda_high = 196./9.+4.*mass[i]*mass[i];
      else if (fermion == F_CLASS_P4)
	(*remez_arg_mc)[i].lambda_high = 50./9.+4.*mass[i]*mass[i];
      break;
    default:
      ERR.General(cname,fname,"RationalBoundsType %d not implemented\n", rat[i].mc_approx.bounds_type);
    }

    (*remez_arg_mc)[i].degree = rat[i].mc_approx.stop_rsd.stop_rsd_len;
    (*remez_arg_mc)[i].field_type = rat[i].field_type;
    (*remez_arg_mc)[i].power_num = rat[i].power_num;
    (*remez_arg_mc)[i].power_den = 2*rat[i].power_den;
    (*remez_arg_mc)[i].precision = rat[i].precision;
    (*remez_arg_mc)[i].valid_approx = 0;
    (*remez_arg_mc)[i].delta_m = 4.0*rat[i].stag_bsn_mass*rat[i].stag_bsn_mass - 
      (*remez_arg_mc)[i].lambda_low;
  }

  //!< Construct approximations
  for (int i=0; i<n_masses; i++) {
    for (int j=0; j<i; j++) {
      //!< Avoid unnecessary reconstruction
//      if (mass[j] == mass[i]) {
    VRB.Result(cname,fname,"i=%d j=%d\n",i,j);
	(*remez_arg_md)[i].valid_approx = 
	  compareApprox((*remez_arg_md)[i], (*remez_arg_md)[j]);

	(*remez_arg_mc)[i].valid_approx = 
	  compareApprox((*remez_arg_mc)[i], (*remez_arg_mc)[j]);
//      }
    }

    //!< First generate the md appropximation
    if (!(*remez_arg_md)[i].valid_approx) {
      AlgRemez remez((*remez_arg_md)[i]);
      remez.generateApprox();
      (*remez_arg_md)[i].valid_approx = 1;
    }

    //!< Now generate the mc approximation
    if (!(*remez_arg_mc)[i].valid_approx) {
      AlgRemez remez((*remez_arg_mc)[i]);
      remez.generateApprox();
      (*remez_arg_mc)[i].valid_approx = 1;
    }
    
  }
  
}

//!< Free up the rational approximations
void AlgActionRational::destroyApprox(RemezArg *remez_arg_md, 
				      RemezArg *remez_arg_mc) 
{
  char *fname = "destroyApprox(RemezArg*,RemezArg*)";

  delete[] remez_arg_md;
  delete[] remez_arg_mc;
}

int AlgActionRational::compareApprox(RemezArg &arg1, RemezArg &arg2) {


#if 0
 printf("field type: %d %d\n",arg1.field_type,arg2.field_type);
 printf("power_num: %d %d\n",arg1.power_num,arg2.power_num);
 printf("power_den: %d %d\n",arg1.power_den,arg2.power_den);
 printf("lambda_low: %e %e\n",arg1.lambda_low,arg2.lambda_low);
 printf("lambda_high: %e %e\n",arg1.lambda_high,arg2.lambda_high);
 printf("precision: %d %d\n",arg1.precision,arg2.precision);
 printf("degree: %d %d\n",arg1.degree,arg2.degree);
#endif

  //!< no need to recalculate approximation if same mass
  if (arg1.field_type == arg2.field_type &&
      arg1.power_num == arg2.power_num &&
      arg1.power_den == arg2.power_den &&
      arg1.lambda_low == arg2.lambda_low &&
      arg1.lambda_high == arg2.lambda_high &&
      arg1.precision == arg2.precision &&
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

void AlgActionRational::generateCgArg(Float *mass,
                                      CgArg ****cg_arg_fg,
                                      CgArg ****cg_arg_md, 
				      CgArg ****cg_arg_mc, const char *label, 
				      RationalDescr *rat)
{

  char *fname = "generateCgArg(F*,Cg****,Cg***,Cg***,char*,RationalDescr*)";
  char fg_label[100], fg_label_i[100], fg_label_ij[100];
  char md_label[100], md_label_i[100], md_label_ij[100];
  char mc_label[100], mc_label_i[100], mc_label_ij[100];

  sprintf(fg_label, "%s_fg", label);
  sprintf(fg_label_i, "%s_fg[i]", label);
  sprintf(fg_label_ij, "%s_fg[i][j]", label);
  sprintf(md_label, "%s_md", label);
  sprintf(md_label_i, "%s_md[i]", label);
  sprintf(md_label_ij, "%s_md[i][j]", label);
  sprintf(mc_label, "%s_mc", label);
  sprintf(mc_label_i, "%s_mc[i]", label);
  sprintf(mc_label_ij, "%s_mc[i][j]", label);

  (*cg_arg_fg) = (CgArg ***) smalloc(n_masses*sizeof(CgArg**),fg_label,fname,cname);
  (*cg_arg_md) = (CgArg ***) smalloc(n_masses*sizeof(CgArg**),md_label,fname,cname);
  (*cg_arg_mc) = (CgArg ***) smalloc(n_masses*sizeof(CgArg**),mc_label,fname,cname);
				  

  for(int i=0; i<n_masses; i++) {
    (*cg_arg_fg)[i] = (CgArg**)smalloc(rat[i].md_approx.stop_rsd.stop_rsd_len*
				       sizeof(CgArg*),fg_label_i,fname,cname);

    (*cg_arg_md)[i] = (CgArg**)smalloc(rat[i].md_approx.stop_rsd.stop_rsd_len*
				       sizeof(CgArg*),md_label_i,fname,cname);
				    
    (*cg_arg_mc)[i] = (CgArg**)smalloc(rat[i].mc_approx.stop_rsd.stop_rsd_len*
				       sizeof(CgArg*),mc_label_i,fname,cname);

    for (int j=0; j<rat[i].md_approx.stop_rsd.stop_rsd_len; j++) {
      // currently force gradient step is using the same CG arg as the
      // normal MD step, except that the stopping condition is different.
      (*cg_arg_fg)[i][j] = (CgArg*)smalloc(sizeof(CgArg),fg_label_ij,fname,cname);
      (*cg_arg_fg)[i][j]->mass = mass[i];
      (*cg_arg_fg)[i][j]->max_num_iter = max_num_iter[i];
      (*cg_arg_fg)[i][j]->stop_rsd = rat[i].stop_rsd_fg_mult *
	rat[i].md_approx.stop_rsd.stop_rsd_val[j];

      (*cg_arg_md)[i][j] = (CgArg*)smalloc(sizeof(CgArg),md_label_ij,fname,cname);
      (*cg_arg_md)[i][j]->mass = mass[i];
      (*cg_arg_md)[i][j]->max_num_iter = max_num_iter[i];
      (*cg_arg_md)[i][j]->stop_rsd = 
	rat[i].md_approx.stop_rsd.stop_rsd_val[j];
    }

    for (int j=0; j<rat[i].mc_approx.stop_rsd.stop_rsd_len; j++) {
      (*cg_arg_mc)[i][j] = (CgArg*)smalloc(sizeof(CgArg),mc_label_ij,fname,cname);
      (*cg_arg_mc)[i][j]->mass = mass[i];
      (*cg_arg_mc)[i][j]->max_num_iter = max_num_iter[i];
      (*cg_arg_mc)[i][j]->stop_rsd = 
	rat[i].mc_approx.stop_rsd.stop_rsd_val[j];
    }

  }
}

void AlgActionRational::destroyCgArg(CgArg ***cg_arg_fg,
                                     CgArg ***cg_arg_md,
                                     CgArg ***cg_arg_mc,
				     const char *label, RemezArg *remez_arg_md,
				     RemezArg *remez_arg_mc) {

  char *fname = "destroyCgArg(Cg***, Cg***,Cg***,char*,RemezArg*,RemezArg*)";
  char fg_label[100], fg_label_i[100], fg_label_ij[100];
  char md_label[100], md_label_i[100], md_label_ij[100];
  char mc_label[100], mc_label_i[100], mc_label_ij[100];

  sprintf(fg_label, "%s_fg", label);
  sprintf(fg_label_i, "%s_fg[i]", label);
  sprintf(fg_label_ij, "%s_fg[i][j]", label);

  sprintf(md_label, "%s_md", label);
  sprintf(md_label_i, "%s_md[i]", label);
  sprintf(md_label_ij, "%s_md[i][j]", label);

  sprintf(mc_label, "%s_mc", label);
  sprintf(mc_label_i, "%s_mc[i]", label);
  sprintf(mc_label_ij, "%s_mc[i][j]", label);

  for (int i=0; i<n_masses; i++) {
    for (int j=0; j<remez_arg_mc[i].degree; j++) {
      sfree(cg_arg_mc[i][j], mc_label_ij, fname, cname);
    }
    for (int j=0; j<remez_arg_md[i].degree; j++) {
      sfree(cg_arg_fg[i][j], fg_label_ij, fname, cname);
      sfree(cg_arg_md[i][j], md_label_ij, fname, cname);
    }
    sfree(cg_arg_mc[i], mc_label_i, fname, cname);
    sfree(cg_arg_md[i], md_label_i, fname, cname);
    sfree(cg_arg_fg[i], fg_label_i, fname, cname);
  }
  sfree(cg_arg_mc, mc_label, fname, cname);
  sfree(cg_arg_md, md_label, fname, cname);
  sfree(cg_arg_fg, fg_label, fname, cname);
}

void AlgActionRational::generateEigArg(EigenDescr eigen) {
  char *fname = "generateEigArg()";

  //!< Setup AlgEig parameters if necessary
  eig_arg.pattern_kind = ARRAY;
  eig_arg.Mass.Mass_len = n_masses;
  eig_arg.Mass.Mass_val = 
    (Float*) smalloc(n_masses*sizeof(Float),"Mass_val", fname, cname);
  eig_arg.N_eig = 1;
  eig_arg.Kalk_Sim = 0;
  eig_arg.MaxCG = eigen.max_num_iter;
  eig_arg.RsdR_a = eigen.stop_rsd;
  eig_arg.RsdR_r = eigen.stop_rsd;
  eig_arg.Rsdlam = eigen.stop_rsd;
  eig_arg.Cv_fact =   0.0;
  eig_arg.N_min = 0;
  eig_arg.N_max = 0;
  eig_arg.N_KS_max = 0;
  eig_arg.n_renorm = 100;
  eig_arg.ProjApsiP = 0;
  eig_arg.print_hsum = 0;
  eig_arg.hsum_dir = 0;
  eig_arg.ncorr = 0;
  eig_arg.fname = 0;
  
  lambda_low = (Float**)smalloc(eig_arg.N_eig*sizeof(Float*),
				"lambda_low",fname,cname);
  lambda_high = (Float**)smalloc(eig_arg.N_eig*sizeof(Float*),
				 "lambda_high",fname,cname);
  for (int i=0; i<eig_arg.N_eig; i++) {
    lambda_low[i] = (Float*)smalloc(n_masses*sizeof(Float),
				    "lambda_low[i]",fname,cname);
    lambda_high[i] = (Float*)smalloc(n_masses*sizeof(Float),
				     "lambda_high[i]",fname,cname);
  }
}

void AlgActionRational::destroyEigArg() {
  char *fname = "destroyEigArg()";

  for (int i=0; i<eig_arg.N_eig; i++) {
    sfree(lambda_high[i], "lambda_low[i]", fname, cname);
    sfree(lambda_low[i], "lambda_high[i]", fname, cname);
  }
  sfree(lambda_low,"lambda_low", fname, cname);
  sfree(lambda_high,"lambda_high",fname, cname);
  
  sfree(eig_arg.Mass.Mass_val, "Mass_val", fname, cname);

}

//!< Set mass i, pole j as being included (used when splitting time scales)
void AlgActionRational::setSplit(int i, int j) {
  char *fname = "setSplit(int, int)";
  if (i<0 || i>n_masses) ERR.General(cname, fname, "Invalid mass\n");
  if (j<0 || j>=remez_arg_md[i].degree)
    ERR.General(cname, fname, "Invalid split parameters: mass %d\n",i);

  if (!splitCheck[i][j]) splitCheck[i][j] = 1;
  else 
    ERR.General(cname, fname, "Mass %d, pole %d has already been included\n",
		i, j);
}

//!< Check that all of the partial fractions have been accounted for
void AlgActionRational::checkSplit() {
  char *fname = "checkSplit()";
  for (int i=0; i<n_masses; i++) {
    for (int j=0; j<remez_arg_md[i].degree; j++) {
      if (!splitCheck[i][j]) 
	ERR.General(cname, fname, "Mass %d, pole %d has not been included\n", 
		    i, j);
    }
  }

}


//!< Check that the approximation bounds are still valid for the mc approx
void AlgActionRational::checkApprox(Float *mass, RemezArg *remez_arg, 
				    EigenDescr eigen) 
{

  char *fname = "checkApprox()";
  
  Lattice &lat = LatticeFactory::Create(fermion, G_CLASS_NONE);
  
  //!< First setup the masses
  for (int i=0; i<n_masses; i++) eig_arg.Mass.Mass_val[i] = mass[i];

  {
    //!< Measure the lowest eigenvalue
    sprintf(eig_file,"%s.%d",eigen.eig_lo_stem,traj);
    eig_arg.fname = eig_file;
    // CJ: rescaling RsdR_a to be in range with rational approximation 
    eig_arg.RsdR_a = eig_arg.RsdR_r * remez_arg[0].lambda_high;
    eig_arg.RitzMatOper = MATPCDAG_MATPC;
    
    AlgEig eig(lat,&ca_eig,&eig_arg);
    eig.run(lambda_low);

    for (int i=0; i<n_masses; i++) {
      if (lambda_low[0][i] < remez_arg[i].lambda_low) {
	ERR.General(cname, fname, 
		    "Lower bound exceeded: mass[%d] = %f, %e < %e\n", 
		    i, mass[i], lambda_low[0][i], remez_arg[i].lambda_low);
      } else {
	VRB.Result(cname,fname,"Lower bound valid: mass[%d] = %f, %e > %e\n", 
		   i, mass[i], lambda_low[0][i], remez_arg[i].lambda_low);
      }
    }

  }
  
  {
    //!< Measure the highest eigenvalue
    sprintf(eig_file,"%s.%d",eigen.eig_hi_stem,traj);
    eig_arg.fname = eig_file;    
    eig_arg.RitzMatOper = NEG_MATPCDAG_MATPC;
    
    AlgEig eig(lat,&ca_eig,&eig_arg);
    eig.run(lambda_high);
    
    for (int i=0; i<n_masses; i++) {
      lambda_high[0][i] *= -1.0;
      if (lambda_high[0][i] > remez_arg[i].lambda_high) {
	ERR.General(cname, fname, 
		    "Upper bound exceeded: mass[%d] = %f, %e > %e\n", 
		    i, mass[i], lambda_high[0][i], remez_arg[i].lambda_high);
      } else {
	VRB.Result(cname,fname,"Upper bound valid: mass[%d] = %f, %e < %e\n", 
		   i, mass[i], lambda_high[0][i], remez_arg[i].lambda_high);
      }
    }
    
  }
  
  LatticeFactory::Destroy();
}

CPS_END_NAMESPACE
