#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of Fmobius class.

*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/f_mobius/f_mobius.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// f_zmobius.C
//
// Fzmobius is derived from FwilsonTypes and is relevant to
// complexified mobius domain wall fermions
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <config.h>
#include <util/lattice.h>
#include <util/dirac_op.h>
#include <util/time_cps.h>
#include <util/enum_func.h>
#include <util/dwf.h>
#include <util/zmobius.h> 
CPS_START_NAMESPACE

using cps::Complex;

#if 0
void zTimesV1PluszTimesV2(Complex *a, Complex b, const Complex *c,
		    cps::Complex d, const cps::Complex *e, int len)
{
//#pragma omp parallel for
    for(int i = 0; i < len/2; ++i) {
    	a[i] = b * c[i] + d*e[i];
    }
}
#endif

Fzmobius::Fzmobius() : FdwfBase(){
  cname = "Fzmobius";

//  ((Zmobus*)f_dirac_op_init_ptr)->zmobius_kappa_b=0;
//  ((Zmobus*)f_dirac_op_init_ptr)->zmobius_kappa_c=0;
//  ((Zmobus*)f_dirac_op_init_ptr)->zmobius_kappa_ratio=0;

  // default preconditioning
  ((Zmobus*)f_dirac_op_init_ptr)->pc_type=ZMOB_PC_SYM2;
  
}
 
Fzmobius::~Fzmobius(){
  Zmobus* dwfarg=(Zmobus*)f_dirac_op_init_ptr;
//  if(dwfarg->zmobius_kappa_b) delete [] dwfarg->zmobius_kappa_b;
//  if(dwfarg->zmobius_kappa_c) delete [] dwfarg->zmobius_kappa_c;
//  if(dwfarg->zmobius_kappa_ratio) delete [] dwfarg->zmobius_kappa_ratio;
  
}

FclassType Fzmobius::Fclass(void) const {
  return F_CLASS_ZMOBIUS;
}

int Fzmobius::FmatInv(Vector *f_out, Vector *f_in, 
		     CgArg *cg_arg, 
		     Float *true_res,
		     CnvFrmType cnv_frm,
		     PreserveType prs_f_in, int dminus)
{
  int iter;
  char *fname = "FmatInv(CgArg*,V*,V*,F*,CnvFrmType)";
  VRB.Func(cname,fname);
  Vector *temp;

  int local_ls = GJP.SnodeSites();
  const int s_node_coor = GJP.SnodeCoor();
  const int ls_stride = 24 * GJP.VolNodeSites()/2;
    
  int size = GJP.VolNodeSites() * local_ls * 2 * Colors() * SpinComponents();
  //if(prs_f_in==PRESERVE_YES){ 
    temp = (Vector*) smalloc(cname,fname, "temp", size * sizeof(Float));
    moveFloat((IFloat*)temp,(IFloat*)f_in, size);
    //  }
  
  Vector *dminus_in = (Vector *)
    smalloc(cname, fname, "dminus_in", size* sizeof(Float) );

    
  DiracOpZMobius dop(*this, f_out, f_in, cg_arg, cnv_frm);

  // First multiply D_- to source
  if(dminus==1) dop. Dminus(dminus_in, f_in);
  else if(dminus==-1){ dop. Dminus(f_out, f_in); return 0;}
  else moveFloat((IFloat*)dminus_in, (IFloat*)f_in, size);

#if 1  
  // Multiply 2*kappa
  // do even / odd 
  for(int ieo=0;ieo<2;++ieo){
    for(int s=0; s<local_ls;++s){
      int glb_s = s + local_ls*s_node_coor;
      const Complex kappa_b =
	1.0 / ( 2. * (GJP.ZMobius_b()[glb_s]
		     *(4. - GJP.DwfHeight()) + GJP.DwfA5Inv()) );
 	VRB.Flow(cname,fname,"s=%d Zmobius_b=%e %e kappa_b=%e %e\n",
	glb_s,GJP.ZMobius_b()[glb_s].real(),GJP.ZMobius_b()[glb_s].imag(),kappa_b.real(),kappa_b.imag());
      int idx = s*ls_stride/2;// "/2" is for complex
      vecTimesEquComplex((Complex*)dminus_in+idx+ieo*size/4,
			 2.0*kappa_b, ls_stride);
    }
  }
#endif
  
  // solve it
  //    Because we use dminus_in keeping f_in, we don't have to preserve src
  iter = dop.MatInv(f_out, dminus_in, true_res, PRESERVE_NO);

  // TIZB check
  //Float norm;
  //norm = f_in->NormSqGlbSum(size);
  //if(!UniqueID()) printf("f_mobius Norm in %.14e\n",norm);



  
//  if(prs_f_in==PRESERVE_YES)
	//moveFloat((IFloat*)f_in,(IFloat*)temp, size);
    

  //norm = f_out->NormSqGlbSum(size);
  //if(!UniqueID()) printf("f_mobius  Norm out %.14e\n",norm);

  //dop.Mat(temp,f_out);  
  //  norm = temp->NormSqGlbSum(size);
  //if(!UniqueID()) printf("f_mobius  Norm Mat*out %.14e\n",norm);
  

  sfree(cname, fname,  "dminus_in",  dminus_in);
  // commented out above, so do here too. TB if(prs_f_in==PRESERVE_YES) 
    sfree(cname, fname,  "temp",  temp);

  // Return the number of iterations
  return iter;
}


// FmatInvMobius: same as FmatInv, except that we use mobius DWF
// formalism to speed up the CG inversion (via constructing initial guess).
// n_restart: How many restarts we perform
int Fzmobius::FmatInv(Vector *f_out,
		      Vector *f_in,
		      MobiusArg *mob_l,
		      MobiusArg *mob_s,
		      Float *true_res,
		      CnvFrmType cnv_frm,
		      PreserveType prs_f_in)

{
  const char *fname = "FmatInv(V*, V*, mdwfArg, mdwfArg, F, Cnv, Preserv, int, F)";
  VRB.Func(fname,cname);
  // n_restart and stop_rsd is taken from mob_l
  int n_restart = mob_l->rsd_vec.rsd_vec_len;
  Float *rsd_vec = mob_l->rsd_vec.rsd_vec_val;

  // this implementation doesn't allow splitting in s direction(yet).
  if(GJP.Snodes() != 1){
    ERR.NotImplemented(cname, fname);
  }
  if(cnv_frm != CNV_FRM_YES){
    ERR.NotImplemented(cname, fname);
  }
  if(n_restart < 2){
    ERR.General(cname, fname, "Value %d is invalid for n_restart.\n", n_restart);
  }

  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fmobius
  // call to mobius_init but is done here again in case the user
  // has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the Fmobius object.
  //----------------------------------------------------------------


  
  int ls = GJP.SnodeSites();
  int mob_l_ls = mob_l->ls;
  int mob_s_ls = mob_s->ls;
  if(mob_l_ls*2 != mob_l->zmobius_b_coeff.zmobius_b_coeff_len
     || mob_l_ls*2 != mob_l->zmobius_c_coeff.zmobius_c_coeff_len)
    ERR.General(cname,fname,"wrong ls for large mobius %d %d %d\n",
		mob_l_ls, mob_l->zmobius_b_coeff.zmobius_b_coeff_len,
		mob_l->zmobius_c_coeff.zmobius_c_coeff_len);
  if(mob_s_ls*2 != mob_s->zmobius_b_coeff.zmobius_b_coeff_len
     || mob_s_ls*2 != mob_s->zmobius_c_coeff.zmobius_c_coeff_len)
    ERR.General(cname,fname,"wrong ls for small mobius %d %d %d\n",
		mob_s_ls, mob_s->zmobius_b_coeff.zmobius_b_coeff_len,
		mob_s->zmobius_c_coeff.zmobius_c_coeff_len);
  Float* mobius_b_l = mob_l->zmobius_b_coeff.zmobius_b_coeff_val;
  Float* mobius_c_l = mob_l->zmobius_c_coeff.zmobius_c_coeff_val;
  Float* mobius_b_s = mob_s->zmobius_b_coeff.zmobius_b_coeff_val;
  Float* mobius_c_s = mob_s->zmobius_c_coeff.zmobius_c_coeff_val;
  int size_4d = GJP.VolNodeSites() * 2 * Colors() * SpinComponents();
  //int size_4d = GJP.VolNodeSites() * 24 ;
  int mob_l_size_5d = size_4d * mob_l_ls;
  int mob_s_size_5d = size_4d * mob_s_ls;

  // input times D_-
  Vector *dminus_in = (Vector *) smalloc(cname, fname, "dminus_in", sizeof(Float) * mob_l_size_5d);
  Vector *tmp_mob_l_5d = (Vector *) smalloc(cname, fname, "tmp_mob_l_5d", sizeof(Float) * mob_l_size_5d);
  Vector *tmp2_mob_l_5d = (Vector *) smalloc(cname, fname, "tmp2_mob_l_5d", sizeof(Float) * mob_l_size_5d);
  Vector *tmp_mob_s_5d = (Vector *) smalloc(cname, fname, "tmp_mob_s_5d", sizeof(Float) * mob_s_size_5d);


  CgArg _cg_arg_l = mob_l->cg;
  CgArg _cg_arg_s = mob_s->cg;

  CgArg *cg_arg_l = &_cg_arg_l;
  CgArg *cg_arg_s = &_cg_arg_s;

  
  // the first time, we solve in Mobius to some degree of accuracy
  {
    Float rsd = cg_arg_l->stop_rsd;
    cg_arg_l->stop_rsd = rsd_vec[0];
    GJP.SnodeSites(mob_l_ls);
    GJP.ZMobius_b(mobius_b_l, mob_l_ls);
    GJP.ZMobius_c(mobius_c_l, mob_l_ls);

    //printf("TIZB first solve %e\n", GJP.ZMobius_b()[0]);
    
    DiracOpZMobius dop(*this, f_out, f_in, cg_arg_l, cnv_frm);
    dop.Dminus(dminus_in,f_in);

#if 1  
  // Multiply 2*kappa
    {

    int local_ls = GJP.SnodeSites();
    const int s_node_coor = GJP.SnodeCoor();
    const int ls_stride = 24 * GJP.VolNodeSites()/2;
    int size = GJP.VolNodeSites() * local_ls * 2 * Colors() * SpinComponents();
    // do even / odd   
    for(int ieo=0;ieo<2;++ieo){
    for(int s=0; s<local_ls;++s){
      int glb_s = s + local_ls*s_node_coor;
      const Complex kappa_b =
	1.0 / ( 2. * (GJP.ZMobius_b()[glb_s]
		     *(4. - GJP.DwfHeight()) + GJP.DwfA5Inv()) );
      int idx = s*ls_stride/2;// "/2" is for complex
      vecTimesEquComplex((Complex*)dminus_in+idx+ieo*size/4,
			 2.0*kappa_b, ls_stride);
    }
  }
  //moveFloat((IFloat*)f_in,(IFloat*)dminus_in, size);
    }
#endif


    int iter = dop.MatInv(f_out, dminus_in, PRESERVE_YES);
#if 0
      for(int i=0;i<GJP.VolNodeSites()*GJP.SnodeSites()*24;i++)
	printf("IN OUT1 %d %e %e\n",i,*((Float*)dminus_in+i),*((Float*)f_out+i) );
#endif
    cg_arg_l->stop_rsd = rsd;
    this->Convert(CANONICAL,dminus_in);
  }

  int restart_cnt;
  for(restart_cnt = 1; restart_cnt < n_restart; ++restart_cnt){
    //constructing the new residue
    tmp2_mob_l_5d->CopyVec(dminus_in, mob_l_size_5d);
    {
      //if(!UniqueID())printf("TIZB 1st MAT\n");
      GJP.SnodeSites(mob_l_ls);
      GJP.ZMobius_b(mobius_b_l, mob_l_ls);
      GJP.ZMobius_c(mobius_c_l, mob_l_ls);
      //printf("TIZB  second %e\n", GJP.ZMobius_b()[0]);
    
      DiracOpZMobius dop(*this, tmp_mob_l_5d, f_out, cg_arg_l, cnv_frm);
      dop.Mat(tmp_mob_l_5d, f_out);
    }
    tmp2_mob_l_5d->VecMinusEquVec(tmp_mob_l_5d, mob_l_size_5d);
    //end constructing the new residue(new residue now in tmp2_mob_l_5d)

#if 0
      for(int i=0;i<GJP.VolNodeSites()*GJP.SnodeSites()*24;i++)
	printf("IN OUT2 %d %e\n",i,*((Float*)tmp2_mob_l_5d+i));
#endif

    cg_arg_s->stop_rsd = rsd_vec[restart_cnt];
    tmp_mob_l_5d->VecZero(mob_l_size_5d);

    // Large PV Solve
    {
      //if(!UniqueID())printf("TIZB Large PV\n");
      Float mass = cg_arg_l->mass;
      Float stop_rsd = cg_arg_l->stop_rsd;
      int max_num_iter = cg_arg_l->max_num_iter;

      cg_arg_l->mass = 1.0;
      cg_arg_l->stop_rsd = 1e-5; //rsd_vec[restart_cnt];
      cg_arg_l->max_num_iter = cg_arg_s->max_num_iter;

      GJP.SnodeSites(mob_l_ls);
      GJP.ZMobius_b(mobius_b_l, mob_l_ls);
      GJP.ZMobius_c(mobius_c_l, mob_l_ls);

      DiracOpZMobius dop(*this, tmp_mob_l_5d, tmp2_mob_l_5d, cg_arg_l, cnv_frm);
      int iter = dop.MatInv(tmp_mob_l_5d, tmp2_mob_l_5d, PRESERVE_YES);

      cg_arg_l->mass = mass;
      cg_arg_l->stop_rsd = stop_rsd;
      cg_arg_l->max_num_iter = max_num_iter;
    }
    SpinProject(tmp2_mob_l_5d, tmp_mob_l_5d, mob_l_ls, 1);
    

    // go to the small Mobius lattice
    tmp_mob_s_5d->VecZero(mob_s_size_5d);
    tmp_mob_s_5d->CopyVec(tmp2_mob_l_5d, size_4d);
    
    SpinProject(tmp_mob_l_5d, tmp_mob_s_5d, mob_s_ls, 0);
    

    {
      //if(!UniqueID())printf("TIZB SMALL PV\n");
      Float mass= cg_arg_s->mass; 
      cg_arg_s->mass = 1.0;
      
      GJP.SnodeSites(mob_s_ls);
      GJP.ZMobius_b(mobius_b_s, mob_s_ls);
      GJP.ZMobius_c(mobius_c_s, mob_s_ls);

      DiracOpZMobius dop(*this, tmp_mob_s_5d, tmp_mob_l_5d, cg_arg_s, cnv_frm);
      dop.Mat(tmp_mob_s_5d, tmp_mob_l_5d);
      
      cg_arg_s->mass = mass;
    }
    


    // SMALL SOLVE
    tmp_mob_l_5d->VecZero(mob_s_size_5d);
    {
      //if(!UniqueID())printf("TIZB SMALL SOLVE ls=%d \n", mob_s_ls);
      GJP.SnodeSites(mob_s_ls);
      GJP.ZMobius_b(mobius_b_s,mob_s_ls);
      GJP.ZMobius_c(mobius_c_s,mob_s_ls);
      //      printf("TIZB  fifth Small Solve %e\n", GJP.ZMobius_b()[0]);

      
      DiracOpZMobius dop(*this, tmp_mob_l_5d, tmp_mob_s_5d, cg_arg_s, cnv_frm);
      dop.MatInv(tmp_mob_l_5d, tmp_mob_s_5d, true_res, PRESERVE_YES);
    }
    SpinProject(tmp_mob_s_5d, tmp_mob_l_5d, mob_s_ls, 1);


    //OV2DWF
    tmp_mob_l_5d->VecNegative(tmp_mob_s_5d, size_4d);
    {
      Vector * tmp_2nd_plane = (Vector *)((Float *)tmp_mob_l_5d + size_4d);
      Vector * src_2nd_plane = (Vector *)((Float *)tmp2_mob_l_5d + size_4d);
      tmp_2nd_plane->CopyVec(src_2nd_plane, mob_l_size_5d - size_4d);
    }
    
    SpinProject(tmp2_mob_l_5d, tmp_mob_l_5d, mob_l_ls, 0);
    
    {
      GJP.SnodeSites(mob_l_ls);
      GJP.ZMobius_b(mobius_b_l,mob_l_ls);
      GJP.ZMobius_c(mobius_c_l,mob_l_ls);

      DiracOpZMobius dop(*this, tmp_mob_l_5d, tmp2_mob_l_5d, cg_arg_l, cnv_frm);
      dop.Mat(tmp_mob_l_5d, tmp2_mob_l_5d);
    }
    
    // LARGE PV Solve
    tmp2_mob_l_5d->VecZero(mob_l_size_5d);
    {
      Float mass = cg_arg_l->mass;
      Float stop_rsd = cg_arg_l->stop_rsd;
      int max_num_iter = cg_arg_l->max_num_iter;

      cg_arg_l->mass = 1.0;
      cg_arg_l->stop_rsd = 1e-5; //cg_arg_s->stop_rsd;
      cg_arg_l->max_num_iter = cg_arg_s->max_num_iter;

      GJP.SnodeSites(mob_l_ls);
      GJP.ZMobius_b(mobius_b_l,mob_l_ls);
      GJP.ZMobius_c(mobius_c_l,mob_l_ls);

      DiracOpZMobius dop(*this, tmp2_mob_l_5d, tmp_mob_l_5d, cg_arg_l, cnv_frm);
      dop.MatInv(tmp2_mob_l_5d, tmp_mob_l_5d, PRESERVE_NO);

      cg_arg_l->mass = mass;
      cg_arg_l->stop_rsd = stop_rsd;
      cg_arg_l->max_num_iter = max_num_iter;
    }

    SpinProject(tmp_mob_l_5d, tmp2_mob_l_5d, mob_l_ls, 1);
    
    tmp_mob_l_5d->CopyVec(tmp_mob_s_5d, size_4d);
    
    SpinProject(tmp2_mob_l_5d, tmp_mob_l_5d, mob_l_ls, 0);

    f_out->VecAddEquVec(tmp2_mob_l_5d, mob_l_size_5d);

  }

  // final solve
  int iter;
  {
    GJP.SnodeSites(mob_l_ls);
    GJP.ZMobius_b(mobius_b_l,mob_l_ls);
    GJP.ZMobius_c(mobius_c_l,mob_l_ls);

    DiracOpZMobius dop(*this, f_out, dminus_in, cg_arg_l, cnv_frm);
    iter = dop.MatInv(f_out, dminus_in, true_res, PRESERVE_YES);
  }

  sfree(cname, fname,  "tmp_mob_l_5d",  tmp_mob_l_5d);
  sfree(cname, fname, "tmp2_mob_l_5d", tmp2_mob_l_5d);
  sfree(cname, fname,  "tmp_mob_s_5d",  tmp_mob_s_5d);
  sfree(cname, fname,  "dminus_in",  dminus_in);

  GJP.SnodeSites(ls);

  return iter;
}

//Copied from base class since it uses Dwf explicitly. 

//------------------------------------------------------------------
// int FeigSolv(Vector **f_eigenv, Float *lambda, int *valid_eig,
//              EigArg *eig_arg, 
//              CnvFrmType cnv_frm = CNV_FRM_YES):
// It solve  A * f_eigenv = lambda * f_eigenv where
// A is the fermion matrix (Dirac operator). The solution
// is done with the Ritz algorithm. eig_arg is the 
// structure that contains all the control parameters, f_eigenv
// is the fermion field eigenvectors, lambda are the
// returned eigenvalues.
// f_eigenv is defined on the whole lattice.
// The function returns the total number of Ritz iterations.
//------------------------------------------------------------------
int Fzmobius::FeigSolv(Vector **f_eigenv, Float *lambda,
		      Float *chirality, int *valid_eig,
		      Float **hsum,
		      EigArg *eig_arg, 
		      CnvFrmType cnv_frm)
{

  int iter;
  char *fname = "FeigSolv(EigArg*,V*,F*,CnvFrmType)";
  VRB.Func(cname,fname);
  VRB.Warn(cname,fname,"Not checked\n");
  CgArg cg_arg;
  cg_arg.mass = eig_arg->mass;
  cg_arg.RitzMatOper = eig_arg->RitzMatOper;
  int N_eig = eig_arg->N_eig;
 
  if(cnv_frm == CNV_FRM_YES)
    for(int i=0; i < N_eig; ++i)
      Fconvert(f_eigenv[i], DWF_4D_EOPREC_EE, StrOrd());
 
  // Call constructor and solve for eigenvectors.
  // Use null pointers to fake out constructor.
  Vector *v1 = (Vector *)0;
  Vector *v2 = (Vector *)0;
 
  DiracOpZMobius dop(*this, v1, v2, &cg_arg, CNV_FRM_NO);
  
  iter = dop.RitzEig(f_eigenv, lambda, valid_eig, eig_arg);
 
  if(cnv_frm == CNV_FRM_YES)
    for(int i=0; i < N_eig; ++i)
      Fconvert(f_eigenv[i], CANONICAL, StrOrd());
 
  // rescale eigenvalues to normal convention
  Float factor = 5. - GJP.DwfHeight();
  int i;
  //for(i=0; i<N_eig; ++i)
  //lambda[i] *= factor;

  // calculate chirality
  int Ncb = NumChkb(cg_arg.RitzMatOper);
  size_t f_size = GJP.VolNodeSites()*2*Colors()*SpinComponents()*Ncb/2;
  Vector *four = (Vector *) smalloc (cname,fname, "four", f_size * sizeof(Float));
  Vector *fourg5 = (Vector *) smalloc (cname,fname, "fourg5", f_size * sizeof(Float));
  Float help;

  for (i=0; i<N_eig; i++) {
    Ffive2four (four, f_eigenv[i], 0, GJP.Snodes()*GJP.SnodeSites()-1,Ncb);

    // normalize four
    factor=four->NormSqNode(f_size);
    glb_sum(&factor);
    factor=1./sqrt(factor);
    four->VecTimesEquFloat(factor,f_size);
    Gamma5(fourg5,four,GJP.VolNodeSites()*Ncb/2);
    chirality[i]= four->ReDotProductNode(fourg5, f_size);
    glb_sum(&chirality[i]);
  }

  // calculate hsum
  // slice sum the eigenvector density to make a 1D vector
  if (eig_arg->print_hsum) {

    if(cnv_frm == CNV_FRM_NO)
      for(int i=0; i < N_eig; ++i)
	Fconvert(f_eigenv[i], CANONICAL, StrOrd());

    Float *f_in = (Float *) smalloc (cname, fname, "f_in", 
		GJP.VolNodeSites()*GJP.SnodeSites()*sizeof(Float));
    
    for(i=0; i < N_eig; ++i) {
      IFloat *fp= (IFloat *) f_eigenv[i];
      int j;
      for (j=0; j<GJP.VolNodeSites()*GJP.SnodeSites(); j++, fp+= 24) 
	f_in[j]= Float (dotProduct (fp,fp,24));

      if (i==0) {
	for (j=0; j<GJP.VolNodeSites()*GJP.SnodeSites(); j++)
	  printf ("%f ", f_in[j]);
      }
      printf ("\n");

      f_eigenv[i]->SliceArraySumFive (hsum[i], f_in, eig_arg->hsum_dir);
    }

    sfree(cname, fname, "f_in", f_in);
  }
  sfree(cname,fname, "four",four);
  sfree(cname,fname, "fourg5",fourg5);
  // Return the number of iterations
  return iter;
}

//------------------------------------------------------------------
// Solve  A * f_eigenv = lambda * f_eigenv where
// A is the fermion matrix (Dirac operator). The solution
// is done with the Lanczos algorithm. eig_arg is the 
// structure that contains all the control parameters, f_eigenv
// is the fermion field eigenvectors, lambda are the
// returned eigenvalues.
// f_eigenv is defined on the whole lattice.
//------------------------------------------------------------------
int Fzmobius::FeigSolv(Vector **f_eigenv, Float *lambda,
		      LanczosArg *eig_arg, 
		      CnvFrmType cnv_frm)
{
  int iter;
  char *fname = "FeigSolv(V*,F*,LanczosArg*,CnvFrmType)";
  VRB.Func(cname,fname);
  CgArg cg_arg;
  cg_arg.mass = eig_arg->mass;
  cg_arg.RitzMatOper = eig_arg->RitzMat_convcheck;
  cg_arg.eigen_shift = eig_arg->eigen_shift;
  // printf("Shift=%e\n",cg_arg.eigen_shift );exit(2);
  //int N_eig = eig_arg->N_eig;
  int nk = eig_arg->nk_lanczos_vectors;
  int np = eig_arg->np_lanczos_vectors;
  int maxiters = eig_arg->maxiters;
  Float stopres = eig_arg->stop_residual;
  MatrixPolynomialArg* cheby_arg = &(eig_arg->matpoly_arg);
  
  if(cnv_frm == CNV_FRM_YES) // convert only nk, not (nk+np)
    for(int i=0; i < nk; ++i) 
      Fconvert(f_eigenv[i], DWF_4D_EOPREC_EE, StrOrd());
 
  // Call constructor and solve for eigenvectors.
  // Use null pointers to fake out constructor.
  Vector *v1 = (Vector *)0;
  Vector *v2 = (Vector *)0;
  DiracOpZMobius dop(*this, v1, v2, &cg_arg, CNV_FRM_NO);
  
  iter = dop.ImpResLanczos(f_eigenv, lambda,  eig_arg); 

  if(cnv_frm == CNV_FRM_YES)
    for(int i=0; i < nk; ++i) // convert only nk, not (nk+np)
      Fconvert(f_eigenv[i], CANONICAL, StrOrd());
  
  return iter;
}


void Fzmobius::Fdslash(Vector *f_out, Vector *f_in, CgArg *cg_arg,
                 CnvFrmType cnv_frm, int dir_flag){
  DiracOpZMobius dop(*this, f_out, f_in, cg_arg, cnv_frm);
  unsigned long offset = GJP.VolNodeSites()*FsiteSize()/ (2*6);
//  dop.Dslash(f_out,f_in,CHKB_ODD,DAG_NO);
  dop.Dslash(f_out,f_in,CHKB_EVEN,DAG_NO);

#if 0 
  // Multiply 2*kappa
  // do even / odd 
{
  int local_ls = GJP.SnodeSites();
  const int s_node_coor = GJP.SnodeCoor();
  const int ls_stride = 24 * GJP.VolNodeSites()/2;
  const Complex half=-0.5;
  int size = GJP.VolNodeSites() * local_ls * 2 * Colors() * SpinComponents();
  for(int ieo=0;ieo<2;++ieo){
    for(int s=0; s<local_ls;++s){
      int glb_s = s + local_ls*s_node_coor;
//      const Complex kappa_b =
//	2.0 / ( 2 * (GJP.ZMobius_b()[glb_s]
//		     *(4 - GJP.DwfHeight()) + GJP.DwfA5Inv()) );
      const Complex kappa_b =
	 1 * (GJP.ZMobius_b()[glb_s]
		     *(4 - GJP.DwfHeight()) + GJP.DwfA5Inv()) ;
      int idx = s*ls_stride/2;// "/2" is for complex
//      vecTimesEquComplex((Complex*)f_out+idx+ieo*size/4,
//			 kappa_b, ls_stride);
//void zTimesV1PluszTimesV2(Complex *a, Complex b, const Complex *c,
//                    Complex d, const Complex *e, int len)
      zTimesV1PluszTimesV2( (Complex*)f_out+idx+ieo*size/4,
			 kappa_b, (Complex*)f_in+idx+ieo*size/4,
			 half, (Complex*)f_out+idx+ieo*size/4, ls_stride);
    }
  }
}
#endif

#if 0
  // Multiply 2*kappa
  // do even / odd 
  for(int ieo=0;ieo<2;++ieo){
    for(int s=0; s<local_ls;++s){
      int glb_s = s + local_ls*s_node_coor;
      const Complex kappa_b =
	1.0 / ( 2 * (GJP.ZMobius_b()[glb_s]
		     *(4 - GJP.DwfHeight()) + GJP.DwfA5Inv()) );
      int idx = s*ls_stride/2;// "/2" is for complex
      vecTimesEquComplex((Complex*)dminus_in+idx+ieo*size/4,
			 2.0*kappa_b, ls_stride);
    }
  }
  //moveFloat((IFloat*)f_in,(IFloat*)dminus_in, size);
#endif

}

    void Fzmobius::MatPc (Vector * f_out, Vector * f_in, Float mass, Float epsilon,
              DagType dag){
	CgArg cg_arg; 
	cg_arg.mass=mass;
	cg_arg.epsilon=epsilon;
	DiracOpZMobius dirac(*this,f_out,f_in,&cg_arg,CNV_FRM_NO);
	if (dag==DAG_NO) dirac.MatPc(f_out,f_in);
	else dirac.MatPcDag(f_out,f_in);
	
    }

CPS_END_NAMESPACE
