#include<config.h>
#ifdef USE_OMP
#include<omp.h>
#endif
CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of FdwfBase class.

  $Id: f_dwf_base.C,v 1.44 2013-06-07 19:26:34 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/f_dwf_base/f_dwf_base.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// f_dwf_base.C
//
// FdwfBase is derived from FwilsonTypes and is relevant to
// domain wall fermions, moved from old Fdwf class
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/qcdio.h>
#include <math.h>
#include <cstdlib>
#include <util/lattice.h>
#include <util/dirac_op.h>
#include <util/dwf.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/vector.h>
#include <util/random.h>
#include <util/error.h>
#include <util/time_cps.h>
#include <util/enum_func.h>
#include <comms/scu.h> // GRF
#include <comms/glb.h>
CPS_START_NAMESPACE

//------------------------------------------------------------------
// Initialize static variables.
//------------------------------------------------------------------

//------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------

// Hack to measure the time for constructor/destructors
static unsigned long called=1;
Float ctor_time=0.;
FdwfBase::FdwfBase()
{
  ctor_time -= dclock();
  cname = "FdwfBase";
  char *fname = "FdwfBase()";
  VRB.Func(cname,fname);

  //----------------------------------------------------------------
  // Check if anisotropy is present and exit since FdwfBase has
  // not been tested for anisotropic lattices.
  //----------------------------------------------------------------
  if(GJP.XiBare() != 1 ||
     GJP.XiV()    != 1 ||
     GJP.XiVXi()  != 1   ){
    ERR.General(cname,fname,
    "XiBare=%g, XiV=%g, XiVXi=%g : FdwfBase has not been tested with anisotropy\n",
                GJP.XiBare(), GJP.XiV(), GJP.XiVXi());
  }

  //----------------------------------------------------------------
  // Do initializations before the dwf library can be used
  //----------------------------------------------------------------
  static Dwf dwf_struct;
  f_dirac_op_init_ptr = &dwf_struct;
  dwf_init((Dwf *) f_dirac_op_init_ptr);
  ctor_time += dclock();
}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
FdwfBase::~FdwfBase()
{
  ctor_time -= dclock();
  char *fname = "~FdwfBase()";
  VRB.Func(cname,fname);

  //----------------------------------------------------------------
  // Un-initialize the dwf library. Memory is set free here.
  //----------------------------------------------------------------
  dwf_end((Dwf *) f_dirac_op_init_ptr);
  ctor_time += dclock();
  if(called%50==0){
    print_time(cname,"(Ctor+Dtor)*50)",ctor_time);
    ctor_time=0.;
  }
  called++;
}


//------------------------------------------------------------------
// FclassType Fclass(void):
// It returns the type of fermion class.
//------------------------------------------------------------------
FclassType FdwfBase::Fclass(void) const {
  return F_CLASS_DWF;
}


#if 0
//------------------------------------------------------------------
// int ExactFlavors() : 
// Returns the number of exact flavors of the matrix that
// is inverted during a molecular dynamics evolution.
//------------------------------------------------------------------
int FdwfBase::ExactFlavors(void) const
{
  return 2;
}


//------------------------------------------------------------------
// int SpinComponents() : 
// Returns the number of spin components.
//------------------------------------------------------------------
int FdwfBase::SpinComponents(void) const
{
  return 4;
}
#endif


//------------------------------------------------------------------
// int FsiteSize() : 
// Returns the number of fermion field components 
// (including real/imaginary) on a site of the 4-D lattice.
//------------------------------------------------------------------
int FdwfBase::FsiteSize(void) const
{
  return 2 * Colors() * SpinComponents() * GJP.SnodeSites();  
  // re/im * colors * spin_components * Ls
}

//------------------------------------------------------------------
// int FchkbEvl() :
// returns 1 => The fermion fields in the evolution
//      or the CG that inverts the evolution matrix
//      are defined on a single checkerboard (half the 
//      lattice).
//------------------------------------------------------------------
int FdwfBase::FchkbEvl(void) const
{
  return 1;
}


//------------------------------------------------------------------
// int FmatEvlInv(Vector *f_out, Vector *f_in, 
//                CgArg *cg_arg, 
//                Float *true_res,
//		  CnvFrmType cnv_frm = CNV_FRM_YES):
// It calculates f_out where A * f_out = f_in and
// A is the preconditioned fermion matrix that appears
// in the HMC evolution (even/odd  preconditioning 
// of [Dirac^dag Dirac]. The inversion is done
// with the conjugate gradient. cg_arg is the structure
// that contains all the control parameters, f_in is the
// fermion field source vector, f_out should be set to be
// the initial guess and on return is the solution.
// f_in and f_out are defined on a checkerboard.
// If true_res !=0 the value of the true residual is returned
// in true_res.
// *true_res = |src - MatPcDagMatPc * sol| / |src|
// The function returns the total number of CG iterations.
//------------------------------------------------------------------
int FdwfBase::FmatEvlInv(Vector *f_out, Vector *f_in, 
		     CgArg *cg_arg, 
		     Float *true_res,
		     CnvFrmType cnv_frm)
{
  int iter;
  char *fname = "FmatEvlInv(CgArg*,V*,V*,F*,CnvFrmType)";
  VRB.Func(cname,fname);

  DiracOpDwf dwf(*this, f_out, f_in, cg_arg, cnv_frm);
  
  iter = dwf.InvCg(&(cg_arg->true_rsd));
  if (true_res) *true_res = cg_arg ->true_rsd;

  // Return the number of iterations
  return iter;
}


//------------------------------------------------------------------
// Overloaded function is same as original but with true_res=0;
//------------------------------------------------------------------
int FdwfBase::FmatEvlInv(Vector *f_out, Vector *f_in, 
		     CgArg *cg_arg, 
		     CnvFrmType cnv_frm)
{ return FmatEvlInv(f_out, f_in, cg_arg, 0, cnv_frm); }


//------------------------------------------------------------------
// int FmatEvlMInv(Vector *f_out, Vector *f_in, 
//                Float shift[], int Nshift, 
//                CgArg **cg_arg, Float *true_res,
//		  CnvFrmType cnv_frm = CNV_FRM_YES):
// It calculates f_out where (A + shift)* f_out = f_in and
// A is the fermion matrix that appears in the HMC 
// evolution ([Dirac^dag Dirac]) and shift is a real shift of the 
// fermion matrix, with Nshift such shifts. The inversion is done 
// with the multishift conjugate gradient. cg_arg is the structure
// that contains all the control parameters, f_in is the
// fermion field source vector, f_out is the array of solution 
// vectors, f_in and f_out are defined on a checkerboard.
// The function returns the total number of CG iterations.
//------------------------------------------------------------------
int FdwfBase::FmatEvlMInv(Vector **f_out, Vector *f_in, Float *shift, 
			  int Nshift, int isz, CgArg **cg_arg, 
			  CnvFrmType cnv_frm, MultiShiftSolveType type, 
			  Float *alpha, Vector **f_out_d)
{
  int iter;
  char *fname = "FmatEvlMInv(V**, V*, .....)";
  VRB.Func(cname,fname);

  int f_size = GJP.VolNodeSites() * FsiteSize() / (FchkbEvl()+1);
  Float dot = f_in -> NormSqGlbSum(f_size);
  VRB.Result(cname,fname,"f_size=%d\n",f_size);

  Float *RsdCG = new Float[Nshift];
  for (int s=0; s<Nshift; s++) RsdCG[s] = cg_arg[s]->stop_rsd;

  Float * out_p = (Float *)f_in;
  VRB.Result(cname,fname,"f_in=%g\n",out_p[0]);
  for(int i =0;i<Nshift;i++){
	Float a_tmp=1.; 
	if (alpha) a_tmp=alpha[i];
  VRB.Result(cname,fname,"%d: shift=%g alpha=%g RsdCG=%g\n",i,shift[i],a_tmp,RsdCG[i]);
  }
  //Fake the constructor
  DiracOpDwf dwf(*this, f_out[0], f_in, cg_arg[0], cnv_frm);

  int return_value= dwf.MInvCG(f_out,f_in,dot,shift,Nshift,isz,RsdCG,type,alpha);  
  out_p = (Float *)f_out[0];
  VRB.Result(cname,fname,"f_out[0]=%g\n",out_p[0]);
  for (int s=0; s<Nshift; s++) cg_arg[s]->true_rsd = RsdCG[s];  

  delete[] RsdCG;
  return return_value;
}

//------------------------------------------------------------------
// Lattice class api to the chronological inverter
//------------------------------------------------------------------
void FdwfBase::FminResExt(Vector *sol, Vector *source, Vector **sol_old, 
			 Vector **vm, int degree, CgArg *cg_arg, CnvFrmType cnv_frm)
{

  char *fname = "FminResExt(V*, V*, V**, V**, int, CgArg *, CnvFrmType)";
  VRB.Func(cname,fname);
  
  DiracOpDwf dwf(*this, sol, source, cg_arg, cnv_frm);
  dwf.MinResExt(sol,source,sol_old,vm,degree);
  
}

struct Results{
  int iter;
  Float time;
};

// ======================================================================
// set_tuning: set MdwfArg to current parameter to be tuned.
//
// collect_results: write tuning state to file, collect tuning result,
// and set up for the next tuning(we only know what to do right after
// a test is finished).
// ======================================================================
static void set_tuning(const CgArg *cg_arg)
{
  const char *cname = "";
  const char *fname = "set_tuning()";

  const MdwfTuning *tuning = GJP.GetMdwfTuning();

  switch(tuning->stage){
  case 0: // test finished
    if(tuning->opti_index < 0){
      GJP.SetMdwfArg(NULL);
    }else{
      MdwfArg *opti_arg = &(tuning->mdwf_arg.mdwf_arg_val[tuning->opti_index]);

      opti_arg->cg_arg = *cg_arg;
      GJP.SetMdwfArg(opti_arg);
    }
    return;
  case 1: // zero time
    GJP.SetMdwfArg(NULL);
    return;
  case 2: // c5
  case 3: // rsd
  case 4: // restart count
    {
      MdwfArg *arg = &(tuning->mdwf_arg.mdwf_arg_val[tuning->index]);
      GJP.SetMdwfArg(arg);

      // modify rsd_vec_len, since we can't change the corresponding value in MdwfTuning.
      arg = GJP.GetMdwfArg();
      arg->cg_arg = *cg_arg;
      arg->rsd_vec.rsd_vec_len = tuning->rc_val;
      return;
    }
  default:
    ERR.General(cname, fname, "Invalid value for tuning->stage (current value is %d).\n", tuning->stage);
    return;
  }
}

// write out one line of the test result.
// if passing a null pointer to the first pointer, write out a header instead.
static void record_results(const struct Results *results)
{
  FILE *fp = Fopen(GJP.GetMdwfTuningRecordFN(), "a+");

  MdwfTuning *tuning = GJP.GetMdwfTuning();

  if(tuning->stage == 1 || tuning->stage == 0 && tuning->opti_index < 0){
    Fprintf(fp, "%5d %17.10e %17.10e %17.10e %5d %17.10e %10d\n", 0, 0.0, 0.0, 0.0, 0, results->time, results->iter);
  }else if(tuning->stage == 0){
    Fprintf(fp, "%5d %17.10e %17.10e %17.10e %5d %17.10e %10d\n",
            tuning->mdwf_arg.mdwf_arg_val[tuning->opti_index].b5.b5_len,
            tuning->mdwf_arg.mdwf_arg_val[tuning->opti_index].b5.b5_val[0],
            tuning->mdwf_arg.mdwf_arg_val[tuning->opti_index].c5.c5_val[0],
            tuning->mdwf_arg.mdwf_arg_val[tuning->opti_index].rsd_vec.rsd_vec_val[1],
            tuning->mdwf_arg.mdwf_arg_val[tuning->opti_index].rsd_vec.rsd_vec_len,
            results->time,
            results->iter);
  }else{
    Fprintf(fp, "%5d %17.10e %17.10e %17.10e %5d %17.10e %10d\n",
            tuning->mdwf_arg.mdwf_arg_val[tuning->index].b5.b5_len,
            tuning->mdwf_arg.mdwf_arg_val[tuning->index].b5.b5_val[0],
            tuning->mdwf_arg.mdwf_arg_val[tuning->index].c5.c5_val[0],
            tuning->mdwf_arg.mdwf_arg_val[tuning->index].rsd_vec.rsd_vec_val[1],
            tuning->rc_val,
            results->time,
            results->iter);
  }
  Fclose(fp);
}

static void collect_results(const struct Results *results)
{
  const char *cname = "";
  const char *fname = "collect_results()";

  MdwfTuning *tuning = GJP.GetMdwfTuning();

  // Don't modify these 2 numbers.
  const Float magic1 = 0.5*(3.0-sqrt(5.0)); //0.382
  const Float magic2 = 1.0 - magic1; //0.618

  switch(tuning->stage){
  case 0: // test finished
    break;
  case 1: { // zero time
    tuning->opti_time = results->time;
    tuning->opti_index = -1;

    tuning->index = 0;
    tuning->stage = 2;

    Float c5_mid = GJP.SnodeSites()/(2.0*tuning->ls_min) - 0.5;
    Float c5_0 = c5_mid - tuning->c5_range;
    Float c5_3 = c5_mid + tuning->c5_range;

    tuning->c5.c5_val[0].val = c5_0;
    tuning->c5.c5_val[1].val = c5_0*magic2 + c5_3*magic1;
    tuning->c5.c5_val[2].val = c5_0*magic1 + c5_3*magic2;
    tuning->c5.c5_val[3].val = c5_3;

    tuning->c5.c5_val[0].dwf_cg = -1;
    tuning->c5.c5_val[1].dwf_cg = -1;
    tuning->c5.c5_val[2].dwf_cg = -1;
    tuning->c5.c5_val[3].dwf_cg = -1;

    tuning->rc_val = 2;
    tuning->mdwf_arg.mdwf_arg_val[0].rsd_vec.rsd_vec_val[0] = 1.0e-2;
    tuning->mdwf_arg.mdwf_arg_val[0].rsd_vec.rsd_vec_val[1] = 1.0e-8;

    int i;
    for(i = 0; i < tuning->ls_min; ++i){
      tuning->mdwf_arg.mdwf_arg_val[0].b5.b5_val[i] = c5_0 + 1.0;
      tuning->mdwf_arg.mdwf_arg_val[0].c5.c5_val[i] = c5_0;
    }
    break;
  }
  case 2: { // c5
    int i = 0;
    while(i<3 && tuning->c5.c5_val[i].dwf_cg != -1) ++i;
    tuning->c5.c5_val[i].dwf_cg = results->iter;

    do ++i; while(i<4 && tuning->c5.c5_val[i].dwf_cg != -1);

    MdwfArg *arg= &(tuning->mdwf_arg.mdwf_arg_val[tuning->index]);
    Float c5;
    if(i == 4) { // we have done all tests and ready to shrink the search range.
      C5State *c5_state = tuning->c5.c5_val;
      if(c5_state[1].dwf_cg <= c5_state[2].dwf_cg && c5_state[2].dwf_cg < c5_state[3].dwf_cg){ // cut [2-3]
        c5_state[3] = c5_state[2];
        c5_state[2] = c5_state[1];
        c5_state[1].val = c5_state[0].val*magic2 + c5_state[3].val*magic1;
        c5_state[1].dwf_cg = -1;

        c5 = c5_state[1].val;
      }else if(c5_state[1].dwf_cg >= c5_state[2].dwf_cg && c5_state[0].dwf_cg > c5_state[1].dwf_cg){ // cut [0-1]
        c5_state[0] = c5_state[1];
        c5_state[1] = c5_state[2];

        c5_state[2].val = c5_state[0].val*magic1 + c5_state[3].val*magic2;
        c5_state[2].dwf_cg = -1;

        c5 = c5_state[2].val;
      }else{ // stop, prepare for next stage
        tuning->stage = 3;

        tuning->rsd_val = 1.0e-8;
        tuning->rsd_time = 1.0e30;

        tuning->rc_val = 2;
        arg->rsd_vec.rsd_vec_val[0] = 1.0e-2;
        arg->rsd_vec.rsd_vec_val[1] = tuning->rsd_val;

        c5 = (c5_state[0].val + c5_state[3].val)*0.5;
      }
    }else{
      c5 = tuning->c5.c5_val[i].val;
    }

    // set c5 values for the MdwfArg copy in MdwfTuning.
    int j;
    int ls = arg->b5.b5_len;
    for(j = 0; j < ls; ++j){
      arg->b5.b5_val[j] = c5 + 1.0;
      arg->c5.c5_val[j] = c5;
    }
    break;
  }
  case 3: { // rsd
    MdwfArg *arg= &(tuning->mdwf_arg.mdwf_arg_val[tuning->index]);
    Float cur_time = results->time;
    if(cur_time > tuning->rsd_time || tuning->rsd_val >= 0.3){ // stop, prepare for next stage
      tuning->stage = 4;

      arg->rsd_vec.rsd_vec_val[1] = tuning->rsd_val/(tuning->rsd_granularity);

      tuning->rc_val = 3;
      arg->rsd_vec.rsd_vec_val[2] = arg->rsd_vec.rsd_vec_val[1];
      tuning->rc_time = tuning->rsd_time;
    }else{ // test a larger stopping condition
      tuning->rsd_val *= tuning->rsd_granularity;
      tuning->rsd_time = cur_time;

      arg->rsd_vec.rsd_vec_val[1] = tuning->rsd_val;
    }
    break;
  }
  case 4: { // restart count
    MdwfArg *arg= &(tuning->mdwf_arg.mdwf_arg_val[tuning->index]);
    Float cur_time = results->time;

    if(cur_time <= tuning->rc_time){ // test a larger restart count
      if(tuning->rc_val < tuning->rc_max){
        tuning->rc_time = cur_time;

        arg->rsd_vec.rsd_vec_val[(tuning->rc_val)++] = arg->rsd_vec.rsd_vec_val[1];
        break;
      }else{
        ++tuning->rc_val;
        tuning->rc_time = cur_time;
        VRB.Warn(cname, fname, "rc_max reached, increase rc_max may yield a better tuning result.\n");
      }
    } // stop, prepare for next ls value or finalize tuning
   
    arg->rsd_vec.rsd_vec_len = --tuning->rc_val;

    // test if this min time is also global min time
    if(tuning->rc_time < tuning->opti_time) {
      tuning->opti_index = tuning->index;
      tuning->opti_time = tuning->rc_time;
    }

    // prepare for next ls
    ++tuning->index;
    int next_ls = tuning->index*2 + tuning->ls_min;
    if(next_ls > tuning->ls_max){ // we have exhausted the search, finalize tuning
      tuning->stage = 0;
    }else{
      tuning->stage = 2;
      int index = tuning->index;

      Float c5_mid = GJP.SnodeSites()/(2.0*next_ls) - 0.5;
      Float c5_0 = c5_mid - tuning->c5_range;
      Float c5_3 = c5_mid + tuning->c5_range;

      tuning->c5.c5_val[0].val = c5_0;
      tuning->c5.c5_val[1].val = c5_0*magic2 + c5_3*magic1;
      tuning->c5.c5_val[2].val = c5_0*magic1 + c5_3*magic2;
      tuning->c5.c5_val[3].val = c5_3;

      tuning->c5.c5_val[0].dwf_cg = -1;
      tuning->c5.c5_val[1].dwf_cg = -1;
      tuning->c5.c5_val[2].dwf_cg = -1;
      tuning->c5.c5_val[3].dwf_cg = -1;

      tuning->rc_val = 2;
      tuning->mdwf_arg.mdwf_arg_val[index].rsd_vec.rsd_vec_val[0] = 1.0e-2;
      tuning->mdwf_arg.mdwf_arg_val[index].rsd_vec.rsd_vec_val[1] = 1.0e-8;

      int i;
      int ls = index*2 + tuning->ls_min;
      for(i = 0; i < ls; ++i) {
        tuning->mdwf_arg.mdwf_arg_val[index].b5.b5_val[i] = c5_0 + 1.0;
        tuning->mdwf_arg.mdwf_arg_val[index].c5.c5_val[i] = c5_0;
      }
    }
    break;
  }
  default:
    ERR.General(cname, fname, "Invalid value for MdwfTuning::stage (current value is %d).\n", tuning->stage);
    break;
  }
  tuning->Encode(GJP.GetMdwfTuningFN(), "mdwf_tuning");
}

//------------------------------------------------------------------
// int FmatInv(Vector *f_out, Vector *f_in, 
//             CgArg *cg_arg, 
//             Float *true_res,
//             CnvFrmType cnv_frm = CNV_FRM_YES,
//             PreserveType prs_f_in = PRESERVE_YES):
// It calculates f_out where A * f_out = f_in and
// A is the fermion matrix (Dirac operator). The inversion
// is done with the conjugate gradient. cg_arg is the 
// structure that contains all the control parameters, f_in 
// is the fermion field source vector, f_out should be set 
// to be the initial guess and on return is the solution.
// f_in and f_out are defined on the whole lattice.
// If true_res !=0 the value of the true residual is returned
// in true_res.
// *true_res = |src - MatPcDagMatPc * sol| / |src|
// cnv_frm is used to specify if f_in should be converted 
// from canonical to fermion order and f_out from fermion 
// to canonical. 
// prs_f_in is used to specify if the source
// f_in should be preserved or not. If not the memory usage
// is less by half the size of a fermion vector.
// The function returns the total number of CG iterations.
//------------------------------------------------------------------
int FdwfBase::FmatInv(Vector *f_out, Vector *f_in, 
		  CgArg *cg_arg, 
		  Float *true_res,
		  CnvFrmType cnv_frm,
		  PreserveType prs_f_in)
{
  const char *fname = "FmatInv(CgArg*,V*,V*,F*,CnvFrmType)";
  VRB.Func(cname, fname);

  if(GJP.GetMdwfTuning() != NULL){
    struct Results results;
    set_tuning(cg_arg);

    results.time = -dclock();
    results.iter = 
      FmatInvMobius(f_out, f_in, cg_arg, GJP.GetMdwfArg(), true_res, cnv_frm, prs_f_in);
    results.time += dclock();

    record_results(&results);
    collect_results(&results);

    return results.iter;

  }else if(GJP.GetMdwfArg() != NULL){
    return FmatInvMobius(f_out, f_in, cg_arg, GJP.GetMdwfArg(), true_res, cnv_frm, prs_f_in);

  }else {
    DiracOpDwf dwf(*this, f_out, f_in, cg_arg, cnv_frm);
    return dwf.MatInv(true_res, prs_f_in);
  }
}

// FmatInvMobius: same as FmatInv, except that we use mobius DWF
// to speed up the CG inversion (via constructing initial guess).
int FdwfBase::FmatInvMobius(Vector *f_out,
                            Vector *f_in,
                            CgArg *cg_arg_dwf,
                            MdwfArg *mdwf_arg,
                            Float *true_res,
                            CnvFrmType cnv_frm,
                            PreserveType prs_f_in)
{
  const char *fname = "FmatInvMobius(V*, V*, ...)";

#ifdef USE_MDWF
  // this implementation doesn't allow splitting in s direction(yet).
  if(GJP.Snodes() != 1) ERR.NotImplemented(cname, fname);
  if(cnv_frm != CNV_FRM_YES) ERR.NotImplemented(cname, fname);

  if(mdwf_arg == NULL) {
    DiracOpDwf dwf(*this, f_out, f_in, cg_arg_dwf, cnv_frm);
    return dwf.MatInv(f_out, f_in, true_res, prs_f_in);
  }
  if(mdwf_arg->rsd_vec.rsd_vec_len < 2){
    ERR.General(cname, fname, "Invalid restart count for MdwfArg.\n");
  }

  if(mdwf_arg->use_mdwf_for_dwf) {
    //======================================================================
    // Construct a Mobius-simulated DWF for use with MDWF library (to
    // supersede the slow CPS DWF solve).
    MdwfArg mobius_dwf;
    {
      mobius_dwf.M5 = GJP.DwfHeight();
      mobius_dwf.use_single_precision = mdwf_arg->use_single_precision;
      //      mobius_dwf.cg_arg = *cg_arg_dwf;
      
      const int dwf_ls = GJP.SnodeSites();
      mobius_dwf.b5.b5_len = dwf_ls;
      mobius_dwf.c5.c5_len = dwf_ls;
      // we need not assign values for rsd_vec, since DiracOpMdwf doesn't interpret it.
      mobius_dwf.b5.b5_val = (Float *)smalloc(cname, fname, "b5_val", sizeof(Float)*dwf_ls);
      mobius_dwf.c5.c5_val = (Float *)smalloc(cname, fname, "c5_val", sizeof(Float)*dwf_ls);
      
      for(int i=0;i<dwf_ls; ++i){
        mobius_dwf.b5.b5_val[i] = 1.0;
        mobius_dwf.c5.c5_val[i] = 0.0;
      }
    }
    //======================================================================
    const Float *rsd_vec = mdwf_arg->rsd_vec.rsd_vec_val;
    
    const int mob_ls = mdwf_arg->b5.b5_len;
    const int dwf_size_5d = GJP.VolNodeSites() * FsiteSize();
    const int mob_size_5d = GJP.VolNodeSites() * mob_ls * SPINOR_SIZE;
    const int size_4d = GJP.VolNodeSites() * SPINOR_SIZE;
    
    Vector *tmp_dwf_5d = (Vector *)smalloc(cname, fname, "tmp_dwf_5d", sizeof(Float) * dwf_size_5d);
    Vector *tmp2_dwf_5d = (Vector *)smalloc(cname, fname, "tmp2_dwf_5d", sizeof(Float) * dwf_size_5d);
    Vector *tmp_mob_5d = (Vector *)smalloc(cname, fname, "tmp_mob_5d", sizeof(Float) * mob_size_5d);
    
    // the first time, we solve in DWF to some degree of accuracy
    {
      mobius_dwf.cg_arg = *cg_arg_dwf;
      mobius_dwf.cg_arg.stop_rsd = rsd_vec[0];
      
      DiracOpMdwf mdwf(*this, &mobius_dwf);
      mdwf.MatInv(f_out, f_in, NULL, PRESERVE_YES);
    }
  
    const int n_restart = mdwf_arg->rsd_vec.rsd_vec_len;
    for(int restart_cnt = 1; restart_cnt < n_restart; ++restart_cnt){
      //constructing the new residue, note: it's not a good idea to use single precision here.
      tmp2_dwf_5d->CopyVec(f_in, dwf_size_5d);
      {
        mobius_dwf.cg_arg = *cg_arg_dwf;
        mobius_dwf.use_single_precision = 0;
        
        DiracOpMdwf mdwf(*this, &mobius_dwf);
        mdwf.Mat(tmp_dwf_5d, f_out);
        mobius_dwf.use_single_precision = mdwf_arg->use_single_precision;
      }
      tmp2_dwf_5d->VecMinusEquVec(tmp_dwf_5d, dwf_size_5d);
      //end constructing the new residue
      
      mdwf_arg->cg_arg.stop_rsd = rsd_vec[restart_cnt];
      
      tmp_dwf_5d->VecZero(dwf_size_5d);
      {
        mobius_dwf.cg_arg.mass = 1.0;
        mobius_dwf.cg_arg.stop_rsd = rsd_vec[restart_cnt];
        mobius_dwf.cg_arg.max_num_iter = mdwf_arg->cg_arg.max_num_iter;
        
        DiracOpMdwf mdwf(*this, &mobius_dwf);
        mdwf.MatInv(tmp_dwf_5d, tmp2_dwf_5d, NULL, PRESERVE_YES);
      }
      SpinProject(tmp2_dwf_5d, tmp_dwf_5d, GJP.SnodeSites(), 1);
      
      tmp_mob_5d->VecZero(mob_size_5d);
      tmp_mob_5d->CopyVec(tmp2_dwf_5d, size_4d);
      
      SpinProject(tmp_dwf_5d, tmp_mob_5d, mob_ls, 0);
      
      {
        Float mass= mdwf_arg->cg_arg.mass; 
        mdwf_arg->cg_arg.mass = 1.0;
        
        DiracOpMdwf mdwf(*this, mdwf_arg);
        mdwf.Mat(tmp_mob_5d, tmp_dwf_5d);
        
        mdwf_arg->cg_arg.mass = mass;
      }
      
      tmp_dwf_5d->VecZero(mob_size_5d);
      
      {
        Float time = -dclock();
        
        DiracOpMdwf mdwf(*this, mdwf_arg);
        mdwf.MatInv(tmp_dwf_5d, tmp_mob_5d, NULL, PRESERVE_YES);
        
        time += dclock();
        print_flops("DiracOpMdwf","MatInv", 0, time);
      }
      
      SpinProject(tmp_mob_5d, tmp_dwf_5d, mob_ls, 1);
      
      //OV2DWF
      tmp_dwf_5d->VecNegative(tmp_mob_5d, size_4d);
      {
        Vector *tmp_2nd_plane = (Vector *)((Float *)tmp_dwf_5d + size_4d);
        Vector *src_2nd_plane = (Vector *)((Float *)tmp2_dwf_5d + size_4d);
        tmp_2nd_plane->CopyVec(src_2nd_plane, dwf_size_5d - size_4d);
      }
      
      SpinProject(tmp2_dwf_5d, tmp_dwf_5d, GJP.SnodeSites(), 0);
      
      {
        mobius_dwf.cg_arg = *cg_arg_dwf;
        DiracOpMdwf mdwf(*this, &mobius_dwf);
        mdwf.Mat(tmp_dwf_5d, tmp2_dwf_5d);
      }
      
      tmp2_dwf_5d->VecZero(dwf_size_5d);
      {
        mobius_dwf.cg_arg.mass = 1.0;
        mobius_dwf.cg_arg.stop_rsd = mdwf_arg->cg_arg.stop_rsd;
        mobius_dwf.cg_arg.max_num_iter = mdwf_arg->cg_arg.max_num_iter;
        
        DiracOpMdwf mdwf(*this, &mobius_dwf);
        mdwf.MatInv(tmp2_dwf_5d, tmp_dwf_5d, NULL, PRESERVE_NO);
      }
      SpinProject(tmp_dwf_5d, tmp2_dwf_5d, GJP.SnodeSites(), 1);
      
      tmp_dwf_5d->CopyVec(tmp_mob_5d, size_4d);
      
      SpinProject(tmp2_dwf_5d, tmp_dwf_5d, GJP.SnodeSites(), 0);
      
      f_out->VecAddEquVec(tmp2_dwf_5d, dwf_size_5d);
    }
    
    sfree(cname, fname,  "tmp_dwf_5d",  tmp_dwf_5d);
    sfree(cname, fname, "tmp2_dwf_5d", tmp2_dwf_5d);
    sfree(cname, fname,  "tmp_mob_5d",  tmp_mob_5d);
    
    int iter;
    // ======================================================================
    // fix up using Mobius, simulating DWF
    {
      mobius_dwf.cg_arg = *cg_arg_dwf;
      mobius_dwf.cg_arg.stop_rsd /= 5.0 - GJP.DwfHeight();
      mobius_dwf.use_single_precision = 0;
      
      DiracOpMdwf mdwf(*this, &mobius_dwf);
      iter = mdwf.MatInv(f_out, f_in, NULL, PRESERVE_YES);
      // renormalize, the MDWF package has a different normalization convention
      f_out->VecTimesEquFloat(5.0-GJP.DwfHeight(), dwf_size_5d);
    }
    // ======================================================================
    sfree(mobius_dwf.b5.b5_val);
    sfree(mobius_dwf.c5.c5_val);
    
    // final solve, using CPS
    DiracOpDwf dwf(*this, f_out, f_in, cg_arg_dwf, cnv_frm);
    dwf.MatInv(f_out, f_in, true_res, prs_f_in);
    return iter;
  }else{
    const Float *rsd_vec = mdwf_arg->rsd_vec.rsd_vec_val;
    
    const int mob_ls = mdwf_arg->b5.b5_len;
    const int dwf_size_5d = GJP.VolNodeSites() * FsiteSize();
    const int mob_size_5d = GJP.VolNodeSites() * mob_ls * SPINOR_SIZE;
    const int size_4d = GJP.VolNodeSites() * SPINOR_SIZE;
    
    Vector *tmp_dwf_5d = (Vector *) smalloc(cname, fname, "tmp_dwf_5d", sizeof(Float) * dwf_size_5d);
    Vector *tmp2_dwf_5d = (Vector *) smalloc(cname, fname, "tmp2_dwf_5d", sizeof(Float) * dwf_size_5d);
    Vector *tmp_mob_5d = (Vector *) smalloc(cname, fname, "tmp_mob_5d", sizeof(Float) * mob_size_5d);
    
    // the first time, we solve in DWF to some degree of accuracy
    {
      Float rsd = cg_arg_dwf->stop_rsd;
      cg_arg_dwf->stop_rsd = rsd_vec[0];
      
      DiracOpDwf dwf(*this, f_out, f_in, cg_arg_dwf, CNV_FRM_YES);
      dwf.MatInv(f_out, f_in, NULL, PRESERVE_YES);
      cg_arg_dwf->stop_rsd = rsd;
    }
    
    const int n_restart = mdwf_arg->rsd_vec.rsd_vec_len;
    for(int restart_cnt = 1; restart_cnt < n_restart; ++restart_cnt){
      //constructing the new residue
      tmp2_dwf_5d->CopyVec(f_in, dwf_size_5d);
      {
        DiracOpDwf dwf(*this, tmp_dwf_5d, f_out, cg_arg_dwf, CNV_FRM_YES);
        dwf.Mat(tmp_dwf_5d, f_out);
      }
      tmp2_dwf_5d->VecMinusEquVec(tmp_dwf_5d, dwf_size_5d);
      //end constructing the new residue
      
      mdwf_arg->cg_arg.stop_rsd = rsd_vec[restart_cnt];
      
      tmp_dwf_5d->VecZero(dwf_size_5d);
      {
        CgArg cg_arg;
        cg_arg.mass = 1.0;
        
        cg_arg.stop_rsd = rsd_vec[restart_cnt];
        cg_arg.max_num_iter = mdwf_arg->cg_arg.max_num_iter;
        
        DiracOpDwf dwf(*this, tmp_dwf_5d, tmp2_dwf_5d, &cg_arg, CNV_FRM_YES);
        dwf.MatInv(tmp_dwf_5d, tmp2_dwf_5d, NULL, PRESERVE_YES);
      }
      SpinProject(tmp2_dwf_5d, tmp_dwf_5d, GJP.SnodeSites(), 1);
      
      tmp_mob_5d->VecZero(mob_size_5d);
      tmp_mob_5d->CopyVec(tmp2_dwf_5d, size_4d);
    
      SpinProject(tmp_dwf_5d, tmp_mob_5d, mob_ls, 0);
      
      {
        Float mass = mdwf_arg->cg_arg.mass; 
        mdwf_arg->cg_arg.mass = 1.0;
        
        DiracOpMdwf mdwf(*this, mdwf_arg);
        mdwf.Mat(tmp_mob_5d, tmp_dwf_5d);
        
        mdwf_arg->cg_arg.mass = mass;
      }
      
      tmp_dwf_5d->VecZero(mob_size_5d);
      
      {
        Float time = -dclock();
        
        DiracOpMdwf mdwf(*this, mdwf_arg);
        mdwf.MatInv(tmp_dwf_5d, tmp_mob_5d, NULL, PRESERVE_YES);
        
        time += dclock();
        print_flops("DiracOpMdwf","MatInv", 0, time);
      }
      
      SpinProject(tmp_mob_5d, tmp_dwf_5d, mob_ls, 1);
    
      //OV2DWF
      tmp_dwf_5d->VecNegative(tmp_mob_5d, size_4d);
      {
        Vector *tmp_2nd_plane = (Vector *)((Float *)tmp_dwf_5d + size_4d);
        Vector *src_2nd_plane = (Vector *)((Float *)tmp2_dwf_5d + size_4d);
        tmp_2nd_plane->CopyVec(src_2nd_plane, dwf_size_5d - size_4d);
      }
      
      SpinProject(tmp2_dwf_5d, tmp_dwf_5d, GJP.SnodeSites(), 0);
      
      {
        DiracOpDwf dwf(*this, tmp_dwf_5d, tmp2_dwf_5d, cg_arg_dwf, CNV_FRM_YES);
        dwf.Mat(tmp_dwf_5d, tmp2_dwf_5d);
      }
      
      tmp2_dwf_5d->VecZero(dwf_size_5d);
      {
        CgArg cg_arg_tmp;
        cg_arg_tmp.mass = 1.0;
        cg_arg_tmp.stop_rsd = mdwf_arg->cg_arg.stop_rsd;
        cg_arg_tmp.max_num_iter = mdwf_arg->cg_arg.max_num_iter;
        DiracOpDwf dwf(*this, tmp2_dwf_5d, tmp_dwf_5d, &cg_arg_tmp, CNV_FRM_YES);
        dwf.MatInv(tmp2_dwf_5d, tmp_dwf_5d, NULL, PRESERVE_NO);
      }
      SpinProject(tmp_dwf_5d, tmp2_dwf_5d, GJP.SnodeSites(), 1);
      
      tmp_dwf_5d->CopyVec(tmp_mob_5d, size_4d);
      
      SpinProject(tmp2_dwf_5d, tmp_dwf_5d, GJP.SnodeSites(), 0);
      
      f_out->VecAddEquVec(tmp2_dwf_5d, dwf_size_5d);
    }
    
    sfree(cname, fname,  "tmp_dwf_5d",  tmp_dwf_5d);
    sfree(cname, fname, "tmp2_dwf_5d", tmp2_dwf_5d);
    sfree(cname, fname,  "tmp_mob_5d",  tmp_mob_5d);
    
    // final solve
    DiracOpDwf dwf(*this, f_out, f_in, cg_arg_dwf, cnv_frm);
    return dwf.MatInv(f_out, f_in, true_res, prs_f_in);
  }
#else
  ERR.NotImplemented(cname, fname);
  return 1;
#endif
}

int FdwfBase::eig_FmatInv(Vector **V, const int vec_len, Float *M, const int nev, const int m, float **U, Rcomplex *invH, const int def_len, const Float *restart, const int restart_len,
		Vector *f_out, Vector *f_in, 
		  CgArg *cg_arg, 
		  Float *true_res,
		  CnvFrmType cnv_frm,
		  PreserveType prs_f_in)
{
  int iter;
  char *fname = "eig_FmatInv(CgArg*,V*,V*,F*,CnvFrmType)";
  VRB.Func(cname,fname);

  DiracOpDwf dwf(*this, f_out, f_in, cg_arg, cnv_frm);
    
  iter = dwf.eig_MatInv(V,vec_len, M, nev, m, U, invH, def_len, restart,restart_len, f_out, f_in, true_res, prs_f_in);

  // Return the number of iterations
  return iter;
}

//------------------------------------------------------------------
// Overloaded function is same as original but with true_res=0;
//------------------------------------------------------------------
int FdwfBase::FmatInv(Vector *f_out, Vector *f_in, 
		  CgArg *cg_arg, 
		  CnvFrmType cnv_frm,
		  PreserveType prs_f_in)
{ return FmatInv(f_out, f_in, cg_arg, 0, cnv_frm, prs_f_in); }


//------------------------------------------------------------------
/*!
  \param five The 5-dimensional field.
  \param four The 4-dimensional field.
  \param s_u The global 5th direction (s) coordinate where the
   upper two components (right chirality) of the 5-dim. field
   take the values of those of the 4-dim. field.
  \param s_l The global 5th direction (s) coordinate where the
   lower two components (left chirality) of the 5-dim. field
   take the values of those of the 4-dim. field.
   \post The 5-dim field is zero everywhere except where the global
   5th direction coordinate (s) is \a s_l or \a s_u, where it takes the values
   explained above.
*/
//------------------------------------------------------------------
void FdwfBase::Ffour2five(Vector *five, Vector *four, int s_u, int s_l, int Ncb)
{
  int x;
  int i;
  Float *field_4D;
  Float *field_5D;
  char *fname = "Ffour2five(V*,V*,i,i)";
  VRB.Func(cname,fname);


//------------------------------------------------------------------
// Initializations
//------------------------------------------------------------------
  int f_size = GJP.VolNodeSites() * FsiteSize()*Ncb/2;
  int ls = GJP.SnodeSites();
  int vol_4d = GJP.VolNodeSites()*Ncb/2;
  int ls_stride = 24 * vol_4d;
  int s_u_local = s_u % GJP.SnodeSites();
  int s_l_local = s_l % GJP.SnodeSites();
  int s_u_node = s_u / GJP.SnodeSites();
  int s_l_node = s_l / GJP.SnodeSites();


//------------------------------------------------------------------
// Set *five using the 4D field *four. 
//------------------------------------------------------------------

  // Set all components of the 5D field to zero.
  //---------------------------------------------------------------
  field_5D  = (Float *) five;
  for(i=0; i<f_size; i++){
    field_5D[i]  = 0.0;
  }

  // Do the two upper spin components if s_u is in the node
  //---------------------------------------------------------------
  if( s_u_node == GJP.SnodeCoor() ){
    field_4D  = (Float *) four;
    field_5D  = (Float *) five;
    field_5D  = field_5D  + s_u_local * ls_stride;
    for(x=0; x<vol_4d; x++){
      for(i=0; i<12; i++){
	field_5D[i]  = field_4D[i];
      }
      field_4D  = field_4D  + 24;
      field_5D  = field_5D  + 24;
    }
  }

  // Do the two lower spin components if s_l is in the node
  //----------------------------------------------------------------
  if( s_l_node == GJP.SnodeCoor() ){
    field_4D  = (Float *) four;
    field_5D  = (Float *) five;
    field_4D  = field_4D  + 12;
    field_5D  = field_5D  + 12 + s_l_local * ls_stride;
    for(x=0; x<vol_4d; x++){
      for(i=0; i<12; i++){
	field_5D[i]  = field_4D[i];
      }
      field_4D  = field_4D  + 24;
      field_5D  = field_5D  + 24;
    }
  }

}


//------------------------------------------------------------------
/*!
  \param four The 4-dimensional field.
  \param five The 5-dimensional field.
  \param s_u The global 5th direction (s) coordinate where 
   the values of the upper two components (right chirality) of the 5-dim. field
   are taken by those of the 4-dim. field.
  \param s_l The global 5th direction coordinate (s) where the values of 
   the lower two components (left chirality) of the 5-dim. field
   are taken by  those of the 4-dim. field.
   \post The 5-dim field is zero everywhere except where the global
   5th direction coordinate is \a s_l or \a s_u, where it takes the values
   explained above.
   \post An identical 4-dim. field is reproduced on all nodes in the s
   direction.
*/
//------------------------------------------------------------------
void FdwfBase::Ffive2four(Vector *four, Vector *five, int s_u, int s_l, int Ncb)
{
  int x;
  int i;
  Float *field_4D;
  Float *field_5D;
  char *fname = "Ffive2four(V*,V*,i,i)";
  VRB.Func(cname,fname);


//------------------------------------------------------------------
// Initializations
//------------------------------------------------------------------
  int ls = GJP.SnodeSites();
  int f_size = GJP.VolNodeSites() * FsiteSize()*Ncb / (ls*2);
  int vol_4d = GJP.VolNodeSites()*Ncb/2;
  int ls_stride = 24 * vol_4d;
  int s_u_local = s_u % GJP.SnodeSites();
  int s_l_local = s_l % GJP.SnodeSites();
  int s_u_node = s_u / GJP.SnodeSites();
  int s_l_node = s_l / GJP.SnodeSites();


//------------------------------------------------------------------
// Set *four using the 5D field *five. 
//------------------------------------------------------------------

  // Set all components of the 4D field to zero.
  //---------------------------------------------------------------
  field_4D  = (Float *) four;
  for(i=0; i<f_size; i++){
    field_4D[i]  = 0.0;
  }

  // Do the two upper spin components if s_u is in the node
  //---------------------------------------------------------------
  if( s_u_node == GJP.SnodeCoor() ){
    field_4D = (Float *) four;
    field_5D = (Float *) five;
    field_5D = field_5D + s_u_local * ls_stride;
    for(x=0; x<vol_4d; x++){
      for(i=0; i<12; i++){
	field_4D[i] = field_5D[i];
      }
      field_4D = field_4D + 24;
      field_5D = field_5D + 24;
    }
  }

  // Do the two lower spin components if s_l is in the node
  //----------------------------------------------------------------
  if( s_l_node == GJP.SnodeCoor() ){
    field_4D = (Float *) four;
    field_5D = (Float *) five;
    field_4D = field_4D + 12;
    field_5D = field_5D + 12 + s_l_local * ls_stride;
    for(x=0; x<vol_4d; x++){
      for(i=0; i<12; i++){
	field_4D[i] = field_5D[i];
      }
      field_4D = field_4D + 24;
      field_5D = field_5D + 24;
    }
  }

  // Sum along s direction to get the same 4D field in all 
  // s node slices.
  //----------------------------------------------------------------
  if( GJP.Snodes() > 1) {
    Float sum;
    field_4D  = (Float *) four;
    for(i=0; i<f_size; i++){
      sum = field_4D[i];
      glb_sum_dir(&sum, 4);
      field_4D[i] = sum;    
    }
  }

}

void FdwfBase::Fsolfour2five(Vector *sol_5d, Vector *sol_4d, Vector *src_5d, CgArg *cg_arg)
{
  const char *fname="Fsolfour2five()";

  const int f_size_4d = GJP.VolNodeSites() * SPINOR_SIZE;
  const int f_size_5d = GJP.VolNodeSites() * FsiteSize();
  const int ls_glb = GJP.Snodes() * GJP.SnodeSites();
  const int s_size = GJP.SnodeSites();
  VRB.Result(cname,fname,"sol_5d=%p sol_4d=%p src_5d=%p cg_arg=%p",sol_5d,sol_4d,src_5d,cg_arg);

  Vector *tmp_5d = (Vector *)smalloc(cname, fname, "tmp_5d", sizeof(Float) * f_size_5d);
  
  tmp_5d->CopyVec(src_5d, f_size_5d);
  
  // construct P^{-1} D(1)^{-1}*src
  sol_5d->VecZero(f_size_5d);
  {
    Float mass = cg_arg->mass;
    cg_arg->mass = 1.;
    
    DiracOpDwf dwf(*this, sol_5d, tmp_5d, cg_arg, CNV_FRM_YES);
    dwf.MatInv(sol_5d, tmp_5d, NULL, PRESERVE_YES);
    cg_arg->mass = mass;
  }
  SpinProject(tmp_5d, sol_5d, s_size, 1);

  // set c_1 = -y_1
  if(GJP.SnodeCoor() == 0)
    tmp_5d->VecNegative(sol_4d, f_size_4d);

  SpinProject(sol_5d, tmp_5d, s_size, 0);
  {
    DiracOpDwf dwf(*this, tmp_5d, sol_5d, cg_arg, CNV_FRM_YES);
    dwf.Mat(tmp_5d, sol_5d);
  }
  {
    Float mass = cg_arg->mass;
    cg_arg->mass = 1.;
    
    DiracOpDwf dwf(*this, sol_5d, tmp_5d, cg_arg, CNV_FRM_YES);
    dwf.MatInv(sol_5d, tmp_5d, NULL, PRESERVE_YES);
    cg_arg->mass = mass;
  }
  SpinProject(tmp_5d, sol_5d, s_size, 1);

  if(GJP.SnodeCoor() == 0)
    tmp_5d->CopyVec(sol_4d, f_size_4d);

  SpinProject(sol_5d, tmp_5d, s_size, 0);

  {
    DiracOpDwf dwf(*this, sol_5d, src_5d, cg_arg, CNV_FRM_YES);
    dwf.MatInv(sol_5d, src_5d, NULL, PRESERVE_YES);
  }

  sfree(cname, fname, "tmp_5d", tmp_5d);
}

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
int FdwfBase::FeigSolv(Vector **f_eigenv, Float *lambda,
		      Float *chirality, int *valid_eig,
		      Float **hsum,
		      EigArg *eig_arg, 
		      CnvFrmType cnv_frm)
{
  int iter;
  char *fname = "FeigSolv(EigArg*,V*,F*,CnvFrmType)";
  VRB.Func(cname,fname);
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
 
  DiracOpDwf dwf(*this, v1, v2, &cg_arg, CNV_FRM_NO);
  
  iter = dwf.RitzEig(f_eigenv, lambda, valid_eig, eig_arg);
 
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
  int f_size = GJP.VolNodeSites()*2*Colors()*SpinComponents()*Ncb/2;
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
int FdwfBase::FeigSolv(Vector **f_eigenv, Float *lambda,
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
  DiracOpDwf dwf(*this, v1, v2, &cg_arg, CNV_FRM_NO);
  
  iter = dwf.ImpResLanczos(f_eigenv, lambda,  eig_arg); 

  if(cnv_frm == CNV_FRM_YES)
    for(int i=0; i < nk; ++i) // convert only nk, not (nk+np)
      Fconvert(f_eigenv[i], CANONICAL, StrOrd());
  
  return iter;
}


//------------------------------------------------------------------
// SetPhi(Vector *phi, Vector *frm1, Vector *frm2, Float mass):
// It sets the pseudofermion field phi from frm1, frm2.
// Note that frm2 is not used.
//------------------------------------------------------------------
Float FdwfBase::SetPhi(Vector *phi, Vector *frm1, Vector *frm2,
		       Float mass, DagType dag){
  char *fname = "SetPhi(V*,V*,V*,F)";
  VRB.Func(cname,fname);
  CgArg cg_arg;
  cg_arg.mass = mass;

  if (phi == 0)
    ERR.Pointer(cname,fname,"phi") ;

  if (frm1 == 0)
    ERR.Pointer(cname,fname,"frm1") ;
  
  DiracOpDwf dwf(*this, frm1, 0, &cg_arg, CNV_FRM_NO) ;
#if 0
  { IFloat *tmp = (IFloat *)frm1;
  printf("frm1[0]=%e\n",*tmp);}
#endif
  if (dag == DAG_YES) dwf.MatPcDag(phi, frm1) ;
  else dwf.MatPc(phi, frm1) ;
#if 0
{ IFloat *tmp = (IFloat *)phi;
  printf("phi[0]=%e\n",*tmp);}
#endif

  return FhamiltonNode(frm1, frm1);
}

#define PROFILE

ForceArg FdwfBase::RHMC_EvolveMomFforce(Matrix *mom, Vector **sol, int degree,
					int isz, Float *alpha, Float mass, 
					Float dt, Vector **sol_d, 
					ForceMeasure force_measure) {
  char *fname = "RHMC_EvolveMomFforce";
  char *force_label=NULL;

  Float L1 = 0.0;
  Float L2 = 0.0;
  Float Linf = 0.0;

  int g_size = GJP.VolNodeSites() * GsiteSize();

  Matrix *mom_tmp;

  if (force_measure == FORCE_MEASURE_YES) {
    mom_tmp = (Matrix*)smalloc(g_size*sizeof(Float),cname, fname, "mom_tmp");
    ((Vector*)mom_tmp) -> VecZero(g_size);
    force_label = new char[100];
  } else {
    mom_tmp = mom;
  }

#ifdef PROFILE
  Float time = -dclock();
  ForceFlops=0;
#endif

  for (int i=0; i<degree; i++){
    ForceArg Fdt = FdwfBase::EvolveMomFforce(mom_tmp,sol[i],mass,alpha[i]*dt);
    Float *out_p = (Float *) sol[i];
    VRB.Result(cname,fname,"sol[%d]=%g\n",i,out_p[0]);
    out_p = (Float *) mom_tmp;
    VRB.Result(cname,fname,"mom_tmp[3]=%g\n",out_p[3]);

    if (force_measure == FORCE_MEASURE_YES) {
      sprintf(force_label, "Rational, mass = %e, pole = %d:", mass, i+isz);
      Fdt.print(dt, force_label);
    }
  }

#ifdef PROFILE
  time += dclock();
  print_flops(cname,fname,ForceFlops,time);
#endif

  // If measuring the force, need to measure and then sum to mom
  if (force_measure == FORCE_MEASURE_YES) {
    for (int i=0; i<g_size/18; i++) {
      Float norm = (mom_tmp+i)->norm();
      Float tmp = sqrt(norm);
      L1 += tmp;
      L2 += norm;
      Linf = (tmp>Linf ? tmp : Linf);
    }
    glb_sum(&L1);
    glb_sum(&L2);
    glb_max(&Linf);

    L1 /= 4.0*GJP.VolSites();
    L2 /= 4.0*GJP.VolSites();

    fTimesV1PlusV2((IFloat*)mom, 1.0, (IFloat*)mom_tmp, (IFloat*)mom, g_size);

    delete[] force_label;
    sfree(mom_tmp, cname, fname, "mom_tmp");
  }

  return ForceArg(L1, sqrt(L2), Linf);
}

//------------------------------------------------------------------
// Float FhamiltonNode(Vector *phi, Vector *chi):
// The fermion Hamiltonian of the node sublattice
// chi must be the solution of Cg with source phi.	       
//------------------------------------------------------------------
Float FdwfBase::FhamiltonNode(Vector *phi, Vector *chi){
  char *fname = "FhamiltonNode(V*,V*)";
  VRB.Func(cname,fname);

  if (phi == 0)
    ERR.Pointer(cname,fname,"phi") ;

  if (chi == 0)
    ERR.Pointer(cname,fname,"chi") ;

  int f_size = GJP.VolNodeSites() * FsiteSize() / 2 ;

  Float ret_val;
  ret_val = phi->ReDotProductNode(chi, f_size ) ;

  // Sum accross s nodes in case Snodes() != 1
  glb_sum_dir(&ret_val, 4) ;

  return ret_val ;

}


//------------------------------------------------------------------
// Float BhamiltonNode(Vector *boson, Float mass):
// The boson Hamiltonian of the node sublattice.
//------------------------------------------------------------------
Float FdwfBase::BhamiltonNode(Vector *boson, Float mass){
  char *fname = "BhamiltonNode(V*,F)";
  VRB.Func(cname,fname);
  CgArg cg_arg;
  cg_arg.mass = mass;

  if (boson == 0)
    ERR.Pointer(cname,fname,"boson");

  int f_size = GJP.VolNodeSites() * FsiteSize() / 2 ;

  char *str_tmp = "bsn_tmp" ;
  Vector *bsn_tmp = (Vector *)
    smalloc(cname,fname,str_tmp,f_size*sizeof(Float));

  DiracOpDwf dwf(*this, boson, bsn_tmp, &cg_arg, CNV_FRM_NO) ;

  dwf.MatPc(bsn_tmp,boson);

  Float ret_val = bsn_tmp->NormSqNode(f_size) ;

  sfree(cname,fname,str_tmp,bsn_tmp);

  // Sum accross s nodes in case Snodes() != 1
  glb_sum_dir(&ret_val, 4) ;

  return ret_val ;
}


//------------------------------------------------------------------
// int FsiteOffsetChkb(const int *x):
// Sets the offsets for the fermion fields on a 
// checkerboard. The fermion field storage order
// is not the canonical one but it is particular
// to the fermion type. x[i] is the 
// ith coordinate where i = {0,1,2,3} = {x,y,z,t}.
//------------------------------------------------------------------
int FdwfBase::FsiteOffsetChkb(const int *x) const {
// ???
//  ERR.NotImplemented(cname, "FsiteOffsetChkb");
  int index = x[4];
  int vol = GJP.NodeSites(4);
  int parity = (x[4]+x[3]+x[2]+x[1]+x[0]+1)%2; //Odd first
  for(int i = 3; i>=0;i--){
    index = index*GJP.NodeSites(i)+x[i];
    vol *= GJP.NodeSites(i);
  }
  index = (index + vol*parity)/2;
  if (0){
	printf("FsiteOffsetChkb:(%d %d %d %d %d) %d\n",
	x[0],x[1],x[2],x[3],x[4],index);
  }

  return index; 
}


//--------------------------------------------------------------------
// void SpinProject():
//
// Does a spin projection along s direction, specifically (the
// matrices appear below apply in s direction, the example are given
// for Ls=5):
// 
// for type == 0:
//
//       [ P- P+ 0  0  0  ]
//       [ 0  P- P+ 0  0  ]
// out = [ 0  0  P- P+ 0  ] in
//       [ 0  0  0  P- P+ ]
//       [ P+ 0  0  0  P- ]
//
// for type == 1:
//
//       [ P- 0  0  0  P+ ]
//       [ P+ P- 0  0  0  ]
// out = [ 0  P+ P- 0  0  ] in
//       [ 0  0  P+ P- 0  ]
//       [ 0  0  0  P+ P- ]
//
// in and out are assumed to be in CANONICAL storage order.
//--------------------------------------------------------------------
void FdwfBase::SpinProject(Vector * out, Vector *in, int s_size, int type)
{
  const char * fname = "SpinProject(V*, V*, int)";

  // this implementation doesn't allow splitting in s direction(yet).
  if(GJP.Snodes() != 1){
    ERR.NotImplemented(cname, fname);
  }

  int lcl_size[5];

  lcl_size[0] = GJP.XnodeSites();
  lcl_size[1] = GJP.YnodeSites();
  lcl_size[2] = GJP.ZnodeSites();
  lcl_size[3] = GJP.TnodeSites();
  lcl_size[4] = s_size;

#define for_each_site_4d(x)                             \
  for(x[0] = 0; x[0] < lcl_size[0]; ++x[0])             \
    for(x[1] = 0; x[1] < lcl_size[1]; ++x[1])           \
      for(x[2] = 0; x[2] < lcl_size[2]; ++x[2])         \
        for(x[3] = 0; x[3] < lcl_size[3]; ++x[3])       

  int lcl_pos[4], color;
  int s;

  int f_size_4d= GJP.VolNodeSites()*2*Colors()*SpinComponents();

  for(s = 0; s < lcl_size[4]; ++s){
    for_each_site_4d(lcl_pos){
      int offset = s*f_size_4d;

      int offset_next;
      if(type == 0){
        offset_next = s == lcl_size[4]-1 ? 0: offset+f_size_4d;
      }else{
        offset_next = s == 0 ? (lcl_size[4]-1)*f_size_4d: offset-f_size_4d;
      }
      
      offset += SPINOR_SIZE * (lcl_pos[0] + lcl_size[0]*(lcl_pos[1] + lcl_size[1]*(lcl_pos[2] + lcl_size[2]*lcl_pos[3])));
      offset_next += SPINOR_SIZE * (lcl_pos[0] + lcl_size[0]*(lcl_pos[1] + lcl_size[1]*(lcl_pos[2] + lcl_size[2]*lcl_pos[3])));
      
      Float * out_site = (Float *)out + offset;
      Float * in_site = (Float *)in + offset;
      Float * in_site_next = (Float *)in + offset_next;
      
      memcpy((void *)out_site,
             (void *)in_site_next,
             SPINOR_SIZE/2*sizeof(Float));
      
      memcpy((void *)(out_site + SPINOR_SIZE/2),
             (void *)(in_site  + SPINOR_SIZE/2),
             SPINOR_SIZE/2*sizeof(Float));
    }
  }
  return;
}
//#endif

//--------------------------------------------------------------------
// void Freflex (Vector *out, Vector *in)
// does the reflexion in s needed for the hermitian D_dwf operator.
//
//--------------------------------------------------------------------
void FdwfBase::Freflex(Vector *out, Vector *in)
{
  char *fname = "Freflex(V*,V*)";
  VRB.Func(cname,fname);

  int i,n,node,f_size_5d,f_size_4d,half_size_5d,half_size_4d;
  Vector *send_buf, *rcv_buf;
  int numblk, blklen,s,s_reflex;

  f_size_5d= GJP.VolNodeSites()*FsiteSize();
  f_size_4d= GJP.VolNodeSites()*2*Colors()*SpinComponents();

  half_size_5d= f_size_5d/2;
  half_size_4d= f_size_4d/2;

  numblk=1;
  blklen=f_size_5d;
  while (blklen > 1023) {
    numblk*=2;
    blklen/=2;
  }

  VRB.Debug (cname, fname,"%d %d %d %d %d %d\n",
	     f_size_5d, f_size_4d, half_size_5d, half_size_4d, 
	     numblk, blklen);

  //reserve space for send and receive buffers
  send_buf= (Vector *) smalloc (f_size_5d*sizeof(IFloat));
  if (send_buf == 0) ERR.Pointer (cname, fname, "send_buf");
  VRB.Smalloc (cname, fname, "send_buf", send_buf, f_size_5d * sizeof(IFloat));
  rcv_buf= (Vector *) smalloc (f_size_5d*sizeof(IFloat));
  if (rcv_buf == 0) ERR.Pointer (cname, fname, "rcv_buf");
  VRB.Smalloc (cname, fname, "rcv_buf", rcv_buf, f_size_5d * sizeof(IFloat));

  for (n=0; n<GJP.Snodes(); n++) {

    VRB.Debug (cname,fname,"n= %d\n", n);
    if (n==0) {
      VRB.Debug (cname,fname,"copy from in to rcv_buf\n");
      // copy from in to rcv_buf
      for (i=0; i<f_size_5d; i++)
	((IFloat *)rcv_buf)[i]= ((IFloat *)in)[i];
    }
    else {
      VRB.Debug (cname,fname,
		"copy from rcv_buf to send_buf, send send_buf to rcv_buf\n");
      // copy from rcv_buf to send_buf
      // send send_buf to rcv_buf
      for (i=0; i<f_size_5d; i++)
	((IFloat *)send_buf)[i]= ((IFloat *)rcv_buf)[i];

      for (i=0; i<numblk;i++)
	getMinusData (((IFloat *) rcv_buf)+i*blklen, 
		      ((IFloat *) send_buf)+i*blklen, blklen, 4);
    }

    // node is the node where the data in rcv_buf comes from originally
    if ((node=GJP.SnodeCoor()-n)<0)
      node+= GJP.Snodes();

    VRB.Debug (cname,fname,"node= %d\n", node);

    if (node==GJP.Snodes()-1-GJP.SnodeCoor()) {
      VRB.Debug (cname,fname,"we have the right data in rcv_buf\n");
      // we have the right data in rcv_buf
      // now do the reflexion

      // in and out are in dwf checkerboard storage order

      for (s=0;s<GJP.SnodeSites();s++) {
	s_reflex= GJP.SnodeSites()-1-s;

	if ((GJP.SnodeSites()%2)==0) {
	  for (i=0; i<half_size_4d; i++) {
	    ((IFloat *)out)[s_reflex*half_size_4d+i]= 
	      ((IFloat *)rcv_buf)[s*half_size_4d+i+half_size_5d];
	    ((IFloat *)out)[s_reflex*half_size_4d+i+half_size_5d]= 
	      ((IFloat *)rcv_buf)[s*half_size_4d+i];
	  }
	}
	else {
	  for (i=0; i<half_size_4d; i++) {
	    ((IFloat *)out)[s_reflex*half_size_4d+i]= 
	      ((IFloat *)rcv_buf)[s*half_size_4d+i];
	    ((IFloat *)out)[s_reflex*half_size_4d+i+half_size_5d]= 
             ((IFloat *)rcv_buf)[s*half_size_4d+i+half_size_5d];
	  }
	}
      }
    }
  }
  
  VRB.Sfree(cname,fname, "send_buf", send_buf);
  sfree(send_buf);
  VRB.Sfree(cname,fname, "rcv_buf", rcv_buf);
  sfree(rcv_buf);

  VRB.FuncEnd (cname,fname);
}


CPS_END_NAMESPACE
