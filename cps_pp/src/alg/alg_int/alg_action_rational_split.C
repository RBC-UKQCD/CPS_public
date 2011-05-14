#include<config.h>
#include<math.h>
CPS_START_NAMESPACE 
//------------------------------------------------------------------
//
// alg_action_rational_split.C
//
// AlgActionRationalSplit represents a subset of a bilinear action
// with a rational approximation kernel.  The rational function is
// split to allow each partial fraction to evolve at independent rates.
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

AlgActionRationalSplit::AlgActionRationalSplit(AlgActionRational &Rat,
					       ActionRationalSplitArg &r_arg)
					       
					       
  : AlgActionRational()
{
  cname = "AlgActionRationalSplit";
  char *fname = "AlgActionRationalSplit(AlgActionRational *, int **)";
  VRB.Func(cname,fname);

  int_type = INT_RATIONAL_SPLIT;
  rat_split_arg = &r_arg;
  rat = &Rat;
  n_masses = 0;
  fermion = rat->getFermion();

  //!< First check n_masses split = n_masses rational
  if (rat_split_arg->fractionSplit.fractionSplit_len != rat->getNmass())
    ERR.General(cname, fname,
		"Inconsistency between RationalSplitArg and AlgActionRational n_masses\n");

  if (rat->getNmass() > 0) {
    fractionSplit = (int**)smalloc(2*sizeof(int*),cname,fname,"fractionSplit");
    fractionSplit[0] = 
      (int*)smalloc(rat->getNmass()*sizeof(int),cname,fname,"fractionSplit[0]");
    fractionSplit[1] = 
      (int*)smalloc(rat->getNmass()*sizeof(int),cname,fname,"fractionSplit[1]");
    
    for (int i=0; i<rat->getNmass(); i++) {
      fractionSplit[0][i] = 
	rat_split_arg->fractionSplit.fractionSplit_val[i].split_low;
      fractionSplit[1][i] = 
	rat_split_arg->fractionSplit.fractionSplit_val[i].split_high;

      //!< Need to check that splits are valid
      for (int j=fractionSplit[0][i];j<fractionSplit[1][i]; j++) 
	rat->setSplit(i,j);
    }
  }

}

AlgActionRationalSplit::~AlgActionRationalSplit() {

  char *fname = "~AlgActionRationalSplit()";
  VRB.Func(cname,fname);
  
  if (rat->getNmass() > 0) {
    sfree(fractionSplit[0], cname, fname, "fractionSplit[0]");
    sfree(fractionSplit[1], cname, fname, "fractionSplit[1]");
    sfree(fractionSplit, cname, fname, "fractionSplit");
  }

}

void AlgActionRationalSplit::heatbath() {
  rat->checkSplit();
  rat->heatbath();
}

Float AlgActionRationalSplit::energy() {
  return rat->energy();
}

void AlgActionRationalSplit::cost(CgStats *cg_stats_global) {
  rat->cost(cg_stats_global);
}

void AlgActionRationalSplit::prepare_fg(Matrix * force, Float dt_ratio)
{
  const char fname[] = "prepare_fg(M*,F)";
  ERR.General(cname, fname, "Force gradient evolution has not been implemented for\n"
              "this class yet. You may delete this stub(the whole prepare_fg() function)\n"
              "and recompile. Then it will use the inherited verion. Hopefully that will\n"
              "work, but no test has been made yet.\n");
}

void AlgActionRationalSplit::evolve(Float dt, int steps) {
  rat->evolve(dt, steps, fractionSplit);
}

int AlgActionRationalSplit::getNmass() {
  return rat->getNmass();
}

Float AlgActionRationalSplit::getMass(int i) {
  return rat->getMass(i);
}

CPS_END_NAMESPACE
