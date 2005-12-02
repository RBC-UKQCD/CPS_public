#include<alg/int_arg.h>
#include<util/error.h>

CPS_START_NAMESPACE

void ActionBilinearArg::resize(int nmass) {

  bilinears.bilinears_len = nmass;
  bilinears.bilinears_val = new BilinearDescr[nmass];

}

void ActionRationalArg::resize(int nmass) {

  bi_arg.resize(nmass);
  rationals.rationals_len = nmass;
  rationals.rationals_val = new RationalDescr[nmass];

}

void ActionRationalArg::resize(int mass, int deg_md, int deg_mc) {

  if (rationals.rationals_len > mass) {

    rationals.rationals_val[mass].md_approx.stop_rsd.stop_rsd_len = deg_md;
    rationals.rationals_val[mass].md_approx.stop_rsd.stop_rsd_val = 
      new Float[deg_md];

    rationals.rationals_val[mass].mc_approx.stop_rsd.stop_rsd_len = deg_mc;
    rationals.rationals_val[mass].mc_approx.stop_rsd.stop_rsd_val = 
      new Float[deg_mc];

  } else {
    char *cname = "ActionRationalArg";
    char *fname = "resize(int mass, int deg_md, int deg_mc)";
    ERR.General(cname, fname, "mass > Nmass");
  }

}

void ActionRationalSplitArg::resize(int nmass) {

  fractionSplit.fractionSplit_len = nmass;
  fractionSplit.fractionSplit_val = new SplitRange[nmass];  

}

void ActionBosonArg::resize(int nmass) {

  bi_arg.resize(nmass);
  bosons.bosons_len = nmass;
  bosons.bosons_val = new BosonDescr[nmass];

}

void ActionFermionArg::resize(int nmass) {

  bi_arg.resize(nmass);
  fermions.fermions_len = nmass;
  fermions.fermions_val = new FermionDescr[nmass];

}

CPS_END_NAMESPACE
