#include<alg/int_arg.h>
#include<util/error.h>

CPS_START_NAMESPACE

void ActionBilinearArg::resize(u_int nmass) {

  bilinears.bilinears_len = nmass;
  bilinears.bilinears_val = new BilinearDescr[nmass];

}

void ActionRationalArg::resize(u_int nmass) {

  bi_arg.resize(nmass);
  rationals.rationals_len = nmass;
  rationals.rationals_val = new RationalDescr[nmass];

}

void ActionRationalArg::resize(u_int mass, int deg_md, int deg_mc) {

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

void ActionRationalSplitArg::resize(u_int nmass) {

  fractionSplit.fractionSplit_len = nmass;
  fractionSplit.fractionSplit_val = new SplitRange[nmass];  

}

void ActionBosonArg::resize(u_int nmass) {

  bi_arg.resize(nmass);
  bosons.bosons_len = nmass;
  bosons.bosons_val = new BosonDescr[nmass];

}

void ActionFermionArg::resize(u_int nmass) {

  bi_arg.resize(nmass);
  fermions.fermions_len = nmass;
  fermions.fermions_val = new FermionDescr[nmass];

}

void ActionQuotientArg::resize(u_int nmass) {

  bi_arg.resize(nmass);
  quotients.quotients_len = nmass;
  quotients.quotients_val = new QuotientDescr[nmass];

}

void ActionRationalQuotientArg::resize(u_int nmass) {

  bi_arg.resize(nmass);
  bsn_mass.bsn_mass_len = nmass;
  bsn_mass.bsn_mass_val = new Float[nmass];
  frm_mass.frm_mass_len = nmass;
  frm_mass.frm_mass_val = new Float[nmass];
  bosons.bosons_len = nmass;
  bosons.bosons_val = new RationalDescr[nmass];
  fermions.fermions_len = nmass;
  fermions.fermions_val = new RationalDescr[nmass];

}

void ActionRationalQuotientArg::resize(u_int mass, int frm_deg_md, int frm_deg_mc, 
			       int bsn_deg_md, int bsn_deg_mc)
{

  if (fermions.fermions_len > mass) {

    fermions.fermions_val[mass].md_approx.stop_rsd.stop_rsd_len = frm_deg_md;
    fermions.fermions_val[mass].md_approx.stop_rsd.stop_rsd_val = 
      new Float[frm_deg_md];

    fermions.fermions_val[mass].mc_approx.stop_rsd.stop_rsd_len = frm_deg_mc;
    fermions.fermions_val[mass].mc_approx.stop_rsd.stop_rsd_val = 
      new Float[frm_deg_mc];

  }

  if (bosons.bosons_len > mass) {
    bosons.bosons_val[mass].md_approx.stop_rsd.stop_rsd_len = bsn_deg_md;
    bosons.bosons_val[mass].md_approx.stop_rsd.stop_rsd_val = 
      new Float[bsn_deg_md];

    bosons.bosons_val[mass].mc_approx.stop_rsd.stop_rsd_len = bsn_deg_mc;
    bosons.bosons_val[mass].mc_approx.stop_rsd.stop_rsd_val = 
      new Float[bsn_deg_mc];

  } else {
    char *cname = "ActionRationalArg";
    char *fname = "resize(i,i,i,i,i)";
    ERR.General(cname, fname, "mass > Nmass");
  }

}

ApproxDescr::ApproxDescr(){
  approx_type=RATIONAL_APPROX_POWER;
  bounds_type=RATIONAL_BOUNDS_AUTOMATIC;
  lambda_low=0.;
  lambda_high=0.;
  stop_rsd.stop_rsd_len=0;
  stop_rsd.stop_rsd_val=NULL;
}

CPS_END_NAMESPACE
