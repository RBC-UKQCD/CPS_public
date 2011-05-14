
/* Required by hmd_Arg and by ActionRationalArg */
//typedef Float FRatArray[MAX_RAT_DEGREE];

/*
class IntArg {

};
*/

class IntABArg {

  IntegratorType type;
  int A_steps;
  int B_steps;
  IntegratorLevel level;
  Float lambda;

};

class ActionArg {

  ForceMeasure force_measure;
  string force_label<>;

};

class BilinearDescr {
  Float mass;
  int max_num_iter;
};

class ActionBilinearArg {

  memfun void resize(u_int nmass);

  FclassType fermion;
  BilinearDescr bilinears<>;
  ActionArg action_arg;

};

class ApproxDescr {

  RationalApproxType approx_type;
  RationalBoundsType bounds_type;
  Float lambda_low;
  Float lambda_high;
  Float stop_rsd<>;
  memfun ApproxDescr(void);

};

class RationalDescr {

  FieldType field_type;
  int power_num;
  int power_den;
  //! The precision used in the Remez algorithm of the RHMC approximation.
  long precision;
  //! The multiplier we used for additional CG inversion in force
  //! gradient integrator.
  Float stop_rsd_fg_mult;
  ApproxDescr md_approx;
  ApproxDescr mc_approx;
  //! The boson mass parameter used for Hasenbusch trick (staggered only)
  Float stag_bsn_mass;
} ;

class EigenDescr {
  //! Whether the approximation bounds are checked
  EigenMeasure eigen_measure;

  //! Tolerance in the eigenvalue measurement
  Float stop_rsd;

  int max_num_iter;

  //! Where the eigenvalues should be output to
  string eig_lo_stem<>;
  string eig_hi_stem<>;
} ;

class ActionRationalArg {

  memfun void resize(u_int nmass);
  memfun void resize(u_int mass, int deg_md, int deg_mc);

  ActionBilinearArg bi_arg;
   
  int remez_generate;
  string rat_poles_file<>;

  RationalDescr rationals<>;

  EigenDescr eigen;

};

class SplitRange {
  int split_low;
  int split_high;
} ;

class ActionRationalSplitArg {

  memfun void resize(u_int nmass);
  SplitRange fractionSplit<>;

};

class BosonDescr {
  //! ~~epsilon parameter for twisted mass wilson fermions
  Float epsilon;
  Float stop_rsd_hb;
} ;

class ActionBosonArg {

  memfun void resize(u_int nmass);
  ActionBilinearArg bi_arg;
  BosonDescr bosons<>;

};

class FermionDescr {
  //! ~~epsilon parameter for twisted mass wilson fermions
  Float epsilon;
  int   chrono;
  //! The multiplier we used for additional CG inversion in force
  //! gradient integrator.
  Float stop_rsd_fg_mult;
  Float stop_rsd_md;
  Float stop_rsd_mc;
} ;

class ActionFermionArg {

  memfun void resize(u_int nmass);

  ActionBilinearArg bi_arg;

  //! Chonological inversion parameter
  FermionDescr fermions<>;

};

class QuotientDescr {
  Float bsn_mass;
  //! ~~epsilon parameter for twisted mass wilson fermions
  Float bsn_mass_epsilon;
  Float frm_mass;
  //! ~~epsilon parameter for twisted mass wilson fermions
  Float frm_mass_epsilon;
  int   chrono;
  Float stop_rsd_hb;
  //! The multiplier we used for additional CG inversion in force
  //! gradient integrator.
  Float stop_rsd_fg_mult;
  Float stop_rsd_md;
  Float stop_rsd_mc;
} ;

class ActionQuotientArg {

  memfun void resize(u_int nmass);

  // Mass parameter here is dummy
  ActionBilinearArg bi_arg;

  //!< Quotient parameters
  QuotientDescr quotients<>;

};

class ActionRationalQuotientArg {

  memfun void resize(u_int nmass);
  memfun void resize(u_int mass, int frm_deg_md, int frm_deg_mc, 
	int bsn_deg_md, int bsn_deg_mc);

  ActionBilinearArg bi_arg;
   
  //! the allowed sprectral spread in the rational approximation
  Float spread;

  int remez_generate;
  string rat_poles_file<>;

  Float bsn_mass<>;
  Float frm_mass<>;

  RationalDescr bosons<>;
  RationalDescr fermions<>;

  EigenDescr eigen;
};

class ActionGaugeArg {

  GclassType gluon;
  ActionArg action_arg;

};
