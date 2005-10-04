
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

class BilinearDescr {
  Float mass;
  int max_num_iter;
};

class ActionBilinearArg {

  memfun void resize(int nmass);

  FclassType fermion;
  BilinearDescr bilinears<>;

};

class ApproxDescr {

  Float stop_rsd;

};

class RationalDescr {

  FieldType field_type;
  Float lambda_low;
  Float lambda_high;
  int power_num;
  int power_den;
  ApproxDescr md_approx<>;
  ApproxDescr mc_approx<>;

} ;

class ActionRationalArg {

  memfun void resize(int nmass);
  memfun void resize(int mass, int deg_md, int deg_mc);

  ActionBilinearArg bi_arg;
   
  //! Whether the RHMC approximation is static or dynamic.
  RatApproxType approx_type;

  //! the allowed sprectral spread in the rational approximation
  Float spread;

  //! The precision used in the Remez algorithm of the RHMC approximation.
  long precision;

  int remez_generate;
  string rat_poles_file<>;

  RationalDescr rationals<>;

};

class SplitRange {
  int split_low;
  int split_high;
} ;

class ActionRationalSplitArg {

  memfun void resize(int nmass);
  SplitRange fractionSplit<>;

};

class BosonDescr {
  Float stop_rsd_hb;
} ;

class ActionBosonArg {

  memfun void resize(int nmass);
  ActionBilinearArg bi_arg;
  BosonDescr bosons<>;

};

class FermionDescr {
  int   chrono;
  Float stop_rsd_md;
  Float stop_rsd_mc;
} ;

class ActionFermionArg {

  memfun void resize(int nmass);

  ActionBilinearArg bi_arg;

  //! Chonological inversion parameter
  FermionDescr fermions<>;

};

class ActionGaugeArg {

  GclassType gluon;

};

