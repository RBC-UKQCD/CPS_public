
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

};

class BilinearDescr {
  Float mass;
  int max_num_iter;
};

class ActionBilinearArg {

  memfun void resize(int nmass);

  FclassType fermion;
  BilinearDescr bilinears<>;
  ActionArg action_arg;

};

class ApproxDescr {

  Float lambda_low;
  Float lambda_high;
  Float stop_rsd<>;

};

class RationalDescr {

  FieldType field_type;
  int power_num;
  int power_den;
  ApproxDescr md_approx;
  ApproxDescr mc_approx;

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

  memfun void resize(int nmass);
  memfun void resize(int mass, int deg_md, int deg_mc);

  ActionBilinearArg bi_arg;
   
  //! The precision used in the Remez algorithm of the RHMC approximation.
  long precision;

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
  ActionArg action_arg;

};

