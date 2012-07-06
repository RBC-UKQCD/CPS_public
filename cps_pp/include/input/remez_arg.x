class RemezArg{

  RationalApproxType approx_type;

  int valid_approx;

  int degree;

  Float lambda_low;
  Float lambda_high;

  int power_num;
  int power_den;

  FieldType field_type;

  Float error;
  Float norm;
  Float residue[MAX_RAT_DEGREE];
  Float pole[MAX_RAT_DEGREE];

  Float norm_inv;
  Float residue_inv[MAX_RAT_DEGREE];
  Float pole_inv[MAX_RAT_DEGREE];

  long precision;

  //!< The shifted mass parameter used when doing staggered Hasenbusch
  Float delta_m;
  memfun RemezArg();
};

class RationalQuotientRemezArg {
      RemezArg bsn_md<>;
      RemezArg bsn_mc<>;
      RemezArg frm_md<>;
      RemezArg frm_mc<>;
};
