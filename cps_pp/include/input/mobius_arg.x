class MobiusArg {

  /*!< parameters specific to Mobius DWF. */
  Float mobius_b_coeff;
  Float mobius_c_coeff;
  Float zmobius_b_coeff<>;
  Float zmobius_c_coeff<>;

  /* 5th dim size */
  int ls;

  /* equivalent of DWF height in Mobius fermions */
  Float M5;

  /* CgArg used by the Mobius library */
  CgArg cg;

  /* stopping condition array for Mobius approximation of DWF */
  /* fermions, the Mobius Dirac operator(DiracOpMdwf) does not */ 
  /* interpret this array. */
  Float rsd_vec<>;

  /* Set to none zero to use single precision in calculation. */
  int use_single_precision;

  memfun   MobiusArg (  ) ;

};
