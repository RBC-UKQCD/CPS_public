
/*! A structure to hold the solver parameters.*/
class CgArg {

  Float mass;			/*!<  The mass parameter. */
  Float epsilon;             /*!< ~~The epsilon parameter for twisted mass Wilson fermions */
  
  int max_num_iter;		/*!<  The maximum number of solver
                                  iterations to do. */
  
  Float stop_rsd;		/*!<  The target residual. */
  Float true_rsd;
  
  enum RitzMatType RitzMatOper; /*!< Which operator to determine eigenvalues
                                  of, if any. */
  enum InverterType Inverter;     /*!< Which solver to use.*/
  
  int bicgstab_n;             /* BiCGstab(n) parameter. */
  
  memfun CgArg();
};

class MdwfArg {
  /*!< parameters specific to Mobius DWF. */
  Float b5<>;
  Float c5<>;
  /* equivalent of DWF height in Mobius fermions */
  Float M5;

  /* CgArg used by the Mobius library */
  CgArg cg_arg;

  /* stopping condition array for Mobius approximation of DWF */
  /* fermions, the Mobius Dirac operator(DiracOpMdwf) does not */ 
  /* interpret this array. */
  Float rsd_vec<>;

  /* Set to none zero to use single precision in calculation. */
  int use_single_precision;
};
