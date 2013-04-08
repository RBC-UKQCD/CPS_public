
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

  string  fname_eigen<>;    /* file name for Low Eigen Modes */

  int neig;         /* number of eigenvectors */

  Float eigen_shift; /* shift for eigen spectrum, place holder. Don't input */
  Float ama_stop_rsd;         /*!<  The target residual for AMA. */
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

  /* Set to none zero to use the Mdwf library for DWF solve */
  int use_mdwf_for_dwf;
};

class C5State {
  Float val;
  int dwf_cg;
};

/* MdwfTuning, used for tuning Mobius preconditioned DWF inverter. */
class MdwfTuning {
  /* stage == 0: tuning finished */
  /* stage == 1: test zero_time */
  /* stage == 2: test c5 */
  /* stage == 3: test rsd */
  /* stage == 4: test rc */
  int stage;

  /* L_s state */
  /* Both ls_min and ls_max must be even. */
  int ls_min;
  int ls_max;

  /* Primary parameter showing which ls we are working on.  */
  /* ls = ls_min + 2*index; */
  int index;

  int opti_index;
  Float opti_time;

  /* optimal MdwfArg and corresponding minimum time for each ls. */
  MdwfArg mdwf_arg<>;

  /* Mdwf B5/C5 value tuning. */
  C5State c5<>;

  /* Mdwf CG inversion residue tuning. */
  Float rsd_val;
  Float rsd_time;

  /* Granularity when searching for an optimal residue. Must be a factor larger than 1. */
  Float rsd_granularity;

  /* controls the size of the region we test for optimal c5 values. */
  Float c5_range;

  /* Mdwf CG restart count tuning. */
  int rc_max;
  int rc_val;
  Float rc_time;
};

/* use this class to initialize the MdwfTuning class. */
class MdwfTuningInitArg {
  /* smallest Ls in Mobius that we wish to test. */
  int ls_min;
  /* largest Ls in Mobius that we wish to test. */
  int ls_max;

  /* maximum allowed restart count. */
  int rc_max;

  /* set this flag so we will use the MDWF inverter for DWF(so we don't */
  /* use the DWF inverter in CPS). */
  int use_mdwf_for_dwf;

  /* Set to none zero to use single precision in calculation. */
  /* Note: there will always be a double precision fixup step in the end. */
  int use_single_precision;

  /* controls the size of the region we test for optimal c5 values. */
  Float c5_range;

  /* Granularity when searching for an optimal residue. Must be a factor larger than 1. */
  Float rsd_granularity;

  /* file to save the tuning state (a vml file that encodes an instance */
  /* of class MdwfTuning).   */
  /* note: if this file exists, we will read the MdwfTuning from the */
  /* file, and the rest parameters in MdwfTuningInitArg will have no */
  /* effect! */
  string tuning_fn<>;

  /* file to record the tuning test result. */
  string tuning_record_fn<>;
};
