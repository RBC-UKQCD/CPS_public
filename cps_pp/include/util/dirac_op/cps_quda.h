#ifndef INCLUDED_CPS_QUDA_H
#define INCLUDED_CPS_QUDA_H

#define MAX(a,b) ((a)>(b)?(a):(b))

CPS_START_NAMESPACE class CPSQuda {
public:
  static constexpr const char *cname = "CPSQuda";

  static QudaPrecision setPrecision (CudaPrecision prec)
  {
    const char *fname = "setPrecision";

    switch (prec) {
    case CUDA_HALF_PRECISION:
      VRB.Debug (cname, fname, "CUDA_HALF_PRECISION\n");
      return QUDA_HALF_PRECISION;
    case CUDA_SINGLE_PRECISION:
      VRB.Debug (cname, fname, "CUDA_SINGLE_PRECISION\n");
      return QUDA_SINGLE_PRECISION;
    case CUDA_DOUBLE_PRECISION:
      VRB.Debug (cname, fname, "CUDA_DOUBLE_PRECISION\n");
      return QUDA_DOUBLE_PRECISION;
    default:
      ERR.General (cname, fname, "Undefined precision %d\n", prec);
    }
  }

  static QudaReconstructType setReconstruct (CudaReconstructType recon)
  {
    const char *fname = "setReconstruct";

    switch (recon) {
    case CUDA_RECONSTRUCT_8:
      return QUDA_RECONSTRUCT_8;
    case CUDA_RECONSTRUCT_12:
      return QUDA_RECONSTRUCT_12;
    case CUDA_RECONSTRUCT_NO:
      return QUDA_RECONSTRUCT_NO;
    default:
      ERR.General (cname, fname, "Undefined reconstruct type %d\n", recon);
    }
  }

  static void ParamSetupGauge (QudaArg & quda_param,
                               QudaGaugeParam & gauge_param, int if_double)
  {

    const char *fname = "ParamSetupGauge()";
    //-----------------------------------
    // Parameter setting for gauge data
    //-----------------------------------

    // Set the CUDA precisions
    gauge_param.reconstruct = setReconstruct (quda_param.reconstruct);
    gauge_param.cuda_prec = setPrecision (quda_param.gauge_prec);

    // Set the CUDA sloppy precisions
    gauge_param.reconstruct_sloppy =
      setReconstruct (quda_param.reconstruct_sloppy);
    gauge_param.cuda_prec_sloppy = setPrecision (quda_param.gauge_prec_sloppy);

//      if(sizeof(Float) == sizeof(double)) {
    VRB.Result (cname, fname, "if_double=%d %d\n", if_double);
    gauge_param.cpu_prec = QUDA_DOUBLE_PRECISION;
//        gauge_param.cpu_prec = QUDA_SINGLE_PRECISION;

    gauge_param.X[0] = GJP.XnodeSites ();
    gauge_param.X[1] = GJP.YnodeSites ();
    gauge_param.X[2] = GJP.ZnodeSites ();
    gauge_param.X[3] = GJP.TnodeSites ();
    gauge_param.anisotropy = GJP.XiBare ();
//      gauge_param.cuda_prec_precondition = QUDA_DOUBLE_PRECISION;
    gauge_param.cuda_prec_precondition = QUDA_HALF_PRECISION;
    gauge_param.reconstruct_precondition =
      setReconstruct (quda_param.reconstruct_sloppy);

    if (GJP.XiDir () != 3) {
      ERR.General (cname, fname, "Anisotropy direction not supported\n");
    }
    //----------------------------------------------------------------------------------------
    // QUDA_FLOAT_GAUGE_ORDER = 1
    // QUDA_FLOAT2_GAUGE_ORDER = 2, // no reconstruct and double precision
    // QUDA_FLOAT4_GAUGE_ORDER = 4, // 8 and 12 reconstruct half and single
    // QUDA_QDP_GAUGE_ORDER,        // expect *gauge[4], even-odd, row-column color
    // QUDA_CPS_WILSON_GAUGE_ORDER, // expect *gauge, even-odd, mu inside, column-row color
    // QUDA_MILC_GAUGE_ORDER,       // expect *gauge, even-odd, mu inside, row-column order
    //
    // Multi-GPU case, we have to use QDP format of gauge data
    //----------------------------------------------------------------------------------------

    gauge_param.gauge_order = QUDA_CPS_WILSON_GAUGE_ORDER;

    gauge_param.gauge_fix = QUDA_GAUGE_FIXED_NO;
    gauge_param.type = QUDA_WILSON_LINKS;

    for (int d = 0; d < 3; d++) {
      if (GJP.Bc (d) != BND_CND_PRD) {
        ERR.General (cname, fname, "Boundary condition not supported\n");
      }
    }
    if (GJP.Tbc () == BND_CND_PRD) {
      gauge_param.t_boundary = QUDA_PERIODIC_T;
    } else {
      gauge_param.t_boundary = QUDA_ANTI_PERIODIC_T;
    }
  }

  static void ParamSetup_dwf (QudaArg & quda_param,
                              QudaGaugeParam & gauge_param,
                              QudaInvertParam & inv_param, int mat_type =
                              1, int if_double = 1) {
    const char *fname =
      "ParamSetup_dwf(QudaArg&, QudaGaugeParam&, QudaInvertParam&, i)";


    ParamSetupGauge (quda_param, gauge_param, if_double);

    if (if_double) {
      inv_param.cpu_prec = QUDA_DOUBLE_PRECISION;
    } else {
      inv_param.cpu_prec = QUDA_SINGLE_PRECISION;
    }
    //-------------------------------------------------
    // Parameter setting for fermion matrix inversion
    //-------------------------------------------------
    inv_param.cuda_prec = setPrecision (quda_param.spinor_prec);
    inv_param.cuda_prec_sloppy = setPrecision (quda_param.spinor_prec_sloppy);
    inv_param.cuda_prec_precondition = QUDA_HALF_PRECISION;

    inv_param.reliable_delta = quda_param.reliable_delta;

    inv_param.Ls = GJP.SnodeSites ();

    //-----------------------------
    // Possible dslash types
    //-----------------------------
    // QUDA_WILSON_DSLASH
    // QUDA_CLOVER_WILSON_DSLASH
    // QUDA_DOMAIN_WALL_DSLASH
    // QUDA_DOMAIN_WALL_4D_DSLASH
    // QUDA_MOBIUS_DWF_DSLASH
    // QUDA_ASQTAD_DSLASH
    // QUDA_TWISTED_MASS_DSLASH
    //-----------------------------
    inv_param.dslash_type = QUDA_DOMAIN_WALL_DSLASH;

    //------------------------------------
    // Possible normalizations
    //------------------------------------
    // QUDA_KAPPA_NORMALIZATION
    // QUDA_MASS_NORMALIZATION
    // QUDA_ASYMMETRIC_MASS_NORMALIZATION
    //------------------------------------
    inv_param.mass_normalization = QUDA_KAPPA_NORMALIZATION;

    inv_param.dagger = QUDA_DAG_NO;

    //------------------------
    // Domain wall parameters
    //------------------------
    inv_param.m5 = -GJP.DwfHeight ();

    switch (mat_type) {
    case 0:
      inv_param.solution_type = QUDA_MATPC_SOLUTION;
      break;
    case 1:
      inv_param.solution_type = QUDA_MATPCDAG_MATPC_SOLUTION;
      break;
    default:
      ERR.General (cname, fname, "Matrix solution type not defined\n");
    }

    inv_param.matpc_type = QUDA_MATPC_ODD_ODD;
    inv_param.preserve_source = QUDA_PRESERVE_SOURCE_NO;
    inv_param.gamma_basis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;

    inv_param.dirac_order = QUDA_DIRAC_ORDER;

    inv_param.input_location = QUDA_CPU_FIELD_LOCATION;
    inv_param.output_location = QUDA_CPU_FIELD_LOCATION;
    inv_param.tune = QUDA_TUNE_YES;
    inv_param.use_init_guess = QUDA_USE_INIT_GUESS_YES;

    //------------------------
    // Possible verbose types
    //------------------------
    // QUDA_SILENT
    // QUDA_SUMMARIZE
    // QUDA_VERBOSE
    // QUDA_DEBUG_VERBOSE
    //------------------------
    inv_param.verbosity = QUDA_VERBOSE;

    // Domain decomposition preconditioner parameters
    inv_param.inv_type_precondition = QUDA_INVALID_INVERTER;
    inv_param.schwarz_type = QUDA_ADDITIVE_SCHWARZ;
    inv_param.precondition_cycle = 1;
    inv_param.tol_precondition = 1.0e-01;
    inv_param.maxiter_precondition = quda_param.maxiter_precondition;
    inv_param.verbosity_precondition = QUDA_VERBOSE;
    inv_param.omega = 1.0;

    gauge_param.ga_pad = 0;
    inv_param.sp_pad = 0;
    inv_param.cl_pad = 0;

#ifdef USE_QMP
    //---------------------------------------------------------------------
    // This part is needed to make buffer memory space for multi-GPU comms
    //---------------------------------------------------------------------
    int x_face_size =
      gauge_param.X[1] * gauge_param.X[2] * gauge_param.X[3] / 2;
    int y_face_size =
      gauge_param.X[0] * gauge_param.X[2] * gauge_param.X[3] / 2;
    int z_face_size =
      gauge_param.X[0] * gauge_param.X[1] * gauge_param.X[3] / 2;
    int t_face_size =
      gauge_param.X[0] * gauge_param.X[1] * gauge_param.X[2] / 2;
    int pad_size = MAX (x_face_size, y_face_size);
    pad_size = MAX (pad_size, z_face_size);
    pad_size = MAX (pad_size, t_face_size);
    gauge_param.ga_pad = pad_size;
#endif
  }

  static void ParamSetup_mdwf (QudaArg & quda_param,
                               QudaGaugeParam & gauge_param,
                               QudaInvertParam & inv_param, int mat_type =
                               1, int if_double = 1) {
    const char *fname =
      "ParamSetup_mdwf(QudaArg&, QudaGaugeParam&, QudaInvertParam&, i)";

    ParamSetupGauge (quda_param, gauge_param, if_double);

    if (if_double) {
      gauge_param.cuda_prec = QUDA_DOUBLE_PRECISION;
      inv_param.cpu_prec = QUDA_DOUBLE_PRECISION;
      inv_param.cuda_prec = QUDA_DOUBLE_PRECISION;
    } else {
      gauge_param.cuda_prec = QUDA_SINGLE_PRECISION;
      inv_param.cpu_prec = QUDA_SINGLE_PRECISION;
      inv_param.cuda_prec = QUDA_SINGLE_PRECISION;
    }

    //-------------------------------------------------
    // Parameter setting for fermion matrix inversion
    //-------------------------------------------------
//      inv_param.cuda_prec = setPrecision(quda_param.spinor_prec);
    inv_param.cuda_prec_sloppy = setPrecision (quda_param.spinor_prec_sloppy);

    inv_param.reliable_delta = quda_param.reliable_delta;

    inv_param.Ls = GJP.SnodeSites ();

    //-----------------------------
    // Possible dslash types
    //-----------------------------
    // QUDA_WILSON_DSLASH
    // QUDA_CLOVER_WILSON_DSLASH
    // QUDA_DOMAIN_WALL_DSLASH
    // QUDA_DOMAIN_WALL_4D_DSLASH
    // QUDA_MOBIUS_DWF_DSLASH
    // QUDA_ASQTAD_DSLASH
    // QUDA_TWISTED_MASS_DSLASH
    //-----------------------------
    inv_param.dslash_type = QUDA_MOBIUS_DWF_DSLASH;

    //------------------------------------
    // Possible normalizations
    //------------------------------------
    // QUDA_KAPPA_NORMALIZATION
    // QUDA_MASS_NORMALIZATION
    // QUDA_ASYMMETRIC_MASS_NORMALIZATION
    //------------------------------------
    inv_param.mass_normalization = QUDA_KAPPA_NORMALIZATION;

    inv_param.dagger = QUDA_DAG_NO;

    //------------------------
    // Domain wall parameters
    //------------------------
    inv_param.m5 = -GJP.DwfHeight ();

    //--------------------------------------------------------------
    // Build MDWF coefficients b_5 & c_5
    //--------------------------------------------------------------
    // Currently, the CPS program uses constant values for b_5, c_5
    // coefficients not array type data.
    //--------------------------------------------------------------
    VRB.Result (cname, fname, "(b_5,c_5) = (%e,%e)\n", GJP.Mobius_b (),
                GJP.Mobius_c ());
    for (int s = 0; s < GJP.SnodeSites (); s++) {
      inv_param.b_5[s] = GJP.Mobius_b ();
      inv_param.c_5[s] = GJP.Mobius_c ();
    }

    switch (mat_type) {
    case 0:
      inv_param.solution_type = QUDA_MATPC_SOLUTION;
      break;
    case 1:
      inv_param.solution_type = QUDA_MATPCDAG_MATPC_SOLUTION;
      break;
    default:
      ERR.General (cname, fname, "Matrix solution type not defined\n");
    }

    switch (GJP.ZMobius_PC_Type ()) {
    case ZMOB_PC_ORIG:
//          inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN_ASYMMETRIC; break;
      inv_param.matpc_type = QUDA_MATPC_ODD_ODD_ASYMMETRIC;
      break;
    case ZMOB_PC_SYM1:
//          inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN; break;
      inv_param.matpc_type = QUDA_MATPC_ODD_ODD;
      break;
    default:
      ERR.General (cname, fname,
                   "Precondition type should be ZMOB_PC_ORIG or ZMOB_PC_SYM1\n");
    }
    inv_param.preserve_source = QUDA_PRESERVE_SOURCE_NO;
    inv_param.gamma_basis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;

    // inv_param.dirac_order = QUDA_CPS_WILSON_DIRAC_ORDER;
    inv_param.dirac_order = QUDA_DIRAC_ORDER;
    // inv_param.siteSubset = QUDA_FULL_SITE_SUBSET;

    inv_param.input_location = QUDA_CPU_FIELD_LOCATION;
    inv_param.output_location = QUDA_CPU_FIELD_LOCATION;
    inv_param.tune = QUDA_TUNE_YES;
    inv_param.use_init_guess = QUDA_USE_INIT_GUESS_YES;
    inv_param.use_alternative_reliable = 1;


    //------------------------
    // Possible verbose types
    //------------------------
    // QUDA_SILENT
    // QUDA_SUMMARIZE
    // QUDA_VERBOSE
    // QUDA_DEBUG_VERBOSE
    //------------------------
    inv_param.verbosity = QUDA_VERBOSE;

    // Domain decomposition preconditioner parameters
    inv_param.inv_type_precondition = QUDA_INVALID_INVERTER;
    inv_param.schwarz_type = QUDA_ADDITIVE_SCHWARZ;
    inv_param.precondition_cycle = 1;
    inv_param.tol_precondition = 1.0e-01;
    inv_param.maxiter_precondition = quda_param.maxiter_precondition;
    inv_param.verbosity_precondition = QUDA_VERBOSE;
    inv_param.omega = 1.0;

    gauge_param.ga_pad = 0;
    inv_param.sp_pad = 0;
    inv_param.cl_pad = 0;

#ifdef USE_QMP
    //---------------------------------------------------------------------
    // This part is needed to make buffer memory space for multi-GPU comms
    //---------------------------------------------------------------------
    int x_face_size =
      gauge_param.X[1] * gauge_param.X[2] * gauge_param.X[3] / 2;
    int y_face_size =
      gauge_param.X[0] * gauge_param.X[2] * gauge_param.X[3] / 2;
    int z_face_size =
      gauge_param.X[0] * gauge_param.X[1] * gauge_param.X[3] / 2;
    int t_face_size =
      gauge_param.X[0] * gauge_param.X[1] * gauge_param.X[2] / 2;
    int pad_size = MAX (x_face_size, y_face_size);
    pad_size = MAX (pad_size, z_face_size);
    pad_size = MAX (pad_size, t_face_size);
    gauge_param.ga_pad = pad_size;
#endif
  }
  static void ParamSetup_mdwf (QudaArg & quda_param,
                               QudaGaugeParam & gauge_param,
                               QudaEigParam & eig_param, int mat_type =
                               1, int if_double = 1) {
    const char *fname =
      "ParamSetup_mdwf(QudaArg&, QudaGaugeParam&, QudaInvertParam&, i)";

    ParamSetupGauge (quda_param, gauge_param, if_double);
    if (if_double) {
      gauge_param.cuda_prec = QUDA_DOUBLE_PRECISION;
      eig_param.cuda_prec_ritz = QUDA_DOUBLE_PRECISION;

    } else {
      gauge_param.cuda_prec = QUDA_SINGLE_PRECISION;
      eig_param.cuda_prec_ritz = QUDA_SINGLE_PRECISION;
    }

    eig_param.location = QUDA_CPU_FIELD_LOCATION;

    gauge_param.ga_pad = 0;

#ifdef USE_QMP
    //---------------------------------------------------------------------
    // This part is needed to make buffer memory space for multi-GPU comms
    //---------------------------------------------------------------------
    int x_face_size =
      gauge_param.X[1] * gauge_param.X[2] * gauge_param.X[3] / 2;
    int y_face_size =
      gauge_param.X[0] * gauge_param.X[2] * gauge_param.X[3] / 2;
    int z_face_size =
      gauge_param.X[0] * gauge_param.X[1] * gauge_param.X[3] / 2;
    int t_face_size =
      gauge_param.X[0] * gauge_param.X[1] * gauge_param.X[2] / 2;
    int pad_size = MAX (x_face_size, y_face_size);
    pad_size = MAX (pad_size, z_face_size);
    pad_size = MAX (pad_size, t_face_size);
    gauge_param.ga_pad = pad_size;
#endif
  }

  typedef struct QudaParams
  {
    QudaGaugeParam gauge_param;
    QudaInvertParam inv_param;
  } QudaParams;
};

CPS_END_NAMESPACE
#endif
