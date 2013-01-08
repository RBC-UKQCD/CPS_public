#include <config.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>


CPS_START_NAMESPACE
//----------------------------------------------------------------
//
// wilson_quda.C
//
// Interface to QUDA Wilson inverter
// 
//------------------------------------------------------------------

CPS_END_NAMESPACE

#include <util/dirac_op.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/error.h>
#include <util/wilson.h>
#include <comms/glb.h>
#include <util/time_cps.h>

#ifdef USE_QUDA
#include <invert_quda.h>
#define MAX(a,b) ((a)>(b)?(a):(b))

CPS_START_NAMESPACE


QudaPrecision setPrecision_wil(CudaPrecision prec) {
  switch (prec) {
    case CUDA_HALF_PRECISION:
      return QUDA_HALF_PRECISION;
    case CUDA_SINGLE_PRECISION:
      return QUDA_SINGLE_PRECISION;
    case CUDA_DOUBLE_PRECISION:
      return QUDA_DOUBLE_PRECISION;
    default:
      ERR.General("", "", "Undefined precision %d\n", prec);
  }
}

QudaReconstructType setReconstruct_wil(CudaReconstructType recon) {
  switch (recon) {
    case CUDA_RECONSTRUCT_8:
      return QUDA_RECONSTRUCT_8;
    case CUDA_RECONSTRUCT_12:
      return QUDA_RECONSTRUCT_12;
    case CUDA_RECONSTRUCT_NO:
      return QUDA_RECONSTRUCT_NO;
  }
}

//--------------------------HJ Kim------------------------------

int DiracOpWilson::QudaInvert(Vector *out, Vector *in, Float *true_res, int mat_type) {

  char *fname = "QudaInvert(V*, V*, F*, int)";

  VRB.ActivateLevel(VERBOSE_FLOW_LEVEL);

  struct timeval start, end;
  gettimeofday(&start,NULL);

  QudaGaugeParam gauge_param = newQudaGaugeParam();
  QudaInvertParam inv_param = newQudaInvertParam();

  int f_size_cb = GJP.VolNodeSites() * lat.FsiteSize() / 2;

  //--------------------------------------
  //  Parameter setting for Gauge Data
  //--------------------------------------

  // set the CUDA precisions
  gauge_param.reconstruct = setReconstruct_wil(QudaParam.reconstruct);
  gauge_param.cuda_prec = setPrecision_wil(QudaParam.gauge_prec);

  // set the CUDA sloppy precisions
  gauge_param.reconstruct_sloppy = setReconstruct_wil(QudaParam.reconstruct_sloppy);
  gauge_param.cuda_prec_sloppy = setPrecision_wil(QudaParam.gauge_prec_sloppy);

  if (sizeof(Float) == sizeof(double)) {
    gauge_param.cpu_prec = QUDA_DOUBLE_PRECISION;
    inv_param.cpu_prec = QUDA_DOUBLE_PRECISION;
  } else {
    gauge_param.cpu_prec = QUDA_SINGLE_PRECISION;
    inv_param.cpu_prec = QUDA_SINGLE_PRECISION;
  }

  gauge_param.X[0] = GJP.XnodeSites();
  gauge_param.X[1] = GJP.YnodeSites();
  gauge_param.X[2] = GJP.ZnodeSites();
  gauge_param.X[3] = GJP.TnodeSites();
  gauge_param.anisotropy = GJP.XiBare();
  gauge_param.cuda_prec_precondition = QUDA_DOUBLE_PRECISION;
  gauge_param.reconstruct_precondition = setReconstruct_wil(QudaParam.reconstruct_sloppy);

  if (GJP.XiDir() != 3) ERR.General(cname, fname, "Anisotropy direction not supported\n");
  
  //---------------------------------------------------
  // QUDA_FLOAT_GAUGE_ORDER = 1
  // QUDA_FLOAT2_GAUGE_ORDER = 2, // no reconstruct and double precision
  // QUDA_FLOAT4_GAUGE_ORDER = 4, // 8 and 12 reconstruct half and single
  // QUDA_QDP_GAUGE_ORDER, // expect *gauge[4], even-odd, row-column color
  // QUDA_CPS_WILSON_GAUGE_ORDER, // expect *gauge, even-odd, mu inside, column-row color
  // QUDA_MILC_GAUGE_ORDER, // expect *gauge, even-odd, mu inside, row-column order
  //
  // MULTI GPU case, we have to use QDP format of gauge data
  //
  //---------------------------------------------------
  
  gauge_param.gauge_order = QUDA_CPS_WILSON_GAUGE_ORDER;

  //---------------------------------------------------
  
  gauge_param.gauge_fix = QUDA_GAUGE_FIXED_NO;
  gauge_param.type = QUDA_WILSON_LINKS;

  for (int d=0; d<3; d++) if (GJP.Bc(d) != BND_CND_PRD) 
    ERR.General(cname, fname, "Boundary condition not supported\n");
  if (GJP.Tbc() == BND_CND_PRD)
    gauge_param.t_boundary = QUDA_PERIODIC_T;
  else
    gauge_param.t_boundary = QUDA_ANTI_PERIODIC_T;
  
  //------------------------------------------
  //  Parameter setting for Matrix invertion
  //------------------------------------------
  inv_param.cuda_prec = setPrecision_wil(QudaParam.spinor_prec);
  inv_param.cuda_prec_sloppy = setPrecision_wil(QudaParam.spinor_prec_sloppy);

  inv_param.maxiter = dirac_arg->max_num_iter;
  inv_param.reliable_delta = QudaParam.reliable_delta;
  
  inv_param.Ls = 1;
  //inv_param.Ls = GJP.SnodeSites();
  
  //--------------------------
  // Possible dslash type
  //--------------------------
  // QUDA_WILSON_DSLASH
  // QUDA_CLOVER_WILSON_DSLASH
  // QUDA_DOMAIN_WALL_DSLASH
  // QUDA_ASQTAD_DSLASH
  // QUDA_TWISTED_MASS_DSLASH
  //--------------------------
  inv_param.dslash_type = QUDA_WILSON_DSLASH;

  //--------------------------------
  // Possible normalization method
  //--------------------------------
  // QUDA_KAPPA_NORMALIZATION
  // QUDA_MASS_NORMALIZATION
  // QUDA_ASYMMETRIC_MASS_NORMALIZATION
  //--------------------------------
  inv_param.kappa = kappa;
  inv_param.mass_normalization = QUDA_KAPPA_NORMALIZATION;
  //inv_param.mass = dirac_arg->mass;
  //inv_param.mass_normalization = QUDA_MASS_NORMALIZATION;
  
  inv_param.dagger = QUDA_DAG_NO;

  switch (mat_type) {
  case 0:
    inv_param.solution_type = QUDA_MATPC_SOLUTION;
    break;
  case 1:
    inv_param.solution_type = QUDA_MATPCDAG_MATPC_SOLUTION;
    break;
  default:
    ERR.General(cname, fname, "Matrix solution type not defined\n");
  }

  inv_param.matpc_type = QUDA_MATPC_ODD_ODD;
  //inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN;
  inv_param.preserve_source = QUDA_PRESERVE_SOURCE_NO;
  inv_param.gamma_basis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;
  //inv_param.gamma_basis = QUDA_UKQCD_GAMMA_BASIS;
  
  inv_param.dirac_order = QUDA_CPS_WILSON_DIRAC_ORDER;
  //inv_param.dirac_order = QUDA_DIRAC_ORDER;

  inv_param.input_location = QUDA_CPU_FIELD_LOCATION;
  inv_param.output_location = QUDA_CPU_FIELD_LOCATION;
  inv_param.tune = QUDA_TUNE_NO;
  inv_param.use_init_guess = QUDA_USE_INIT_GUESS_YES;

  //--------------------------
  // Possible verbose type
  //--------------------------
  // QUDA_SILENT
  // QUDA_SUMMARIZE
  // QUDA_VERBOSE
  // QUDA_DEBUG_VERBOSE
  //--------------------------
  inv_param.verbosity = QUDA_VERBOSE;

  switch (dirac_arg->Inverter) {
  case CG:
    inv_param.inv_type = QUDA_CG_INVERTER;
    inv_param.solve_type = QUDA_NORMEQ_PC_SOLVE;
    break;
  case BICGSTAB:
    inv_param.inv_type = QUDA_BICGSTAB_INVERTER;
    inv_param.solve_type = QUDA_DIRECT_PC_SOLVE;
    break;
  default:
    inv_param.inv_type = QUDA_CG_INVERTER;
    inv_param.solve_type = QUDA_NORMEQ_PC_SOLVE;
    break;
  }
  
  // domain decomposition preconditioner parameters
  inv_param.inv_type_precondition = QUDA_INVALID_INVERTER;
  inv_param.schwarz_type = QUDA_ADDITIVE_SCHWARZ;
  inv_param.precondition_cycle = 1;
  inv_param.tol_precondition = 1e-1;
  inv_param.maxiter_precondition = 10;
  inv_param.verbosity_precondition = QUDA_VERBOSE;
  inv_param.prec_precondition = QUDA_HALF_PRECISION;
  inv_param.omega = 1.0;

  
  gauge_param.ga_pad = 0; // 24*24*24/2;
  inv_param.sp_pad = 0; // 24*24*24/2;
  inv_param.cl_pad = 0; // 24*24*24/2;

#ifdef USE_QMP
  //------------------------------------------
  // This part is needed to make buffer memory
  // space for multi GPU Comm.
  //------------------------------------------
  int x_face_size = gauge_param.X[1]*gauge_param.X[2]*gauge_param.X[3]/2;
  int y_face_size = gauge_param.X[0]*gauge_param.X[2]*gauge_param.X[3]/2;
  int z_face_size = gauge_param.X[0]*gauge_param.X[1]*gauge_param.X[3]/2;
  int t_face_size = gauge_param.X[0]*gauge_param.X[1]*gauge_param.X[2]/2;
  int pad_size =MAX(x_face_size, y_face_size);
  pad_size = MAX(pad_size, z_face_size);
  pad_size = MAX(pad_size, t_face_size);
  gauge_param.ga_pad = pad_size;
#endif

  loadGaugeQuda((void*)gauge_field, &gauge_param);

  Vector *x = (Vector*)smalloc(f_size_cb * sizeof(Float));
  Vector *r = (Vector*)smalloc(f_size_cb * sizeof(Float));
  x->VecZero(f_size_cb); 
  r->VecZero(f_size_cb); 
  
  //----------------------------------------------
  //  Calculate Flops value 
  //----------------------------------------------
  Float flops = 0.0;
  Float matvec_flops = (2*1320+48)*GJP.VolNodeSites()/2;
  if (mat_type == 1) matvec_flops *= 2; // double flops since normal equations
  //----------------------------------------------

  //----------------------------------------------
  // Calculate Stop condition
  //----------------------------------------------
  Float in_norm2 = in->NormSqGlbSum(f_size_cb);
  Float stop = dirac_arg->stop_rsd * dirac_arg->stop_rsd * in_norm2;
  
  int total_iter = 0, k = 0;
  
  // Initial residual
  if (mat_type == 0)
  {
    MatPc(r,out);
  }
  else
  {
    MatPcDagMatPc(r,out);
  }

  r->FTimesV1MinusV2(1.0,in,r,f_size_cb);
  Float r2 = r->NormSqGlbSum(f_size_cb);
  flops += 4*f_size_cb + matvec_flops;

  VRB.Flow(cname, fname, "0 iterations, res^2 = %1.15e, restart = 0\n", r2);

  while (r2 > stop && k < QudaParam.max_restart) {
    inv_param.tol = dirac_arg->stop_rsd;
    if(sqrt(stop/r2)>inv_param.tol) 
    {
      inv_param.tol = sqrt(stop/r2);
    }

    x->VecZero(f_size_cb);
    //---------------------------------
    //  Inversion sequence start
    //---------------------------------
    invertQuda(x, r, &inv_param);
    
    // Update solution
    out->VecAddEquVec(x, f_size_cb);
    
    //------------------------------------
    // Calculate new residual
    if (mat_type == 0)
      MatPc(r, out);
    else
      MatPcDagMatPc(r, out);

    r->FTimesV1MinusV2(1.0,in,r,f_size_cb);
    r2 = r->NormSqGlbSum(f_size_cb);
    //------------------------------------

    k++;
    total_iter += inv_param.iter + 1;
    flops += 1e9*inv_param.gflops + 8*f_size_cb + matvec_flops;

    VRB.Flow(cname, fname, "Gflops = %e, Seconds = %e, Gflops/s = %f\n", 
             inv_param.gflops, inv_param.secs, inv_param.gflops / inv_param.secs);

    VRB.Flow(cname, fname, "True |res| / |src| = %1.15e, iter = %d, restart = %d\n",  
             sqrt(r2)/sqrt(in_norm2), total_iter, k);
  }
  
  gettimeofday(&end,NULL);
  print_flops(cname,fname,flops,&start,&end);
    
  VRB.Flow(cname, fname, "Cuda Space Required. Spinor:%f + Gauge:%f GiB\n", 
	   inv_param.spinorGiB, gauge_param.gaugeGiB);

  VRB.Flow(cname, fname, "True |res| / |src| = %1.15e, iter = %d, restart = %d\n", 
	     sqrt(r2)/sqrt(in_norm2), total_iter, k);

  if (true_res) *true_res = sqrt(r2);

  //----------------------------------------
  //  Finalize QUDA memory and API
  //----------------------------------------
  freeGaugeQuda();
  //----------------------------------------

  sfree(x);
  sfree(r);

  //VRB.DeactivateLevel(VERBOSE_FLOW_LEVEL);
  return total_iter;
}

CPS_END_NAMESPACE
#endif
