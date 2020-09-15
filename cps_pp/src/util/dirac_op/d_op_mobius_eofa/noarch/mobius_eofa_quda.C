#include <config.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>

CPS_START_NAMESPACE
//----------------------------------------------------------------
//
// mobius_eofa_quda.C
//
// Interface to QUDA mobius Domain Wall inverter for EOFA
// 
//------------------------------------------------------------------
CPS_END_NAMESPACE

#include <comms/glb.h>
#include <util/dirac_op.h>
#include <util/dwf.h>
#include <util/error.h>
#include <util/gjp.h>
#include <util/lattice.h>
#include <util/site.h>
#include <util/time_cps.h>
#include <util/verbose.h>
#include <util/wilson.h>

#ifdef USE_QUDA

#include <invert_quda.h>
#include <util/dirac_op/cps_quda.h>

//#define QUDA_NEW_INTERFACE

CPS_START_NAMESPACE

CPSQuda::QudaParams DiracOpMobiusEOFA::GetQudaParams(int mat_type)
{
  const char* fname = "GetQudaParams(i)";

  CPSQuda::QudaParams qp = { newQudaGaugeParam(), newQudaInvertParam() };

  CPSQuda::ParamSetup_mdwf(QudaParam, qp.gauge_param, qp.inv_param, mat_type);

  qp.inv_param.maxiter = dirac_arg->max_num_iter;
  qp.inv_param.mass    = dirac_arg->mass;
  switch(dirac_arg->Inverter){
  case BICGSTAB:
    ERR.General(cname, fname, "BiCGSTAB inverter not implemented\n");
    break;
  case CG:
  default:
    qp.inv_param.inv_type = QUDA_CG_INVERTER;
    qp.inv_param.solve_type = QUDA_NORMOP_PC_SOLVE;
    break;
  }

#ifndef USE_QUDA_LANCZOS
  // Additional EOFA-specific setup courtesy of J. Tu
  qp.inv_param.dslash_type                    = QUDA_MOBIUS_DWF_EOFA_DSLASH;
  qp.inv_param.mass                           = this->m1;
  qp.inv_param.mq1                            = this->m1;
  qp.inv_param.mq2                            = this->m2;
  qp.inv_param.mq3                            = this->m3;
  qp.inv_param.eofa_shift                     = this->a;
  qp.inv_param.eofa_pm                        = (this->pm == 1) ? 1 : 0;
  qp.inv_param.kappa                          = 1.0 / ( 2.0 * ( 4.0 + this->m1 ) );
  qp.inv_param.solver_normalization           = QUDA_DEFAULT_NORMALIZATION;
  qp.inv_param.preserve_source                = QUDA_PRESERVE_SOURCE_YES;
  qp.inv_param.use_sloppy_partial_accumulator = 0;
  qp.inv_param.solution_accumulator_pipeline  = 1;
  qp.inv_param.max_res_increase               = 20000;
  qp.inv_param.tol                            = this->dirac_arg->stop_rsd;
  qp.inv_param.tol_restart                    = 1.0e-03;
  qp.inv_param.maxiter                        = this->dirac_arg->max_num_iter;
  qp.inv_param.solution_type                  = QUDA_MAT_SOLUTION;
#else
  ERR.General(cname,fname,"You may be using QUDA branch without EOFA");
#endif


  return qp;
}

void DiracOpMobiusEOFA::Mat(Vector* out, Vector* in)
{
  const char* fname = "Mat(V*,V*)";

  CPSQuda::QudaParams qp = DiracOpMobiusEOFA::GetQudaParams(1);

#ifndef USE_QUDA_LANCZOS
#if 0
  dslashQuda_mobius_eofa(out, in, &qp.inv_param, QUDA_ODD_PARITY, 2);
#else
  dslashQuda(out, in, &qp.inv_param, QUDA_ODD_PARITY);
#endif
#endif
}

void DiracOpMobiusEOFA::MatHerm(Vector* out, Vector* in)
{
  const char* fname = "MatHerm(V*,V*)";

  const size_t f_size = GJP.VolNodeSites() * lat.FsiteSize() / ( lat.FchkbEvl() + 1 );
  
  Vector* tmp = static_cast<Vector*>( smalloc(cname, fname, "tmp", f_size*sizeof(Float)) );

  CPSQuda::QudaParams qp = DiracOpMobiusEOFA::GetQudaParams(1);

  loadGaugeQuda(static_cast<void*>(gauge_field), &qp.gauge_param);

  // Check that data has been re-ordered for QUDA by constructor
  if(lat.StrOrd() != DWF_4D_EOPREC_EE){ 
    ERR.General(cname, fname, "DiracOpMobiusEOFA expects DWF_4D_EOPREC_EE order for QUDA inverter\n");
  }

#ifndef USE_QUDA_LANCZOS
//  dslashQuda_mobius_eofa(tmp, in, &qp.inv_param, QUDA_ODD_PARITY, 2);
  dslashQuda(out, in, &qp.inv_param, QUDA_ODD_PARITY);
#endif
  lat.Convert(CANONICAL, tmp);
  lat.g5R5(out, tmp);
  lat.Convert(DWF_4D_EOPREC_EE, out);

  sfree(cname, fname, "tmp", tmp);
}

int DiracOpMobiusEOFA::QudaInvert(Vector* out, Vector* in, Float* true_res, int mat_type)
{
  const char* fname = "QudaInvert(V*,V*,F*,int)";
  VRB.Func(cname, fname);

  struct timeval start, end;
  gettimeofday(&start, nullptr);

  CPSQuda::QudaParams qp = DiracOpMobiusEOFA::GetQudaParams(mat_type);

  loadGaugeQuda(static_cast<void*>(gauge_field), &qp.gauge_param);

  // Check that data has been re-ordered for QUDA by constructor
  if(lat.StrOrd() != DWF_4D_EOPREC_EE){ 
    ERR.General(cname, fname, "DiracOpMobiusEOFA expects DWF_4D_EOPREC_EE order for QUDA inverter\n");
  }
  
  // Do the inversion
  invertQuda(out, in, &qp.inv_param);
  
  // Normalize the solution to match CPS conventions
  const size_t f_size = GJP.VolNodeSites() * lat.FsiteSize() / ( lat.FchkbEvl() + 1 );
  out -> VecTimesEquFloat(1.0 / ( GJP.Mobius_b() * ( 4.0 - GJP.DwfHeight() ) + 1.0 ), f_size);

  //  Finalize QUDA memory and API
  freeGaugeQuda();
  
  VRB.FuncEnd(cname, fname);

  return qp.inv_param.iter;
}

CPS_END_NAMESPACE

#endif
