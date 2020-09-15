#include <config.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string>

CPS_START_NAMESPACE
//----------------------------------------------------------------
//
// dwf_quda.C
//
// Interface to QUDA Domain Wall inverter
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
#include <util/dwf.h>


#ifdef USE_QUDA
#include <invert_quda.h>
#include <util/dirac_op/cps_quda.h>

#define QUDA_NEW_INTERFACE

CPS_START_NAMESPACE

int DiracOpDwf::QudaInvert(Vector* out, Vector* in, Float* true_res, int mat_type) 
{
  const char* fname = "QudaInvert(V*, V*, F*, int)";
  VRB.Func(cname, fname);

  // VRB.ActivateLevel(VERBOSE_FLOW_LEVEL);

  struct timeval start, end;
  gettimeofday(&start, NULL);

  QudaGaugeParam gauge_param = newQudaGaugeParam();
  QudaInvertParam inv_param = newQudaInvertParam();
  CPSQuda::ParamSetup_dwf(QudaParam, gauge_param, inv_param, mat_type);

  inv_param.maxiter = dirac_arg -> max_num_iter;
  inv_param.mass = dirac_arg -> mass;
  switch(dirac_arg -> Inverter)
  {
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

  loadGaugeQuda((void*) gauge_field, &gauge_param);

  // F_size_cb is defined by "GJP.VolNodeSites()/2*lat.FsiteSize()".
  // And 'lat.FsiteSize()' has the information of 5th dimension.
  // So, for the total lattice volume, we don`t need inv_param.Ls.
  size_t f_size_cb = static_cast<size_t>( GJP.VolNodeSites() * lat.FsiteSize() / 2 );

  Vector* x = (Vector*) smalloc(f_size_cb * sizeof(Float));
  Vector* r = (Vector*) smalloc(f_size_cb * sizeof(Float));
  x -> VecZero(f_size_cb);
  r -> VecZero(f_size_cb);
  
  //----------------------------------------------
  //  Calculate Flops value 
  //----------------------------------------------
  Float flops = 0.0;
  Float matvec_flops = (2*1320+48)*GJP.VolNodeSites()/2;
  if(mat_type == 1){ matvec_flops *= 2; }  // double flops since normal equations
  //----------------------------------------------

  //----------------------------------------------
  // Calculate Stop condition
  //----------------------------------------------
  Float in_norm2 = in -> NormSqGlbSum(f_size_cb);
  Float stop = (dirac_arg -> stop_rsd) * (dirac_arg -> stop_rsd) * in_norm2;
  
  int total_iter = 0, k = 0;
  
  // Initial residual
  if (mat_type == 0)
  {
    //MatPc(r, out); // CPU version of Dslash operator
    MatQuda(r, out, &inv_param); // GPU version of Dslash operator
  }
  else
  {
    //MatPcDagMatPc(r, out); // CPU version of Dslash operator
    MatDagMatQuda(r, out, &inv_param); // GPU version of Dslash operator
  }
  r -> FTimesV1MinusV2(1.0,in,r,f_size_cb);
  Float r2 = r -> NormSqGlbSum(f_size_cb);

  flops += 4*f_size_cb + matvec_flops;

  VRB.Flow(cname, fname, "0 iterations, res^2 = %1.15e, restart = 0 stop=%e max_restart=%d\n", r2, stop, QudaParam.max_restart);

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
    {
      //MatPc(r, out); // CPU version of Dslash operator
      MatQuda(r, out, &inv_param); // GPU version of Dslash operator
    }
    else
    {
      //MatPcDagMatPc(r, out); // CPU version of Dslash operator
      MatDagMatQuda(r, out, &inv_param); // GPU version of Dslash operator
    }

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
  CGflops += flops;
    
//  VRB.Flow(cname, fname, "Cuda Space Required. Spinor:%f + Gauge:%f GiB\n", inv_param.spinorGiB, gauge_param.gaugeGiB);

  VRB.Flow(cname, fname, "True |res| / |src| = %1.15e, iter = %d, restart = %d\n", 
	     sqrt(r2)/sqrt(in_norm2), total_iter, k);

  if (true_res) *true_res = sqrt(r2);

  //----------------------------------------
  //  Finalize QUDA memory and API
  //----------------------------------------
  freeGaugeQuda();
  //----------------------------------------

  if(x){ sfree(x); x = NULL; }
  if(r){ sfree(r); r = NULL; }
  VRB.FuncEnd(cname, fname);

  return total_iter;
}

CPS_END_NAMESPACE

#endif

