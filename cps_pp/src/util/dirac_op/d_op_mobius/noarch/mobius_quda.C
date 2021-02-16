#include <config.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <string>

//----------------------------------------------------------------
//
// mobius_quda.C
//
// Interface to QUDA mobius Domain Wall inverter
// 
//------------------------------------------------------------------
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
//#define USE_QUDA_SPLIT_GRID
#ifdef USE_QUDA_SPLIT_GRID
#include <split_grid.h>
#endif
CPS_START_NAMESPACE 
int DiracOpMobius::QudaInvert (std::vector < Vector * >out,
			   std::vector < Vector * >in,
			   std::vector < Float > &true_res,
			   std::vector < int >&iters, int mat_type)
{

  const char *fname = "QudaInvert(v<V*>,v<V*>,v<F>,v<int>,i)";
  VRB.Func (cname, fname);
  int len = out.size ();
  assert (in.size () == len);
  iters.resize (len, 0);
  true_res.resize (len, 0.);

  struct timeval start, end;
  gettimeofday (&start, NULL);


  QudaGaugeParam gauge_param = newQudaGaugeParam ();
  QudaInvertParam inv_param = newQudaInvertParam ();
  CPSQuda::ParamSetup_mdwf (QudaParam, gauge_param, inv_param, mat_type);

#ifdef USE_QUDA_MADWF_ML
  inv_param.Ls_cheap = QudaParam.Ls_cheap;
  inv_param.perform_mspcg_madwf_ml_training =
    QudaParam.perform_mspcg_madwf_ml_training;
  if (inv_param.perform_mspcg_madwf_ml_training)
    VRB.Result (cname, fname,
		"inv_param.perform_mspcg_madwf_ml_training=true\n");
  inv_param.use_mspcg_madwf_ml_training =
    QudaParam.use_mspcg_madwf_ml_training;
  if (inv_param.use_mspcg_madwf_ml_training)
    VRB.Result (cname, fname, "inv_param.use_mspcg_madwf_ml_training=true\n");
  inv_param.maxiter_precondition = QudaParam.maxiter_precondition;
  inv_param.mu = QudaParam.mu;
#endif

  inv_param.maxiter = dirac_arg->max_num_iter;
  inv_param.mass = dirac_arg->mass;
  switch (dirac_arg->Inverter)
    {
    case BICGSTAB:
      inv_param.inv_type = QUDA_BICGSTAB_INVERTER;
      inv_param.solve_type = QUDA_DIRECT_PC_SOLVE;
      break;
    case CG:
    case CG_FIXED_ITER:
    default:
      inv_param.inv_type = QUDA_CG_INVERTER;
      inv_param.solve_type = QUDA_NORMOP_PC_SOLVE;
#ifdef USE_QUDA_MSPCG
      if ((GJP.ZMobius_PC_Type () == ZMOB_PC_SYM1) && (mat_type == 1))
#ifndef QUDA_NEW_INTERFACE
	inv_param.inv_type = QUDA_MSPCG_INVERTER;
#else
	inv_param.inv_type = QUDA_PCG_INVERTER;
      inv_param.inv_type_precondition = QUDA_CG_INVERTER;
      inv_param.schwarz_type = QUDA_ADDITIVE_SCHWARZ;
#endif
#endif
      break;
    }

  loadGaugeQuda ((void *) gauge_field, &gauge_param);

//  size_t f_size_cb = (size_t) GJP.VolNodeSites () * lat.FsiteSize () / 2;
//  assert (lat.half_size == f_size_cb);
  size_t f_size_cb = lat.half_size;

  std::vector < Vector * >x, r;
  for (int i = 0; i < len; i++)
    {
      x.push_back ((Vector *) smalloc (f_size_cb * sizeof (Float)));
      r.push_back ((Vector *) smalloc (f_size_cb * sizeof (Float)));
      x[i]->VecZero (f_size_cb);
      r[i]->VecZero (f_size_cb);
    }

  //----------------------------------------------
  //  Calculate Flops value 
  //----------------------------------------------
  Float flops = 0.0;
  Float matvec_flops = (2 * 1320 + 48) * GJP.VolNodeSites () / 2;
  if (mat_type == 1)
    matvec_flops *= 2;		// double flops since normal equations
  //----------------------------------------------

  //----------------------------------------------
  // Calculate Stop condition
  //----------------------------------------------
  std::vector < Float > in_norm2, r2, stop;
  for (int i = 0; i < len; i++) {
      in_norm2.push_back (in[i]->NormSqGlbSum (f_size_cb));
      VRB.Flow (cname, fname, "src %d: in_norm2=%e\n", in_norm2[i]);
      stop.push_back (dirac_arg->stop_rsd * dirac_arg->stop_rsd *
		      in_norm2[i]);
    }

  int total_iter = 0, k = 0;
  int MatDMat_test = 0;
  int cg_test = 1;

  //-----------------------------------temp-------------------------------------------

  Float total_time = -dclock ();
  Float inv_time = 0.;
  Float qudamat_time = 0.;
  if (cg_test == 1) {
      qudamat_time -= dclock ();
      // Initial residual
      for (int i = 0; i < len; i++)
	if (mat_type == 0) {
	    MatQuda (r[i], out[i], &inv_param);
	  }
	else {
	    MatDagMatQuda (r[i], out[i], &inv_param);
	  }
      qudamat_time += dclock ();
      for (int i = 0; i < len; i++) {
	  r[i]->FTimesV1MinusV2 (1.0, in[i], r[i], f_size_cb);
	  r2.push_back (r[i]->NormSqGlbSum (f_size_cb));
	}

      flops += 4 * f_size_cb + matvec_flops;

      for (int i = 0; i < len; i++)
	VRB.Result (cname, fname,
		    "src %d: 0 iterations, res^2 = %1.15e, restart = 0\n", i,
		    r2[i]);

      int if_continue = 0;
      for (int i = 0; i < len; i++)
	if (r2[i] > stop[i])
	  if_continue = 1;
      if (k > QudaParam.max_restart)
	if_continue = 0;
      if (total_iter >= dirac_arg->max_num_iter * len)
	if_continue = 0;

      while (if_continue) {

	  if ((k > 0) && dirac_arg->Inverter == CG_LOWMODE_DEFL)
	    InvLowModeApprox (x, r, dirac_arg->fname_eigen, dirac_arg->neig);
	  else
	    for (int i = 0; i < len; i++)
	      x[i]->VecZero (f_size_cb);
	  //---------------------------------
	  //  Inversion sequence start
	  //---------------------------------
	  inv_time -= dclock ();

//         inv_param.verbosity = QUDA_DEBUG_VERBOSE;
#ifdef USE_QUDA_SPLIT_GRID
          int Nsplit=1;
	  for (int i = 0; i < 4; i++) 
            Nsplit *= QudaParam.split_grid[i];
          if ((Nsplit>1)&&(len%Nsplit==0)){
            // Run split grid
	  for (int i = 0; i < 4; i++) 
            inv_param.split_grid[i] = QudaParam.split_grid[i];
            inv_param.num_src = Nsplit;

	  inv_param.tol = 1e10;
	  for (int i = 0; i < len; i++) {
	      VRB.Result (cname, fname,
			  " src %d: r2 stop = %0.7e %0.7e, k = %d max_restart=%d \n",
			  i, r2[i], stop[i], k, QudaParam.max_restart);
	      if ((sqrt (stop[i] / r2[i])*0.99) < inv_param.tol) {
		  inv_param.tol = sqrt (stop[i] / r2[i]) * 0.9;
	      }
          }
	  if (inv_param.tol < dirac_arg->inner_rsd)
	  inv_param.tol = dirac_arg->inner_rsd;
	  VRB.Result (cname, fname,"Nsplit=%d inv_param.tol=%e\n",Nsplit,inv_param.tol);
//            std::vector < void * >xsplit(Nsplit), rsplit(Nsplit);
          void **xsplit = (void**)malloc(Nsplit*sizeof(void*));
          void **rsplit = (void**)malloc(Nsplit*sizeof(void*));
          for(int i=0 ;i<len;i+=Nsplit){
            for(int j=0 ;j<Nsplit;j++){
              xsplit[j]=x[i+j];
              rsplit[j]=r[i+j];
            }
	  VRB.Result (cname, fname,"i=%d Nsplit=%d inv_param.tol=%e\n",i,Nsplit,inv_param.tol);
            invertSplitGridQuda(xsplit,rsplit, &inv_param, (void *)gauge_field, &gauge_param);
	      VRB.Result (cname, fname,
			  "Total Gflops = %e, Seconds = %e, Gflops/s = %f\n",
			  inv_param.gflops, inv_param.secs,
			  inv_param.gflops / inv_param.secs);
            }
          free(xsplit);free(rsplit);
          } else 
#endif
{
	  for (int i = 0; i < len; i++) {
	      VRB.Result (cname, fname,
			  " src %d: r2 stop = %0.7e %0.7e, k = %d max_restart=%d \n",
			  i, r2[i], stop[i], k, QudaParam.max_restart);
	      inv_param.tol = dirac_arg->stop_rsd;
	      if (sqrt (stop[i] / r2[i]) > inv_param.tol) {
		  inv_param.tol = sqrt (stop[i] / r2[i]) * 0.99;
		}
	      if (inv_param.tol < dirac_arg->inner_rsd)
		inv_param.tol = dirac_arg->inner_rsd;


	      invertQuda (x[i], r[i], &inv_param);
	      VRB.Result (cname, fname,
			  "Total Gflops = %e, Seconds = %e, Gflops/s = %f\n",
			  inv_param.gflops, inv_param.secs,
			  inv_param.gflops / inv_param.secs);

	      iters[i] += inv_param.iter + 1;
	      total_iter += inv_param.iter + 1;
	  }
}
	  VRB.Debug (cname, fname, "invertQuda() done\n");
	  inv_time += dclock ();

	  // Update solution
	  for (int i = 0; i < len; i++)
	    out[i]->VecAddEquVec (x[i], f_size_cb);

	  qudamat_time -= dclock ();
	  //------------------------------------
	  // Calculate new residual
	  for (int i = 0; i < len; i++)
	    {
	      if (mat_type == 0)
		MatQuda (r[i], out[i], &inv_param);
	      else
		MatDagMatQuda (r[i], out[i], &inv_param);
//                MatPcDagMatPc (r[i], out[i]);
	    }
	  qudamat_time += dclock ();

	  for (int i = 0; i < len; i++)
	    {
	      r[i]->FTimesV1MinusV2 (1.0, in[i], r[i], f_size_cb);
	      r2[i] = r[i]->NormSqGlbSum (f_size_cb);
	    }
	  //------------------------------------
	  k++;
	  //dividing by the number of nodes
	  flops +=
	    1e9 * inv_param.gflops / (Float) GJP.TotalNodes () +
	    8 * f_size_cb + matvec_flops * len;

	  for (int i = 0; i < len; i++)
	    VRB.Result (cname, fname,
			"True |res| / |src| = %1.15e, iter = %d, restart = %d\n",
			sqrt (r2[i]) / sqrt (in_norm2[i]), iters[i], k);

	  if_continue = 0;
	  for (int i = 0; i < len; i++)
	    if (r2[i] > stop[i])
	      if_continue = 1;
	  if (k > QudaParam.max_restart)
	    if_continue = 0;
	  if (total_iter >= dirac_arg->max_num_iter * len)
	    if_continue = 0;
	  VRB.Result (cname, fname, "k=%d total_iter=%d if_continue=%d\n", k,
		      total_iter, if_continue);
	}

      if (total_iter >= (len * dirac_arg->max_num_iter))
	{
	  if (dirac_arg->Inverter == CG)
	    {
	      VRB.Result (cname, fname,
			  "Not converged! Trying CG one last time\n");
	      inv_time -= dclock ();
	      inv_param.inv_type = QUDA_CG_INVERTER;
	      for (int i = 0; i < len; i++)
		invertQuda (x[i], r[i], &inv_param);
	      VRB.Debug (cname, fname, "invertQuda() done\n");
	      inv_time += dclock ();
	      if ((inv_param.iter + 1) >= dirac_arg->max_num_iter)
		ERR.General (cname, fname, "CG not converged\n");
	    }
	  else
	    {
// CG_FIXED_ITER or CG_LOWMMODE_DEFL    
	      VRB.Result (cname, fname,
			  "Not converged! hope it is intentional?\n");
	      VRB.Warn (cname, fname,
			"Not converged! hope it is intentional?\n");
	    }
	}

      gettimeofday (&end, NULL);
      print_flops (cname, fname, flops, &start, &end);
//    VRB.Result(cname, fname, "Cuda Space Required. Spinor:%f + Gauge:%f GiB\n", 
//        inv_param.spinorGiB, gauge_param.gaugeGiB);
      for (int i = 0; i < len; i++)
	{
	  VRB.Result (cname, fname,
		      "True |res| / |src| = %1.15e, iter = %d, restart = %d\n",
		      sqrt (r2[i]) / sqrt (in_norm2[i]), iters[i], k);
//    if (true_res)
	  true_res[i] = sqrt (r2[i]);
	}

    }
  total_time += dclock ();
  VRB.Result (cname, fname, "inv_time=%g qudamat_time=%g total_time=%g\n",
	      inv_time, qudamat_time, total_time);
  //----------------------------------------
  //  Finalize QUDA memory and API
  //----------------------------------------
  freeGaugeQuda ();
  //----------------------------------------
  for (int i = 0; i < len; i++)
    {
      if (x[i])
	{
	  sfree (x[i]);
	  x[i] = NULL;
	}
      if (r[i])
	{
	  sfree (r[i]);
	  r[i] = NULL;
	}
    }
  VRB.FuncEnd (cname, fname);

  return total_iter;
}

int
DiracOpMobius::QudaInvert (Vector * out, Vector * in, Float * true_res,
			   int mat_type)
{
  const char *fname = "QudaInvert(V*,V*,F*,int)";
  VRB.Func (cname, fname);


  struct timeval start, end;
  gettimeofday (&start, NULL);


  QudaGaugeParam gauge_param = newQudaGaugeParam ();
  QudaInvertParam inv_param = newQudaInvertParam ();
  CPSQuda::ParamSetup_mdwf (QudaParam, gauge_param, inv_param, mat_type);

#ifdef USE_QUDA_MADWF_ML
  inv_param.Ls_cheap = QudaParam.Ls_cheap;
  inv_param.perform_mspcg_madwf_ml_training =
    QudaParam.perform_mspcg_madwf_ml_training;
  if (inv_param.perform_mspcg_madwf_ml_training)
    VRB.Result (cname, fname,
		"inv_param.perform_mspcg_madwf_ml_training=true\n");
  inv_param.use_mspcg_madwf_ml_training =
    QudaParam.use_mspcg_madwf_ml_training;
  if (inv_param.use_mspcg_madwf_ml_training)
    VRB.Result (cname, fname, "inv_param.use_mspcg_madwf_ml_training=true\n");
  inv_param.maxiter_precondition = QudaParam.maxiter_precondition;
  inv_param.mu = QudaParam.mu;
#endif

  inv_param.maxiter = dirac_arg->max_num_iter;
  inv_param.mass = dirac_arg->mass;
  switch (dirac_arg->Inverter)
    {
    case BICGSTAB:
      inv_param.inv_type = QUDA_BICGSTAB_INVERTER;
      inv_param.solve_type = QUDA_DIRECT_PC_SOLVE;
      break;
    case CG:
    default:
      inv_param.inv_type = QUDA_CG_INVERTER;
//    inv_param.solve_type = QUDA_NORMEQ_PC_SOLVE;
      inv_param.solve_type = QUDA_NORMOP_PC_SOLVE;
#ifdef USE_QUDA_MSPCG
      if ((GJP.ZMobius_PC_Type () == ZMOB_PC_SYM1) && (mat_type == 1))
#ifndef QUDA_NEW_INTERFACE
	inv_param.inv_type = QUDA_MSPCG_INVERTER;
#else
	inv_param.inv_type = QUDA_PCG_INVERTER;
      inv_param.inv_type_precondition = QUDA_CG_INVERTER;
      inv_param.schwarz_type = QUDA_ADDITIVE_SCHWARZ;
#endif
#endif
      break;
    }

  loadGaugeQuda ((void *) gauge_field, &gauge_param);

  // F_size_cb is defined by "GJP.VolNodeSites()/2*lat.FsiteSize()".
  // And 'lat.FsiteSize()' has the information of 5th dimension.
  // So, for the total lattice volume, we don`t need inv_param.Ls.
//  size_t f_size_cb = (size_t) GJP.VolNodeSites () * lat.FsiteSize () / 2;
//  assert (lat.half_size == f_size_cb);
  size_t f_size_cb = lat.half_size;

  Vector *x = (Vector *) smalloc (f_size_cb * sizeof (Float));
  Vector *r = (Vector *) smalloc (f_size_cb * sizeof (Float));
  x->VecZero (f_size_cb);
  r->VecZero (f_size_cb);

  //----------------------------------------------
  //  Calculate Flops value 
  //----------------------------------------------
  Float flops = 0.0;
  Float matvec_flops = (2 * 1320 + 48) * GJP.VolNodeSites () / 2;
  if (mat_type == 1)
    matvec_flops *= 2;		// double flops since normal equations
  //----------------------------------------------

  //----------------------------------------------
  // Calculate Stop condition
  //----------------------------------------------
  Float in_norm2 = in->NormSqGlbSum (f_size_cb);
  Float stop = dirac_arg->stop_rsd * dirac_arg->stop_rsd * in_norm2;

  int total_iter = 0, k = 0;
  int MatDMat_test = 0;
  int cg_test = 1;
//----------------------Debug code------------------------
  if (MatDMat_test)
    {
      Vector *in_tmp = (Vector *) smalloc (f_size_cb * sizeof (Float));
      in_tmp->VecZero (f_size_cb);
      Vector *r_tmp = (Vector *) smalloc (f_size_cb * sizeof (Float));
      r_tmp->VecZero (f_size_cb);
      Float *in_vec = (Float *) in_tmp;
      for (size_t tmp = 0; tmp < f_size_cb; tmp++)
	in_vec[tmp] = 0.0;
      in_vec[24] = 1.0;
      in_vec[1560] = 1.0;
      in_vec[3048] = 1.0;
      in_vec[3072] = 1.0;
      in_vec[4584] = 1.0;
      in_vec[4608] = 1.0;
      in_vec[6144] = 1.0;
      // in_vec[96] = 1.0; in_vec[384] = 1.0; in_vec[1536] = 1.0; in_vec[1512] = 1.0; 
      Float r2 = in_tmp->NormSqGlbSum (f_size_cb);
      VRB.Flow (cname, fname, "CPU input res^2 = %1.15e\n", r2);

      MatPc (r, in);
      MatQuda (r_tmp, in, &inv_param);
      Float *cpu_vec = (Float *) r;
      Float *gpu_vec = (Float *) r_tmp;

      int err_count = 0;
      int xid, yid, zid, tid, sid, spin, clr, idx;

      for (sid = 0; sid < GJP.SnodeSites (); sid++)
	for (tid = 0; tid < GJP.TnodeSites (); tid++)
	  for (zid = 0; zid < GJP.ZnodeSites (); zid++)
	    for (yid = 0; yid < GJP.YnodeSites (); yid++)
	      for (xid = 0; xid < GJP.XnodeSites (); xid++)
		{
		  if ((tid + zid + yid + xid) % 2 == 0)
		    {
		      for (spin = 0; spin < 4; spin++)
			for (clr = 0; clr < 3; clr++)
			  {
			    idx =
			      sid * GJP.TnodeSites () * GJP.ZnodeSites () *
			      GJP.YnodeSites () * GJP.XnodeSites () * 12 +
			      ((tid * GJP.ZnodeSites () * GJP.YnodeSites () *
				GJP.XnodeSites ()
				+ zid * GJP.YnodeSites () * GJP.XnodeSites ()
				+ yid * GJP.XnodeSites () + xid) / 2) * 24 +
			      6 * spin + 2 * clr;
			    Float diff1 = fabs (gpu_vec[idx] - cpu_vec[idx]);
			    if (diff1 > 10e-12)
			      {
				printf
				  ("(%d,%d,%d,%d,%d,s %d,c %d) real: GPU : %e\t CPU : %e\t GPU-CPU = %e\n",
				   sid, tid, zid, yid, xid, spin, clr,
				   gpu_vec[idx], cpu_vec[idx], diff1);
				err_count++;
				printf ("error!!!!!\n");
			      }
			    Float diff2 =
			      fabs (gpu_vec[idx + 1] - cpu_vec[idx + 1]);
			    if (diff2 > 10e-12)
			      {
				printf
				  ("(%d,%d,%d,%d,%d,s %d,c %d) img: GPU : %e\t CPU : %e\t GPU-CPU = %e\n",
				   sid, tid, zid, yid, xid, spin, clr,
				   gpu_vec[idx + 1], cpu_vec[idx + 1], diff2);
				err_count++;
				printf ("error!!!!!\n");
			      }
			  }
		    }
//      printf("\n");
		}
//    for(idx = 0; idx < f_size_cb; idx++)
//    {
//      Float diff1 = fabs(gpu_vec[idx]-cpu_vec[idx]);
//      printf("%d : GPU : %e\t CPU : %e\t GPU-CPU = %e\n",idx,gpu_vec[idx], cpu_vec[idx], diff1);
//      if(diff1 > 10e-12)
//        printf("error!!!!!\n");
//    }
      if (err_count)
	printf ("Node %d: error count = %d\n", UniqueID (), err_count);
      r2 = r_tmp->NormSqGlbSum (f_size_cb);
      VRB.Flow (cname, fname, "GPU out res^2 = %1.15e\n", r2);
      r2 = r->NormSqGlbSum (f_size_cb);
      VRB.Flow (cname, fname, "CPU out res^2 = %1.15e\n", r2);
      if (r_tmp != NULL)
	{
	  sfree (r_tmp);
	  r_tmp = NULL;
	}
      if (in_tmp != NULL)
	{
	  sfree (in_tmp);
	  in_tmp = NULL;
	}
    }
  //-----------------------------------temp-------------------------------------------

  Float total_time = -dclock ();
  Float inv_time = 0.;
  Float qudamat_time = 0.;
  if (cg_test == 1)
    {
      qudamat_time -= dclock ();
      // Initial residual
      if (mat_type == 0)
	{
	  MatQuda (r, out, &inv_param);
	}
      else
	{
	  MatDagMatQuda (r, out, &inv_param);
//      MatPcDagMatPc(r, out);
	}
      qudamat_time += dclock ();
      r->FTimesV1MinusV2 (1.0, in, r, f_size_cb);
      Float r2 = r->NormSqGlbSum (f_size_cb);

      flops += 4 * f_size_cb + matvec_flops;

      VRB.Result (cname, fname, "0 iterations, res^2 = %1.15e, restart = 0\n",
		  r2);

      while (r2 > stop && k < QudaParam.max_restart)
	{
	  VRB.Result (cname, fname,
		      " r2 stop = %0.7e %0.7e, k = %d max_restart=%d \n", r2,
		      stop, k, QudaParam.max_restart);
	  inv_param.tol = dirac_arg->stop_rsd;
	  if (sqrt (stop / r2) > inv_param.tol)
	    {
	      inv_param.tol = sqrt (stop / r2) * 0.5;
	    }

	  if ((k % 2 ) && dirac_arg->Inverter == CG_LOWMODE_DEFL)
	    InvLowModeApprox (x, r, dirac_arg->fname_eigen, dirac_arg->neig,
			      true_res);
	  else
	    x->VecZero (f_size_cb);
	  //---------------------------------
	  //  Inversion sequence start
	  //---------------------------------
	  inv_time -= dclock ();
	  invertQuda (x, r, &inv_param);
	  VRB.Debug (cname, fname, "invertQuda() done\n");
	  inv_time += dclock ();

	  // Update solution
	  out->VecAddEquVec (x, f_size_cb);

	  qudamat_time -= dclock ();
	  //------------------------------------
	  // Calculate new residual
	  if (mat_type == 0)
	    MatQuda (r, out, &inv_param);
	  else
	    MatDagMatQuda (r, out, &inv_param);
//        MatPcDagMatPc (r, out);
	  qudamat_time += dclock ();

	  r->FTimesV1MinusV2 (1.0, in, r, f_size_cb);
	  r2 = r->NormSqGlbSum (f_size_cb);
	  //------------------------------------
	  k++;
	  total_iter += inv_param.iter + 1;
	  //dividing by the number of nodes
	  flops +=
	    1e9 * inv_param.gflops / (Float) GJP.TotalNodes () +
	    8 * f_size_cb + matvec_flops;

	  VRB.Result (cname, fname,
		      "Total Gflops = %e, Seconds = %e, Gflops/s = %f\n",
		      inv_param.gflops, inv_param.secs,
		      inv_param.gflops / inv_param.secs);
	  VRB.Result (cname, fname,
		      "True |res| / |src| = %1.15e, iter = %d, restart = %d\n",
		      sqrt (r2) / sqrt (in_norm2), total_iter, k);
	}

      if (total_iter >= dirac_arg->max_num_iter)
	{
	  VRB.Result (cname, fname,
		      "Not converged! Trying CG one last time\n");
	  inv_time -= dclock ();
	  inv_param.inv_type = QUDA_CG_INVERTER;
	  invertQuda (x, r, &inv_param);
	  VRB.Debug (cname, fname, "invertQuda() done\n");
	  inv_time += dclock ();
	  if ((inv_param.iter + 1) >= dirac_arg->max_num_iter)
	    ERR.General (cname, fname, "Not converged!\n");
	}

      gettimeofday (&end, NULL);
      print_flops (cname, fname, flops, &start, &end);
//    VRB.Result(cname, fname, "Cuda Space Required. Spinor:%f + Gauge:%f GiB\n", 
//        inv_param.spinorGiB, gauge_param.gaugeGiB);
//    VRB.Result (cname, fname, "True |res| / |src| = %1.15e, iter = %d, restart = %d\n", sqrt (r2) / sqrt (in_norm2), total_iter, k);
      if (true_res)
	*true_res = sqrt (r2);
    }
  total_time += dclock ();
  VRB.Result (cname, fname, "inv_time=%g qudamat_time=%g total_time=%g\n",
	      inv_time, qudamat_time, total_time);
  //----------------------------------------
  //  Finalize QUDA memory and API
  //----------------------------------------
  freeGaugeQuda ();
  //----------------------------------------
  if (x)
    {
      sfree (x);
      x = NULL;
    }
  if (r)
    {
      sfree (r);
      r = NULL;
    }
  VRB.FuncEnd (cname, fname);

  return total_iter;
}

int
DiracOpMobius::QudaMInvCG (Vector ** out, Vector * in, Float in_norm,
			   Float * shift, int Nshift, int isz,
			   Float * RsdCG, MultiShiftSolveType type,
			   Float * alpha)
{
  const char *fname = "QudaMInvCG(V*,V*,F*,int)";


  VRB.Result (cname, fname, "Started\n");
  struct timeval start, end;
  gettimeofday (&start, NULL);


  QudaGaugeParam gauge_param = newQudaGaugeParam ();
  QudaInvertParam inv_param = newQudaInvertParam ();
  CPSQuda::ParamSetup_mdwf (QudaParam, gauge_param, inv_param, 1);
// overwrite to single, as half appears to be really poor!
  inv_param.cuda_prec_sloppy = QUDA_HALF_PRECISION;

  inv_param.inv_type = QUDA_CG_INVERTER;
  inv_param.solution_type = QUDA_MATPCDAG_MATPC_SOLUTION;
  inv_param.solve_type = QUDA_NORMOP_PC_SOLVE;
// for master and/or old QUDA branches
//  inv_param.solve_type = QUDA_NORMEQ_PC_SOLVE;

  inv_param.maxiter = dirac_arg->max_num_iter;
  inv_param.mass = dirac_arg->mass;
  inv_param.num_offset = Nshift;

  loadGaugeQuda ((void *) gauge_field, &gauge_param);

  // F_size_cb is defined by "GJP.VolNodeSites()/2*lat.FsiteSize()".
  // And 'lat.FsiteSize()' has the information of 5th dimension.
  // So, for the total lattice volume, we don`t need inv_param.Ls.
  size_t f_size_cb = (size_t) GJP.VolNodeSites () * lat.FsiteSize () / 2;
  assert (lat.half_size == f_size_cb);

  Vector *x = (Vector *) smalloc (f_size_cb * sizeof (Float));
  Vector *r = (Vector *) smalloc (f_size_cb * sizeof (Float));
  x->VecZero (f_size_cb);
  r->VecZero (f_size_cb);
  Float r2 = r->NormSqGlbSum (f_size_cb);

  //----------------------------------------------
  //  Calculate Flops value 
  //----------------------------------------------
  Float flops = 0.0;
  Float matvec_flops = (2 * 1320 + 48) * GJP.VolNodeSites () / 2;
  matvec_flops *= 2;		// double flops since normal equations
  //----------------------------------------------

  //----------------------------------------------
  // Calculate Stop condition
  //----------------------------------------------
  Float in_norm2 = in->NormSqGlbSum (f_size_cb);
  Float stop = dirac_arg->stop_rsd * dirac_arg->stop_rsd * in_norm2;

  int total_iter = 0, k = 0;
//  int MatDMat_test = 0;
  int cg_test = 1;

  Float total_time = -dclock ();
  Float inv_time = 0.;
//  Float qudamat_time=0.;
//  if(cg_test==1) {
//    while (r2 > stop && k < QudaParam.max_restart) {
  inv_param.tol = dirac_arg->stop_rsd;
  if (sqrt (stop / r2) > inv_param.tol)
    {
      inv_param.tol = sqrt (stop / r2);
    }

  x->VecZero (f_size_cb);
  //---------------------------------
  //  Inversion sequence start
  //---------------------------------
  inv_time -= dclock ();
  void *OutMulti[inv_param.num_offset];
  for (int i = 0; i < inv_param.num_offset; i++)
    {
      if (type == MULTI)
	OutMulti[i] = static_cast < Vector * >(out[i]);
      else
	{
	  OutMulti[i] =
	    (Vector *) smalloc (cname, fname, "OutMulti[i]",
				f_size_cb * sizeof (Float));
	}
      inv_param.offset[i] = shift[i];
      inv_param.tol_offset[i] = RsdCG[i];
      VRB.Result (cname, fname, "%d: OutMulti %p offset %g tol_offset %g\n",
		  i, OutMulti[i], inv_param.offset[i],
		  inv_param.tol_offset[i]);
    }
  invertMultiShiftQuda (OutMulti, in, &inv_param);
  inv_time += dclock ();
  x->VecZero (f_size_cb);
  r->VecZero (f_size_cb);
  for (int i = 0; i < inv_param.num_offset; i++)
    {
      Float d;
#if 1
      MatDagMatQuda (r, (Vector *) OutMulti[i], &inv_param);
#else
      MatPcDagMatPc (r, (Vector *) OutMulti[i]);
#endif
      r->FTimesV1PlusV2 (shift[i], (Vector *) OutMulti[i], r, f_size_cb);
      r->FTimesV1PlusV2 (-1., in, r, f_size_cb);
      d = r->NormSqGlbSum (f_size_cb);
      Float *f_p = (Float *) OutMulti[i];
      VRB.Result (cname, fname,
		  "Pole %d : shift %g (%g %g) True |res| / |src| = %1.15e\n",
		  i, shift[i], *f_p, *(f_p + 1), sqrt (d) / sqrt (in_norm2));
    }

  total_iter += inv_param.iter + 1;
  flops += 1e9 * inv_param.gflops + 8 * f_size_cb + matvec_flops;

  VRB.Result (cname, fname, "Gflops = %e, Seconds = %e, Gflops/s = %f\n",
	      inv_param.gflops, inv_param.secs,
	      inv_param.gflops / inv_param.secs);
  VRB.Result (cname, fname,
	      "True |res| / |src| = %1.15e, iter = %d, restart = %d\n",
	      sqrt (r2) / sqrt (in_norm2), total_iter, k);
//    }
  if (type == SINGLE)
    {
//    out[0]->VecZero (f_size_cb); REALLY??
      for (int i = 0; i < inv_param.num_offset; i++)
	{
	  Vector *v_p = (Vector *) OutMulti[i];
	  out[0]->FTimesV1PlusV2 (alpha[i], v_p, out[0], f_size_cb);
	  sfree (cname, fname, "OutMulti[i]", OutMulti[i]);
	}

    }

  gettimeofday (&end, NULL);
  print_flops (cname, fname, flops, &start, &end);
//    VRB.Flow(cname, fname, "Cuda Space Required. Spinor:%f + Gauge:%f GiB\n", 
//        inv_param.spinorGiB, gauge_param.gaugeGiB);
  VRB.Flow (cname, fname,
	    "True |res| / |src| = %1.15e, iter = %d, restart = %d\n",
	    sqrt (r2) / sqrt (in_norm2), total_iter, k);
//     if (true_res) *true_res = sqrt(r2);
//  }
  total_time += dclock ();
  VRB.Result (cname, fname, "inv_time=%g total_time=%g\n", inv_time,
	      total_time);
  //----------------------------------------
  //  Finalize QUDA memory and API
  //----------------------------------------
  freeGaugeQuda ();
  //----------------------------------------
  if (x)
    {
      sfree (x);
      x = NULL;
    }
  if (r)
    {
      sfree (r);
      r = NULL;
    }
//      free(OutMulti);OutMulti=NULL;

  VRB.FuncEnd (cname, fname);
  return total_iter;
}



#if (defined USE_QUDA_LANCZOS ) || (defined QUDA_NEW_INTERFACE)
int
DiracOpMobius::QudaLanczos (Vector ** V,	//Lanczos vectors, eigenvectors of RitzMat on return
			    Float * alpha,	// eigenvalues
			    LanczosArg * lanczos_arg)
{
  const char *fname = "QudaLanczos(V*,V*,F*,int)";


  VRB.Func (cname, fname);
  struct timeval start, end;
  gettimeofday (&start, NULL);

  int nk = lanczos_arg->nk_lanczos_vectors;	// number of wanted +extra eigenvectors
  int nt = lanczos_arg->nt_lanczos_vectors;	// number of wanted eigenvectors
  int np = lanczos_arg->np_lanczos_vectors;	// extra, unwanted eigenvectors
  int MaxIters = lanczos_arg->maxiters;	// maximum number of restarting
  Float StopRes = lanczos_arg->stop_residual;	// stop when residual is smaller than this
  PrecType prec = lanczos_arg->precision;

  MatrixPolynomialArg *cheby_arg = &lanczos_arg->matpoly_arg;

  QudaGaugeParam gauge_param = newQudaGaugeParam ();
  QudaInvertParam inv_param = newQudaInvertParam ();
  QudaEigParam eig_param = newQudaEigParam ();

  int if_double = 1;
  if (prec == PREC_SINGLE)
    if_double = 0;
  CPSQuda::ParamSetup_mdwf (QudaParam, gauge_param, inv_param, 1, if_double);
  CPSQuda::ParamSetup_mdwf (QudaParam, gauge_param, eig_param, 1, if_double);
  eig_param.invert_param = &inv_param;



  inv_param.cuda_prec_sloppy = QUDA_SINGLE_PRECISION;
  inv_param.inv_type = QUDA_CG_INVERTER;
  inv_param.solution_type = QUDA_MATPCDAG_MATPC_SOLUTION;
  inv_param.solve_type = QUDA_NORMOP_PC_SOLVE;
  inv_param.tol = StopRes;
//  inv_param.verbosity = QUDA_DEBUG_VERBOSE;
  inv_param.verbosity = QUDA_VERBOSE;
  inv_param.gamma_basis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;
  inv_param.cuda_prec_precondition = QUDA_SINGLE_PRECISION;

// for master and/or old QUDA branches
//  eig_param.solve_type = QUDA_NORMEQ_PC_SOLVE;
  eig_param.spectrum = QUDA_SPECTRUM_SR_EIG;
  eig_param.use_poly_acc = QUDA_BOOLEAN_YES;
//  eig_param.use_poly_acc = QUDA_BOOLEAN_NO;
  eig_param.use_norm_op = QUDA_BOOLEAN_YES;
  eig_param.use_dagger = QUDA_BOOLEAN_NO;
  eig_param.poly_deg = cheby_arg->Npol + 1;
  eig_param.require_convergence = QUDA_BOOLEAN_NO;	// wants to return and at least write tunecache
  eig_param.tol = StopRes;
  eig_param.a_min = pow (cheby_arg->params.params_val[1], 2);
  eig_param.a_max = pow (cheby_arg->params.params_val[0], 2);
#ifndef QUDA_NEW_INTERFACE
  eig_param.eig_type = QUDA_EIG_TR_LANCZOS;
  eig_param.nEv = nk;
  eig_param.nKr = nk + np;
  eig_param.nConv = nt;
#else
  eig_param.eig_type = QUDA_EIG_TR_LANCZOS;
  int bsize = 1;
#if 0
  bsize = 16;
  while (nk % bsize)
    {
      bsize = nk % bsize;
      VRB.Result (cname, fname, "nk=%d bsize=%d\n", nk, bsize);
    }
  if (bsize < 1)
    bsize = 1;
  while (np % bsize)
    {
      bsize = np % bsize;
      VRB.Result (cname, fname, "np=%d bsize=%d\n", np, bsize);
    }
  if (bsize < 1)
    bsize = 1;
  if (bsize > 1)
    eig_param.eig_type = QUDA_EIG_BLK_TR_LANCZOS;
#endif
  eig_param.n_ev = nk;
  eig_param.n_kr = nk + np;
  if (bsize > 1)
    VRB.Result (cname, fname,
		"Running block lanczos: Nk=%d Np=%d block_size=%d\n", nk, np,
		bsize);
  eig_param.block_size = bsize;
  eig_param.n_conv = nt;
#endif
  eig_param.tol = StopRes;
  eig_param.check_interval = 1;
  eig_param.max_restarts = MaxIters;
  strcpy (eig_param.vec_infile, "");
  strcpy (eig_param.vec_outfile, "");
//  eig_param.eig_type = QUDA_LANCZOS;

  inv_param.maxiter = dirac_arg->max_num_iter;
  inv_param.mass = dirac_arg->mass;
  VRB.Result (cname, fname, "maxiter=%d mass=%e\n", inv_param.maxiter,
	      inv_param.mass);

  loadGaugeQuda ((void *) gauge_field, &gauge_param);

//  size_t f_size_cb = (size_t) GJP.VolNodeSites () * lat.FsiteSize () / 2;
//  assert (lat.half_size == f_size_cb);
  size_t f_size_cb = lat.half_size;

  void **host_evecs =
    (void **) smalloc (cname, fname, "host_evecs", nt * sizeof (void *));
  for (int i = 0; i < nt; i++)
    {
      host_evecs[i] = (void *) V[i];
    }
  double _Complex *host_evals =
    (double _Complex *) smalloc (nk * sizeof (double _Complex));
  VRB.Result (cname, fname, "before eigensolveQuda() done\n");
  VRB.Result (cname, fname, "lanczos_arg->conv_check=%d nt=%d\n",
	      lanczos_arg->conv_check, nt);
  eigensolveQuda (host_evecs, host_evals, &eig_param);
  VRB.Result (cname, fname, "eigensolveQuda() done\n");
//  exit(-42);
//  cps::sync();

  Vector *r = (Vector *) smalloc (f_size_cb * sizeof (Float));
  r->VecZero (f_size_cb);
  Vector *V_tmp = (Vector *) smalloc (f_size_cb * sizeof (Float));
  V_tmp->VecZero (f_size_cb);

  for (int i = 0; i < nt; i++)
    alpha[i] = creal (host_evals[i]);
  if (lanczos_arg->conv_check)
    {
      Vector *V_test = (Vector *) smalloc (f_size_cb * sizeof (Float));
      for (int i = 0; i < nt; i += lanczos_arg->conv_check)
	{
	  Vector *V_d = V[i];
	  VRB.Result (cname, fname, "eval[%d]=%e %e\n", i,
		      creal (host_evals[i]), cimag (host_evals[i]));
	  if (prec == PREC_SINGLE)
	    {
//        inv_param.cuda_prec = QUDA_SINGLE_PRECISION;
	      inv_param.cpu_prec = QUDA_DOUBLE_PRECISION;
	      Float *F_d = (Float *) V_test;
	      float *F_f = (float *) V[i];
#pragma omp parallel for
	      for (size_t j = 0; j < f_size_cb; j++)
		F_d[j] = F_f[j];
	      V_d = V_test;
	    }
#if 1
	  MatDagMatQuda (r, V_d, &inv_param);
#else
	  MatPcDagMatPc (r, V[i]);
#endif
	  printf ("Node %d: V_d\n", UniqueID ());
	  Float V2 = V_d->NormSqGlbSum (f_size_cb);
	  printf ("Node %d: r2\n", UniqueID ());
	  Float r2 = r->NormSqGlbSum (f_size_cb);
	  printf ("Node %d: r\n", UniqueID ());
	  r->FTimesV1MinusV2 (alpha[i], V_d, r, f_size_cb);
	  VRB.Result (cname, fname,
		      "norm(V) norm(r) norm(r-alpha.V) = %e %e %e\n", V2, r2,
		      r->NormSqGlbSum (f_size_cb));
	  exit (-42);
	}
      sfree (V_test);
    }

  //----------------------------------------------
  //  Calculate Flops value 
  //----------------------------------------------
  Float flops = 0.0;
  Float matvec_flops = (2 * 1320 + 48) * GJP.VolNodeSites () / 2;
  matvec_flops *= 2;		// double flops since normal equations
  //----------------------------------------------

  //----------------------------------------
  //  Finalize QUDA memory and API
  //----------------------------------------
  freeGaugeQuda ();
  //----------------------------------------
  if (r)
    {
      sfree (r);
      r = NULL;
    }
  if (V_tmp)
    {
      sfree (V_tmp);
      V_tmp = NULL;
    }

  sfree (host_evecs);
  sfree (host_evals);
  VRB.Result (cname, fname, "Done\n");
//  exit(-42);
  return 0;			//QUDA does not return the number of ops or iters currently
}
#endif

CPS_END_NAMESPACE
#endif
