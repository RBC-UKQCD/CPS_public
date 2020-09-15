#ifndef _TWOPOINT_FUNCTION_GENERIC_H
#define _TWOPOINT_FUNCTION_GENERIC_H

#include "prop_sitematrix_getter.h"
#include "spin_flav_op.h"

CPS_START_NAMESPACE

//Assume form
// coeff * \sum_y e^{i(p_psibar + p_psi)y} Tr{  ( prop1(y,-p_psi;tsrc) )^dag SinkOp prop2(y,p_psibar;tsrc) SrcOp }
//use_opposite_sink_mom optionally flips the sign of the sink momentum to the 'wrong' value - used in testing the flavor projection in the paper
//User is required to provide the propagators with the right momentum
template<typename MatrixType>
void twoPointFunctionGeneric(fMatrix<Rcomplex> &into, const int tsrc, const Complex &coeff,
			     const SrcSnkOp<MatrixType> &sink_op, const SrcSnkOp<MatrixType> &src_op,
			     const ThreeMomentum &p_psibar, const ThreeMomentum &p_psi,
			     const PropSiteMatrixGetter &prop_dag, const PropSiteMatrixGetter &prop_undag,
			     const PropSplane splane = SPLANE_BOUNDARY,
			     bool use_opposite_sink_mom = false){
  ThreeMomentum p_tot_src = p_psibar + p_psi;
  ThreeMomentum p_tot_snk = -p_tot_src; //mom_phase computes exp(-p.x)
  if(use_opposite_sink_mom) p_tot_snk = -p_tot_snk;
  
  //if(!UniqueID()) printf("Computing 2pt LW with src momentum %s and snk momentum %s\n",p_tot_src.str().c_str(),p_tot_snk.str().c_str());

  const int Lt = GJP.TnodeSites()*GJP.Tnodes();
  //#define TWOPT_TEST

#ifndef TWOPT_TEST
  const int nthread = omp_get_max_threads();
  basicComplexArray<Rcomplex> tmp(Lt,nthread); //defaults to zero for all elements
#else
  basicComplexArray<Rcomplex> tmp(Lt,1); 
#endif

  int vol3d = GJP.VolNodeSites()/GJP.TnodeSites();

#pragma omp parallel for
  for(int x=0;x<GJP.VolNodeSites();x++){
    int pos[4];
    int rem = x;
    for(int i=0;i<4;i++){ pos[i] = rem % GJP.NodeSites(i); rem /= GJP.NodeSites(i); }

    int x3d_lcl = x % vol3d;
    int t_glb = pos[3] + GJP.TnodeCoor() * GJP.TnodeSites();
    int tdis_glb = (t_glb - tsrc + Lt) % Lt; //t_glb = 0 .. Lt-1 -> tdis_glb = Lt-tsrc .. Lt-1, 0 .. Lt-tsrc-1

    //Actually getting the prop with the chosen tdis_glb depends on the periodicity of the propagator in question: cf propwrapper.h
    MatrixType prop1_site;
    prop_dag.siteMatrix(prop1_site,x3d_lcl,tdis_glb,splane);
    prop1_site.hconj();
    sink_op.rightMultiply(prop1_site);
    
    MatrixType prop2_site;
    prop_undag.siteMatrix(prop2_site,x3d_lcl,tdis_glb,splane);
    src_op.rightMultiply(prop2_site);

    std::complex<double> phase = coeff * mom_phase(p_tot_snk, pos);

#ifdef TWOPT_TEST
# pragma omp critical
    {
      tmp[tdis_glb] += phase * Trace(prop1_site, prop2_site);
    }
#else
    tmp(tdis_glb, omp_get_thread_num()) += phase * Trace(prop1_site, prop2_site);
#endif
  }
#ifndef TWOPT_TEST
  tmp.threadSum();
#endif
  tmp.nodeSum();

  for(int tdis=0;tdis<Lt;tdis++)
    into(tsrc, tdis) = tmp[tdis];
}

template<typename MatrixType>
void twoPointFunctionWallSinkGeneric(fMatrix<Rcomplex> &into, const int tsrc, const Complex &coeff,
				     const SrcSnkOp<MatrixType> &sink_op, const SrcSnkOp<MatrixType> &src_op,
				     const WallSinkPropSiteMatrixGetter<MatrixType> &prop1W, const WallSinkPropSiteMatrixGetter<MatrixType> &prop2W){
  const int Lt = GJP.TnodeSites()*GJP.Tnodes();
  basicComplexArray<Rcomplex> tmp(Lt,1); 

  //WallSinkProp are available for all times on every node, so no need to nodeSum

#pragma omp_parallel for
  for(int t_dis=0;t_dis<Lt;t_dis++){
    MatrixType prop1_t;
    prop1W.siteMatrix(prop1_t,t_dis);
    prop1_t.hconj();
    sink_op.rightMultiply(prop1_t);

    MatrixType prop2_t;
    prop2W.siteMatrix(prop2_t,t_dis);
    src_op.rightMultiply(prop2_t);

    tmp[t_dis] = coeff * Trace(prop1_t, prop2_t);
  }

  for(int tdis=0;tdis<Lt;tdis++)
    into(tsrc, tdis) = tmp[tdis];
}


CPS_END_NAMESPACE

#endif
