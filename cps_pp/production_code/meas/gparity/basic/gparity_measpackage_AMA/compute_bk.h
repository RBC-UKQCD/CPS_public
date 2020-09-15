#ifndef _COMPUTE_BK_AMA_H
#define _COMPUTE_BK_AMA_H

CPS_START_NAMESPACE

//We compute B_K for a given first kaon timeslice (t0) and a fixed K->K separation for each operator insertion time. The sink kaon timeslice is t1 = (t0 + tsep) % Lt
//The matrix is indexed as [t0][(top-t0+Lt)%Lt]
//Source momenta of the strange quark props are needed for the flavor projection.
//It is assumed that the total kaon momentum is zero, and we project onto zero momentum at the operator insertion
void GparityBK(fMatrix<Rcomplex> &into, const int t0, const int t1,
	       const PropSiteMatrixGetter &prop_h_t0, const PropSiteMatrixGetter &prop_l_t0, const ThreeMomentum &p_psi_h_t0,
	       const PropSiteMatrixGetter &prop_h_t1, const PropSiteMatrixGetter &prop_l_t1, const ThreeMomentum &p_psi_h_t1,
	       const bool do_flav_project = true
	       ){
  const int Lt = GJP.TnodeSites()*GJP.Tnodes();

  const int nthread = omp_get_max_threads();
  basicComplexArray<Rcomplex> tmp(Lt,nthread); //defaults to zero for all elements

  FlavorMatrix kaon_proj_t0 = getProjector(p_psi_h_t0);
  FlavorMatrix kaon_proj_t1 = getProjector(p_psi_h_t1);

  int vol3d = GJP.VolNodeSites()/GJP.TnodeSites();

#pragma omp parallel for
  for(int x=0;x<GJP.VolNodeSites();x++){
    int pos[4];
    int rem = x;
    for(int i=0;i<4;i++){ pos[i] = rem % GJP.NodeSites(i); rem /= GJP.NodeSites(i); }

    int x3d_lcl = x % vol3d;
    int t_glb = pos[3] + GJP.TnodeCoor() * GJP.TnodeSites(); //operator insertion time
    int tdis0_glb = t_glb - t0; //linear time coordinate
    int tdis1_glb = t_glb - t1;
    
    int tdis_into = (tdis0_glb + Lt)% Lt; //output time coordinate modulo Lt

    SpinColorFlavorMatrix prop_l_t0_site;
    prop_l_t0.siteMatrix(prop_l_t0_site,x3d_lcl,tdis0_glb);
    if(do_flav_project) prop_l_t0_site *= kaon_proj_t0;

    SpinColorFlavorMatrix prop_h_dag_t0_site;
    prop_h_t0.siteMatrix(prop_h_dag_t0_site,x3d_lcl,tdis0_glb);
    prop_h_dag_t0_site.hconj();
    
    SpinColorFlavorMatrix prop_prod_t0 = prop_l_t0_site * prop_h_dag_t0_site;

    SpinColorFlavorMatrix prop_l_t1_site;
    prop_l_t1.siteMatrix(prop_l_t1_site,x3d_lcl,tdis1_glb);
    if(do_flav_project) prop_l_t1_site *= kaon_proj_t1;

    SpinColorFlavorMatrix prop_h_dag_t1_site;
    prop_h_t1.siteMatrix(prop_h_dag_t1_site,x3d_lcl,tdis1_glb);
    prop_h_dag_t1_site.hconj();

    SpinColorFlavorMatrix prop_prod_t1 = prop_l_t1_site * prop_h_dag_t1_site;

    for(int mu=0;mu<4;mu++){
      for(int Gamma = 0; Gamma < 2; Gamma++){  //\gamma^\mu and \gamma^\mu\gamma^5
	SpinColorFlavorMatrix part1 = prop_prod_t0;
	if(Gamma == 1) part1.gl(-5);
	part1.gl(mu);
	part1.pr(F0);

	SpinColorFlavorMatrix part2 = prop_prod_t1;
	if(Gamma == 1) part2.gl(-5);
	part2.gl(mu);
	part2.pr(F0);

	tmp(tdis_into, omp_get_thread_num()) += 2.0*Trace(part1)*Trace(part2);
	tmp(tdis_into, omp_get_thread_num()) += -2.0*Trace(part1, part2);
      }
    }
  }
  tmp.threadSum();
  tmp.nodeSum();

  for(int tdis=0;tdis<Lt;tdis++)
    into(t0, tdis) = tmp[tdis];
}




void StandardBK(fMatrix<Rcomplex> &into, const int t0, const int t1,
	       const PropSiteMatrixGetter &prop_h_t0, const PropSiteMatrixGetter &prop_l_t0,
	       const PropSiteMatrixGetter &prop_h_t1, const PropSiteMatrixGetter &prop_l_t1){
  const int Lt = GJP.TnodeSites()*GJP.Tnodes();

  const int nthread = omp_get_max_threads();
  basicComplexArray<Rcomplex> tmp(Lt,nthread); //defaults to zero for all elements

  int vol3d = GJP.VolNodeSites()/GJP.TnodeSites();

#pragma omp_parallel for
  for(int x=0;x<GJP.VolNodeSites();x++){
    int pos[4];
    int rem = x;
    for(int i=0;i<4;i++){ pos[i] = rem % GJP.NodeSites(i); rem /= GJP.NodeSites(i); }

    int x3d_lcl = x % vol3d;
    int t_glb = pos[3] + GJP.TnodeCoor() * GJP.TnodeSites(); //operator insertion time
    int tdis0_glb = t_glb - t0; //linear time coordinate
    int tdis1_glb = t_glb - t1;
    
    int tdis_into = (tdis0_glb +Lt)% Lt; //output time coordinate modulo Lt

    WilsonMatrix prop_l_t0_site;
    prop_l_t0.siteMatrix(prop_l_t0_site,x3d_lcl,tdis0_glb);

    WilsonMatrix prop_h_dag_t0_site;
    prop_h_t0.siteMatrix(prop_h_dag_t0_site,x3d_lcl,tdis0_glb);
    prop_h_dag_t0_site.hconj();
    
    WilsonMatrix prop_prod_t0 = prop_l_t0_site * prop_h_dag_t0_site;

    WilsonMatrix prop_l_t1_site;
    prop_l_t1.siteMatrix(prop_l_t1_site,x3d_lcl,tdis1_glb);

    WilsonMatrix prop_h_dag_t1_site;
    prop_h_t1.siteMatrix(prop_h_dag_t1_site,x3d_lcl,tdis1_glb);
    prop_h_dag_t1_site.hconj();

    WilsonMatrix prop_prod_t1 = prop_l_t1_site * prop_h_dag_t1_site;

    for(int mu=0;mu<4;mu++){
      for(int Gamma = 0; Gamma < 2; Gamma++){  //\gamma^\mu and \gamma^\mu\gamma^5
	WilsonMatrix part1 = prop_prod_t0;
	if(Gamma == 1) part1.gl(-5);
	part1.gl(mu);

	WilsonMatrix part2 = prop_prod_t1;
	if(Gamma == 1) part2.gl(-5);
	part2.gl(mu);

	tmp(tdis_into, omp_get_thread_num()) += 2.0*Trace(part1)*Trace(part2);
	tmp(tdis_into, omp_get_thread_num()) += -2.0*Trace(part1, part2);
      }
    }
  }
  tmp.threadSum();
  tmp.nodeSum();

  for(int tdis=0;tdis<Lt;tdis++)
    into(t0, tdis) = tmp[tdis];
}




inline void getBKsnkPropBcAndWrapperTsnk(TbcStatus &time_bc_t1, int &t1, const TbcStatus &time_bc_t0, const int t0, const int tsep){
  const int Lt = GJP.Tnodes()*GJP.TnodeSites();
  time_bc_t1 = time_bc_t0;
  t1 = t0 + tsep;
  if(t1 >= Lt){
    if(time_bc_t0.isCombinedType()){ //Use F(t+Lt) = B(t) and B(t+Lt) = F(t)
      time_bc_t1.swapTbcCombination();
      t1 -= Lt;
    }else if(time_bc_t0.getSingleBc() == BND_CND_PRD){
      t1 -= Lt;
    }else if(time_bc_t0.getSingleBc() == BND_CND_APRD){
      ERR.General("","getBKsnkPropBcAndWrapperTsnk","- sign from tsnk prop crossing boundary not implemented yet\n"); //G(t-Lt) = -G(t), need to pass the minus sign into the function
    }
  }
  assert(t1>=0 && t1<Lt);
}


CPS_END_NAMESPACE

#endif
