//Pion 2pt LW functions pseudoscalar and axial sinks
void measurePion2ptLWGparity(const PropMomContainer &props, const PropPrecision status, const TbcStatus &time_bc, const std::vector<int> &tslices, const MesonMomenta &ll_meson_momenta,
		      const std::string &results_dir, const int conf){
  //Loop over light-light meson momenta
  const int Lt = GJP.Tnodes()*GJP.TnodeSites();

  for(int pidx=0;pidx<ll_meson_momenta.nMom();pidx++){
    ThreeMomentum p_psibar = ll_meson_momenta.getMomentum(SrcPsiBar,pidx);
    ThreeMomentum p_psi = ll_meson_momenta.getMomentum(SrcPsi,pidx);
	  
    ThreeMomentum p_prop_dag = ll_meson_momenta.getMomentum(DaggeredProp,pidx);
    ThreeMomentum p_prop_undag = ll_meson_momenta.getMomentum(UndaggeredProp,pidx);

    Pion2PtSinkOp sink_ops[5] = { AX, AY, AZ, AT, P };//Generate a flavor 'f' gauge fixed wall momentum propagator from given timeslice. Momenta are in units of pi/2L
    std::string sink_op_stub[5] = { "AX", "AY", "AZ", "AT", "P" };
	
    for(int op=0;op<5;op++){

      fMatrix<Rcomplex> results(Lt,Lt); //[tsrc][tsnk-tsrc]
      fMatrix<Rcomplex> results_wrongsinkmom(Lt,Lt); //[tsrc][tsnk-tsrc]  wrong sink momentum (should give zero within statistics)
      fMatrix<Rcomplex> results_wrongproj(Lt,Lt); //[tsrc][tsnk-tsrc]  opposite projection op (optional, used for paper)
      fMatrix<Rcomplex> results_wrongproj_wrongsinkmom(Lt,Lt); //wrong sink mom and wrong projector - non-zero within statistics as discussed in paper

      for(int s=0;s<tslices.size();s++){
	const int tsrc = tslices[s];
	    
	std::auto_ptr<PropSiteMatrixGetter> prop_dag(PropSiteMatrixGetterFactory::get(props,Light,status,tsrc,p_prop_dag,time_bc)); //prop that is daggered
	std::auto_ptr<PropSiteMatrixGetter> prop_undag(PropSiteMatrixGetterFactory::get(props,Light,status,tsrc,p_prop_undag,time_bc));

	pionTwoPointLWGparity(results,tsrc,sink_ops[op],p_psibar,p_psi,*prop_dag,*prop_undag,SPLANE_BOUNDARY,false,false); //right proj op, right sink mom

	if(op == 4){ //only pseudoscalar sink op
	  pionTwoPointLWGparity(results_wrongsinkmom,tsrc,sink_ops[op],p_psibar,p_psi,*prop_dag,*prop_undag,SPLANE_BOUNDARY,false,true); //right proj op, wrong sink mom
	  pionTwoPointLWGparity(results_wrongproj,tsrc,sink_ops[op],p_psibar,p_psi,*prop_dag,*prop_undag,SPLANE_BOUNDARY,true,false); //wrong proj op, right sink mom
	  pionTwoPointLWGparity(results_wrongproj_wrongsinkmom,tsrc,sink_ops[op],p_psibar,p_psi,*prop_dag,*prop_undag,SPLANE_BOUNDARY,true,true); //wrong proj op, wrong sink mom
	}
      }
      writePion2ptLW(results, results_dir, sink_op_stub[op], p_psibar, p_psi, status, time_bc, conf, "");
      if(op == 4){
	writePion2ptLW(results_wrongsinkmom, results_dir, sink_op_stub[op], p_psibar, p_psi, status, time_bc, conf, "_wrongsinkmom");
	writePion2ptLW(results_wrongproj, results_dir, sink_op_stub[op], p_psibar, p_psi, status, time_bc, conf, "_wrongproj");
	writePion2ptLW(results_wrongproj_wrongsinkmom, results_dir, sink_op_stub[op], p_psibar, p_psi, status, time_bc, conf, "_wrongproj_wrongsinkmom");
      }
    }

    //Also do A4 A4 
    fMatrix<Rcomplex> results(Lt,Lt);
    for(int s=0;s<tslices.size();s++){
      const int tsrc = tslices[s];
      
      std::auto_ptr<PropSiteMatrixGetter> prop_dag(PropSiteMatrixGetterFactory::get(props,Light,status,tsrc,p_prop_dag,time_bc)); //prop that is daggered
      std::auto_ptr<PropSiteMatrixGetter> prop_undag(PropSiteMatrixGetterFactory::get(props,Light,status,tsrc,p_prop_undag,time_bc));

      pionTwoPointA4A4LWGparity(results,tsrc,p_psibar,p_psi,*prop_dag,*prop_undag);
    }
    writeBasic2ptLW(results,results_dir,"pion_AT_AT_LW",p_psibar,p_psi,status,time_bc,conf);
  }
}

//Pion 2pt WW functions pseudoscalar sink
void measurePion2ptPPWWGparity(const PropMomContainer &props, const PropPrecision status, const TbcStatus &time_bc, const std::vector<int> &tslices, const MesonMomenta &ll_meson_momenta, Lattice &lat,
			const std::string &results_dir, const int conf){
  const int Lt = GJP.Tnodes()*GJP.TnodeSites();

  //Loop over light-light meson momenta
  for(int pidx=0;pidx<ll_meson_momenta.nMom();pidx++){
    ThreeMomentum p1 = ll_meson_momenta.getMomentum(SrcPsiBar,pidx); //label momentum of psibar at source as p1
    ThreeMomentum p2 = ll_meson_momenta.getMomentum(SrcPsi,pidx);    //label momentum of psi at source as p2
	  
    ThreeMomentum p_prop_dag = -p2; //cf below
    ThreeMomentum p_prop_undag = p1; 
	  
    const ThreeMomentum &p_psibar_src = p1;
    const ThreeMomentum &p_psi_src = p2; //always the same

    //Consider two scenarios with the same total sink momentum

    //1) The quarks each have the same sink momenta as they do at the source (up to the necessary - sign in the phase)
    // \sum_{x1,x2,y1,y2} 
    //<
    //  [[ exp(-i[-p2].x1)\bar\psi(x1,t) A exp(-i[-p1].x2)\psi(x2,t) ]]
    //  *
    //  [[ exp(-i p1.y1)\bar\psi(y1,0) B exp(-i p2.y2)\psi(y2,0) ]]
    //>
    // = 
    //Tr( 
    //   [\sum_{x1,y2} exp(-i[-p2].x1) exp(-i p2.y2) G(y2,0;x1,t)] A
    //   *
    //   [\sum_{x2,y1} exp(-i[-p1].x2) exp(-i p1.y1) G(x2,t;y1,0)] B
    //  )
    // = 
    //Tr( 
    //   g5 [\sum_{x1,y2} exp(-ip2.x1) exp(-i[-p2].y2) G(x1,t;y2,0)]^\dagger g5 A
    //   *
    //   [\sum_{x2,y1} exp(-i[-p1].x2) exp(-i p1.y1) G(x2,t;y1,0)] B
    //  )

    ThreeMomentum p_prop_dag_snk_keep = p2; 
    ThreeMomentum p_prop_undag_snk_keep = -p1;
    ThreeMomentum p_psi_snk_keep = -p1;

    //2) The quarks have their sink momenta exchanged
    // \sum_{x1,x2,y1,y2} 
    //<
    //  [[ exp(-i[-p1].x1)\bar\psi(x1,t) A exp(-i[-p2].x2)\psi(x2,t) ]]
    //  *
    //  [[ exp(-i p1.y1)\bar\psi(y1,0) B exp(-i p2.y2)\psi(y2,0) ]]
    //>
    // = 
    //Tr( 
    //   [\sum_{x1,y2} exp(-i[-p1].x1) exp(-i p2.y2) G(y2,0;x1,t)] A
    //   *
    //   [\sum_{x2,y1} exp(-i[-p2].x2) exp(-i p1.y1) G(x2,t;y1,0)] B
    //  )
    // = 
    //Tr( 
    //   g5 [\sum_{x1,y2} exp(-ip1.x1) exp(-i[-p2].y2) G(x1,t;y2,0)]^\dagger g5 A
    //   *
    //   [\sum_{x2,y1} exp(-i[-p2].x2) exp(-i p1.y1) G(x2,t;y1,0)] B
    //  )

    ThreeMomentum p_prop_dag_snk_exch = p1; 
    ThreeMomentum p_prop_undag_snk_exch = -p2;
    ThreeMomentum p_psi_snk_exch = -p2;

    fMatrix<Rcomplex> results_momkeep(Lt,Lt); //[tsrc][tsnk-tsrc]
    fMatrix<Rcomplex> results_momexch(Lt,Lt);

    for(int s=0;s<tslices.size();s++){
      const int tsrc = tslices[s];
      
      //Prop1
      std::auto_ptr<WallSinkPropSiteMatrixGetter<SpinColorFlavorMatrix> > prop_dag_FT_keep(WallSinkPropSiteMatrixGetterFactory<SpinColorFlavorMatrix>::get(props,Light,status,tsrc,p_prop_dag,p_prop_dag_snk_keep,time_bc,lat));
      std::auto_ptr<WallSinkPropSiteMatrixGetter<SpinColorFlavorMatrix> > prop_dag_FT_exch(WallSinkPropSiteMatrixGetterFactory<SpinColorFlavorMatrix>::get(props,Light,status,tsrc,p_prop_dag,p_prop_dag_snk_exch,time_bc,lat));

      //Prop2
      std::auto_ptr<WallSinkPropSiteMatrixGetter<SpinColorFlavorMatrix> > prop_undag_FT_keep(WallSinkPropSiteMatrixGetterFactory<SpinColorFlavorMatrix>::get(props,Light,status,tsrc,p_prop_undag,p_prop_undag_snk_keep,time_bc,lat));
      std::auto_ptr<WallSinkPropSiteMatrixGetter<SpinColorFlavorMatrix> > prop_undag_FT_exch(WallSinkPropSiteMatrixGetterFactory<SpinColorFlavorMatrix>::get(props,Light,status,tsrc,p_prop_undag,p_prop_undag_snk_exch,time_bc,lat));
 
      pionTwoPointPPWWGparity(results_momkeep, tsrc, p_psi_snk_keep, p_psi_src, *prop_dag_FT_keep, *prop_undag_FT_keep);
      pionTwoPointPPWWGparity(results_momexch, tsrc, p_psi_snk_exch, p_psi_src, *prop_dag_FT_exch, *prop_undag_FT_exch);
    }
    writeBasic2ptLW(results_momkeep,results_dir,"pion_P_P_WW_momkeep",p_psibar_src,p_psi_src,status,time_bc,conf);
    writeBasic2ptLW(results_momexch,results_dir,"pion_P_P_WW_momexch",p_psibar_src,p_psi_src,status,time_bc,conf);
  }
}

//SU(2) singlet 2pt LW functions (only for G-parity)
void measureLightFlavorSingletLW(const PropMomContainer &props, const PropPrecision status, const TbcStatus &time_bc, const std::vector<int> &tslices, const MesonMomenta &meson_momenta,
				 const std::string results_dir, const int conf){
  assert(GJP.Gparity());

  const int Lt = GJP.Tnodes()*GJP.TnodeSites();

  //Loop over light-light meson momenta
  for(int pidx=0;pidx<meson_momenta.nMom();pidx++){
    ThreeMomentum p_psibar = meson_momenta.getMomentum(SrcPsiBar,pidx);
    ThreeMomentum p_psi = meson_momenta.getMomentum(SrcPsi,pidx);
	  
    ThreeMomentum p_prop_dag = meson_momenta.getMomentum(DaggeredProp,pidx);
    ThreeMomentum p_prop_undag = meson_momenta.getMomentum(UndaggeredProp,pidx);
	  
    fMatrix<Rcomplex> results(Lt,Lt); //[tsrc][tsnk-tsrc]
    for(int s=0;s<tslices.size();s++){
      const int tsrc = tslices[s];
	    
      std::auto_ptr<PropSiteMatrixGetter> prop_dag(PropSiteMatrixGetterFactory::get(props,Light,status,tsrc,p_prop_dag,time_bc)); //prop that is daggered
      std::auto_ptr<PropSiteMatrixGetter> prop_undag(PropSiteMatrixGetterFactory::get(props,Light,status,tsrc,p_prop_undag,time_bc));

      lightFlavorSingletLWGparity(results,tsrc,p_psibar,p_psi,*prop_dag,*prop_undag);
    }
    writeBasic2ptLW(results,results_dir,"light_pseudo_flav_singlet_LW",p_psibar,p_psi,status,time_bc,conf);
  }
}


//Kaon 2pt LW
void measureKaon2ptLWGparity(const PropMomContainer &props, const PropPrecision status, const TbcStatus &time_bc, const std::vector<int> &tslices, const MesonMomenta &meson_momenta,
		      const std::string results_dir, const int conf){

  const int Lt = GJP.Tnodes()*GJP.TnodeSites();

  for(int pidx=0;pidx<meson_momenta.nMom();pidx++){
    assert(meson_momenta.getQuarkType(SrcPsi,pidx) == Heavy); 
    assert(meson_momenta.getQuarkType(SrcPsiBar,pidx) == Light);

    ThreeMomentum p_psibar_l = meson_momenta.getMomentum(SrcPsiBar,pidx);
    ThreeMomentum p_psi_h = meson_momenta.getMomentum(SrcPsi,pidx);
	  
    ThreeMomentum p_prop_dag_h = meson_momenta.getMomentum(DaggeredProp,pidx); //daggered prop is heavy here
    ThreeMomentum p_prop_undag_l = meson_momenta.getMomentum(UndaggeredProp,pidx);

    fMatrix<Rcomplex> results_PP(Lt,Lt); //[tsrc][tsnk-tsrc]
    fMatrix<Rcomplex> results_A4physP(Lt,Lt);
    fMatrix<Rcomplex> results_A4unphysP(Lt,Lt);
    fMatrix<Rcomplex> results_A4combA4comb(Lt,Lt);

    for(int s=0;s<tslices.size();s++){
      const int tsrc = tslices[s];
	    
      std::auto_ptr<PropSiteMatrixGetter> prop_dag_h(PropSiteMatrixGetterFactory::get(props,Heavy,status,tsrc,p_prop_dag_h,time_bc)); 
      std::auto_ptr<PropSiteMatrixGetter> prop_undag_l(PropSiteMatrixGetterFactory::get(props,Light,status,tsrc,p_prop_undag_l,time_bc));

      kaonTwoPointPPLWGparity(results_PP, tsrc, p_psibar_l, p_psi_h, *prop_dag_h, *prop_undag_l);
      kaonTwoPointA4PhysPLWGparity(results_A4physP, tsrc, p_psibar_l, p_psi_h, *prop_dag_h, *prop_undag_l);
      kaonTwoPointA4UnphysPLWGparity(results_A4unphysP, tsrc, p_psibar_l, p_psi_h, *prop_dag_h, *prop_undag_l);
      kaonTwoPointA4combA4combLWGparity(results_A4combA4comb, tsrc, p_psibar_l, p_psi_h, *prop_dag_h, *prop_undag_l);
    }
    writeBasic2ptLW(results_PP,results_dir,"kaon_P_P_LW",p_psibar_l,p_psi_h,status,time_bc,conf);
    writeBasic2ptLW(results_A4physP,results_dir,"kaon_ATphys_P_LW",p_psibar_l,p_psi_h,status,time_bc,conf);
    writeBasic2ptLW(results_A4unphysP,results_dir,"kaon_ATunphys_P_LW",p_psibar_l,p_psi_h,status,time_bc,conf);
    writeBasic2ptLW(results_A4combA4comb,results_dir,"kaon_ATcomb_ATcomb_LW",p_psibar_l,p_psi_h,status,time_bc,conf);
  }
}

//Kaon 2pt LW functions pseudoscalar sink (cf pion version for comments)
void measureKaon2ptPPWWGparity(const PropMomContainer &props, const PropPrecision status, const TbcStatus &time_bc, const std::vector<int> &tslices, const MesonMomenta &meson_momenta, Lattice &lat,
			const std::string &results_dir, const int conf){
  const int Lt = GJP.Tnodes()*GJP.TnodeSites();

  //Loop over light-light meson momenta
  for(int pidx=0;pidx<meson_momenta.nMom();pidx++){
    assert(meson_momenta.getQuarkType(SrcPsi,pidx) == Heavy); 
    assert(meson_momenta.getQuarkType(SrcPsiBar,pidx) == Light);

    ThreeMomentum p1 = meson_momenta.getMomentum(SrcPsiBar,pidx); //label momentum of psibar_l at source as p1
    ThreeMomentum p2 = meson_momenta.getMomentum(SrcPsi,pidx);    //label momentum of psi_h at source as p2
	  
    ThreeMomentum p_prop_h_dag = -p2; //cf below
    ThreeMomentum p_prop_l_undag = p1; 
	  
    const ThreeMomentum &p_psibar_l_src = p1;
    const ThreeMomentum &p_psi_h_src = p2; //always the same

    ThreeMomentum p_prop_h_dag_snk_keep = p2; 
    ThreeMomentum p_prop_l_undag_snk_keep = -p1;
    ThreeMomentum p_psi_l_snk_keep = -p1;

    ThreeMomentum p_prop_h_dag_snk_exch = p1; 
    ThreeMomentum p_prop_l_undag_snk_exch = -p2;
    ThreeMomentum p_psi_l_snk_exch = -p2;

    fMatrix<Rcomplex> results_momkeep(Lt,Lt); //[tsrc][tsnk-tsrc]
    fMatrix<Rcomplex> results_momexch(Lt,Lt);

    for(int s=0;s<tslices.size();s++){
      const int tsrc = tslices[s];
      
      //Heavy prop
      std::auto_ptr<WallSinkPropSiteMatrixGetter<SpinColorFlavorMatrix> > prop_h_dag_FT_keep(WallSinkPropSiteMatrixGetterFactory<SpinColorFlavorMatrix>::get(props,Heavy,status,tsrc,p_prop_h_dag,p_prop_h_dag_snk_keep,time_bc,lat));
      std::auto_ptr<WallSinkPropSiteMatrixGetter<SpinColorFlavorMatrix> > prop_h_dag_FT_exch(WallSinkPropSiteMatrixGetterFactory<SpinColorFlavorMatrix>::get(props,Heavy,status,tsrc,p_prop_h_dag,p_prop_h_dag_snk_exch,time_bc,lat));

      //Light prop
      std::auto_ptr<WallSinkPropSiteMatrixGetter<SpinColorFlavorMatrix> > prop_l_undag_FT_keep(WallSinkPropSiteMatrixGetterFactory<SpinColorFlavorMatrix>::get(props,Light,status,tsrc,p_prop_l_undag,p_prop_l_undag_snk_keep,time_bc,lat));
      std::auto_ptr<WallSinkPropSiteMatrixGetter<SpinColorFlavorMatrix> > prop_l_undag_FT_exch(WallSinkPropSiteMatrixGetterFactory<SpinColorFlavorMatrix>::get(props,Light,status,tsrc,p_prop_l_undag,p_prop_l_undag_snk_exch,time_bc,lat));

      kaonTwoPointPPWWGparity(results_momkeep, tsrc, p_psi_l_snk_keep, p_psi_h_src, *prop_h_dag_FT_keep, *prop_l_undag_FT_keep);
      kaonTwoPointPPWWGparity(results_momexch, tsrc, p_psi_l_snk_exch, p_psi_h_src, *prop_h_dag_FT_exch, *prop_l_undag_FT_exch);
    }
    writeBasic2ptLW(results_momkeep,results_dir,"kaon_P_P_WW_momkeep",p_psibar_l_src,p_psi_h_src,status,time_bc,conf);
    writeBasic2ptLW(results_momexch,results_dir,"kaon_P_P_WW_momexch",p_psibar_l_src,p_psi_h_src,status,time_bc,conf);
  }
}

//Measure BK with source kaons on each of the timeslices t0 in prop_tsources and K->K time separations tseps
//Can use standard P or A time BCs but you will need to use closer-together kaon sources to avoid round-the-world effects. These can be eliminated by using the F=P+A and B=P-A combinations
//Either can be specified using the appropriate time_bc parameter below
//Can optionally choose to disable the source/sink flavor projection
void measureBKGparity(const PropMomContainer &props, const PropPrecision status, const std::vector<int> &prop_tsources, const std::vector<int> &tseps, const MesonMomenta &meson_momenta,
	       const TbcStatus &time_bc_t0,
	       const std::string &results_dir, const int conf, const bool do_flavor_project = true){
  const int Lt = GJP.Tnodes()*GJP.TnodeSites();
  
  //<  \bar\psi_l s3 g5 \psi_h   O_VV+AA   \bar\psi_l s3 g5 \psi_h >

  //Do all combinations of source and sink kaon momenta that have the same total momentum. This allows us to look at alternate quark momentum combinations
  //In the meson_momenta, prop index 0 is the strange quark (as in the standard 2pt function case), and is the propagator that is daggered. Prop index 1 is the light quark.
  for(int p0idx=0;p0idx<meson_momenta.nMom();p0idx++){
    assert(meson_momenta.getQuarkType(SrcPsi,p0idx) == Heavy);
    assert(meson_momenta.getQuarkType(SrcPsiBar,p0idx) == Light);

    ThreeMomentum p_psibar_l_t0 = meson_momenta.getMomentum(SrcPsiBar,p0idx);
    ThreeMomentum p_psi_h_t0 = meson_momenta.getMomentum(SrcPsi,p0idx);

    ThreeMomentum p_prop_h_t0 = meson_momenta.getMomentum(DaggeredProp,p0idx);
    ThreeMomentum p_prop_l_t0 = meson_momenta.getMomentum(UndaggeredProp,p0idx);

    for(int p1idx=p0idx;p1idx<meson_momenta.nMom();p1idx++){
      if(meson_momenta.getMomentum(Total,p1idx) != meson_momenta.getMomentum(Total,p0idx)) continue;

      ThreeMomentum p_psibar_l_t1 = meson_momenta.getMomentum(SrcPsiBar,p1idx);
      ThreeMomentum p_psi_h_t1 = meson_momenta.getMomentum(SrcPsi,p1idx);
      
      ThreeMomentum p_prop_h_t1 = meson_momenta.getMomentum(DaggeredProp,p1idx);
      ThreeMomentum p_prop_l_t1 = meson_momenta.getMomentum(UndaggeredProp,p1idx);

      for(int tspi=0;tspi<tseps.size();tspi++){
	fMatrix<Rcomplex> results(Lt,Lt);

	for(int t0i=0;t0i<prop_tsources.size();t0i++){
	  int t0 = prop_tsources[t0i];
	  int t1; TbcStatus time_bc_t1(time_bc_t0);
	  getBKsnkPropBcAndWrapperTsnk(time_bc_t1,t1, time_bc_t0, t0, tseps[tspi]);

	  //check t1 is in the vector, if not skip
	  if(std::find(prop_tsources.begin(), prop_tsources.end(), t1) == prop_tsources.end()) continue;

	  std::auto_ptr<PropSiteMatrixGetter> prop_h_t0(PropSiteMatrixGetterFactory::get(props,Heavy,status,t0,p_prop_h_t0,time_bc_t0));
	  std::auto_ptr<PropSiteMatrixGetter> prop_l_t0(PropSiteMatrixGetterFactory::get(props,Light,status,t0,p_prop_l_t0,time_bc_t0));

	  std::auto_ptr<PropSiteMatrixGetter> prop_h_t1(PropSiteMatrixGetterFactory::get(props,Heavy,status,t1,p_prop_h_t1,time_bc_t1));
	  std::auto_ptr<PropSiteMatrixGetter> prop_l_t1(PropSiteMatrixGetterFactory::get(props,Light,status,t1,p_prop_l_t1,time_bc_t1));

	  GparityBK(results, t0, t1,
		    *prop_h_t0, *prop_l_t0, p_psi_h_t0,
		    *prop_h_t1, *prop_l_t1, p_psi_h_t1,
		    do_flavor_project);
	}
	{
	  std::ostringstream os;
	  os << results_dir << "/BK_srcMom" << p_psibar_l_t0.file_str() << "_plus" << p_psi_h_t0.file_str() 
	     << "snkMom" << p_psibar_l_t1.file_str() << "_plus" << p_psi_h_t1.file_str()
	     << "_srcTbc" << time_bc_t0.getTag()
	     << "_tsep" << tseps[tspi]
	     << (status == Sloppy ? "_sloppy" : "_exact") 	    
	     << '.' << conf;
	  results.write(os.str());
	}
      }
    }
  }
}

//Note: Mres is only properly defined with APRD time BCs. A runtime check is performed
//Can optionally choose to disable the source/sink flavor projection 
void measureMresGparity(const PropMomContainer &props, const PropPrecision status, const TbcStatus &time_bc, const std::vector<int> &tslices, const MesonMomenta &meson_momenta,
		 const std::string &results_dir, const int conf, const bool do_flavor_project = true){
  const int Lt = GJP.Tnodes()*GJP.TnodeSites();

  for(int pidx=0;pidx<meson_momenta.nMom();pidx++){
    ThreeMomentum p_psibar = meson_momenta.getMomentum(SrcPsiBar,pidx);
    ThreeMomentum p_psi = meson_momenta.getMomentum(SrcPsi,pidx);
	  
    ThreeMomentum p_prop_dag = meson_momenta.getMomentum(DaggeredProp,pidx);
    ThreeMomentum p_prop_undag = meson_momenta.getMomentum(UndaggeredProp,pidx);
	  
    fMatrix<Rcomplex> pion(Lt,Lt); //[tsrc][tsnk-tsrc]
    fMatrix<Rcomplex> j5q(Lt,Lt);
    
    for(int s=0;s<tslices.size();s++){
      const int tsrc = tslices[s];
	    
      std::auto_ptr<PropSiteMatrixGetter> prop_dag(PropSiteMatrixGetterFactory::get(props,Light,status,tsrc,p_prop_dag,time_bc)); //prop that is daggered
      std::auto_ptr<PropSiteMatrixGetter> prop_undag(PropSiteMatrixGetterFactory::get(props,Light,status,tsrc,p_prop_undag,time_bc));

      J5Gparity(pion,tsrc,p_psibar,p_psi,*prop_dag,*prop_undag,SPLANE_BOUNDARY,do_flavor_project); 
      J5Gparity(j5q,tsrc,p_psibar,p_psi,*prop_dag,*prop_undag,SPLANE_MIDPOINT,do_flavor_project);
    }
    writeBasic2ptLW(pion,results_dir,"J5_LW",p_psibar,p_psi,status,time_bc,conf);
    writeBasic2ptLW(j5q,results_dir,"J5q_LW",p_psibar,p_psi,status,time_bc,conf);
  }
}
