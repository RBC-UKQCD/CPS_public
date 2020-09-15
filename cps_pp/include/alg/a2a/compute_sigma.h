#ifndef _COMPUTE_SIGMA_H
#define _COMPUTE_SIGMA_H

#include<alg/a2a/required_momenta.h>
#include<alg/a2a/mesonfield_computemany.h>

//Compute stationary sigma meson two-point function with and without GPBC
CPS_START_NAMESPACE

//Policy for stationary sigma with and without GPBC
class StationarySigmaMomentaPolicy{
public:
  void setupMomenta(const int &ngp){
    RequiredMomentum<StationarySigmaMomentaPolicy> *tt = static_cast<RequiredMomentum<StationarySigmaMomentaPolicy>*>(this);
    if(ngp == 0){
      tt->addP("(0,0,0) + (0,0,0)");
    }else if(ngp == 1){
      tt->addPandMinusP("(-1,0,0) + (1,0,0)");
      tt->addPandMinusP("(-3,0,0) + (3,0,0)");
    }else if(ngp == 2){
      tt->addPandMinusP("(-1,-1,0) + (1,1,0)");
      tt->addPandMinusP("(3,-1,0) + (-3,1,0)");
      tt->addPandMinusP("(-1,3,0) + (1,-3,0)");
    }else if(ngp == 3){
      tt->addPandMinusP("(-1,-1,-1) + (1,1,1)");
      tt->addPandMinusP("(3,-1,-1) + (-3,1,1)");
      tt->addPandMinusP("(-1,3,-1) + (1,-3,1)");
      tt->addPandMinusP("(-1,-1,3) + (1,1,-3)");
    }else{
      ERR.General("StationarySigmaMomentaPolicy","setupMomenta","ngp cannot be >3\n");
    }
  }
};


template<typename mf_Policies>
class ComputeSigma{
 public:
  typedef typename A2Asource<typename mf_Policies::SourcePolicies::ComplexType, typename mf_Policies::SourcePolicies::DimensionPolicy, typename mf_Policies::SourcePolicies::AllocPolicy>::FieldType::InputParamType FieldParamType;

#ifdef USE_DESTRUCTIVE_FFT
  typedef A2AvectorW<mf_Policies> Wtype;
  typedef A2AvectorV<mf_Policies> Vtype;
#else
  typedef const A2AvectorW<mf_Policies> Wtype;
  typedef const A2AvectorV<mf_Policies> Vtype;
#endif
  
  //Computes sigma meson fields and saves to disk
  static void computeAndWrite(const std::string &work_dir, const int traj,
			      Wtype &W, Vtype &V, const Float &rad, Lattice &lattice,
			      const FieldParamType &src_setup_params = NullObject()){

    typedef typename mf_Policies::ComplexType ComplexType;
    typedef typename mf_Policies::ScalarComplexType ScalarComplexType;
    typedef typename mf_Policies::SourcePolicies SourcePolicies;
    typedef A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorVfftw> MesonFieldType;
    typedef std::vector<MesonFieldType> MesonFieldVectorType;
    
    int Lt = GJP.Tnodes()*GJP.TnodeSites();

    RequiredMomentum<StationarySigmaMomentaPolicy> momenta;

    std::vector<Wtype*> Wspecies(1,&W);
    std::vector<Vtype*> Vspecies(1,&V);
    
    if(GJP.Gparity()){
      typedef A2AflavorProjectedExpSource<SourcePolicies> ExpSrcType;
      typedef A2AflavorProjectedHydrogenSource<SourcePolicies> HydSrcType;

      int pbase[3]; //we reset the momentum for each computation so we technically don't need this - however the code demands a valid momentum
      GparityBaseMomentum(pbase,+1);

#ifdef SIGMA_DO_SOURCES_SEPARATELY
      //Multi-src multi-mom strategy consumes a lot of memory - too much for a 64-node job on Cori I. This option does the two sources separately, reducing the memory usage by a factor of 2 at the loss of computational efficiency
      typedef GparitySourceShiftInnerProduct<ComplexType,ExpSrcType, flavorMatrixSpinColorContract<0,ComplexType,true,false> > ExpInnerType;
      typedef GparityFlavorProjectedShiftSourceStorage<mf_Policies, ExpInnerType> ExpStorageType;

      typedef GparitySourceShiftInnerProduct<ComplexType,HydSrcType, flavorMatrixSpinColorContract<0,ComplexType,true,false> > HydInnerType;
      typedef GparityFlavorProjectedShiftSourceStorage<mf_Policies, HydInnerType> HydStorageType;

      ExpSrcType exp_src(rad,pbase,src_setup_params); //1s
      HydSrcType hyd_src(2,0,0,rad,pbase,src_setup_params); //2s
      
      ExpInnerType exp_gunit_s0_inner(sigma0, exp_src);
      ExpStorageType exp_mf_store(exp_gunit_s0_inner,exp_src);

      HydInnerType hyd_gunit_s0_inner(sigma0, hyd_src);
      HydStorageType hyd_mf_store(hyd_gunit_s0_inner,hyd_src);

      for(int pidx=0;pidx<momenta.nMom();pidx++){
	ThreeMomentum p_w = momenta.getWmom(pidx,false);
	ThreeMomentum p_v = momenta.getVmom(pidx,false);
	exp_mf_store.addCompute(0,0, p_w,p_v);	
	hyd_mf_store.addCompute(0,0, p_w,p_v);	
      }
      if(!UniqueID()) printf("Computing sigma meson fields with 1s source\n");

      ComputeMesonFields<mf_Policies,ExpStorageType>::compute(exp_mf_store,Wspecies,Vspecies,lattice
#  ifdef NODE_DISTRIBUTE_MESONFIELDS
							      ,true
#  endif
							      );
      
      if(!UniqueID()) printf("Computing sigma meson fields with 2s source\n");
      ComputeMesonFields<mf_Policies,HydStorageType>::compute(hyd_mf_store,Wspecies,Vspecies,lattice
#  ifdef NODE_DISTRIBUTE_MESONFIELDS
							      ,true
#  endif
							      );
#else      
      typedef Elem<ExpSrcType, Elem<HydSrcType,ListEnd > > SrcList;
      typedef A2AmultiSource<SrcList> MultiSrcType;      
      //typedef SCFspinflavorInnerProduct<0,ComplexType,MultiSrcType,true,false> MultiInnerType; //unit matrix spin structure
      //typedef GparityFlavorProjectedMultiSourceStorage<mf_Policies, MultiInnerType> StorageType;

      //Allows for more memory efficient computation algorithm
      typedef GparitySourceShiftInnerProduct<ComplexType,MultiSrcType, flavorMatrixSpinColorContract<0,ComplexType,true,false> > MultiInnerType;
      typedef GparityFlavorProjectedShiftSourceStorage<mf_Policies, MultiInnerType> StorageType;
      
      MultiSrcType src;
      src.template getSource<0>().setup(rad,pbase,src_setup_params); //1s
      src.template getSource<1>().setup(2,0,0,rad,pbase,src_setup_params); //2s
      
      MultiInnerType gunit_s0_inner(sigma0, src);
      StorageType mf_store(gunit_s0_inner,src);

      for(int pidx=0;pidx<momenta.nMom();pidx++){
	ThreeMomentum p_w = momenta.getWmom(pidx,false);
	ThreeMomentum p_v = momenta.getVmom(pidx,false);
	mf_store.addCompute(0,0, p_w,p_v);	
      }
      if(!UniqueID()) printf("Computing sigma meson fields with 1s/2s sources\n");

      ComputeMesonFields<mf_Policies,StorageType>::compute(mf_store,Wspecies,Vspecies,lattice
#  ifdef NODE_DISTRIBUTE_MESONFIELDS
							   ,true
#  endif
							   );

#endif

      std::string src_names[2] = {"1s","2s"};
      if(!UniqueID()) printf("Writing sigma meson fields to disk\n");
      for(int pidx=0;pidx<momenta.nMom();pidx++){
	ThreeMomentum p_wdag = -momenta.getWmom(pidx,false);
	ThreeMomentum p_v = momenta.getVmom(pidx,false);
	
	for(int s=0;s<2;s++){
	  std::ostringstream os; //momenta in units of pi/2L
	  os << work_dir << "/traj_" << traj << "_sigma_mfwv_mom" << p_wdag.file_str() << "_plus" << p_v.file_str() << "_hyd" << src_names[s] << "_rad" << rad << ".dat";
#ifdef SIGMA_DO_SOURCES_SEPARATELY
	  MesonFieldVectorType &mf_q = (s == 0 ? exp_mf_store[pidx] : hyd_mf_store[pidx] );
#else
	  MesonFieldVectorType &mf_q = mf_store(s,pidx);
#endif

#ifdef NODE_DISTRIBUTE_MESONFIELDS
	  nodeGetMany(1,&mf_q);
#endif
	  MesonFieldType::write(os.str(),mf_q);
	  for(int t=0;t<Lt;t++) mf_q[t].free_mem(); //no longer needed
	}
      } 

    }else{
      typedef A2AexpSource<SourcePolicies> SrcType;
      typedef SCspinInnerProduct<0,ComplexType,SrcType> InnerType;
      typedef BasicSourceStorage<mf_Policies,InnerType> StorageType;
      
      SrcType src(rad,src_setup_params);
      InnerType gunit_inner(src);

      StorageType mf_store(gunit_inner);

      for(int pidx=0;pidx<momenta.nMom();pidx++){
	ThreeMomentum p_w = momenta.getWmom(pidx,false);
	ThreeMomentum p_v = momenta.getVmom(pidx,false);
	mf_store.addCompute(0,0, p_w,p_v);	
      }
      ComputeMesonFields<mf_Policies,StorageType>::compute(mf_store,Wspecies,Vspecies,lattice);
      
      for(int pidx=0;pidx<momenta.nMom();pidx++){
	ThreeMomentum p_wdag = -momenta.getWmom(pidx,false);
	ThreeMomentum p_v = momenta.getVmom(pidx,false);
	
	std::ostringstream os; //momenta in units of pi/2L
	os << work_dir << "/traj_" << traj << "_sigma_mfwv_mom" << p_wdag.file_str() << "_plus" << p_v.file_str() << "_hyd1s_rad" << rad << ".dat";
	MesonFieldType::write(os.str(),mf_store[pidx]);
      } 
    }
  }

};

CPS_END_NAMESPACE

#endif

