#ifndef _COMPUTE_KAON_H
#define _COMPUTE_KAON_H

#include<alg/a2a/required_momenta.h>
#include<alg/a2a/mesonfield_computemany.h>

//Compute stationary kaon two-point function with and without GPBC
CPS_START_NAMESPACE

//Policy for stationary kaon with and without GPBC
class StationaryKaonMomentaPolicy{
public:
  void setupMomenta(const int &ngp){
    RequiredMomentum<StationaryKaonMomentaPolicy> *tt = static_cast<RequiredMomentum<StationaryKaonMomentaPolicy>*>(this);
    if(ngp == 0){
      tt->addP("(0,0,0) + (0,0,0)");
    }else if(ngp == 1){
      tt->addP("(-1,0,0) + (1,0,0)");
    }else if(ngp == 2){
      tt->addP("(-1,-1,0) + (1,1,0)");
    }else if(ngp == 3){
      tt->addP("(-1,-1,-1) + (1,1,1)");
    }else{
      ERR.General("StationaryKaonMomentaPolicy","setupMomenta","ngp cannot be >3\n");
    }
  }
};


template<typename mf_Policies>
class ComputeKaon{
 public:
  typedef typename A2Asource<typename mf_Policies::SourcePolicies::ComplexType, typename mf_Policies::SourcePolicies::DimensionPolicy, typename mf_Policies::SourcePolicies::AllocPolicy>::FieldType::InputParamType FieldParamType;
#ifdef USE_DESTRUCTIVE_FFT
  typedef A2AvectorW<mf_Policies> Wtype;
  typedef A2AvectorV<mf_Policies> Vtype;
#else
  typedef const A2AvectorW<mf_Policies> Wtype;
  typedef const A2AvectorV<mf_Policies> Vtype;
#endif

  //Compute the two-point function using a hydrogen-wavefunction source of radius 'rad'
  //result is indexed by (tsrc, tsep)  where tsep is the source-sink separation
  static void compute(fMatrix<typename mf_Policies::ScalarComplexType> &into,
		      Wtype &W, Vtype &V, 
		      Wtype &W_s, Vtype &V_s,
		      const Float &rad, Lattice &lattice,
		      const FieldParamType &src_setup_params = NullObject()){
    typedef typename mf_Policies::ComplexType ComplexType;
    typedef typename mf_Policies::ScalarComplexType ScalarComplexType;
    typedef typename mf_Policies::SourcePolicies SourcePolicies;

    int Lt = GJP.Tnodes()*GJP.TnodeSites();
    into.resize(Lt,Lt);

    RequiredMomentum<StationaryKaonMomentaPolicy> kaon_momentum;
    ThreeMomentum p_w_src = kaon_momentum.getWmom(0);
    ThreeMomentum p_v_src = kaon_momentum.getVmom(0);
    
    ThreeMomentum p_w_snk = -p_w_src; //sink momentum is opposite source
    ThreeMomentum p_v_snk = -p_v_src;

    std::vector<Wtype*> Wspecies(2); Wspecies[0] = &W; Wspecies[1] = &W_s;
    std::vector<Vtype*> Vspecies(2); Vspecies[0] = &V; Vspecies[1] = &V_s;
    
    //Construct the meson fields
    std::vector<A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorVfftw> > mf_ls(Lt);
    std::vector<A2AmesonField<mf_Policies,A2AvectorWfftw,A2AvectorVfftw> > mf_sl(Lt);
    
    if(!GJP.Gparity()){
      typedef A2AexpSource<SourcePolicies> SourceType;
      typedef SCspinInnerProduct<15,ComplexType, SourceType> InnerType;
      typedef BasicSourceStorage<mf_Policies,InnerType> StorageType;
      
      SourceType src(rad,src_setup_params);
      InnerType g5_inner(src);
      StorageType mf_store(g5_inner);
      
      mf_store.addCompute(0,1, p_w_src,p_v_src);	
      mf_store.addCompute(1,0, p_w_snk,p_v_snk);

      ComputeMesonFields<mf_Policies,StorageType>::compute(mf_store,Wspecies,Vspecies,lattice);
      mf_ls = mf_store[0];
      mf_sl = mf_store[1];
      
    }else{ //For GPBC we need a different smearing function for source and sink because the flavor structure depends on the momentum of the V field, which is opposite between source and sink
      typedef A2AflavorProjectedExpSource<SourcePolicies> SourceType;
      typedef SCFspinflavorInnerProduct<15, ComplexType, SourceType > InnerType;
      typedef GparityFlavorProjectedBasicSourceStorage<mf_Policies,InnerType> StorageType;

      int pbase[3]; //we reset the momentum for each computation so we technically don't need this - however the code demands a valid momentum
      GparityBaseMomentum(pbase,+1);
      
      SourceType src(rad,pbase,src_setup_params);
      InnerType g5_s0_inner(sigma0,src);
      StorageType mf_store(g5_s0_inner,src);

      mf_store.addCompute(0,1, p_w_src,p_v_src);	
      mf_store.addCompute(1,0, p_w_snk,p_v_snk);

#ifndef NODE_DISTRIBUTE_MESONFIELDS
      ComputeMesonFields<mf_Policies,StorageType>::compute(mf_store,Wspecies,Vspecies,lattice);
      mf_ls = mf_store[0];
      mf_sl = mf_store[1];
#else
      ComputeMesonFields<mf_Policies,StorageType>::compute(mf_store,Wspecies,Vspecies,lattice, true);
      for(int t=0;t<Lt;t++){
	mf_store[0][t].nodeGet();
	mf_ls[t].move(mf_store[0][t]);
	mf_store[1][t].nodeGet();
	mf_sl[t].move(mf_store[1][t]);
      }
#endif
    }

    //Compute the two-point function
    trace(into,mf_sl,mf_ls);
    into *= ScalarComplexType(0.5,0);
    rearrangeTsrcTsep(into); //rearrange temporal ordering
  }

};

CPS_END_NAMESPACE

#endif

