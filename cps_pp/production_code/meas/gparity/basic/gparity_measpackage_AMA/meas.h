#ifndef _MEAS_GP_H  
#define _MEAS_GP_H  

#include <algorithm>
#include <memory>
#include "pion_twopoint.h"
#include "kaon_twopoint.h"
#include "compute_bk.h"
#include "enums.h"
#include "prop_tag.h"
#include <alg/eigen/Krylov_5d.h>

CPS_START_NAMESPACE

QPropWMomSrc* randomSolutionPropagator(const bool store_midprop, Lattice &latt){
  if(!UniqueID()) printf("Generating random propagator solution vector\n");
  CommonArg c_arg;
  QPropWMomSrc* ret = new QPropWMomSrc(latt,&c_arg); //this constructor does nothing
  const int msize = 12*12*2;
  ret->Allocate(PROP);
  for(int f=0;f<GJP.Gparity()+1;f++){
    for(int x=0;x<GJP.VolNodeSites();x++){
      LRG.AssignGenerator(x,f);
      Float* off = (Float*)&(ret->SiteMatrix(x,f));
      for(int i=0;i<msize;i++) off[i] = LRG.Urand(0.5,-0.5,FOUR_D);
    }
  }
  if(store_midprop){
    ret->Allocate(MIDPROP);
    for(int f=0;f<GJP.Gparity()+1;f++){
      for(int x=0;x<GJP.VolNodeSites();x++){
	LRG.AssignGenerator(x,f);
	Float* off = (Float*)&(ret->MidPlaneSiteMatrix(x,f));
	for(int i=0;i<msize;i++) off[i] = LRG.Urand(0.5,-0.5,FOUR_D);
      }
    }
  }
  return ret;
}


//Generate a flavor 'f' gauge fixed wall momentum propagator from given timeslice. Momenta are in units of pi/2L
//Eigenvectors must be those appropriate to the choice of temporal BC, 'time_bc'

//Note: If using Fbfm or FGrid, the current temporal BC listed in GJP.Tbc() must be applied to the bfm/Grid internal gauge field (i.e. minuses on t-links at boundard for APRD) prior to using this method. Internally
//it changes the bc to 'time_bc' but it changes it back at the end.
QPropWMomSrc* computePropagator(const double mass, const double stop_prec, const int t, const int flav, const int p[3], const BndCndType time_bc, const bool store_midprop, 
				Lattice &latt,  BFM_Krylov::Lanczos_5d<double> *deflate = NULL, const bool random_solution = false){ 
  if(random_solution) return randomSolutionPropagator(store_midprop,latt);

  multi1d<float> *eval_conv = NULL;

  if(deflate != NULL){
    if(latt.Fclass() != F_CLASS_BFM && latt.Fclass() != F_CLASS_BFM_TYPE2)
      ERR.General("","computePropagator","Deflation only implemented for Fbfm\n");
    if(Fbfm::use_mixed_solver){
      //Have to convert evals to single prec
      eval_conv = new multi1d<float>(deflate->bl.size());
      for(int i=0;i<eval_conv->size();i++) eval_conv->operator[](i) = deflate->bl[i];
      dynamic_cast<Fbfm&>(latt).set_deflation<float>(&deflate->bq,eval_conv,0); //last argument is really obscure - it's the number of eigenvectors subtracted from the solution to produce a high-mode inverse - we want zero here
    }else dynamic_cast<Fbfm&>(latt).set_deflation(&deflate->bq,&deflate->bl,0);
  }

  CommonArg c_arg;
  
  CgArg cg;
  cg.mass = mass;
  cg.max_num_iter = 10000;
  cg.stop_rsd = stop_prec;
  cg.true_rsd = stop_prec;
  cg.RitzMatOper = NONE;
  cg.Inverter = CG;
  cg.bicgstab_n = 0;

  QPropWArg qpropw_arg;
  qpropw_arg.cg = cg;
  qpropw_arg.x = 0;
  qpropw_arg.y = 0;
  qpropw_arg.z = 0;
  qpropw_arg.t = t;
  qpropw_arg.flavor = flav; 
  qpropw_arg.ensemble_label = "ens";
  qpropw_arg.ensemble_id = "ens_id";
  qpropw_arg.StartSrcSpin = 0;
  qpropw_arg.EndSrcSpin = 4;
  qpropw_arg.StartSrcColor = 0;
  qpropw_arg.EndSrcColor = 3;
  qpropw_arg.gauge_fix_src = 1;
  qpropw_arg.gauge_fix_snk = 0;
  qpropw_arg.store_midprop = store_midprop ? 1 : 0; //for mres

  //Switching boundary conditions is poorly implemented in Fbfm (and probably FGrid)
  //For traditional lattice types the BC is applied by modififying the gauge field when the Dirac operator is created and reverting when destroyed. This only ever happens internally - no global instance of the Dirac operator exists
  //On the other hand, Fbfm does all its inversion internally and doesn't instantiate a CPS Dirac operator. We therefore have to manually force Fbfm to change its internal gauge field by applying BondCond
  bool is_wrapper_type = ( latt.Fclass() == F_CLASS_BFM || latt.Fclass() == F_CLASS_BFM_TYPE2 ); //I hate this!

  BndCndType init_tbc = GJP.Tbc();
  BndCndType target_tbc = time_bc;

  GJP.Bc(3,target_tbc);
  if(is_wrapper_type) latt.BondCond();  //Apply new BC to internal gauge fields

  QPropWMomSrc* ret = new QPropWMomSrc(latt,&qpropw_arg,const_cast<int*>(p),&c_arg);

  //Restore the BCs
  if(is_wrapper_type) latt.BondCond();  //unapply existing BC
  GJP.Bc(3,init_tbc);
  
  if(deflate != NULL) dynamic_cast<Fbfm&>(latt).unset_deflation();
  if(eval_conv !=NULL) delete eval_conv;

  return ret;
}

QPropWMomSrc* computePropagator(const double mass, const double stop_prec, const int t, const int flav, const ThreeMomentum &p, const BndCndType time_bc, const bool store_midprop, 
				Lattice &latt, BFM_Krylov::Lanczos_5d<double> *deflate = NULL, const bool random_solution = false){ 
  return computePropagator(mass,stop_prec,t,flav,p.ptr(),time_bc,store_midprop,latt,deflate,random_solution);
}

void quarkInvert(PropMomContainer &props, const QuarkType qtype, const PropPrecision pp, const double stop_prec, const double mass, const BndCndType time_bc,
		 const std::vector<int> &tslices, const QuarkMomenta &quark_momenta, const bool store_midprop, 
		 Lattice &lattice, BFM_Krylov::Lanczos_5d<double> *lanc = NULL, const bool random_solution = false){
  if(!UniqueID()) printf("Computing %s %s quark propagators\n", pp == Sloppy ? "sloppy":"exact", qtype==Light ? "light" : "heavy");
  double time = -dclock();
  
  for(int s=0;s<tslices.size();s++){
    const int tsrc = tslices[s];
    
    for(int pidx=0;pidx<quark_momenta.nMom();pidx++){
      const ThreeMomentum &p = quark_momenta.getMom(pidx);
      if(!UniqueID()) std::cout << "Starting inversion for prop on timeslice " << tsrc << " with momentum phase " << p.str() << '\n';  

      if(GJP.Gparity()){
	QPropWMomSrc* prop_f0 = computePropagator(mass,stop_prec,tsrc,0,p.ptr(),time_bc,store_midprop,lattice,lanc,random_solution);
	QPropWMomSrc* prop_f1 = computePropagator(mass,stop_prec,tsrc,1,p.ptr(),time_bc,store_midprop,lattice,lanc,random_solution);
	
	//Add both + and - source momentum  (PropMomContainer manages prop memory)
	PropWrapper prop_pplus(prop_f0,prop_f1,false);
	props.insert(prop_pplus, propTag(qtype,pp,tsrc,p,time_bc));
	
	PropWrapper prop_pminus(prop_f0,prop_f1,true);
	props.insert(prop_pminus, propTag(qtype,pp,tsrc,-p,time_bc));
      }else{
	QPropWMomSrc* prop = computePropagator(mass,stop_prec,tsrc,0,p.ptr(),time_bc,store_midprop,lattice,lanc,random_solution);
	PropWrapper propw(prop);
	props.insert(propw, propTag(qtype,pp,tsrc,p,time_bc));
      }
    }
  }
  print_time("main","Inversions",time + dclock());
}

//Combine quarks with P and A Tbcs into F=P+A and B=P-A types which are added to the PropMomContainer with appropriate tags
static void quarkCombine(PropMomContainer &props, const QuarkType qtype, const PropPrecision pp, const std::vector<int> &tslices, const QuarkMomenta &quark_momenta){
  if(!UniqueID()) printf("Combining %s %s quark propagators with different Tbcs\n", pp == Sloppy ? "sloppy":"exact", qtype==Light ? "light" : "heavy");
  double time = -dclock();
  
  for(int s=0;s<tslices.size();s++){
    const int tsrc = tslices[s];
    
    for(int pidx=0;pidx<quark_momenta.nMom();pidx++){
      const ThreeMomentum &p = quark_momenta.getMom(pidx);
      if(!UniqueID()) std::cout << "Starting combination of props on timeslice " << tsrc << " with momentum phase " << p.str() << '\n';  

      PropWrapper prop_P = props.get(propTag(qtype,pp,tsrc,p,BND_CND_PRD));
      PropWrapper prop_A = props.get(propTag(qtype,pp,tsrc,p,BND_CND_APRD));

      PropWrapper combF = PropWrapper::combinePA(prop_P,prop_A,CombinationF);
      props.insert(combF, propTag(qtype,pp,tsrc,p,CombinationF));

      PropWrapper combB = PropWrapper::combinePA(prop_P,prop_A,CombinationB);
      props.insert(combB, propTag(qtype,pp,tsrc,p,CombinationB));

      if(GJP.Gparity()){
	//We can just change the flip flag and the momentum for the -ve mom counterparts without using up extra memory
	combF.setFlip(true);
	props.insert(combF, propTag(qtype,pp,tsrc,-p,CombinationF));
	
	combB.setFlip(true);
	props.insert(combB, propTag(qtype,pp,tsrc,-p,CombinationB));
      }
    }
  }
  print_time("main","Combinations",time + dclock());
}

inline std::auto_ptr<BFM_Krylov::Lanczos_5d<double> > doLanczos(GnoneFbfm &lattice, const LancArg lanc_arg, const BndCndType time_bc){
  if(lanc_arg.N_get == 0) return std::auto_ptr<BFM_Krylov::Lanczos_5d<double> >(NULL);

  BndCndType init_tbc = GJP.Tbc();
  BndCndType target_tbc = time_bc;

  GJP.Bc(3,target_tbc);
  lattice.BondCond();  //Apply BC to internal gauge fields

  bfm_evo<double> &dwf_d = static_cast<Fbfm&>(lattice).bd;
  std::auto_ptr<BFM_Krylov::Lanczos_5d<double> > ret(new BFM_Krylov::Lanczos_5d<double>(dwf_d, const_cast<LancArg&>(lanc_arg)));
  ret->Run();
  if(Fbfm::use_mixed_solver){
    //Convert eigenvectors to single precision
    ret->toSingle();
  }

  //Restore the BCs
  lattice.BondCond();  //unapply BC
  GJP.Bc(3,init_tbc);

  return ret;
}




void writeBasic2ptLW(fMatrix<Rcomplex> &results, const std::string &results_dir, const std::string &descr, const ThreeMomentum &p_psibar, const ThreeMomentum &p_psi, 
		     const PropPrecision status, const TbcStatus &time_bc, const int conf, const std::string &extra_descr = ""){
  std::ostringstream os; 
  os << results_dir << '/' << descr << "_mom" << p_psibar.file_str() << "_plus" << p_psi.file_str() << (status == Sloppy ? "_sloppy" : "_exact") << "_tbc" << time_bc.getTag() << time_bc.getTag() << extra_descr << '.' << conf;
  results.write(os.str());
}


void writePion2ptLW(fMatrix<Rcomplex> &results, const std::string &results_dir, const std::string &snk_op, const ThreeMomentum &p_psibar, const ThreeMomentum &p_psi, 
		    const PropPrecision status, const TbcStatus &time_bc, const int conf, const std::string &extra_descr){
  std::ostringstream os; 
  os << results_dir << "/pion_" << snk_op << "_P_LW_mom" << p_psibar.file_str() << "_plus" << p_psi.file_str() << (status == Sloppy ? "_sloppy" : "_exact") << "_tbc" << time_bc.getTag() << time_bc.getTag() << extra_descr << '.' << conf;
  results.write(os.str());
}

#include "meas_standard.h"
#include "meas_gparity.h"


void measurePion2ptLW(const PropMomContainer &props, const PropPrecision status, const TbcStatus &time_bc, const std::vector<int> &tslices, const MesonMomenta &ll_meson_momenta,
		      const std::string &results_dir, const int conf){
  if(!UniqueID()) printf("Computing pion 2pt LW with %s props\n", status == Sloppy ? "sloppy" : "exact");
  double time = -dclock();

  if(GJP.Gparity()) measurePion2ptLWGparity(props,status,time_bc, tslices,ll_meson_momenta,results_dir, conf);
  else measurePion2ptLWStandard(props,status,time_bc,tslices,ll_meson_momenta,results_dir, conf);

  print_time("main","Pion 2pt LW",time + dclock());
}

void measurePion2ptPPWW(const PropMomContainer &props, const PropPrecision status, const TbcStatus &time_bc, const std::vector<int> &tslices, const MesonMomenta &ll_meson_momenta, Lattice &lat,
			const std::string &results_dir, const int conf){
  if(!UniqueID()) printf("Computing pion 2pt WW with %s props\n", status == Sloppy ? "sloppy" : "exact");
  double time = -dclock();

  if(GJP.Gparity()) measurePion2ptPPWWGparity(props,status,time_bc,tslices,ll_meson_momenta,lat,results_dir,conf);
  else measurePion2ptPPWWStandard(props,status,time_bc,tslices,ll_meson_momenta,lat,results_dir,conf);

  print_time("main","Pion 2pt WW",time + dclock());
}

void measureKaon2ptLW(const PropMomContainer &props, const PropPrecision status, const TbcStatus &time_bc, const std::vector<int> &tslices, const MesonMomenta &meson_momenta,
		      const std::string results_dir, const int conf){
  if(!UniqueID()) printf("Computing kaon 2pt LW with %s props\n", status == Sloppy ? "sloppy" : "exact");
  double time = -dclock();

  if(GJP.Gparity()) measureKaon2ptLWGparity(props,status,time_bc,tslices,meson_momenta,results_dir,conf);
  else measureKaon2ptLWStandard(props,status,time_bc,tslices,meson_momenta,results_dir,conf);

  print_time("main","Kaon 2pt LW",time + dclock());
}

//Kaon 2pt LW functions pseudoscalar sink (cf pion version for comments)
void measureKaon2ptPPWW(const PropMomContainer &props, const PropPrecision status, const TbcStatus &time_bc, const std::vector<int> &tslices, const MesonMomenta &meson_momenta, Lattice &lat,
			const std::string &results_dir, const int conf){

  if(!UniqueID()) printf("Computing pion 2pt WW with %s props\n", status == Sloppy ? "sloppy" : "exact");
  double time = -dclock();

  if(GJP.Gparity()) measureKaon2ptPPWWGparity(props,status,time_bc,tslices,meson_momenta,lat,results_dir,conf);
  else measureKaon2ptPPWWStandard(props,status,time_bc,tslices,meson_momenta,lat,results_dir,conf);

  print_time("main","Kaon 2pt WW",time + dclock());
}

//Measure BK with source kaons on each of the timeslices t0 in prop_tsources and K->K time separations tseps
//Can use standard P or A time BCs but you will need to use closer-together kaon sources to avoid round-the-world effects. These can be eliminated by using the F=P+A and B=P-A combinations
//Either can be specified using the appropriate time_bc parameter below
//For G-parity can optionally choose to disable the source/sink flavor projection (ignored for standard BCs)
void measureBK(const PropMomContainer &props, const PropPrecision status, const std::vector<int> &prop_tsources, const std::vector<int> &tseps, const MesonMomenta &meson_momenta,
	       const TbcStatus &time_bc_t0,
	       const std::string &results_dir, const int conf, const bool do_flavor_project = true){
  if(!UniqueID()) printf("Computing BK with %s props\n", status == Sloppy ? "sloppy" : "exact");
  double time = -dclock();

  if(GJP.Gparity()) measureBKGparity(props,status,prop_tsources,tseps,meson_momenta,time_bc_t0,results_dir,conf,do_flavor_project);
  else measureBKStandard(props,status,prop_tsources,tseps,meson_momenta,time_bc_t0,results_dir,conf);

  print_time("main","BK",time + dclock());
}


//Note: Mres is only properly defined with APRD time BCs. A runtime check is *not* performed
//For G-parity can optionally choose to disable the source/sink flavor projection (ignored for standard BCs)
void measureMres(const PropMomContainer &props, const PropPrecision status, const TbcStatus &time_bc, const std::vector<int> &tslices, const MesonMomenta &meson_momenta,
		 const std::string &results_dir, const int conf, const bool do_flavor_project = true){
  if(!UniqueID()) printf("Computing J5 and J5q (mres) with %s props\n", status == Sloppy ? "sloppy" : "exact");
  double time = -dclock();

  if(GJP.Gparity()) measureMresGparity(props,status,time_bc,tslices,meson_momenta,results_dir,conf,do_flavor_project);
  else measureMresStandard(props,status,time_bc,tslices,meson_momenta,results_dir,conf);

  print_time("main","J5 and J5q",time + dclock());
}







CPS_END_NAMESPACE


#endif
