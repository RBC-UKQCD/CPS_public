#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <util/qcdio.h>
#ifdef PARALLEL
#include <comms/sysfunc_cps.h>
#endif
#include <comms/scu.h>
#include <comms/glb.h>

#include <util/lattice.h>
#include <util/time_cps.h>
#include <util/smalloc.h>

#include <util/command_line.h>

#include<unistd.h>
#include<config.h>

#include <alg/lanc_arg.h>
#include <util/spincolorflavormatrix.h>
#include <alg/propmanager.h>
#include <alg/fix_gauge_arg.h>
#include <alg/alg_gparitycontract.h>
#include <alg/prop_dft.h>

#include <alg/alg_smear.h>
#include <alg/alg_tcharge.h>
#include <alg/alg_wilsonflow.h>
#include <alg/alg_actiondensity.h>

#include <string>
#include <sstream>
#include <bitset>

#include <util/lat_cont.h>

#ifdef USE_BFM 

//CK: these are redefined by BFM (to the same values)
#undef ND
#undef SPINOR_SIZE
#undef HALF_SPINOR_SIZE
#undef GAUGE_SIZE
#endif

#include <util/omp_wrapper.h>

CPS_START_NAMESPACE

void QuarkMomCombination::reset(){ 
  quark_mom.clear(); 
  props.clear();
  allowed_combinations.clear();
  chosen_momcomb = -1;
  allowed_comb_calculated = false;
}
void QuarkMomCombination::add_prop(const QPropWcontainer &prop, const bool &conj){
  std::vector<std::vector<int> > qmom = prop.convert<const QPropWcontainer>().get_allowed_momenta();
  if(qmom.size()==0) ERR.General(cname,"add_prop","Propagator has no allowed momenta!\n");
  if(conj) //complex conjugation of prop inverts the allowed sink momenta
    for(int i=0;i<qmom.size();i++)
      for(int j=0;j<3;j++) qmom[i][j]*=-1;
  quark_mom.push_back(qmom);

  props.push_back(std::pair<QPropWcontainer const*,bool>(&prop,conj));
}
void QuarkMomCombination::calc_allowed_combinations(){
  if(allowed_comb_calculated) return;
  allowed_combinations.clear();
  std::vector<int> start(3,0);
  mom_comb_recurse(allowed_combinations,start,0);
  if(allowed_combinations.size()==0){ ERR.General(cname,"calc_allowed_combinations","Found no momentum combinations!\n"); }
  allowed_comb_calculated = true;
}
void QuarkMomCombination::mom_comb_recurse(std::vector<std::vector<int> > &into, std::vector<int> &cur_branch, const int &depth) const{
  if(depth == quark_mom.size()){
    into.push_back(cur_branch);
  }else{
    for(int poss=0;poss<quark_mom[depth].size();poss++){ //loop over possible momenta for this quark
      const std::vector<int> &choice = quark_mom[depth][poss];
      std::vector<int> new_branch(cur_branch);
      for(int i=0;i<3;i++) new_branch[i]+=choice[i];
      mom_comb_recurse(into,new_branch,depth+1);
    }
  }
}

bool QuarkMomCombination::contains(const std::vector<int> &what){
  calc_allowed_combinations();
  for(int i=0;i<allowed_combinations.size();i++){
    if(what == allowed_combinations[i]) return true;
  }
  return false;
}
bool QuarkMomCombination::contains(const int *what){
  std::vector<int> w(3); for(int i=0;i<3;i++) w[i] = what[i];
  return contains(w);
}
const std::vector<std::vector<int> > &QuarkMomCombination::get_allowed_momenta() const{
  if(!allowed_comb_calculated) ERR.General(cname,"get_allowed_momenta","Momentum combinations have not yet been calculated");
  return allowed_combinations;
}

static void MomCombError(const char* cname, const char* fname, const std::vector<int> &desired_mom, const QuarkMomCombination &momcomb){
  if(!UniqueID()){
    printf("%s::%s Desired momentum (%d,%d,%d) not available. Available combinations are:\n",cname,fname,desired_mom[0],desired_mom[1],desired_mom[2]);
    const std::vector<std::vector<int> > &mc = momcomb.get_allowed_momenta();
    for(int i=0;i<mc.size();i++) printf("%d: (%d,%d,%d)\n",i,mc[i][0],mc[i][1],mc[i][2]);
  }
  fflush(stdout);
  sync();
  exit(-1);
}

void QuarkMomCombination::set_desired_momentum(const std::vector<int> &what){
  calc_allowed_combinations();
  int pos = -1;
  for(int i=0;i<allowed_combinations.size();i++){
    if(what == allowed_combinations[i]){ pos = i; break; }
  }
  if(pos == -1) MomCombError(cname,"set_desired_momentum",what,*this);

  chosen_momcomb = pos;
}
void QuarkMomCombination::set_desired_momentum(const int *what){
  std::vector<int> w(3); for(int i=0;i<3;i++) w[i] = what[i];
  return set_desired_momentum(w);
}

//Get the momentum in lattice units
std::vector<Float> QuarkMomCombination::get_p() const{
  if(chosen_momcomb == -1) ERR.General(cname,"get_p","Must set the desired momentum combination before calling this function");
  
  const static Float Pi_const = M_PI;
  const std::vector<int> &momphase = allowed_combinations[chosen_momcomb];

  std::vector<Float> out(3);
  
  for(int d=0;d<3;d++){
    //NOTE: In G-parity directions, momentum is discretised in odd units of \pi/2L rather than even/odd units of 2\pi/L (periodic/antiperiodic).
    int ntwists = 0; //number of times the angle offset for twisted BCs is applied; depends on quark combination, e.g. for G*G^\dag, ntwists = 0 as the momenta cancel
    if(GJP.Bc(d) == BND_CND_TWISTED || GJP.Bc(d) == BND_CND_GPARITY_TWISTED){
      for(int propidx=0;propidx<props.size();propidx++){
	if(props[propidx].second) --ntwists;
	else ++ntwists;
      }
    }

    Float L = GJP.Nodes(d)*GJP.NodeSites(d);
    Float mom_unit;
    Float twist = 0;
    if(GJP.Bc(d) == BND_CND_GPARITY) mom_unit = Pi_const/(2.0*L);
    else if(GJP.Bc(d) == BND_CND_PRD) mom_unit = 2.0*Pi_const/L;
    else if(GJP.Bc(d) == BND_CND_APRD) mom_unit = Pi_const/L;
    else if(GJP.Bc(d) == BND_CND_TWISTED){ 
      mom_unit = 2.0*Pi_const/L;
      twist = ntwists*GJP.TwistAngle(d)*Pi_const/L;
    }
    else if(GJP.Bc(d) == BND_CND_GPARITY_TWISTED){
      //To reproduce regular G-parity, twist angle = pi
      //In that case BCs only allow quark momenta to be odd-integer multiples of pi/2L
      //The integer momenta given in the regular G-parity case are therefore 1,3,5,7... = pi/2L, 3pi/2L, 5pi/2L, 7pi/2L ...
      //Here the pi/2L offset is fixed by the twisting, so the corresponding integer momenta to those above are 0,1,2,3... (in unit of pi/L) = pi/2L, 3pi/2L, 5pi/2L, 7pi/2L...

      mom_unit = Pi_const/L;
      twist = ntwists*GJP.TwistAngle(d)*Pi_const/(2.0*L);
    }
    else ERR.General(cname,"sink_phasefac(int *,const int &)","Unknown boundary condition\n");
    
    out[d] = momphase[d]*mom_unit + twist;
  }
  return out;
}


Rcomplex QuarkMomCombination::phase_factor(const int *pos) const{
  std::vector<Float> p(get_p());
  Float pdotx = 0.0;
  for(int d=0;d<3;d++) pdotx += p[d]*pos[d];
  return Rcomplex(cos(pdotx),sin(pdotx));
}
Rcomplex QuarkMomCombination::phase_factor(const int &site) const{
  int pos[4];
  int rem = site;
  for(int i=0;i<4;i++){
    pos[i] = rem % GJP.NodeSites(i) + GJP.NodeCoor(i)*GJP.NodeSites(i);
    rem /= GJP.NodeSites(i);
  }
  return phase_factor(pos);
}




void ContractionQuarkMomCombination::add_contraction(const int &contraction, const QuarkMomCombination &cmomcomb){
  momcomb.push_back(cmomcomb);
  contraction_map[contraction] = momcomb.size()-1;
}
void ContractionQuarkMomCombination::same(const int &contraction, const int &same_as_contraction){
  if(!contraction_map.count(same_as_contraction))
    ERR.General("ContractionQuarkMomCombination","same(const int &contraction, const int &same_as_contraction)",
		"same_as_contraction index does not correspond to an existing contraction");
  
  contraction_map[contraction] = contraction_map[same_as_contraction];
}

void ContractionQuarkMomCombination::set_desired_momentum(const std::vector<int> &what){
  for(int i=0;i<momcomb.size();i++) momcomb[i].set_desired_momentum(what);
}
void ContractionQuarkMomCombination::set_desired_momentum(const int *what){
  for(int i=0;i<momcomb.size();i++) momcomb[i].set_desired_momentum(what);
}

Complex ContractionQuarkMomCombination::phase_factor(const int &contraction, const int* pos) const{ 
  if(!contraction_map.count(contraction)) ERR.General("ContractionQuarkMomCombination","phase_factor(const int &contraction, const int* pos)","Unknown contraction idx");
  return  momcomb[ contraction_map.at(contraction) ].phase_factor(pos); 
}
Complex ContractionQuarkMomCombination::phase_factor(const int &contraction, const int& site) const{ 
  if(!contraction_map.count(contraction)) ERR.General("ContractionQuarkMomCombination","phase_factor(const int &contraction, const int& site)","Unknown contraction idx");
  return  momcomb[ contraction_map.at(contraction) ].phase_factor(site); 
}
std::vector<Float> ContractionQuarkMomCombination::get_p(const int &contraction) const{
  if(!contraction_map.count(contraction)) ERR.General("ContractionQuarkMomCombination","get_p(..)","Unknown contraction idx");
  return momcomb[ contraction_map.at(contraction) ].get_p();
}


AlgGparityContract::AlgGparityContract(Lattice & latt, CommonArg& c_arg, GparityContractArg& arg): Alg(latt,&c_arg), args(&arg), binary_write(false){ 
  cname = "AlgGparityContract"; 
}
AlgGparityContract::AlgGparityContract(Lattice & latt, CommonArg& c_arg): Alg(latt,&c_arg), args(NULL), binary_write(false){
  cname = "AlgGparityContract"; 
}

void AlgGparityContract::global_coord(const int &site, int *into_vec){
  int rem = site;
  for(int i=0;i<4;i++){
    into_vec[i] = rem % GJP.NodeSites(i) + GJP.NodeCoor(i)*GJP.NodeSites(i);
    rem /= GJP.NodeSites(i);
  }
}

//Momentum phase. Momenta are in units of 2pi/L
Rcomplex AlgGparityContract::phase_factor(const Float *p, const int* global_pos){
  const static Float pi = M_PI;
  Float pdotx = 0.0;
  for(int d=0;d<3;d++) pdotx += global_pos[d]*p[d]*2*pi/(GJP.NodeSites(d)*GJP.Nodes(d));
  return Rcomplex(cos(pdotx),sin(pdotx));
}

Rcomplex AlgGparityContract::phase_factor(const Float *p, const int &site){
  int pos[4]; global_coord(site,pos);
  return phase_factor(p,pos);
}


void AlgGparityContract::run(const int &conf_idx){
  if(args == NULL) ERR.General(cname,"run(const int)","args pointer has not been set\n");
  run(conf_idx,*args);
}

void AlgGparityContract::run(const int &conf_idx, const GparityContractArg& job){
  //Calculate propagators first. When contracting on only a single thread
  //this is not strictly necessary as the QPropWcontainer will calculate
  //the prop if it has not already been done. However in a multi-threaded
  //inversion, all the threads try to calculate the prop independently, and it will crash.
  PropManager::calcProps(AlgLattice());

  for(int i=0;i<job.meas.meas_len;i++){
    spectrum(job.meas.meas_val[i],conf_idx);   
  }
}

//Left multiply by a gamma matrix structure in QDP-style conventions:
//\Gamma(n) = \gamma_1^n1 \gamma_2^n2  \gamma_3^n3 \gamma_4^n4    where ni are bit fields: n4 n3 n2 n1
void AlgGparityContract::qdp_gl(WilsonMatrix &wmat, const int &gidx){
  std::bitset<4> gmask(gidx);
  if(gmask == std::bitset<4>(15)){ wmat.gl(-5); return; }

  //CPS conventions \gamma^4 = gl(3) etc
  for(int i=3;i>=0;i--) if(gmask[i]) wmat.gl(i);
}
//Right multiply by a gamma matrix structure in QDP-style conventions:
//\Gamma(n) = \gamma_1^n1 \gamma_2^n2  \gamma_3^n3 \gamma_4^n4    where ni are bit fields: n4 n3 n2 n1 
void AlgGparityContract::qdp_gr(WilsonMatrix &wmat, const int &gidx){
  std::bitset<4> gmask(gidx);
  if(gmask == std::bitset<4>(15)){ wmat.gr(-5); return; }

  //CPS conventions \gamma^4 = gl(3) etc
  for(int i=0;i<4;i++) if(gmask[i]) wmat.gr(i);
}

void AlgGparityContract::qdp_gl(SpinColorFlavorMatrix &wmat, const int &gidx){
  std::bitset<4> gmask(gidx);
  if(gmask == std::bitset<4>(15)){ wmat.gl(-5); return; }

  //CPS conventions \gamma^4 = gl(3) etc
  for(int i=3;i>=0;i--) if(gmask[i]) wmat.gl(i);
}
void AlgGparityContract::qdp_gr(SpinColorFlavorMatrix &wmat, const int &gidx){
  std::bitset<4> gmask(gidx);
  if(gmask == std::bitset<4>(15)){ wmat.gr(-5); return; }

  //CPS conventions \gamma^4 = gl(3) etc
  for(int i=0;i<4;i++) if(gmask[i]) wmat.gr(i);
}
//Coefficient when matrix is transposed or conjugated
Float AlgGparityContract::qdp_gcoeff(const int &gidx, const bool &transpose, const bool &conj){
  if(gidx == 0 || gidx==2 || gidx==8 || gidx==15){ return 1.0; } //unit matrix, gamma^2, gamma^4 and gamma^5 are hermitian and real
  else if(gidx==1 || gidx==4){ //gamma^1 and gamma^3 are hermitian and imaginary
    if(transpose && conj) return 1.0;
    else return -1.0;
  }

  std::bitset<4> gmask(gidx);
  Float out(1);
  if(transpose) out *= -1.0; //- sign for reordering 2 or 3 gamma matrices
  
  std::bitset<4> mask(1);
  for(int i=0;i<4;i++, mask<<=1){
    if(gmask[i]){
      Float mult = qdp_gcoeff((int)mask.to_ulong(),transpose,conj);
      out*=mult;
    }
  }
  return out;
}

Float AlgGparityContract::pauli_coeff(const int &pidx, const bool &transpose, const bool &conj){
  if( (transpose && conj) || (!transpose && !conj) ) return 1.0; //all Hermitian
  
  if(pidx == 2 && (transpose || conj)) return -1.0; //sigma_2
  else return 1.0; //all others are real and invariant under transpose
}

void AlgGparityContract::meson_LL_std(QPropWcontainer &prop, const int* sink_mom, const int &gamma_idx_1, const int &gamma_idx_2, FILE *fp){  
  /*Mesons comprising $ \bar u $ and $ d$*/
  std::ostringstream os; os << "LL_MESON " << gamma_idx_1 << " " << gamma_idx_2;
  CorrelationFunction corrfunc(os.str().c_str(),CorrelationFunction::THREADED);

  if(UniqueID()==0) printf("Doing LL_MESON %d %d contraction\n",gamma_idx_1,gamma_idx_2);

  /*Mesons comprising $ \bar u $ and $ d$*/
  /*Require a "CorrelationFunction &corrfunc"*/
  /*Require propagator "QPropWcontainer &prop_src_y_0_pcon corresponding to \mathcal{G}^{(0)}_{x,y}*/
  QPropWcontainer &prop_src_y_0_pcon = prop;

  /*Fourier transform on sink index x*/
  ContractionQuarkMomCombination cmomenta;  
  {
    /*[-1 ]*[{\rm tr}_{sc,0}\left\{S_2 \gamma^5 \mathcal{G}^{(0) \dagger}_{x,y} \gamma^5 S_1 \mathcal{G}^{(0)}_{x,y}\right\}_{0}  ]*/
    QuarkMomCombination momcomb;
    momcomb.add_prop(prop_src_y_0_pcon, false);
    momcomb.add_prop(prop_src_y_0_pcon, true);
    
    cmomenta.add_contraction(0,momcomb);
  }
  /*Require a 3-component array 'desired_mom_x' representing the required momentum at this sink position*/
  const int *desired_mom_x = sink_mom;
  cmomenta.set_desired_momentum(desired_mom_x);


  corrfunc.setNcontractions(1);
#pragma omp parallel for default(shared)
  for(int x=0;x<GJP.VolNodeSites();x++){
    int x_pos_vec[4];
    global_coord(x,x_pos_vec);
  
    /*Get all WilsonMatrices needed*/
    WilsonMatrix prop_snk_x_0_src_y_0_hconj_wmat(prop_src_y_0_pcon.getProp(AlgLattice()).SiteMatrix(x,0));
    prop_snk_x_0_src_y_0_hconj_wmat.hconj();
  
    WilsonMatrix& prop_snk_x_0_src_y_0_wmat = prop_src_y_0_pcon.getProp(AlgLattice()).SiteMatrix(x,0);
  
    /*Starting contraction 0*/
    /*[-1 ]*[{\rm tr}_{sc,0}\left\{S_2 \gamma^5 \mathcal{G}^{(0) \dagger}_{x,y} \gamma^5 S_1 \mathcal{G}^{(0)}_{x,y}\right\}_{0}  ]*/
  
    {
      Rcomplex contraction(-1 , 0);
      contraction *= cmomenta.phase_factor(0,x_pos_vec);
//sink_phasefac(desired_mom_x,x_pos_vec,false);
    
      Rcomplex result_subdiag1(1.0);
      {
	WilsonMatrix sdiag1_trset0_wmat_prod(prop_snk_x_0_src_y_0_hconj_wmat);
	sdiag1_trset0_wmat_prod.gl(-5);
	//sdiag1_trset0_wmat_prod = S_2 * sdiag1_trset0_wmat_prod;
	qdp_gl(sdiag1_trset0_wmat_prod,gamma_idx_2);

	sdiag1_trset0_wmat_prod.gr(-5);
	//sdiag1_trset0_wmat_prod *= S_1;
	qdp_gr(sdiag1_trset0_wmat_prod,gamma_idx_1);

	sdiag1_trset0_wmat_prod *= prop_snk_x_0_src_y_0_wmat;
	Rcomplex sdiag1_trset0_cmplx;
	sdiag1_trset0_cmplx = sdiag1_trset0_wmat_prod.Trace();
	result_subdiag1 *= sdiag1_trset0_cmplx;
      }
      contraction *= result_subdiag1;
    
    
      corrfunc(omp_get_thread_num(),0,x_pos_vec[3]) += contraction;
    }
  }
  corrfunc.write(fp);
}

void AlgGparityContract::meson_LL_gparity(QPropWcontainer &prop, const int* sink_mom, const int &gamma_idx_1, const int &gamma_idx_2, FILE *fp){
  /*Mesons comprising $ \bar u $ and $ d$*/
  /*Require a "CorrelationFunction &corrfunc"*/
  std::ostringstream os; os << "LL_MESON " << gamma_idx_1 << " " << gamma_idx_2;
  CorrelationFunction corrfunc(os.str().c_str(),CorrelationFunction::THREADED);

  if(UniqueID()==0) printf("Doing LL_MESON %d %d contraction\n",gamma_idx_1,gamma_idx_2);

  /*Mesons comprising $ \bar u $ and $ d$*/
  /*Require a "CorrelationFunction &corrfunc"*/
  /*Require propagator "QPropWcontainer &prop_src_y_u_d_eitherflav_pcon corresponding to \mathcal{G}^{[u/d] }_{x,y} with source of either flavour (full prop matrix is generated using single flavour source). Source must be real.*/
  QPropWcontainer &prop_src_y_u_d_eitherflav_pcon = prop;

  /*Fourier transform on sink index x*/
  /*Require a 3-component array 'desired_mom_x' representing the required momentum at this sink position*/

  ContractionQuarkMomCombination cmomenta;  
  {
    QuarkMomCombination momcomb;
    momcomb.add_prop(prop_src_y_u_d_eitherflav_pcon, false);
    momcomb.add_prop(prop_src_y_u_d_eitherflav_pcon, false);

    cmomenta.add_contraction(0,momcomb);
    cmomenta.same(1,0);
  }
  const int *desired_mom_x = sink_mom;
  cmomenta.set_desired_momentum(desired_mom_x);

  corrfunc.setNcontractions(2);
#pragma omp parallel for default(shared)
  for(int x=0;x<GJP.VolNodeSites();x++){
    int x_pos_vec[4];
    global_coord(x,x_pos_vec);
  
    /*Get all SpinColorFlavorMatrices needed*/
    SpinColorFlavorMatrix prop_ud_snk_x_src_y_trans_scfmat(prop_src_y_u_d_eitherflav_pcon , AlgLattice(), x);
    prop_ud_snk_x_src_y_trans_scfmat.transpose();
  
    SpinColorFlavorMatrix prop_ud_snk_x_src_y_scfmat(prop_src_y_u_d_eitherflav_pcon , AlgLattice(), x);
  
    /*Starting contraction 0*/
    /*[1 ]*[{\rm tr}_{scf,0}\left\{C S_2^T F_1 F_\updownarrow \mathcal{G}^{[u/d] T}_{x,y} F_1 F_\updownarrow C S_1 \mathcal{G}^{[u/d] }_{x,y}\right\}_{0}  ]*/
  
    {
      Rcomplex contraction(1 , 0);
      contraction *= cmomenta.phase_factor(0,x_pos_vec);
    
      Rcomplex result_subdiag1(1.0);
      {
	SpinColorFlavorMatrix sdiag1_trset0_scfmat_prod(prop_ud_snk_x_src_y_trans_scfmat);
	sdiag1_trset0_scfmat_prod.pl(Fud);
	sdiag1_trset0_scfmat_prod.pl(F1);
	//sdiag1_trset0_scfmat_prod = S_2.trans() * sdiag1_trset0_scfmat_prod;
	qdp_gl(sdiag1_trset0_scfmat_prod,gamma_idx_2);
	result_subdiag1 *= qdp_gcoeff(gamma_idx_2,true,false);

	sdiag1_trset0_scfmat_prod.ccl(-1);
	sdiag1_trset0_scfmat_prod.pr(F1);
	sdiag1_trset0_scfmat_prod.pr(Fud);
	sdiag1_trset0_scfmat_prod.ccr(1);
	//sdiag1_trset0_scfmat_prod *= S_1;
	qdp_gr(sdiag1_trset0_scfmat_prod,gamma_idx_1);

	sdiag1_trset0_scfmat_prod *= prop_ud_snk_x_src_y_scfmat;
	Rcomplex sdiag1_trset0_cmplx;
	sdiag1_trset0_cmplx = sdiag1_trset0_scfmat_prod.Trace();
	result_subdiag1 *= sdiag1_trset0_cmplx;
      }
      contraction *= result_subdiag1;
    
    
      corrfunc(omp_get_thread_num(),0,x_pos_vec[3]) += contraction;
    }
    /*Starting contraction 1*/
    /*[{\rm tr}_{scf,0}\left\{S_2 C F_0 F_\updownarrow \mathcal{G}^{[u/d] T}_{x,y} F_1 F_\updownarrow C S_1 \mathcal{G}^{[u/d] }_{x,y}\right\}_{0}  ]*[1 ]*/
  
    {
      Rcomplex contraction(1 , 0);
      contraction *= cmomenta.phase_factor(1,x_pos_vec);
    
      Rcomplex result_subdiag0(1.0);
      {
	SpinColorFlavorMatrix sdiag0_trset0_scfmat_prod(prop_ud_snk_x_src_y_trans_scfmat);
	sdiag0_trset0_scfmat_prod.pl(Fud);
	sdiag0_trset0_scfmat_prod.pl(F0);
	sdiag0_trset0_scfmat_prod.ccl(-1);
	//sdiag0_trset0_scfmat_prod = S_2 * sdiag0_trset0_scfmat_prod;
	qdp_gl(sdiag0_trset0_scfmat_prod,gamma_idx_2);

	sdiag0_trset0_scfmat_prod.pr(F1);
	sdiag0_trset0_scfmat_prod.pr(Fud);
	sdiag0_trset0_scfmat_prod.ccr(1);
	//sdiag0_trset0_scfmat_prod *= S_1;
	qdp_gr(sdiag0_trset0_scfmat_prod,gamma_idx_1);

	sdiag0_trset0_scfmat_prod *= prop_ud_snk_x_src_y_scfmat;
	Rcomplex sdiag0_trset0_cmplx;
	sdiag0_trset0_cmplx = sdiag0_trset0_scfmat_prod.Trace();
	result_subdiag0 *= sdiag0_trset0_cmplx;
      }
      contraction *= result_subdiag0;
    
    
      corrfunc(omp_get_thread_num(),1,x_pos_vec[3]) += contraction;
    }
  }
  corrfunc.write(fp);
}

void AlgGparityContract::contract_LL_mesons(const ContractionTypeLLMesons &args, const int &conf_idx){
  std::ostringstream file; file << args.file << "." << conf_idx;

  FILE *fp;
  if ((fp = Fopen(file.str().c_str(), "w")) == NULL) {
    ERR.FileW("CorrelationFunction","write(const char *file)",file.str().c_str());
  }
  PropagatorContainer &prop = PropManager::getProp(args.prop_L);

  if(prop.type() == QPROPW_TYPE){
    //loop through LL meson correlation functions
    if(GJP.Gparity()){
      for(int g1=0;g1<16;g1++){
	for(int g2=0;g2<16;g2++){
	  meson_LL_gparity(prop.convert<QPropWcontainer>(), args.sink_mom, g1, g2, fp);
	}
      }
    }else{
      for(int g1=0;g1<16;g1++){
	for(int g2=0;g2<16;g2++){
	  meson_LL_std(prop.convert<QPropWcontainer>(), args.sink_mom, g1, g2, fp);
	}
      }
    }
  }else ERR.General("AlgGparityContract","contract_LL_mesons(const ContractionTypeLLMesons &args, const int &conf_idx)","Not implemented for types other than QPROPW_TYPE\n");

  Fclose(fp);
}


void AlgGparityContract::meson_HL_gparity(QPropWcontainer &prop_H,QPropWcontainer &prop_L, const int* sink_mom, const int &gamma_idx_1, const int &gamma_idx_2, FILE *fp){
  QPropWcontainer &prop_src_y_u_d_eitherflav_pcon = prop_L;
  QPropWcontainer &prop_src_y_sprime_s_eitherflav_pcon = prop_H;
  const int *desired_mom_x = sink_mom;

  {
    /*<<(\bar s',s)*(\bar s,s')>>*/
    /*Require a "CorrelationFunction &corrfunc" with option "CorrelationFunction::THREADED"*/
    std::ostringstream os; os << "HL_MESON_SPRIME_S_S_SPRIME " << gamma_idx_1 << " " << gamma_idx_2;
    CorrelationFunction corrfunc(os.str().c_str(),CorrelationFunction::THREADED);

    if(UniqueID()==0) printf("Doing HL_MESON_SPRIME_S_S_SPRIME %d %d contraction\n",gamma_idx_1,gamma_idx_2);

    /*Require propagator "QPropWcontainer &prop_src_y_sprime_s_eitherflav_pcon corresponding to \mathcal{G}^{[s^\prime/s] }_{x,y} with source of either flavour (full prop matrix is generated using single flavour source). Source must be real.*/

    /*Fourier transform on sink index x*/
    /*Require a 3-component array 'desired_mom_x' representing the required momentum at this sink position*/
    ContractionQuarkMomCombination cmomenta;      
    {
      QuarkMomCombination momcomb;
      momcomb.add_prop(prop_src_y_sprime_s_eitherflav_pcon, false);
      momcomb.add_prop(prop_src_y_sprime_s_eitherflav_pcon, false);

      cmomenta.add_contraction(0,momcomb);
      cmomenta.same(1,0);
    }
    cmomenta.set_desired_momentum(desired_mom_x);

    corrfunc.setNcontractions(2);
#pragma omp parallel for default(shared)
    for(int x=0;x<GJP.VolNodeSites();x++){
      int x_pos_vec[4];
      global_coord(x,x_pos_vec);
  
      /*Get all SpinColorFlavorMatrices needed*/
      SpinColorFlavorMatrix prop_sprimes_snk_x_src_y_scfmat(prop_src_y_sprime_s_eitherflav_pcon , AlgLattice(), x);
  
      SpinColorFlavorMatrix prop_sprimes_snk_x_src_y_trans_scfmat(prop_src_y_sprime_s_eitherflav_pcon , AlgLattice(), x);
      prop_sprimes_snk_x_src_y_trans_scfmat.transpose();
  
      /*Starting contraction 0*/
      /*[{\rm tr}_{scf,0}\left\{C \Gamma[g2] F_1 F_\updownarrow \mathcal{G}^{[s^\prime/s] T}_{x,y} F_1 F_\updownarrow C \Gamma[g1] \mathcal{G}^{[s^\prime/s] }_{x,y}\right\}_{0}  ]*[f_\Gamma(g2,T) ]*/
  
      {
	Rcomplex contraction(1 , 0);
	contraction *= qdp_gcoeff(gamma_idx_2,true,false);
	contraction *= cmomenta.phase_factor(0,x_pos_vec);
    
	Rcomplex result_subdiag0(1.0);
	{
	  SpinColorFlavorMatrix sdiag0_trset0_scfmat_prod(prop_sprimes_snk_x_src_y_trans_scfmat);
	  sdiag0_trset0_scfmat_prod.pl(Fud);
	  sdiag0_trset0_scfmat_prod.pl(F1);
	  qdp_gl(sdiag0_trset0_scfmat_prod,gamma_idx_2);
	  sdiag0_trset0_scfmat_prod.ccl(-1);
	  sdiag0_trset0_scfmat_prod.pr(F1);
	  sdiag0_trset0_scfmat_prod.pr(Fud);
	  sdiag0_trset0_scfmat_prod.ccr(1);
	  qdp_gr(sdiag0_trset0_scfmat_prod,gamma_idx_1);
	  sdiag0_trset0_scfmat_prod *= prop_sprimes_snk_x_src_y_scfmat;
	  Rcomplex sdiag0_trset0_cmplx;
	  sdiag0_trset0_cmplx = sdiag0_trset0_scfmat_prod.Trace();
	  result_subdiag0 *= sdiag0_trset0_cmplx;
	}
	contraction *= result_subdiag0;
    
    
	corrfunc(omp_get_thread_num(),0,x_pos_vec[3]) += contraction;
      }
      /*Starting contraction 1*/
      /*[{\rm tr}_{scf,0}\left\{\Gamma[g2] C F_0 F_\updownarrow \mathcal{G}^{[s^\prime/s] T}_{x,y} F_1 F_\updownarrow C \Gamma[g1] \mathcal{G}^{[s^\prime/s] }_{x,y}\right\}_{0}  ]*/
  
      {
	Rcomplex contraction(1 , 0);
	contraction *= cmomenta.phase_factor(1,x_pos_vec);
    
	Rcomplex result_subdiag0(1.0);
	{
	  SpinColorFlavorMatrix sdiag0_trset0_scfmat_prod(prop_sprimes_snk_x_src_y_trans_scfmat);
	  sdiag0_trset0_scfmat_prod.pl(Fud);
	  sdiag0_trset0_scfmat_prod.pl(F0);
	  sdiag0_trset0_scfmat_prod.ccl(-1);
	  qdp_gl(sdiag0_trset0_scfmat_prod,gamma_idx_2);
	  sdiag0_trset0_scfmat_prod.pr(F1);
	  sdiag0_trset0_scfmat_prod.pr(Fud);
	  sdiag0_trset0_scfmat_prod.ccr(1);
	  qdp_gr(sdiag0_trset0_scfmat_prod,gamma_idx_1);
	  sdiag0_trset0_scfmat_prod *= prop_sprimes_snk_x_src_y_scfmat;
	  Rcomplex sdiag0_trset0_cmplx;
	  sdiag0_trset0_cmplx = sdiag0_trset0_scfmat_prod.Trace();
	  result_subdiag0 *= sdiag0_trset0_cmplx;
	}
	contraction *= result_subdiag0;
    
    
	corrfunc(omp_get_thread_num(),1,x_pos_vec[3]) += contraction;
      }
    }

    corrfunc.write(fp);
  }


  {
    /*<<(\bar d,s)*(\bar s,d)>>*/
    /*Require a "CorrelationFunction &corrfunc" with option "CorrelationFunction::THREADED"*/
    std::ostringstream os; os << "HL_MESON_D_S_S_D " << gamma_idx_1 << " " << gamma_idx_2;
    CorrelationFunction corrfunc(os.str().c_str(),CorrelationFunction::THREADED);

    if(UniqueID()==0) printf("Doing HL_MESON_D_S_S_D %d %d contraction\n",gamma_idx_1,gamma_idx_2);

    /*Require propagator "QPropWcontainer &prop_src_y_u_d_eitherflav_pcon corresponding to \mathcal{G}^{[u/d] }_{x,y} with source of either flavour (full prop matrix is generated using single flavour source). Source must be real.*/
    /*Require propagator "QPropWcontainer &prop_src_y_sprime_s_eitherflav_pcon corresponding to \mathcal{G}^{[s^\prime/s] }_{x,y} with source of either flavour (full prop matrix is generated using single flavour source). Source must be real.*/

    /*Fourier transform on sink index x*/
    /*Require a 3-component array 'desired_mom_x' representing the required momentum at this sink position*/
    ContractionQuarkMomCombination cmomenta;
    {
      QuarkMomCombination momcomb;
      momcomb.add_prop(prop_src_y_u_d_eitherflav_pcon, true);
      momcomb.add_prop(prop_src_y_sprime_s_eitherflav_pcon, false);
      cmomenta.add_contraction(0,momcomb);
    }
    cmomenta.set_desired_momentum(desired_mom_x);

    corrfunc.setNcontractions(1);
#pragma omp parallel for default(shared)
    for(int x=0;x<GJP.VolNodeSites();x++){
      int x_pos_vec[4];
      global_coord(x,x_pos_vec);
  
      /*Get all SpinColorFlavorMatrices needed*/
      SpinColorFlavorMatrix prop_ud_snk_x_src_y_hconj_scfmat(prop_src_y_u_d_eitherflav_pcon , AlgLattice(), x);
      prop_ud_snk_x_src_y_hconj_scfmat.hconj();
  
      SpinColorFlavorMatrix prop_sprimes_snk_x_src_y_scfmat(prop_src_y_sprime_s_eitherflav_pcon , AlgLattice(), x);
  
      /*Starting contraction 0*/
      /*[-1 ]*[{\rm tr}_{scf,0}\left\{\gamma^5 \Gamma[g1] F_0 \mathcal{G}^{[s^\prime/s] }_{x,y} F_0 \Gamma[g2] \gamma^5 \mathcal{G}^{[u/d] \dagger}_{x,y}\right\}_{0}  ]*/
  
      {
	Rcomplex contraction(-1 , 0);
	contraction *= cmomenta.phase_factor(0,x_pos_vec);
    
	Rcomplex result_subdiag1(1.0);
	{
	  SpinColorFlavorMatrix sdiag1_trset0_scfmat_prod(prop_sprimes_snk_x_src_y_scfmat);
	  sdiag1_trset0_scfmat_prod.pl(F0);
	  qdp_gl(sdiag1_trset0_scfmat_prod,gamma_idx_1);
	  sdiag1_trset0_scfmat_prod.gl(-5);
	  sdiag1_trset0_scfmat_prod.pr(F0);
	  qdp_gr(sdiag1_trset0_scfmat_prod,gamma_idx_2);
	  sdiag1_trset0_scfmat_prod.gr(-5);
	  sdiag1_trset0_scfmat_prod *= prop_ud_snk_x_src_y_hconj_scfmat;
	  Rcomplex sdiag1_trset0_cmplx;
	  sdiag1_trset0_cmplx = sdiag1_trset0_scfmat_prod.Trace();
	  result_subdiag1 *= sdiag1_trset0_cmplx;
	}
	contraction *= result_subdiag1;
    
    
	corrfunc(omp_get_thread_num(),0,x_pos_vec[3]) += contraction;
      }
    }
    corrfunc.write(fp);
  }

  {
    /*<<(\bar d,s)*(\bar u,s')>>*/
    /*Require a "CorrelationFunction &corrfunc" with option "CorrelationFunction::THREADED"*/
    std::ostringstream os; os << "HL_MESON_D_S_U_SPRIME " << gamma_idx_1 << " " << gamma_idx_2;
    CorrelationFunction corrfunc(os.str().c_str(),CorrelationFunction::THREADED);

    if(UniqueID()==0) printf("Doing HL_MESON_D_S_U_SPRIME %d %d contraction\n",gamma_idx_1,gamma_idx_2);

    /*Require propagator "QPropWcontainer &prop_src_y_u_d_eitherflav_pcon corresponding to \mathcal{G}^{[u/d] }_{x,y} with source of either flavour (full prop matrix is generated using single flavour source). Source must be real.*/
    /*Require propagator "QPropWcontainer &prop_src_y_sprime_s_eitherflav_pcon corresponding to \mathcal{G}^{[s^\prime/s] }_{x,y} with source of either flavour (full prop matrix is generated using single flavour source). Source must be real.*/

    /*Fourier transform on sink index x*/
    /*Require a 3-component array 'desired_mom_x' representing the required momentum at this sink position*/
    ContractionQuarkMomCombination cmomenta;
    {
      QuarkMomCombination momcomb;
      momcomb.add_prop(prop_src_y_u_d_eitherflav_pcon, true);
      momcomb.add_prop(prop_src_y_sprime_s_eitherflav_pcon, false);
      cmomenta.add_contraction(0,momcomb);
    }
    cmomenta.set_desired_momentum(desired_mom_x);

    corrfunc.setNcontractions(1);
#pragma omp parallel for default(shared)
    for(int x=0;x<GJP.VolNodeSites();x++){
      int x_pos_vec[4];
      global_coord(x,x_pos_vec);
  
      /*Get all SpinColorFlavorMatrices needed*/
      SpinColorFlavorMatrix prop_ud_snk_x_src_y_hconj_scfmat(prop_src_y_u_d_eitherflav_pcon , AlgLattice(), x);
      prop_ud_snk_x_src_y_hconj_scfmat.hconj();
  
      SpinColorFlavorMatrix prop_sprimes_snk_x_src_y_scfmat(prop_src_y_sprime_s_eitherflav_pcon , AlgLattice(), x);
  
      /*Starting contraction 0*/
      /*[{\rm tr}_{scf,0}\left\{\gamma^5 \Gamma[g1] F_0 \mathcal{G}^{[s^\prime/s] }_{x,y} F_1 C \Gamma[g2] \gamma^5 C \mathcal{G}^{[u/d] \dagger}_{x,y}\right\}_{0}  ]*[f_\Gamma(g2,T) ]*/
  
      {
	Rcomplex contraction(1 , 0);
	contraction *= qdp_gcoeff(gamma_idx_2,true,false);
	contraction *= cmomenta.phase_factor(0,x_pos_vec);
    
	Rcomplex result_subdiag0(1.0);
	{
	  SpinColorFlavorMatrix sdiag0_trset0_scfmat_prod(prop_sprimes_snk_x_src_y_scfmat);
	  sdiag0_trset0_scfmat_prod.pl(F0);
	  qdp_gl(sdiag0_trset0_scfmat_prod,gamma_idx_1);
	  sdiag0_trset0_scfmat_prod.gl(-5);
	  sdiag0_trset0_scfmat_prod.pr(F1);
	  sdiag0_trset0_scfmat_prod.ccr(1);
	  qdp_gr(sdiag0_trset0_scfmat_prod,gamma_idx_2);
	  sdiag0_trset0_scfmat_prod.gr(-5);
	  sdiag0_trset0_scfmat_prod.ccr(1);
	  sdiag0_trset0_scfmat_prod *= prop_ud_snk_x_src_y_hconj_scfmat;
	  Rcomplex sdiag0_trset0_cmplx;
	  sdiag0_trset0_cmplx = sdiag0_trset0_scfmat_prod.Trace();
	  result_subdiag0 *= sdiag0_trset0_cmplx;
	}
	contraction *= result_subdiag0;
    
    
	corrfunc(omp_get_thread_num(),0,x_pos_vec[3]) += contraction;
      }
    }
    corrfunc.write(fp);

  }

  {
    /*<<(\bar d,s')*(\bar s',d)>>*/
    /*Require a "CorrelationFunction &corrfunc" with option "CorrelationFunction::THREADED"*/
    std::ostringstream os; os << "HL_MESON_D_SPRIME_SPRIME_D " << gamma_idx_1 << " " << gamma_idx_2;
    CorrelationFunction corrfunc(os.str().c_str(),CorrelationFunction::THREADED);

    if(UniqueID()==0) printf("Doing HL_MESON_D_SPRIME_SPRIME_D %d %d contraction\n",gamma_idx_1,gamma_idx_2);

    /*Require propagator "QPropWcontainer &prop_src_y_u_d_eitherflav_pcon corresponding to \mathcal{G}^{[u/d] }_{x,y} with source of either flavour (full prop matrix is generated using single flavour source). Source must be real.*/
    /*Require propagator "QPropWcontainer &prop_src_y_sprime_s_eitherflav_pcon corresponding to \mathcal{G}^{[s^\prime/s] }_{x,y} with source of either flavour (full prop matrix is generated using single flavour source). Source must be real.*/

    /*Fourier transform on sink index x*/
    /*Require a 3-component array 'desired_mom_x' representing the required momentum at this sink position*/
    ContractionQuarkMomCombination cmomenta;
    {
      QuarkMomCombination momcomb;
      momcomb.add_prop(prop_src_y_u_d_eitherflav_pcon, true);
      momcomb.add_prop(prop_src_y_sprime_s_eitherflav_pcon, true);
      cmomenta.add_contraction(0,momcomb);
    }
    cmomenta.set_desired_momentum(desired_mom_x);

    corrfunc.setNcontractions(1);
#pragma omp parallel for default(shared)
    for(int x=0;x<GJP.VolNodeSites();x++){
      int x_pos_vec[4];
      global_coord(x,x_pos_vec);
  
      /*Get all SpinColorFlavorMatrices needed*/
      SpinColorFlavorMatrix prop_sprimes_snk_x_src_y_cconj_scfmat(prop_src_y_sprime_s_eitherflav_pcon , AlgLattice(), x);
      prop_sprimes_snk_x_src_y_cconj_scfmat.cconj();
  
      SpinColorFlavorMatrix prop_ud_snk_x_src_y_hconj_scfmat(prop_src_y_u_d_eitherflav_pcon , AlgLattice(), x);
      prop_ud_snk_x_src_y_hconj_scfmat.hconj();
  
      /*Starting contraction 0*/
      /*[{\rm tr}_{scf,0}\left\{\gamma^5 \Gamma[g1] \gamma^5 C F_0 F_\updownarrow \mathcal{G}^{[s^\prime/s] *}_{x,y} F_1 F_\updownarrow \gamma^5 C \Gamma[g2] \gamma^5 \mathcal{G}^{[u/d] \dagger}_{x,y}\right\}_{0}  ]*/
  
      {
	Rcomplex contraction(1 , 0);
	contraction *= cmomenta.phase_factor(0,x_pos_vec);
    
	Rcomplex result_subdiag0(1.0);
	{
	  SpinColorFlavorMatrix sdiag0_trset0_scfmat_prod(prop_sprimes_snk_x_src_y_cconj_scfmat);
	  sdiag0_trset0_scfmat_prod.pl(Fud);
	  sdiag0_trset0_scfmat_prod.pl(F0);
	  sdiag0_trset0_scfmat_prod.ccl(-1);
	  sdiag0_trset0_scfmat_prod.gl(-5);
	  qdp_gl(sdiag0_trset0_scfmat_prod,gamma_idx_1);
	  sdiag0_trset0_scfmat_prod.gl(-5);
	  sdiag0_trset0_scfmat_prod.pr(F1);
	  sdiag0_trset0_scfmat_prod.pr(Fud);
	  sdiag0_trset0_scfmat_prod.gr(-5);
	  sdiag0_trset0_scfmat_prod.ccr(1);
	  qdp_gr(sdiag0_trset0_scfmat_prod,gamma_idx_2);
	  sdiag0_trset0_scfmat_prod.gr(-5);
	  sdiag0_trset0_scfmat_prod *= prop_ud_snk_x_src_y_hconj_scfmat;
	  Rcomplex sdiag0_trset0_cmplx;
	  sdiag0_trset0_cmplx = sdiag0_trset0_scfmat_prod.Trace();
	  result_subdiag0 *= sdiag0_trset0_cmplx;
	}
	contraction *= result_subdiag0;
    
    
	corrfunc(omp_get_thread_num(),0,x_pos_vec[3]) += contraction;
      }
    }
    
    corrfunc.write(fp);
  }

  {
    /*<<(\bar d,s')*(\bar u,s)>>*/
    /*Require a "CorrelationFunction &corrfunc" with option "CorrelationFunction::THREADED"*/
    /*Require propagator "QPropWcontainer &prop_src_y_u_d_eitherflav_pcon corresponding to \mathcal{G}^{[u/d] }_{x,y} with source of either flavour (full prop matrix is generated using single flavour source). Source must be real.*/    
    std::ostringstream os; os << "HL_MESON_D_SPRIME_U_S " << gamma_idx_1 << " " << gamma_idx_2;
    CorrelationFunction corrfunc(os.str().c_str(),CorrelationFunction::THREADED);

    if(UniqueID()==0) printf("Doing HL_MESON_D_SPRIME_U_S %d %d contraction\n",gamma_idx_1,gamma_idx_2);

    /*Require propagator "QPropWcontainer &prop_src_y_sprime_s_eitherflav_pcon corresponding to \mathcal{G}^{[s^\prime/s] }_{x,y} with source of either flavour (full prop matrix is generated using single flavour source). Source must be real.*/

    /*Fourier transform on sink index x*/
    /*Require a 3-component array 'desired_mom_x' representing the required momentum at this sink position*/
    ContractionQuarkMomCombination cmomenta;
    {
      QuarkMomCombination momcomb;
      momcomb.add_prop(prop_src_y_u_d_eitherflav_pcon, true);
      momcomb.add_prop(prop_src_y_sprime_s_eitherflav_pcon, true);
      cmomenta.add_contraction(0,momcomb);
    }
    cmomenta.set_desired_momentum(desired_mom_x);

    corrfunc.setNcontractions(1);
#pragma omp parallel for default(shared)
    for(int x=0;x<GJP.VolNodeSites();x++){
      int x_pos_vec[4];
      global_coord(x,x_pos_vec);
  
      /*Get all SpinColorFlavorMatrices needed*/
      SpinColorFlavorMatrix prop_ud_snk_x_src_y_hconj_scfmat(prop_src_y_u_d_eitherflav_pcon , AlgLattice(), x);
      prop_ud_snk_x_src_y_hconj_scfmat.hconj();
  
      SpinColorFlavorMatrix prop_sprimes_snk_x_src_y_cconj_scfmat(prop_src_y_sprime_s_eitherflav_pcon , AlgLattice(), x);
      prop_sprimes_snk_x_src_y_cconj_scfmat.cconj();
  
      /*Starting contraction 0*/
      /*[{\rm tr}_{scf,0}\left\{\gamma^5 \Gamma[g1] \gamma^5 C F_0 F_\updownarrow \mathcal{G}^{[s^\prime/s] *}_{x,y} F_0 F_\updownarrow \gamma^5 \Gamma[g2] \gamma^5 C \mathcal{G}^{[u/d] \dagger}_{x,y}\right\}_{0}  ]*[f_\Gamma(g2,T) ]*/
  
      {
	Rcomplex contraction(1 , 0);
	contraction *= qdp_gcoeff(gamma_idx_2,true,false);
	contraction *= cmomenta.phase_factor(0,x_pos_vec);
    
	Rcomplex result_subdiag0(1.0);
	{
	  SpinColorFlavorMatrix sdiag0_trset0_scfmat_prod(prop_sprimes_snk_x_src_y_cconj_scfmat);
	  sdiag0_trset0_scfmat_prod.pl(Fud);
	  sdiag0_trset0_scfmat_prod.pl(F0);
	  sdiag0_trset0_scfmat_prod.ccl(-1);
	  sdiag0_trset0_scfmat_prod.gl(-5);
	  qdp_gl(sdiag0_trset0_scfmat_prod,gamma_idx_1);
	  sdiag0_trset0_scfmat_prod.gl(-5);
	  sdiag0_trset0_scfmat_prod.pr(F0);
	  sdiag0_trset0_scfmat_prod.pr(Fud);
	  sdiag0_trset0_scfmat_prod.gr(-5);
	  qdp_gr(sdiag0_trset0_scfmat_prod,gamma_idx_2);
	  sdiag0_trset0_scfmat_prod.gr(-5);
	  sdiag0_trset0_scfmat_prod.ccr(1);
	  sdiag0_trset0_scfmat_prod *= prop_ud_snk_x_src_y_hconj_scfmat;
	  Rcomplex sdiag0_trset0_cmplx;
	  sdiag0_trset0_cmplx = sdiag0_trset0_scfmat_prod.Trace();
	  result_subdiag0 *= sdiag0_trset0_cmplx;
	}
	contraction *= result_subdiag0;
    
    
	corrfunc(omp_get_thread_num(),0,x_pos_vec[3]) += contraction;
      }
    }
    corrfunc.write(fp);
  }

  {
    /*<<(\bar u,s)*(\bar s,u)>>*/
    /*Require a "CorrelationFunction &corrfunc" with option "CorrelationFunction::THREADED"*/
    std::ostringstream os; os << "HL_MESON_U_S_S_U " << gamma_idx_1 << " " << gamma_idx_2;
    CorrelationFunction corrfunc(os.str().c_str(),CorrelationFunction::THREADED);

    if(UniqueID()==0) printf("Doing HL_MESON_U_S_S_U %d %d contraction\n",gamma_idx_1,gamma_idx_2);

    /*Require propagator "QPropWcontainer &prop_src_y_u_d_eitherflav_pcon corresponding to \mathcal{G}^{[u/d] }_{x,y} with source of either flavour (full prop matrix is generated using single flavour source). Source must be real.*/
    /*Require propagator "QPropWcontainer &prop_src_y_sprime_s_eitherflav_pcon corresponding to \mathcal{G}^{[s^\prime/s] }_{x,y} with source of either flavour (full prop matrix is generated using single flavour source). Source must be real.*/

    /*Fourier transform on sink index x*/
    /*Require a 3-component array 'desired_mom_x' representing the required momentum at this sink position*/
    ContractionQuarkMomCombination cmomenta;
    {
      QuarkMomCombination momcomb;
      momcomb.add_prop(prop_src_y_u_d_eitherflav_pcon, false);
      momcomb.add_prop(prop_src_y_sprime_s_eitherflav_pcon, false);
      cmomenta.add_contraction(0,momcomb);
    }
    cmomenta.set_desired_momentum(desired_mom_x);

    corrfunc.setNcontractions(1);
#pragma omp parallel for default(shared)
    for(int x=0;x<GJP.VolNodeSites();x++){
      int x_pos_vec[4];
      global_coord(x,x_pos_vec);
  
      /*Get all SpinColorFlavorMatrices needed*/
      SpinColorFlavorMatrix prop_ud_snk_x_src_y_scfmat(prop_src_y_u_d_eitherflav_pcon , AlgLattice(), x);
  
      SpinColorFlavorMatrix prop_sprimes_snk_x_src_y_trans_scfmat(prop_src_y_sprime_s_eitherflav_pcon , AlgLattice(), x);
      prop_sprimes_snk_x_src_y_trans_scfmat.transpose();
  
      /*Starting contraction 0*/
      /*[{\rm tr}_{scf,0}\left\{C \Gamma[g2] F_1 F_\updownarrow \mathcal{G}^{[s^\prime/s] T}_{x,y} F_0 F_\updownarrow \Gamma[g1] C \mathcal{G}^{[u/d] }_{x,y}\right\}_{0}  ]*[f_\Gamma(g2,T) ]*[f_\Gamma(g1,T) ]*/
  
      {
	Rcomplex contraction(1 , 0);
	contraction *= qdp_gcoeff(gamma_idx_2,true,false)*qdp_gcoeff(gamma_idx_1,true,false);
	contraction *= cmomenta.phase_factor(0,x_pos_vec);
    
	Rcomplex result_subdiag0(1.0);
	{
	  SpinColorFlavorMatrix sdiag0_trset0_scfmat_prod(prop_sprimes_snk_x_src_y_trans_scfmat);
	  sdiag0_trset0_scfmat_prod.pl(Fud);
	  sdiag0_trset0_scfmat_prod.pl(F1);
	  qdp_gl(sdiag0_trset0_scfmat_prod,gamma_idx_2);
	  sdiag0_trset0_scfmat_prod.ccl(-1);
	  sdiag0_trset0_scfmat_prod.pr(F0);
	  sdiag0_trset0_scfmat_prod.pr(Fud);
	  qdp_gr(sdiag0_trset0_scfmat_prod,gamma_idx_1);
	  sdiag0_trset0_scfmat_prod.ccr(1);
	  sdiag0_trset0_scfmat_prod *= prop_ud_snk_x_src_y_scfmat;
	  Rcomplex sdiag0_trset0_cmplx;
	  sdiag0_trset0_cmplx = sdiag0_trset0_scfmat_prod.Trace();
	  result_subdiag0 *= sdiag0_trset0_cmplx;
	}
	contraction *= result_subdiag0;
    
    
	corrfunc(omp_get_thread_num(),0,x_pos_vec[3]) += contraction;
      }
    }
    corrfunc.write(fp);
  }


  {
    /*<<(\bar u,s)*(\bar d,s')>>*/
    /*Require a "CorrelationFunction &corrfunc" with option "CorrelationFunction::THREADED"*/
    std::ostringstream os; os << "HL_MESON_U_S_D_SPRIME " << gamma_idx_1 << " " << gamma_idx_2;
    CorrelationFunction corrfunc(os.str().c_str(),CorrelationFunction::THREADED);

    if(UniqueID()==0) printf("Doing HL_MESON_U_S_D_SPRIME %d %d contraction\n",gamma_idx_1,gamma_idx_2);

    /*Require propagator "QPropWcontainer &prop_src_y_u_d_eitherflav_pcon corresponding to \mathcal{G}^{[u/d] }_{x,y} with source of either flavour (full prop matrix is generated using single flavour source). Source must be real.*/
    /*Require propagator "QPropWcontainer &prop_src_y_sprime_s_eitherflav_pcon corresponding to \mathcal{G}^{[s^\prime/s] }_{x,y} with source of either flavour (full prop matrix is generated using single flavour source). Source must be real.*/

    /*Fourier transform on sink index x*/
    /*Require a 3-component array 'desired_mom_x' representing the required momentum at this sink position*/
    ContractionQuarkMomCombination cmomenta;
    {
      QuarkMomCombination momcomb;
      momcomb.add_prop(prop_src_y_u_d_eitherflav_pcon, false);
      momcomb.add_prop(prop_src_y_sprime_s_eitherflav_pcon, false);
      cmomenta.add_contraction(0,momcomb);
    }
    cmomenta.set_desired_momentum(desired_mom_x);

    corrfunc.setNcontractions(1);
#pragma omp parallel for default(shared)
    for(int x=0;x<GJP.VolNodeSites();x++){
      int x_pos_vec[4];
      global_coord(x,x_pos_vec);
  
      /*Get all SpinColorFlavorMatrices needed*/
      SpinColorFlavorMatrix prop_ud_snk_x_src_y_scfmat(prop_src_y_u_d_eitherflav_pcon , AlgLattice(), x);
  
      SpinColorFlavorMatrix prop_sprimes_snk_x_src_y_trans_scfmat(prop_src_y_sprime_s_eitherflav_pcon , AlgLattice(), x);
      prop_sprimes_snk_x_src_y_trans_scfmat.transpose();
  
      /*Starting contraction 0*/
      /*[{\rm tr}_{scf,0}\left\{\Gamma[g2] C F_0 F_\updownarrow \mathcal{G}^{[s^\prime/s] T}_{x,y} F_0 F_\updownarrow \Gamma[g1] C \mathcal{G}^{[u/d] }_{x,y}\right\}_{0}  ]*[f_\Gamma(g1,T) ]*/
  
      {
	Rcomplex contraction(1 , 0);
	contraction *= qdp_gcoeff(gamma_idx_1,true,false);
	contraction *= cmomenta.phase_factor(0,x_pos_vec);
    
	Rcomplex result_subdiag0(1.0);
	{
	  SpinColorFlavorMatrix sdiag0_trset0_scfmat_prod(prop_sprimes_snk_x_src_y_trans_scfmat);
	  sdiag0_trset0_scfmat_prod.pl(Fud);
	  sdiag0_trset0_scfmat_prod.pl(F0);
	  sdiag0_trset0_scfmat_prod.ccl(-1);
	  qdp_gl(sdiag0_trset0_scfmat_prod,gamma_idx_2);
	  sdiag0_trset0_scfmat_prod.pr(F0);
	  sdiag0_trset0_scfmat_prod.pr(Fud);
	  qdp_gr(sdiag0_trset0_scfmat_prod,gamma_idx_1);
	  sdiag0_trset0_scfmat_prod.ccr(1);
	  sdiag0_trset0_scfmat_prod *= prop_ud_snk_x_src_y_scfmat;
	  Rcomplex sdiag0_trset0_cmplx;
	  sdiag0_trset0_cmplx = sdiag0_trset0_scfmat_prod.Trace();
	  result_subdiag0 *= sdiag0_trset0_cmplx;
	}
	contraction *= result_subdiag0;
    
    
	corrfunc(omp_get_thread_num(),0,x_pos_vec[3]) += contraction;
      }
    }


    corrfunc.write(fp);
  }


  {
    /*<<(\bar u,s')*(\bar s',u)>>*/
    /*Require a "CorrelationFunction &corrfunc" with option "CorrelationFunction::THREADED"*/
    std::ostringstream os; os << "HL_MESON_U_SPRIME_SPRIME_U " << gamma_idx_1 << " " << gamma_idx_2;
    CorrelationFunction corrfunc(os.str().c_str(),CorrelationFunction::THREADED);

    if(UniqueID()==0) printf("Doing HL_MESON_U_SPRIME_SPRIME_U %d %d contraction\n",gamma_idx_1,gamma_idx_2);

    /*Require propagator "QPropWcontainer &prop_src_y_u_d_eitherflav_pcon corresponding to \mathcal{G}^{[u/d] }_{x,y} with source of either flavour (full prop matrix is generated using single flavour source). Source must be real.*/
    /*Require propagator "QPropWcontainer &prop_src_y_sprime_s_eitherflav_pcon corresponding to \mathcal{G}^{[s^\prime/s] }_{x,y} with source of either flavour (full prop matrix is generated using single flavour source). Source must be real.*/

    /*Fourier transform on sink index x*/
    /*Require a 3-component array 'desired_mom_x' representing the required momentum at this sink position*/
    ContractionQuarkMomCombination cmomenta;
    {
      QuarkMomCombination momcomb;
      momcomb.add_prop(prop_src_y_u_d_eitherflav_pcon, false);
      momcomb.add_prop(prop_src_y_sprime_s_eitherflav_pcon, true);
      cmomenta.add_contraction(0,momcomb);
    }
    cmomenta.set_desired_momentum(desired_mom_x);

    corrfunc.setNcontractions(1);
#pragma omp parallel for default(shared)
    for(int x=0;x<GJP.VolNodeSites();x++){
      int x_pos_vec[4];
      global_coord(x,x_pos_vec);
  
      /*Get all SpinColorFlavorMatrices needed*/
      SpinColorFlavorMatrix prop_ud_snk_x_src_y_scfmat(prop_src_y_u_d_eitherflav_pcon , AlgLattice(), x);
  
      SpinColorFlavorMatrix prop_sprimes_snk_x_src_y_hconj_scfmat(prop_src_y_sprime_s_eitherflav_pcon , AlgLattice(), x);
      prop_sprimes_snk_x_src_y_hconj_scfmat.hconj();
  
      /*Starting contraction 0*/
      /*[-1 ]*[{\rm tr}_{scf,0}\left\{C \Gamma[g2] \gamma^5 C F_1 \mathcal{G}^{[s^\prime/s] \dagger}_{x,y} F_1 \gamma^5 C \Gamma[g1] C \mathcal{G}^{[u/d] }_{x,y}\right\}_{0}  ]*[f_\Gamma(g2,T) ]*[f_\Gamma(g1,T) ]*/
  
      {
	Rcomplex contraction(-1 , 0);
	contraction *= qdp_gcoeff(gamma_idx_2,true,false)*qdp_gcoeff(gamma_idx_1,true,false);
	contraction *= cmomenta.phase_factor(0,x_pos_vec);
    
	Rcomplex result_subdiag1(1.0);
	{
	  SpinColorFlavorMatrix sdiag1_trset0_scfmat_prod(prop_sprimes_snk_x_src_y_hconj_scfmat);
	  sdiag1_trset0_scfmat_prod.pl(F1);
	  sdiag1_trset0_scfmat_prod.ccl(-1);
	  sdiag1_trset0_scfmat_prod.gl(-5);
	  qdp_gl(sdiag1_trset0_scfmat_prod,gamma_idx_2);
	  sdiag1_trset0_scfmat_prod.ccl(-1);
	  sdiag1_trset0_scfmat_prod.pr(F1);
	  sdiag1_trset0_scfmat_prod.gr(-5);
	  sdiag1_trset0_scfmat_prod.ccr(1);
	  qdp_gr(sdiag1_trset0_scfmat_prod,gamma_idx_1);
	  sdiag1_trset0_scfmat_prod.ccr(1);
	  sdiag1_trset0_scfmat_prod *= prop_ud_snk_x_src_y_scfmat;
	  Rcomplex sdiag1_trset0_cmplx;
	  sdiag1_trset0_cmplx = sdiag1_trset0_scfmat_prod.Trace();
	  result_subdiag1 *= sdiag1_trset0_cmplx;
	}
	contraction *= result_subdiag1;
    
    
	corrfunc(omp_get_thread_num(),0,x_pos_vec[3]) += contraction;
      }
    }

    corrfunc.write(fp);
  }



  {
    /*<<(\bar u,s')*(\bar d,s)>>*/
    /*Require a "CorrelationFunction &corrfunc" with option "CorrelationFunction::THREADED"*/
    std::ostringstream os; os << "HL_MESON_U_SPRIME_D_S " << gamma_idx_1 << " " << gamma_idx_2;
    CorrelationFunction corrfunc(os.str().c_str(),CorrelationFunction::THREADED);

    if(UniqueID()==0) printf("Doing HL_MESON_U_SPRIME_D_S %d %d contraction\n",gamma_idx_1,gamma_idx_2);
    /*Require propagator "QPropWcontainer &prop_src_y_u_d_eitherflav_pcon corresponding to \mathcal{G}^{[u/d] }_{x,y} with source of either flavour (full prop matrix is generated using single flavour source). Source must be real.*/
    /*Require propagator "QPropWcontainer &prop_src_y_sprime_s_eitherflav_pcon corresponding to \mathcal{G}^{[s^\prime/s] }_{x,y} with source of either flavour (full prop matrix is generated using single flavour source). Source must be real.*/

    /*Fourier transform on sink index x*/
    /*Require a 3-component array 'desired_mom_x' representing the required momentum at this sink position*/
    ContractionQuarkMomCombination cmomenta;
    {
      QuarkMomCombination momcomb;
      momcomb.add_prop(prop_src_y_u_d_eitherflav_pcon, false);
      momcomb.add_prop(prop_src_y_sprime_s_eitherflav_pcon, true);
      cmomenta.add_contraction(0,momcomb);
    }
    cmomenta.set_desired_momentum(desired_mom_x);

    corrfunc.setNcontractions(1);
#pragma omp parallel for default(shared)
    for(int x=0;x<GJP.VolNodeSites();x++){
      int x_pos_vec[4];
      global_coord(x,x_pos_vec);
  
      /*Get all SpinColorFlavorMatrices needed*/
      SpinColorFlavorMatrix prop_sprimes_snk_x_src_y_hconj_scfmat(prop_src_y_sprime_s_eitherflav_pcon , AlgLattice(), x);
      prop_sprimes_snk_x_src_y_hconj_scfmat.hconj();
  
      SpinColorFlavorMatrix prop_ud_snk_x_src_y_scfmat(prop_src_y_u_d_eitherflav_pcon , AlgLattice(), x);
  
      /*Starting contraction 0*/
      /*[{\rm tr}_{scf,0}\left\{\Gamma[g2] \gamma^5 F_0 \mathcal{G}^{[s^\prime/s] \dagger}_{x,y} F_1 \gamma^5 C \Gamma[g1] C \mathcal{G}^{[u/d] }_{x,y}\right\}_{0}  ]*[f_\Gamma(g1,T) ]*/
  
      {
	Rcomplex contraction(1 , 0);
	contraction *= qdp_gcoeff(gamma_idx_1,true,false);
	contraction *= cmomenta.phase_factor(0,x_pos_vec);
    
	Rcomplex result_subdiag0(1.0);
	{
	  SpinColorFlavorMatrix sdiag0_trset0_scfmat_prod(prop_sprimes_snk_x_src_y_hconj_scfmat);
	  sdiag0_trset0_scfmat_prod.pl(F0);
	  sdiag0_trset0_scfmat_prod.gl(-5);
	  qdp_gl(sdiag0_trset0_scfmat_prod,gamma_idx_2);
	  sdiag0_trset0_scfmat_prod.pr(F1);
	  sdiag0_trset0_scfmat_prod.gr(-5);
	  sdiag0_trset0_scfmat_prod.ccr(1);
	  qdp_gr(sdiag0_trset0_scfmat_prod,gamma_idx_1);
	  sdiag0_trset0_scfmat_prod.ccr(1);
	  sdiag0_trset0_scfmat_prod *= prop_ud_snk_x_src_y_scfmat;
	  Rcomplex sdiag0_trset0_cmplx;
	  sdiag0_trset0_cmplx = sdiag0_trset0_scfmat_prod.Trace();
	  result_subdiag0 *= sdiag0_trset0_cmplx;
	}
	contraction *= result_subdiag0;
    
    
	corrfunc(omp_get_thread_num(),0,x_pos_vec[3]) += contraction;
      }
    }


    corrfunc.write(fp);
  }
}


void AlgGparityContract::meson_HL_std(QPropWcontainer &prop_H,QPropWcontainer &prop_L, const int* sink_mom, const int &gamma_idx_1, const int &gamma_idx_2, FILE *fp){
  /*Mesons comprising $ \bar u $ and $ s$*/
  /*Require a "CorrelationFunction &corrfunc"*/
  std::ostringstream os; os << "HL_MESON_L_H " << gamma_idx_1 << " " << gamma_idx_2;
  CorrelationFunction corrfunc(os.str().c_str(),CorrelationFunction::THREADED);

  if(UniqueID()==0) printf("Doing HL_MESON_L_H %d %d contraction\n",gamma_idx_1,gamma_idx_2);

  /*Require propagator "QPropWcontainer &prop_src_y_0_pcon corresponding to \mathcal{G}^{(0)}_{x,y}*/
  /*Require propagator "QPropWcontainer &prop_src_y_s_pcon corresponding to \mathcal{G}^{(s)}_{x,y}*/
  QPropWcontainer &prop_src_y_0_pcon = prop_L;
  QPropWcontainer &prop_src_y_s_pcon = prop_H;

  /*Fourier transform on sink index x*/
  /*Require a 3-component array 'desired_mom_x' representing the required momentum at this sink position*/
  const int *desired_mom_x = sink_mom;
  ContractionQuarkMomCombination cmomenta;
  {
    QuarkMomCombination momcomb;
    momcomb.add_prop(prop_src_y_0_pcon, true);
    momcomb.add_prop(prop_src_y_s_pcon, false);
    cmomenta.add_contraction(0,momcomb);
  }
  cmomenta.set_desired_momentum(desired_mom_x);

  corrfunc.setNcontractions(1);
#pragma omp parallel for default(shared)
  for(int x=0;x<GJP.VolNodeSites();x++){
    int x_pos_vec[4];
    global_coord(x,x_pos_vec);
  
    /*Get all WilsonMatrices needed*/
    WilsonMatrix& prop_snk_x_s_src_y_s_wmat = prop_src_y_s_pcon.getProp(AlgLattice()).SiteMatrix(x,0);
  
    WilsonMatrix prop_snk_x_0_src_y_0_hconj_wmat(prop_src_y_0_pcon.getProp(AlgLattice()).SiteMatrix(x,0));
    prop_snk_x_0_src_y_0_hconj_wmat.hconj();
  
    /*Starting contraction 0*/
    /*[-1 ]*[{\rm tr}_{sc,0}\left\{\gamma^5 S_1 \mathcal{G}^{(s)}_{x,y} S_2 \gamma^5 \mathcal{G}^{(0) \dagger}_{x,y}\right\}_{0}  ]*/
  
    {
      Rcomplex contraction(-1 , 0);
      contraction *= cmomenta.phase_factor(0,x_pos_vec);
    
      Rcomplex result_subdiag1(1.0);
      {
	WilsonMatrix sdiag1_trset0_wmat_prod(prop_snk_x_s_src_y_s_wmat);
	qdp_gl(sdiag1_trset0_wmat_prod,gamma_idx_1);

	sdiag1_trset0_wmat_prod.gl(-5);
	qdp_gr(sdiag1_trset0_wmat_prod,gamma_idx_2);

	sdiag1_trset0_wmat_prod.gr(-5);
	sdiag1_trset0_wmat_prod *= prop_snk_x_0_src_y_0_hconj_wmat;
	Rcomplex sdiag1_trset0_cmplx;
	sdiag1_trset0_cmplx = sdiag1_trset0_wmat_prod.Trace();
	result_subdiag1 *= sdiag1_trset0_cmplx;
      }
      contraction *= result_subdiag1;
    
    
      corrfunc(omp_get_thread_num(),0,x_pos_vec[3]) += contraction;
    }
  }
  corrfunc.write(fp);
}

void AlgGparityContract::contract_HL_mesons(const ContractionTypeHLMesons &args, const int &conf_idx){
  std::ostringstream file; file << args.file << "." << conf_idx;

  FILE *fp;
  if ((fp = Fopen(file.str().c_str(), "w")) == NULL) {
    ERR.FileW("CorrelationFunction","write(const char *file)",file.str().c_str());
  }

  PropagatorContainer &prop_H = PropManager::getProp(args.prop_H);
  PropagatorContainer &prop_L = PropManager::getProp(args.prop_L);

  if(prop_H.type() == QPROPW_TYPE && prop_L.type() == QPROPW_TYPE){
  //loop through HL meson correlation functions
  if(GJP.Gparity()){
    for(int g1=0;g1<16;g1++){
      for(int g2=0;g2<16;g2++){
	meson_HL_gparity(prop_H.convert<QPropWcontainer>(),prop_L.convert<QPropWcontainer>(), args.sink_mom, g1, g2, fp);
      }
    }
  }else{
    for(int g1=0;g1<16;g1++){
      for(int g2=0;g2<16;g2++){
  	meson_HL_std(prop_H.convert<QPropWcontainer>(),prop_L.convert<QPropWcontainer>(), args.sink_mom, g1, g2, fp);
      }
    }
  }
  }else ERR.General("AlgGparityContract","contract_HL_mesons(const ContractionTypeHLMesons &args, const int &conf_idx)","Not implemented for types other than QPROPW_TYPE\n");

  Fclose(fp);
}

void AlgGparityContract::contract_OVVpAA_gparity(const ContractionTypeOVVpAA &args, const int &conf_idx){
  std::ostringstream os; os << "O_VV_P_AA";
  CorrelationFunction corrfunc(os.str().c_str(),CorrelationFunction::THREADED);

  contract_OVVpAA_gparity(corrfunc,args);

  std::ostringstream file; file << args.file << "." << conf_idx;
  corrfunc.write(file.str().c_str());
}  

void AlgGparityContract::contract_OVVpAA_gparity(CorrelationFunction &corrfunc, const ContractionTypeOVVpAA &args){
  /*Require a "CorrelationFunction &corrfunc"*/
  /*Require propagator "QPropWcontainer &prop_src_z_u_d_eitherflav_pcon corresponding to \mathcal{G}^{[u/d] }_{y,z} with source of either flavour (full prop matrix is generated using single flavour source). Source must be real.*/
  /*Require propagator "QPropWcontainer &prop_src_z_sprime_s_eitherflav_pcon corresponding to \mathcal{G}^{[s^\prime/s] }_{y,z} with source of either flavour (full prop matrix is generated using single flavour source). Source must be real.*/
  /*Require propagator "QPropWcontainer &prop_src_x_u_d_eitherflav_pcon corresponding to \mathcal{G}^{[u/d] }_{y,x} with source of either flavour (full prop matrix is generated using single flavour source). Source must be real.*/
  /*Require propagator "QPropWcontainer &prop_src_x_sprime_s_eitherflav_pcon corresponding to \mathcal{G}^{[s^\prime/s] }_{y,x} with source of either flavour (full prop matrix is generated using single flavour source). Source must be real.*/

  if(UniqueID()==0) printf("Doing OVVpAA contractions with G-parity BCs\n");
  const char* fname = "contract_OVVpAA_gparity(CorrelationFunction &corrfunc, const ContractionTypeOVVpAA &args)";

  QPropWcontainer &prop_src_z_u_d_eitherflav_pcon = QPropWcontainer::verify_convert(PropManager::getProp(args.prop_L_t1),cname,fname);
  QPropWcontainer &prop_src_z_sprime_s_eitherflav_pcon = QPropWcontainer::verify_convert(PropManager::getProp(args.prop_H_t1),cname,fname);
  QPropWcontainer &prop_src_x_u_d_eitherflav_pcon = QPropWcontainer::verify_convert(PropManager::getProp(args.prop_L_t0),cname,fname);
  QPropWcontainer &prop_src_x_sprime_s_eitherflav_pcon = QPropWcontainer::verify_convert(PropManager::getProp(args.prop_H_t0),cname,fname);

  /*Fourier transform on sink index y*/
  /*Require a 3-component array 'desired_mom_y' representing the required momentum at this sink position*/
  int desired_mom_y[] = {0,0,0};

  ContractionQuarkMomCombination cmomenta;
  {
    QuarkMomCombination momcomb;
    momcomb.add_prop(prop_src_z_u_d_eitherflav_pcon, false);
    momcomb.add_prop(prop_src_z_sprime_s_eitherflav_pcon, true);
    momcomb.add_prop(prop_src_x_u_d_eitherflav_pcon, false);
    momcomb.add_prop(prop_src_x_sprime_s_eitherflav_pcon, true);
    cmomenta.add_contraction(0,momcomb);
    cmomenta.same(1,0);
    cmomenta.same(2,0);
    cmomenta.same(3,0);
  }
  cmomenta.set_desired_momentum(desired_mom_y);


  corrfunc.setNcontractions(4);

#pragma omp parallel for default(shared)
  for(int y=0;y<GJP.VolNodeSites();y++){
    int y_pos_vec[4];
    global_coord(y,y_pos_vec);
  
    /*Get all SpinColorFlavorMatrices needed*/
    SpinColorFlavorMatrix prop_sprimes_snk_y_src_x_hconj_scfmat(prop_src_x_sprime_s_eitherflav_pcon , AlgLattice(), y);
    prop_sprimes_snk_y_src_x_hconj_scfmat.hconj();
  
    SpinColorFlavorMatrix prop_sprimes_snk_y_src_z_hconj_scfmat(prop_src_z_sprime_s_eitherflav_pcon , AlgLattice(), y);
    prop_sprimes_snk_y_src_z_hconj_scfmat.hconj();
  
    SpinColorFlavorMatrix prop_ud_snk_y_src_z_scfmat(prop_src_z_u_d_eitherflav_pcon , AlgLattice(), y);
  
    SpinColorFlavorMatrix prop_ud_snk_y_src_x_scfmat(prop_src_x_u_d_eitherflav_pcon , AlgLattice(), y);
  
    /*Starting contraction 0*/
    /*[{\rm tr}_{scf,0}\left\{\mathcal{G}^{[s^\prime/s] \dagger}_{y,z} F_0 \gamma^\mu \gamma^5 \mathcal{G}^{[u/d] }_{y,z}\right\}_{0}  ]*[{\rm tr}_{scf,0}\left\{\mathcal{G}^{[s^\prime/s] \dagger}_{y,x} F_0 \gamma^\mu \gamma^5 \mathcal{G}^{[u/d] }_{y,x}\right\}_{0}  ]*[0.5 ]*/
  
    {
      Rcomplex contraction(0.5 , 0);
      contraction *= cmomenta.phase_factor(0,y_pos_vec);
    
      /*Contraction contains implicit sum over 1 free gamma matrix indices*/
      Rcomplex gsum_result(0.0);
      for(int gidx0=0;gidx0<4;gidx0++){
	Rcomplex g_prod(1.0);
      
	Rcomplex result_subdiag0(1.0);
	{
	  SpinColorFlavorMatrix sdiag0_trset0_scfmat_prod(prop_sprimes_snk_y_src_z_hconj_scfmat);
	  sdiag0_trset0_scfmat_prod.pr(F0);
	  sdiag0_trset0_scfmat_prod.gr(gidx0);
	  sdiag0_trset0_scfmat_prod.gr(-5);
	  sdiag0_trset0_scfmat_prod *= prop_ud_snk_y_src_z_scfmat;
	  Rcomplex sdiag0_trset0_cmplx;
	  sdiag0_trset0_cmplx = sdiag0_trset0_scfmat_prod.Trace();
	  result_subdiag0 *= sdiag0_trset0_cmplx;
	}
	g_prod *= result_subdiag0;
      
	Rcomplex result_subdiag1(1.0);
	{
	  SpinColorFlavorMatrix sdiag1_trset0_scfmat_prod(prop_sprimes_snk_y_src_x_hconj_scfmat);
	  sdiag1_trset0_scfmat_prod.pr(F0);
	  sdiag1_trset0_scfmat_prod.gr(gidx0);
	  sdiag1_trset0_scfmat_prod.gr(-5);
	  sdiag1_trset0_scfmat_prod *= prop_ud_snk_y_src_x_scfmat;
	  Rcomplex sdiag1_trset0_cmplx;
	  sdiag1_trset0_cmplx = sdiag1_trset0_scfmat_prod.Trace();
	  result_subdiag1 *= sdiag1_trset0_cmplx;
	}
	g_prod *= result_subdiag1;
      
	gsum_result += g_prod;
      }
      contraction *= gsum_result;
    
      corrfunc(omp_get_thread_num(),0,y_pos_vec[3]) += contraction;
    }
    /*Starting contraction 1*/
    /*[-0.5 ]*[{\rm tr}_{scf,0}\left\{\mathcal{G}^{[s^\prime/s] \dagger}_{y,x} F_0 \gamma^\mu \gamma^5 \mathcal{G}^{[u/d] }_{y,z} \mathcal{G}^{[s^\prime/s] \dagger}_{y,z} F_0 \gamma^\mu \gamma^5 \mathcal{G}^{[u/d] }_{y,x}\right\}_{0}  ]*/
  
    {
      Rcomplex contraction(-0.5 , 0);
      contraction *= cmomenta.phase_factor(1,y_pos_vec);
    
      /*Contraction contains implicit sum over 1 free gamma matrix indices*/
      Rcomplex gsum_result(0.0);
      for(int gidx0=0;gidx0<4;gidx0++){
	Rcomplex g_prod(1.0);
      
	Rcomplex result_subdiag1(1.0);
	{
	  SpinColorFlavorMatrix sdiag1_trset0_scfmat_prod(prop_sprimes_snk_y_src_x_hconj_scfmat);
	  sdiag1_trset0_scfmat_prod.pr(F0);
	  sdiag1_trset0_scfmat_prod.gr(gidx0);
	  sdiag1_trset0_scfmat_prod.gr(-5);
	  sdiag1_trset0_scfmat_prod *= prop_ud_snk_y_src_z_scfmat;
	  sdiag1_trset0_scfmat_prod *= prop_sprimes_snk_y_src_z_hconj_scfmat;
	  sdiag1_trset0_scfmat_prod.pr(F0);
	  sdiag1_trset0_scfmat_prod.gr(gidx0);
	  sdiag1_trset0_scfmat_prod.gr(-5);
	  sdiag1_trset0_scfmat_prod *= prop_ud_snk_y_src_x_scfmat;
	  Rcomplex sdiag1_trset0_cmplx;
	  sdiag1_trset0_cmplx = sdiag1_trset0_scfmat_prod.Trace();
	  result_subdiag1 *= sdiag1_trset0_cmplx;
	}
	g_prod *= result_subdiag1;
      
	gsum_result += g_prod;
      }
      contraction *= gsum_result;
    
      corrfunc(omp_get_thread_num(),1,y_pos_vec[3]) += contraction;
    }
    /*Starting contraction 2*/
    /*[{\rm tr}_{scf,0}\left\{\mathcal{G}^{[s^\prime/s] \dagger}_{y,z} F_0 \gamma^\mu \mathcal{G}^{[u/d] }_{y,z}\right\}_{0}  ]*[{\rm tr}_{scf,0}\left\{\mathcal{G}^{[s^\prime/s] \dagger}_{y,x} F_0 \gamma^\mu \mathcal{G}^{[u/d] }_{y,x}\right\}_{0}  ]*[0.5 ]*/
  
    {
      Rcomplex contraction(0.5 , 0);
      contraction *= cmomenta.phase_factor(2,y_pos_vec);
    
      /*Contraction contains implicit sum over 1 free gamma matrix indices*/
      Rcomplex gsum_result(0.0);
      for(int gidx0=0;gidx0<4;gidx0++){
	Rcomplex g_prod(1.0);
      
	Rcomplex result_subdiag0(1.0);
	{
	  SpinColorFlavorMatrix sdiag0_trset0_scfmat_prod(prop_sprimes_snk_y_src_z_hconj_scfmat);
	  sdiag0_trset0_scfmat_prod.pr(F0);
	  sdiag0_trset0_scfmat_prod.gr(gidx0);
	  sdiag0_trset0_scfmat_prod *= prop_ud_snk_y_src_z_scfmat;
	  Rcomplex sdiag0_trset0_cmplx;
	  sdiag0_trset0_cmplx = sdiag0_trset0_scfmat_prod.Trace();
	  result_subdiag0 *= sdiag0_trset0_cmplx;
	}
	g_prod *= result_subdiag0;
      
	Rcomplex result_subdiag1(1.0);
	{
	  SpinColorFlavorMatrix sdiag1_trset0_scfmat_prod(prop_sprimes_snk_y_src_x_hconj_scfmat);
	  sdiag1_trset0_scfmat_prod.pr(F0);
	  sdiag1_trset0_scfmat_prod.gr(gidx0);
	  sdiag1_trset0_scfmat_prod *= prop_ud_snk_y_src_x_scfmat;
	  Rcomplex sdiag1_trset0_cmplx;
	  sdiag1_trset0_cmplx = sdiag1_trset0_scfmat_prod.Trace();
	  result_subdiag1 *= sdiag1_trset0_cmplx;
	}
	g_prod *= result_subdiag1;
      
	gsum_result += g_prod;
      }
      contraction *= gsum_result;
    
      corrfunc(omp_get_thread_num(),2,y_pos_vec[3]) += contraction;
    }
    /*Starting contraction 3*/
    /*[-0.5 ]*[{\rm tr}_{scf,0}\left\{\mathcal{G}^{[s^\prime/s] \dagger}_{y,x} F_0 \gamma^\mu \mathcal{G}^{[u/d] }_{y,z} \mathcal{G}^{[s^\prime/s] \dagger}_{y,z} F_0 \gamma^\mu \mathcal{G}^{[u/d] }_{y,x}\right\}_{0}  ]*/
  
    {
      Rcomplex contraction(-0.5 , 0);
      contraction *= cmomenta.phase_factor(3,y_pos_vec);
    
      /*Contraction contains implicit sum over 1 free gamma matrix indices*/
      Rcomplex gsum_result(0.0);
      for(int gidx0=0;gidx0<4;gidx0++){
	Rcomplex g_prod(1.0);
      
	Rcomplex result_subdiag1(1.0);
	{
	  SpinColorFlavorMatrix sdiag1_trset0_scfmat_prod(prop_sprimes_snk_y_src_x_hconj_scfmat);
	  sdiag1_trset0_scfmat_prod.pr(F0);
	  sdiag1_trset0_scfmat_prod.gr(gidx0);
	  sdiag1_trset0_scfmat_prod *= prop_ud_snk_y_src_z_scfmat;
	  sdiag1_trset0_scfmat_prod *= prop_sprimes_snk_y_src_z_hconj_scfmat;
	  sdiag1_trset0_scfmat_prod.pr(F0);
	  sdiag1_trset0_scfmat_prod.gr(gidx0);
	  sdiag1_trset0_scfmat_prod *= prop_ud_snk_y_src_x_scfmat;
	  Rcomplex sdiag1_trset0_cmplx;
	  sdiag1_trset0_cmplx = sdiag1_trset0_scfmat_prod.Trace();
	  result_subdiag1 *= sdiag1_trset0_cmplx;
	}
	g_prod *= result_subdiag1;
      
	gsum_result += g_prod;
      }
      contraction *= gsum_result;
    
      corrfunc(omp_get_thread_num(),3,y_pos_vec[3]) += contraction;
    }
  }
  corrfunc.sumLattice();
}

void AlgGparityContract::contract_OVVpAA_std(const ContractionTypeOVVpAA &args, const int &conf_idx){

  /*Require a "CorrelationFunction &corrfunc"*/
  std::ostringstream os; os << "O_VV_P_AA";
  CorrelationFunction corrfunc(os.str().c_str(),CorrelationFunction::THREADED);

  if(UniqueID()==0) printf("Doing OVVpAA contractions\n");

  /*Require propagator "QPropWcontainer &prop_src_z_0_pcon corresponding to \mathcal{G}^{(0)}_{y,z}*/
  /*Require propagator "QPropWcontainer &prop_src_z_s_pcon corresponding to \mathcal{G}^{(s)}_{y,z}*/
  /*Require propagator "QPropWcontainer &prop_src_x_0_pcon corresponding to \mathcal{G}^{(0)}_{y,x}*/
  /*Require propagator "QPropWcontainer &prop_src_x_s_pcon corresponding to \mathcal{G}^{(s)}_{y,x}*/
  const char* fname = "contract_OVVpAA_std(const ContractionTypeOVVpAA &args, const int &conf_idx)";

  QPropWcontainer &prop_src_z_0_pcon = QPropWcontainer::verify_convert(PropManager::getProp(args.prop_L_t1),cname,fname);
  QPropWcontainer &prop_src_z_s_pcon = QPropWcontainer::verify_convert(PropManager::getProp(args.prop_H_t1),cname,fname);
  QPropWcontainer &prop_src_x_0_pcon = QPropWcontainer::verify_convert(PropManager::getProp(args.prop_L_t0),cname,fname);
  QPropWcontainer &prop_src_x_s_pcon = QPropWcontainer::verify_convert(PropManager::getProp(args.prop_H_t0),cname,fname);

  /*Fourier transform on sink index y*/
  /*Require a 3-component array 'desired_mom_y' representing the required momentum at this sink position*/
  int desired_mom_y[] = {0,0,0};

  ContractionQuarkMomCombination cmomenta;
  {
    QuarkMomCombination momcomb;
    momcomb.add_prop(prop_src_z_0_pcon, false);
    momcomb.add_prop(prop_src_z_s_pcon, true);
    momcomb.add_prop(prop_src_x_0_pcon, false);
    momcomb.add_prop(prop_src_x_s_pcon, true);
    cmomenta.add_contraction(0,momcomb);
    cmomenta.same(1,0);
    cmomenta.same(2,0);
    cmomenta.same(3,0);
  }
  cmomenta.set_desired_momentum(desired_mom_y);

  corrfunc.setNcontractions(4);
#pragma omp parallel for default(shared)
  for(int y=0;y<GJP.VolNodeSites();y++){
    int y_pos_vec[4];
    global_coord(y,y_pos_vec);
  
    /*Get all WilsonMatrices needed*/
    WilsonMatrix prop_snk_y_s_src_z_s_hconj_wmat(prop_src_z_s_pcon.getProp(AlgLattice()).SiteMatrix(y,0));
    prop_snk_y_s_src_z_s_hconj_wmat.hconj();
  
    WilsonMatrix& prop_snk_y_0_src_z_0_wmat = prop_src_z_0_pcon.getProp(AlgLattice()).SiteMatrix(y,0);
  
    WilsonMatrix& prop_snk_y_0_src_x_0_wmat = prop_src_x_0_pcon.getProp(AlgLattice()).SiteMatrix(y,0);
  
    WilsonMatrix prop_snk_y_s_src_x_s_hconj_wmat(prop_src_x_s_pcon.getProp(AlgLattice()).SiteMatrix(y,0));
    prop_snk_y_s_src_x_s_hconj_wmat.hconj();
  
    /*Starting contraction 0*/
    /*[{\rm tr}_{sc,0}\left\{\mathcal{G}^{(s) \dagger}_{y,z} \gamma^\mu \gamma^5 \mathcal{G}^{(0)}_{y,z}\right\}_{0}  ]*[{\rm tr}_{sc,0}\left\{\mathcal{G}^{(s) \dagger}_{y,x} \gamma^\mu \gamma^5 \mathcal{G}^{(0)}_{y,x}\right\}_{0}  ]*[2 ]*/
  
    {
      Rcomplex contraction(2 , 0);
      contraction *= cmomenta.phase_factor(0,y_pos_vec);
    
      /*Contraction contains implicit sum over 1 free gamma matrix indices*/
      Rcomplex gsum_result(0.0);
      for(int gidx0=0;gidx0<4;gidx0++){
	Rcomplex g_prod(1.0);
      
	Rcomplex result_subdiag0(1.0);
	{
	  WilsonMatrix sdiag0_trset0_wmat_prod(prop_snk_y_s_src_z_s_hconj_wmat);
	  sdiag0_trset0_wmat_prod.gr(gidx0);
	  sdiag0_trset0_wmat_prod.gr(-5);
	  sdiag0_trset0_wmat_prod *= prop_snk_y_0_src_z_0_wmat;
	  Rcomplex sdiag0_trset0_cmplx;
	  sdiag0_trset0_cmplx = sdiag0_trset0_wmat_prod.Trace();
	  result_subdiag0 *= sdiag0_trset0_cmplx;
	}
	g_prod *= result_subdiag0;
      
	Rcomplex result_subdiag1(1.0);
	{
	  WilsonMatrix sdiag1_trset0_wmat_prod(prop_snk_y_s_src_x_s_hconj_wmat);
	  sdiag1_trset0_wmat_prod.gr(gidx0);
	  sdiag1_trset0_wmat_prod.gr(-5);
	  sdiag1_trset0_wmat_prod *= prop_snk_y_0_src_x_0_wmat;
	  Rcomplex sdiag1_trset0_cmplx;
	  sdiag1_trset0_cmplx = sdiag1_trset0_wmat_prod.Trace();
	  result_subdiag1 *= sdiag1_trset0_cmplx;
	}
	g_prod *= result_subdiag1;
      
	gsum_result += g_prod;
      }
      contraction *= gsum_result;
    
      corrfunc(omp_get_thread_num(),0,y_pos_vec[3]) += contraction;
    }
    /*Starting contraction 1*/
    /*[-2 ]*[{\rm tr}_{sc,0}\left\{\mathcal{G}^{(s) \dagger}_{y,x} \gamma^\mu \gamma^5 \mathcal{G}^{(0)}_{y,z} \mathcal{G}^{(s) \dagger}_{y,z} \gamma^\mu \gamma^5 \mathcal{G}^{(0)}_{y,x}\right\}_{0}  ]*/
  
    {
      Rcomplex contraction(-2 , 0);
      contraction *= cmomenta.phase_factor(1,y_pos_vec);
    
      /*Contraction contains implicit sum over 1 free gamma matrix indices*/
      Rcomplex gsum_result(0.0);
      for(int gidx0=0;gidx0<4;gidx0++){
	Rcomplex g_prod(1.0);
      
	Rcomplex result_subdiag1(1.0);
	{
	  WilsonMatrix sdiag1_trset0_wmat_prod(prop_snk_y_s_src_x_s_hconj_wmat);
	  sdiag1_trset0_wmat_prod.gr(gidx0);
	  sdiag1_trset0_wmat_prod.gr(-5);
	  sdiag1_trset0_wmat_prod *= prop_snk_y_0_src_z_0_wmat;
	  sdiag1_trset0_wmat_prod *= prop_snk_y_s_src_z_s_hconj_wmat;
	  sdiag1_trset0_wmat_prod.gr(gidx0);
	  sdiag1_trset0_wmat_prod.gr(-5);
	  sdiag1_trset0_wmat_prod *= prop_snk_y_0_src_x_0_wmat;
	  Rcomplex sdiag1_trset0_cmplx;
	  sdiag1_trset0_cmplx = sdiag1_trset0_wmat_prod.Trace();
	  result_subdiag1 *= sdiag1_trset0_cmplx;
	}
	g_prod *= result_subdiag1;
      
	gsum_result += g_prod;
      }
      contraction *= gsum_result;
    
      corrfunc(omp_get_thread_num(),1,y_pos_vec[3]) += contraction;
    }
    /*Starting contraction 2*/
    /*[{\rm tr}_{sc,0}\left\{\mathcal{G}^{(s) \dagger}_{y,z} \gamma^\mu \mathcal{G}^{(0)}_{y,z}\right\}_{0}  ]*[{\rm tr}_{sc,0}\left\{\mathcal{G}^{(s) \dagger}_{y,x} \gamma^\mu \mathcal{G}^{(0)}_{y,x}\right\}_{0}  ]*[2 ]*/
  
    {
      Rcomplex contraction(2 , 0);
      contraction *= cmomenta.phase_factor(2,y_pos_vec);
    
      /*Contraction contains implicit sum over 1 free gamma matrix indices*/
      Rcomplex gsum_result(0.0);
      for(int gidx0=0;gidx0<4;gidx0++){
	Rcomplex g_prod(1.0);
      
	Rcomplex result_subdiag0(1.0);
	{
	  WilsonMatrix sdiag0_trset0_wmat_prod(prop_snk_y_s_src_z_s_hconj_wmat);
	  sdiag0_trset0_wmat_prod.gr(gidx0);
	  sdiag0_trset0_wmat_prod *= prop_snk_y_0_src_z_0_wmat;
	  Rcomplex sdiag0_trset0_cmplx;
	  sdiag0_trset0_cmplx = sdiag0_trset0_wmat_prod.Trace();
	  result_subdiag0 *= sdiag0_trset0_cmplx;
	}
	g_prod *= result_subdiag0;
      
	Rcomplex result_subdiag1(1.0);
	{
	  WilsonMatrix sdiag1_trset0_wmat_prod(prop_snk_y_s_src_x_s_hconj_wmat);
	  sdiag1_trset0_wmat_prod.gr(gidx0);
	  sdiag1_trset0_wmat_prod *= prop_snk_y_0_src_x_0_wmat;
	  Rcomplex sdiag1_trset0_cmplx;
	  sdiag1_trset0_cmplx = sdiag1_trset0_wmat_prod.Trace();
	  result_subdiag1 *= sdiag1_trset0_cmplx;
	}
	g_prod *= result_subdiag1;
      
	gsum_result += g_prod;
      }
      contraction *= gsum_result;
    
      corrfunc(omp_get_thread_num(),2,y_pos_vec[3]) += contraction;
    }
    /*Starting contraction 3*/
    /*[-2 ]*[{\rm tr}_{sc,0}\left\{\mathcal{G}^{(s) \dagger}_{y,x} \gamma^\mu \mathcal{G}^{(0)}_{y,z} \mathcal{G}^{(s) \dagger}_{y,z} \gamma^\mu \mathcal{G}^{(0)}_{y,x}\right\}_{0}  ]*/
  
    {
      Rcomplex contraction(-2 , 0);
      contraction *= cmomenta.phase_factor(3,y_pos_vec);
    
      /*Contraction contains implicit sum over 1 free gamma matrix indices*/
      Rcomplex gsum_result(0.0);
      for(int gidx0=0;gidx0<4;gidx0++){
	Rcomplex g_prod(1.0);
      
	Rcomplex result_subdiag1(1.0);
	{
	  WilsonMatrix sdiag1_trset0_wmat_prod(prop_snk_y_s_src_x_s_hconj_wmat);
	  sdiag1_trset0_wmat_prod.gr(gidx0);
	  sdiag1_trset0_wmat_prod *= prop_snk_y_0_src_z_0_wmat;
	  sdiag1_trset0_wmat_prod *= prop_snk_y_s_src_z_s_hconj_wmat;
	  sdiag1_trset0_wmat_prod.gr(gidx0);
	  sdiag1_trset0_wmat_prod *= prop_snk_y_0_src_x_0_wmat;
	  Rcomplex sdiag1_trset0_cmplx;
	  sdiag1_trset0_cmplx = sdiag1_trset0_wmat_prod.Trace();
	  result_subdiag1 *= sdiag1_trset0_cmplx;
	}
	g_prod *= result_subdiag1;
      
	gsum_result += g_prod;
      }
      contraction *= gsum_result;
    
      corrfunc(omp_get_thread_num(),3,y_pos_vec[3]) += contraction;
    }
  }

  std::ostringstream file; file << args.file << "." << conf_idx;
  corrfunc.write(file.str().c_str());
}


void AlgGparityContract::contract_OVVpAA(const ContractionTypeOVVpAA &args, const int &conf_idx){
  if(GJP.Gparity()){
    return contract_OVVpAA_gparity(args,conf_idx);
  }else{
    return contract_OVVpAA_std(args,conf_idx);
  }
}

void AlgGparityContract::measure_topological_charge(const ContractionTypeTopologicalCharge &args, const int &conf_idx){
  const char *fname = "measure_topological_charge";
  
  // calculate the topological charge. Need to copy the lattice since
  // we need to smear it first. Use Chulwoo's "lattice container"
  LatticeContainer lat_cont;
  lat_cont.Get(AlgLattice()); //backup the lattice

  std::ostringstream file; file << args.file << "." << conf_idx;
  std::string filename = file.str();

  CommonArg common_arg;
  common_arg.filename = const_cast<char*>(filename.c_str()); //Oh my god I hate C-strings

  ApeSmearArg ape_arg;
  ape_arg.tolerance = args.ape_su3_proj_tolerance;
  ape_arg.orthog = args.ape_orthog;
  ape_arg.coef = args.ape_coef;

  AlgApeSmear ape(AlgLattice(), &common_arg, &ape_arg, args.ape_smear_su3_project);
  AlgTcharge  tcharge(AlgLattice(), &common_arg);
  for (int i = 0; i < args.n_ape_smearing_cycles; ++i) {
    VRB.Result(cname,fname,"%i\n",i);
    VRB.Result(cname,fname,"   running tcharge\n"); tcharge.run();
    VRB.Result(cname,fname,"   running ape\n"); ape.run();
    VRB.Result(cname,fname,"   running ape\n"); ape.run();
    VRB.Result(cname,fname,"   running ape\n"); ape.run();
  }
  tcharge.run();
  // restore the lattice
  lat_cont.Set(AlgLattice());
  common_arg.filename = NULL; //OH MY GOD I HATE C-STRINGS!
}

void AlgGparityContract::measure_wilson_flow(const ContractionTypeWilsonFlow &args, const int &conf_idx){
  const char *fname = "measure_wilson_flow";
  
  std::ostringstream file; file << args.file << "." << conf_idx;
  std::string filename = file.str();

  LatticeContainer lat_cont;
  lat_cont.Get(AlgLattice()); //backup the lattice

  CommonArg carg_null; //write no output using the nasty CommonArg system. Instead write in a more sensible format below

  for(int i=0;i<=args.n_steps;++i){
    AlgActionDensity ad(AlgLattice(),&carg_null);
    Float action_density;
    ad.smartrun(&action_density);

    Float t = i*args.time_step;
    
    FILE *fp;
    if( (fp = Fopen(filename.c_str(), "a")) == NULL ) ERR.FileA(cname,fname,filename.c_str());
    //Prints mean value of (1/2)tr(F_mu_nu F_mu_nu), in units of a^4
    Fprintf(fp, "%e %15e\n", t, action_density);
    Fclose(fp);

    if(i!=args.n_steps){
      AlgWilsonFlow wf(AlgLattice(),&carg_null, args.time_step);
      wf.smartrun();
    }
  }

  //Restore the lattice
  lat_cont.Set(AlgLattice());
}


void AlgGparityContract::measure_mres(const ContractionTypeMres &args, const int &conf_idx){
  std::ostringstream file; file << args.file << "." << conf_idx;

  FILE *fp;
  if ((fp = Fopen(file.str().c_str(), "w")) == NULL) {
    ERR.FileW(cname,"measure_mres",file.str().c_str());
  }
  CorrelationFunction pion("pion",1,CorrelationFunction::THREADED);
  CorrelationFunction j5_q("j5q",1,CorrelationFunction::THREADED);

  if(GJP.Gparity())
    measure_mres_gparity(args,pion,j5_q);
  else
    measure_mres(args,pion,j5_q);

  pion.write(fp);
  j5_q.write(fp);
  
  Fclose(fp);
}


void AlgGparityContract::measure_mres(const ContractionTypeMres &args, CorrelationFunction &pion, CorrelationFunction &j5_q){
  if(GJP.Gparity()) ERR.General(cname,"measure_mres(...)","This is not a G-parity measurement\n");
  if(GJP.Snodes()!=1) ERR.General(cname,"measure_mres(...)","Assumes only 1 node in s-direction\n");
  if(pion.threadType() != CorrelationFunction::THREADED || j5_q.threadType() != CorrelationFunction::THREADED) ERR.General(cname,"measure_mres(...)","Assumes multi-thread CorrelationFunctions\n");
  if(pion.nContractions() !=1 || j5_q.nContractions() != 1) ERR.General(cname,"measure_mres(...)","CorrelationFunctions must have space for only one contraction\n");

  PropagatorContainer *pc = &PropManager::getProp(args.prop);
  if(pc->type()!=QPROPW_TYPE) ERR.General(cname,"measure_mres(...)","Propagator must be QPropW type\n");
  QPropWcontainer &prop_pcon = pc->convert<QPropWcontainer>();
  if(!prop_pcon.hasAttr<StoreMidpropAttrArg>()) ERR.General(cname,"measure_mres(...)","Propagator must have midprop stored to form mres\n");

  QPropW & qp = prop_pcon.getProp(AlgLattice());

#pragma omp parallel for default(shared)
  for(int i = 0; i < GJP.VolNodeSites(); ++i) {
    int x[4];
    global_coord(i,x);
    
    // J5 contraction (pion)
    WilsonMatrix p[2]  = {qp[i], qp[i]};
    p[1].hconj();
    // J5q contraction (midplane)
    WilsonMatrix q[2]  = {qp(i), qp(i)};
    q[1].hconj();
    
    Rcomplex pion_incr = Trace(p[0], p[1]);
    Rcomplex j5q_incr = Trace(q[0], q[1]);
    
    pion(omp_get_thread_num(),0,x[3]) += pion_incr;
    j5_q(omp_get_thread_num(),0,x[3]) += j5q_incr;
  }
  
  pion.sumLattice();
  j5_q.sumLattice();
}

void AlgGparityContract::measure_mres_gparity(const ContractionTypeMres &args, CorrelationFunction &pion, CorrelationFunction &j5_q){
  if(!GJP.Gparity()) ERR.General(cname,"measure_mres_gparity(...)","This is a G-parity measurement\n");
  if(GJP.Snodes()!=1) ERR.General(cname,"measure_mres_gparity(...)","Assumes only 1 node in s-direction\n");
  if(pion.threadType() != CorrelationFunction::THREADED || j5_q.threadType() != CorrelationFunction::THREADED) ERR.General(cname,"measure_mres_gparity(...)","Assumes multi-thread CorrelationFunctions\n");
  if(pion.nContractions() !=1 || j5_q.nContractions() != 1) ERR.General(cname,"measure_mres_gparity(...)","CorrelationFunctions must have space for only one contraction\n");

  PropagatorContainer *pc = &PropManager::getProp(args.prop);
  if(pc->type()!=QPROPW_TYPE) ERR.General(cname,"measure_mres_gparity(...)","Propagator must be QPropW type\n");
  QPropWcontainer &prop_pcon = pc->convert<QPropWcontainer>();
  if(!prop_pcon.hasAttr<StoreMidpropAttrArg>()) ERR.General(cname,"measure_mres_gparity(...)","Propagator must have midprop stored to form mres\n");

#define DEFINITION_ONE

#ifdef DEFINITION_ONE
  int src_p[3];
  prop_pcon.momentum(src_p); //in units of pi/2L

  //In this version we explicitly create pion operators with the appropriate G-parity momentum projection.
  Float mom[3] = {0,0,0};
  for(int d=0;d<3;d++) mom[d] = Float(2*src_p[d]) / 4.0; //units of 2pi/L

#pragma omp parallel for default(shared)
  for(int i = 0; i < GJP.VolNodeSites(); ++i) {
    int x[4];
    global_coord(i,x);
    Rcomplex phase = phase_factor(mom,x);
    
    // J5 contraction (pion)  (factor of 1/2 not included)
    SpinColorFlavorMatrix p[2];
    p[0].generate(prop_pcon,AlgLattice(),i,SPLANE_BOUNDARY);
    p[1] = p[0];
    p[1].flipSourceMomentum(); //NOTE: This assumes the source matrix structure \eta (not including the phase factor) obeys  C\gamma^5 \sigma_2 \eta^* = \eta C\gamma^5 \sigma_2
                               //      which is a property of pretty much all standard source types (wall, cosine wall, point, etc, also with gauge fixing matrix)

    p[1].hconj();
    // J5q contraction (midplane)  (factor of 1/2 not included)
    SpinColorFlavorMatrix q[2];
    q[0].generate(prop_pcon,AlgLattice(),i,SPLANE_MIDPOINT);
    q[1] = q[0];
    q[1].flipSourceMomentum();
    q[1].hconj();

    p[0].pr(sigma3);
    p[1].pr(sigma3);
    Rcomplex pion_incr = Trace(p[0], p[1]);

    q[0].pr(sigma3);
    q[1].pr(sigma3);
    Rcomplex j5q_incr = Trace(q[0], q[1]);

    pion(omp_get_thread_num(),0,x[3]) += phase*pion_incr;
    j5_q(omp_get_thread_num(),0,x[3]) += phase*j5q_incr;
  }
#else  
  //A more naive definition is just to do exactly what we normally do; take    \bar\psi \gamma^5 \psi  as the source and sink operators and forget about G-parity
  //just summing over the sink flavour index as if it was part of the spatial position index, just like we would in the single flavour approach

  //Testing suggests both definitions agree very precisely. For a calculation of mres on its own, this method is probably better because you can use a naive wall
  //source. Using a cosine source in the above gives ~40% larger errors, whereas using wall sources gives errors that are almost exactly the same as with this method.
  //However with the above method you need to generate 2 wall sources, one of each flavour, in order to construct the SCF matrix. Might be able to get away without
  //the momentum projection and using one wall source - not yet tested.

#pragma omp parallel for default(shared)
  for(int i = 0; i < GJP.VolNodeSites(); ++i) {
    int x[4];
    global_coord(i,x);

    for(int f=0;f<2;++f){
      // J5 contraction (pion)
      WilsonMatrix p[2];
      p[0] = prop_pcon.getProp(AlgLattice()).SiteMatrix(i,f);
      p[1] = p[0];
      p[1].hconj();
      // J5q contraction (midplane)
      WilsonMatrix q[2];
      q[0] = prop_pcon.getProp(AlgLattice()).MidPlaneSiteMatrix(i,f);
      q[1] = q[0];
      q[1].hconj();
    
      Rcomplex pion_incr = Trace(p[0], p[1]);
      Rcomplex j5q_incr = Trace(q[0], q[1]);
    
      pion(omp_get_thread_num(),0,x[3]) += pion_incr;
      j5_q(omp_get_thread_num(),0,x[3]) += j5q_incr;
    }
  }


#endif


  pion.sumLattice();
  j5_q.sumLattice();
}


inline static void set_matrix(Float into[], const MatIdxAndCoeff * from, const int &size){
  for(int i=0;i<size;i++) into[ from[i].idx ] = from[i].coeff;
}

void AlgGparityContract::contract_a2a_bilinear(const ContractionTypeA2ABilinear &args, const int &conf){
  CorrelationFunction corr_a2a("A2A Bilinear",1, CorrelationFunction::THREADED);
  contract_a2a_bilinear(args,corr_a2a);
  std::ostringstream file; file << args.file << "." << conf;
  corr_a2a.write(file.str().c_str());
}

void AlgGparityContract::contract_a2a_bilinear(const ContractionTypeA2ABilinear &args, CorrelationFunction &corr_a2a){
#if 0
  Float src_gamma_matrix_linear_comb[16];
  Float snk_gamma_matrix_linear_comb[16];
  for(int i=0;i<16;i++){
    src_gamma_matrix_linear_comb[i] = 0.0;
    snk_gamma_matrix_linear_comb[i] = 0.0;
  }
  set_matrix(src_gamma_matrix_linear_comb, args.source_spin_matrix.source_spin_matrix_val, args.source_spin_matrix.source_spin_matrix_len);
  set_matrix(snk_gamma_matrix_linear_comb, args.sink_spin_matrix.sink_spin_matrix_val, args.sink_spin_matrix.sink_spin_matrix_len);

  Float src_flav_matrix_linear_comb[4];
  Float snk_flav_matrix_linear_comb[4];
  for(int i=0;i<4;i++){
    src_flav_matrix_linear_comb[i] = 0.0;
    snk_flav_matrix_linear_comb[i] = 0.0;
  }
  set_matrix(src_flav_matrix_linear_comb, args.source_flavor_matrix.source_flavor_matrix_val, args.source_flavor_matrix.source_flavor_matrix_len);
  set_matrix(snk_flav_matrix_linear_comb, args.sink_flavor_matrix.sink_flavor_matrix_val, args.sink_flavor_matrix.sink_flavor_matrix_len);

  MFqdpMatrix source_wdagv(MFstructure::W, MFstructure::V, true, false);  //W* V
  source_wdagv.set_matrix(src_gamma_matrix_linear_comb, src_flav_matrix_linear_comb);

  MFqdpMatrix sink_wdagv(MFstructure::W, MFstructure::V, true, false);  //W* V
  sink_wdagv.set_matrix(snk_gamma_matrix_linear_comb, snk_flav_matrix_linear_comb);
  
  MFBasicSource src_smearing;
  MFBasicSource::set_smearing(src_smearing, args.source_smearing);
  src_smearing.fft_src();
  
  MFBasicSource snk_smearing;
  MFBasicSource::set_smearing(snk_smearing, args.sink_smearing);
  snk_smearing.fft_src();
  
  A2APropbfm & prop_snk_src = PropManager::getProp(args.prop_snk_src).convert<A2ApropContainer>().getProp(AlgLattice()); //prop from snk -> src
  A2APropbfm & prop_src_snk = PropManager::getProp(args.prop_src_snk).convert<A2ApropContainer>().getProp(AlgLattice()); //prop from src -> snk
  
  MesonField2 mf_src(prop_src_snk,prop_snk_src, source_wdagv, src_smearing); // W*(t_src) V(t_src; t_snk)  where V has origin on t_snk and terminates at t_src and hence belongs to prop_snk_src , and W belongs to prop_src_snk
  MesonField2 mf_snk(prop_snk_src,prop_src_snk, sink_wdagv, snk_smearing); // W*(t_snk) V(t_snk; t_src)  where V has origin on t_src and terminates at t_snk and hence belongs to prop_src_snk , and W belongs to prop_snk_src
  
  //Calculate correlation function
  MesonField2::contract(mf_src,mf_snk,0,corr_a2a);
#endif
}

void AlgGparityContract::spectrum(const GparityMeasurement &measargs,const int &conf_idx){
  if(measargs.type == CONTRACTION_TYPE_LL_MESONS) contract_LL_mesons(measargs.GparityMeasurement_u.contraction_type_ll_mesons, conf_idx);
  else if(measargs.type == CONTRACTION_TYPE_HL_MESONS) contract_HL_mesons(measargs.GparityMeasurement_u.contraction_type_hl_mesons, conf_idx);
  else if(measargs.type == CONTRACTION_TYPE_O_VV_P_AA) contract_OVVpAA(measargs.GparityMeasurement_u.contraction_type_o_vv_p_aa, conf_idx);
  else if(measargs.type == CONTRACTION_TYPE_ALL_BILINEARS) contract_all_bilinears(measargs.GparityMeasurement_u.contraction_type_all_bilinears, conf_idx);
  else if(measargs.type == CONTRACTION_TYPE_ALL_WALLSINK_BILINEARS_SPECIFIC_MOMENTUM) contract_all_wallsink_bilinears_specific_momentum(measargs.GparityMeasurement_u.contraction_type_all_wallsink_bilinears_specific_momentum,
																	conf_idx);
  else if(measargs.type == CONTRACTION_TYPE_FOURIER_PROP) contract_fourier_prop(measargs.GparityMeasurement_u.contraction_type_fourier_prop, conf_idx);
  else if(measargs.type == CONTRACTION_TYPE_BILINEAR_VERTEX) contract_bilinear_vertex(measargs.GparityMeasurement_u.contraction_type_bilinear_vertex, conf_idx);
  else if(measargs.type == CONTRACTION_TYPE_QUADRILINEAR_VERTEX) contract_quadrilinear_vertex(measargs.GparityMeasurement_u.contraction_type_quadrilinear_vertex, conf_idx);
  else if(measargs.type == CONTRACTION_TYPE_TOPOLOGICAL_CHARGE) measure_topological_charge(measargs.GparityMeasurement_u.contraction_type_topological_charge, conf_idx);
  else if(measargs.type == CONTRACTION_TYPE_MRES) measure_mres(measargs.GparityMeasurement_u.contraction_type_mres, conf_idx);

  else if(measargs.type == CONTRACTION_TYPE_A2A_BILINEAR) contract_a2a_bilinear(measargs.GparityMeasurement_u.contraction_type_a2a_bilinear, conf_idx);

  else if(measargs.type == CONTRACTION_TYPE_WILSON_FLOW) measure_wilson_flow(measargs.GparityMeasurement_u.contraction_type_wilson_flow, conf_idx);

  else ERR.General("AlgGparityContract","spectrum(...)","Invalid contraction type");
}

//NOTE: This code computes  \sum_x tr( M1 G_1^dag(x) M2 G_2(x) ) for all M1, M2  (spin without GPBC and spin/flavor otherwise)
//      Which is the contraction for the usual 2pt function    \sum_xy < \bar\psi_1(x) MSNK \psi_2(x) \bar\psi_2(y) MSRC \psi_1(y) > 
//                                                           = \sum_xy tr( P_1(y,x) MSNK P_2(x,y) MSRC ) =  \sum_xy tr( g5 P_1^dag(x,y) g5 MSNK P_2(x,y) MSRC )
//                                                           = \sum_xy tr( MSRC g5 P_1^dag(x,y) g5 MSNK P_2(x,y) )
//                                                           = \sum_x tr( MSRC g5 G_1^dag(x,y) g5 MSNK G_2(x,y) )
//      THUS  M1 = MSRC g5   and  M2 = g5 MSNK    (note the ordering).
//These same conventions apply to the wallsink bilinears below

void AlgGparityContract::contract_all_bilinears(const ContractionTypeAllBilinears &args, const int &conf_idx){
  std::ostringstream filestr; filestr << args.file << "." << conf_idx;
  std::string file = filestr.str();
  if(!UniqueID()) printf("Contracting all bilinears comprising propagators %s and %s and saving to file %s\n",args.prop_1, args.prop_2, file.c_str());
  
  if(GJP.Gparity()){
    ContractedBilinear<SpinColorFlavorMatrix> conbil;
    _multimom_helper<ContractedBilinear<SpinColorFlavorMatrix> >::add_momenta(conbil,args.momenta.momenta_val, args.momenta.momenta_len);
    conbil.write(args.prop_1, args.op1, args.prop_2, args.op2, file, AlgLattice(), binary_write);
  }else{
    ContractedBilinear<WilsonMatrix> conbil;
    _multimom_helper<ContractedBilinear<WilsonMatrix> >::add_momenta(conbil,args.momenta.momenta_val, args.momenta.momenta_len);
    conbil.write(args.prop_1, args.op1, args.prop_2, args.op2, file, AlgLattice(), binary_write);   
  }
}

void AlgGparityContract::contract_all_wallsink_bilinears_specific_momentum(const ContractionTypeAllWallSinkBilinearsSpecificMomentum &args, const int &conf_idx){
  std::ostringstream filestr; filestr << args.file << "." << conf_idx;
  std::string file = filestr.str();
  if(!UniqueID()) printf("Contracting all wallsink bilinears with chosen momenta comprising propagators %s and %s and saving to file %s\n",args.prop_1, args.prop_2, file.c_str());
  
  if(GJP.Gparity()){
    ContractedWallSinkBilinearSpecMomentum<SpinColorFlavorMatrix> conbil;
    if(args.cosine_sink==1) conbil.enableCosineSink();
    _multimom_helper<ContractedWallSinkBilinearSpecMomentum<SpinColorFlavorMatrix> >::add_momenta(conbil,args.momenta.momenta_val, args.momenta.momenta_len);
    conbil.write(args.prop_1, args.op1, args.prop_2, args.op2, file, AlgLattice(), binary_write);
  }else{
    ContractedWallSinkBilinearSpecMomentum<WilsonMatrix> conbil;
    if(args.cosine_sink==1) conbil.enableCosineSink();
    _multimom_helper<ContractedWallSinkBilinearSpecMomentum<WilsonMatrix> >::add_momenta(conbil,args.momenta.momenta_val, args.momenta.momenta_len);
    conbil.write(args.prop_1, args.op1, args.prop_2, args.op2, file, AlgLattice(), binary_write);
  }
}

void AlgGparityContract::contract_fourier_prop(const ContractionTypeFourierProp &args, const int &conf_idx){
  std::ostringstream filestr; filestr << args.file << "." << conf_idx;
  const char* file = filestr.str().c_str();

  if(GJP.Gparity()){
    FourierProp<SpinColorFlavorMatrix> ftprop;
    if(args.gauge_fix) ftprop.gaugeFixSink(true);
    else ftprop.gaugeFixSink(false);

    _multimom_helper<FourierProp<SpinColorFlavorMatrix> >::add_momenta(ftprop,args.momenta.momenta_val, args.momenta.momenta_len);
    ftprop.write(args.prop, file, AlgLattice());
  }else{
    FourierProp<WilsonMatrix> ftprop;
    if(args.gauge_fix) ftprop.gaugeFixSink(true);
    else ftprop.gaugeFixSink(false);

    _multimom_helper<FourierProp<WilsonMatrix> >::add_momenta(ftprop,args.momenta.momenta_val, args.momenta.momenta_len);
    ftprop.write(args.prop, file, AlgLattice());   
  }
}

void AlgGparityContract::contract_bilinear_vertex(const ContractionTypeBilinearVertex &args, const int &conf_idx){
  std::ostringstream filestr; filestr << args.file << "." << conf_idx;
  const char* file = filestr.str().c_str();

  if(GJP.Gparity()){
    PropagatorBilinear<SpinColorFlavorMatrix> bil;
    _multimom_helper<PropagatorBilinear<SpinColorFlavorMatrix> >::add_momenta(bil,args.momenta.momenta_val, args.momenta.momenta_len);
    bil.write(args.prop_1,OpDagger, args.prop_2, OpNone, file, AlgLattice());
  }else{
    PropagatorBilinear<WilsonMatrix> bil;
    _multimom_helper<PropagatorBilinear<WilsonMatrix> >::add_momenta(bil,args.momenta.momenta_val, args.momenta.momenta_len);
    bil.write(args.prop_1,OpDagger, args.prop_2, OpNone, file, AlgLattice());   
  }
}

void AlgGparityContract::contract_quadrilinear_vertex(const ContractionTypeQuadrilinearVertex &args, const int &conf_idx){
  std::ostringstream filestr; filestr << args.file << "." << conf_idx;
  const char* file = filestr.str().c_str();

  if(GJP.Gparity()){
    PropagatorQuadrilinear<SpinColorFlavorMatrix> quad;
    _multimom_helper<PropagatorQuadrilinear<SpinColorFlavorMatrix> >::add_momenta(quad,args.momenta.momenta_val, args.momenta.momenta_len);
    if(args.spin_structs.spin_structs_len == 0)
      quad.write(args.prop_1,OpDagger, args.prop_2, OpNone, 
		 args.prop_3,OpDagger, args.prop_4, OpNone,
		 file, AlgLattice());
    else{
      FILE *fp;
      if ((fp = Fopen(file, "w")) == NULL) ERR.FileW("AlgGparityContract","contract_quadrilinear_vertex %s",file);

      for(int i=0;i<args.spin_structs.spin_structs_len;i++){
	const QuadrilinearSpinStructure & s = args.spin_structs.spin_structs_val[i];
	for(int g1=0;g1< s.Gamma1.Gamma1_len; g1++) for(int s1=0;s1< s.Sigma1.Sigma1_len; s1++)
	for(int g2=0;g2< s.Gamma2.Gamma2_len; g2++) for(int s2=0;s2< s.Sigma2.Sigma2_len; s2++)
	   quad.write(args.prop_1,OpDagger, args.prop_2, OpNone, 
		      args.prop_3,OpDagger, args.prop_4, OpNone,
		      s.Gamma1.Gamma1_val[g1], s.Sigma1.Sigma1_val[s1],
		      s.Gamma2.Gamma2_val[g2], s.Sigma2.Sigma2_val[s2],
		      fp, AlgLattice());	  
      }
      Fclose(fp);
    }
  }else{
    PropagatorQuadrilinear<WilsonMatrix> quad;
    _multimom_helper<PropagatorQuadrilinear<WilsonMatrix> >::add_momenta(quad,args.momenta.momenta_val, args.momenta.momenta_len);
    quad.write(args.prop_1,OpDagger, args.prop_2, OpNone, 
	       args.prop_3,OpDagger, args.prop_4, OpNone,
	       file, AlgLattice());   
  }
}




CPS_END_NAMESPACE
