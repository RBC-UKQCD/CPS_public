#include<bitset>
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
#include <alg/do_arg.h>
#include <alg/no_arg.h>
#include <alg/common_arg.h>
#include <alg/hmd_arg.h>
#include <alg/alg_plaq.h>
#include <alg/alg_rnd_gauge.h>
#include <alg/threept_arg.h>
#include <alg/threept_prop_arg.h>
#include <alg/alg_threept.h>
#include <util/smalloc.h>
#if(0==1)
 #include <ReadLattice.h>
 #include <WriteLattice.h>
#endif

#include <util/ReadLatticePar.h>
#include <util/WriteLatticePar.h>

#include <util/command_line.h>
#include <sstream>

#include<unistd.h>
#include<config.h>

#include <alg/qpropw.h>
#include <alg/qpropw_arg.h>

#include <alg/alg_fix_gauge.h>
#include <alg/fix_gauge_arg.h>

#include <util/data_shift.h>

#include <alg/prop_attribute_arg.h>
#include <alg/gparity_contract_arg.h>
#include <alg/propmanager.h>
#include <alg/alg_gparitycontract.h>
#include <alg/prop_dft.h>

#include <util/spincolorflavormatrix.h>

#ifdef USE_BFM

//CK: these are redefined by BFM (to the same values)
#undef ND
#undef SPINOR_SIZE
#undef HALF_SPINOR_SIZE
#undef GAUGE_SIZE

#include<util/lattice/fbfm.h>
#include <chroma.h>
#include <pthread.h>
#endif


#ifdef USE_OMP
#include <omp.h>
#endif

using namespace std;
USING_NAMESPACE_CPS

#define SETUP_ARRAY(OBJ,ARRAYNAME,TYPE,SIZE)	\
  OBJ . ARRAYNAME . ARRAYNAME##_len = SIZE; \
  OBJ . ARRAYNAME . ARRAYNAME##_val = new TYPE [SIZE]

#define ELEM(OBJ,ARRAYNAME,IDX) OBJ . ARRAYNAME . ARRAYNAME##_val[IDX]


void print(const WilsonMatrix &w){
  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
      cps::Complex c = w(i,0,j,0);
      printf("(%.4f %.4f) ",c.real(),c.imag());
    }
    printf("\n");
  }
  printf("\n");
}

bool test_equals(const WilsonMatrix &a, const WilsonMatrix &b, const double &eps){
  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
      for(int aa=0;aa<3;aa++){
	for(int bb=0;bb<3;bb++){
	  cps::Complex ca = a(i,aa,j,bb);
	  cps::Complex cb = b(i,aa,j,bb);
	  if( fabs(ca.real()-cb.real()) > eps || fabs(ca.imag()-cb.imag()) > eps ) return false;
	}
      }
    }
  }
  return true;
}
bool test_equals(const SpinColorFlavorMatrix &a, const SpinColorFlavorMatrix &b, const double &eps){
  for(int i=0;i<2;i++)
    for(int j=0;j<2;j++)
      if(!test_equals( a(i,j), b(i,j) ,eps) ) return false;

  return true;
}
void print(const SpinColorFlavorMatrix &w){
  for(int i=0;i<2;i++){
    for(int j=0;j<2;j++){
      printf("Flav idx %d %d\n",i,j);
      print(w(i,j));      
    }
  }
}

void global_coord(const int &site, int *into_vec){
  int shift_x = GJP.XnodeCoor()*GJP.XnodeSites();
  int shift_y = GJP.YnodeCoor()*GJP.YnodeSites();
  int shift_z = GJP.ZnodeCoor()*GJP.ZnodeSites();
  int shift_t = GJP.TnodeCoor()*GJP.TnodeSites();
  //Local lattice dimensions:
  int size_x = GJP.XnodeSites();
  int size_y = GJP.YnodeSites();
  int size_z = GJP.ZnodeSites();
  int size_t = GJP.TnodeSites();
  int size_xy = size_x*size_y;
  int spatial_vol = (GJP.VolNodeSites()/GJP.TnodeSites()); // =size_x*size_y_size_z

  into_vec[3] = site/spatial_vol + shift_t;
  into_vec[2] = (site%spatial_vol)/size_xy + shift_z;
  into_vec[1] = (site%size_xy)/size_x + shift_y;
  into_vec[0] = site%size_x + shift_x;
}
const Float Pi_const(3.141592654);

Rcomplex sink_phasefac(const int *momphase, const int *pos,const bool &is_cconj=false){
  //momphase is the sum of the phase factors from the propagators forming the contraction
  //NOTE: In G-parity directions, momentum is discretised in odd units of \pi/2L rather than even/odd units of 2\pi/L (periodic/antiperiodic).
  Float pdotx = 0.0;
  Float sgn = 1.0;
  if(is_cconj) sgn = -1.0;

  for(int d=0;d<3;d++){
    Float mom_unit;
    if(GJP.Bc(d) == BND_CND_GPARITY) mom_unit = Pi_const/( (Float) 2*GJP.Nodes(d)*GJP.NodeSites(d));
    else if(GJP.Bc(d) == BND_CND_PRD) mom_unit = 2.0*Pi_const/( (Float) GJP.Nodes(d)*GJP.NodeSites(d));
    else if(GJP.Bc(d) == BND_CND_APRD) mom_unit = Pi_const/( (Float) GJP.Nodes(d)*GJP.NodeSites(d));
    else ERR.General("","sink_phasefac(int *,const int &)","Unknown boundary condition\n");
    
    pdotx += sgn*momphase[d]*pos[d]*mom_unit;
  }
  return Rcomplex(cos(pdotx),sin(pdotx));
}


Rcomplex sink_phasefac(int *momphase,const int &site){
  int pos[4]; global_coord(site,pos);
  return sink_phasefac(momphase,pos);
}
void sum_momphase(int *into, PropagatorContainer &prop, const bool &is_cconj){
  int propmom[3]; prop.momentum(propmom);
  if(is_cconj){
    for(int i=0;i<3;i++) into[i]-=propmom[i];
  }else{
    for(int i=0;i<3;i++) into[i]+=propmom[i];
  }
}

//Left multiply by a gamma matrix structure in QDP-style conventions:
//\Gamma(n) = \gamma_1^n1 \gamma_2^n2  \gamma_3^n3 \gamma_4^n4    where ni are bit fields: n4 n3 n2 n1
void qdp_gl(WilsonMatrix &wmat, const int &gidx){
  std::bitset<4> gmask(gidx);
  if(gmask == std::bitset<4>(15)){ wmat.gl(-5); return; }

  //CPS conventions \gamma^4 = gl(3) etc
  for(int i=3;i>=0;i--) if(gmask[i]) wmat.gl(i);
}
//Right multiply by a gamma matrix structure in QDP-style conventions:
//\Gamma(n) = \gamma_1^n1 \gamma_2^n2  \gamma_3^n3 \gamma_4^n4    where ni are bit fields: n4 n3 n2 n1 
void qdp_gr(WilsonMatrix &wmat, const int &gidx){
  std::bitset<4> gmask(gidx);
  if(gmask == std::bitset<4>(15)){ wmat.gr(-5); return; }

  //CPS conventions \gamma^4 = gl(3) etc
  for(int i=0;i<4;i++) if(gmask[i]) wmat.gr(i);
}

void qdp_gl(SpinColorFlavorMatrix &wmat, const int &gidx){
  std::bitset<4> gmask(gidx);
  if(gmask == std::bitset<4>(15)){ wmat.gl(-5); return; }

  //CPS conventions \gamma^4 = gl(3) etc
  for(int i=3;i>=0;i--) if(gmask[i]) wmat.gl(i);
}
void qdp_gr(SpinColorFlavorMatrix &wmat, const int &gidx){
  std::bitset<4> gmask(gidx);
  if(gmask == std::bitset<4>(15)){ wmat.gr(-5); return; }

  //CPS conventions \gamma^4 = gl(3) etc
  for(int i=0;i<4;i++) if(gmask[i]) wmat.gr(i);
}
//Coefficient when matrix is transposed or conjugated
Float qdp_gcoeff(const int &gidx, const bool &transpose, const bool &conj){
  if(gidx==2 || gidx==8 || gidx==15){ return 1.0; } //gamma^2, gamma^4 and gamma^5 are hermitian and real
  else if(gidx==1 || gidx==4){ //gamma^1 and gamma^3 are hermitian and imaginary
    if(transpose && conj) return 1.0;
    else return -1.0;
  } 
  std::bitset<4> gmask(gidx);
  Float out(1);
  if(transpose) out *= -1.0; //- sign for reordering 2 or 3 gamma matrices
    
  std::bitset<4> mask(1);
  for(int i=0;i<4;i++, mask<<1){
    if( (mask &gmask).any() ) out *= qdp_gcoeff((int)mask.to_ulong(),transpose,conj);
  }
  return out;
}

void calc_pi_minus(CorrelationFunction &corrfunc, const char *prop, Lattice &lattice, const int *desired_mom_x){
  int g1(15), g2(15);
  //pi^-
  /*<<(\bar u,d)*(\bar d,u)>>*/
  /*Require a "CorrelationFunction &corrfunc"*/
  /*Require propagator "PropagatorContainer &prop_src_y_u_d_eitherflav_pcon corresponding to \mathcal{G}^{[u/d] }_{x,y} with source of either flavour (full prop matrix is generated using single flavour source). Source must be real.*/
  PropagatorContainer &prop_src_y_u_d_eitherflav_pcon = PropManager::getProp(prop);

  /*Fourier transform on sink index x*/
  /*Require a 3-component array 'desired_mom_x' representing the required momentum at this sink position*/
  {
    QuarkMomCombination momcomb;
    momcomb.add_prop(prop_src_y_u_d_eitherflav_pcon, false);
    momcomb.add_prop(prop_src_y_u_d_eitherflav_pcon, false);
    bool desired_mom_available(momcomb.contains(desired_mom_x));
    /*Create an appropriate error message if !desired_mom_available*/
    if(!desired_mom_available){ std::cout << "pion correlation functions with cosine sources test momentum error\n"; exit(-1); }
  }

  corrfunc.setNcontractions(2);
  for(int x=0;x<GJP.VolNodeSites();x++){
    int x_pos_vec[4];
    global_coord(x,x_pos_vec);
  
    /*Get all SpinColorFlavorMatrices needed*/
    SpinColorFlavorMatrix prop_ud_snk_x_src_y_scfmat(prop_src_y_u_d_eitherflav_pcon , lattice, x);
  
    SpinColorFlavorMatrix prop_ud_snk_x_src_y_trans_scfmat(prop_src_y_u_d_eitherflav_pcon , lattice, x);
    prop_ud_snk_x_src_y_trans_scfmat.transpose();
  
    /*Starting contraction 0*/
    /*[{\rm tr}_{scf,0}\left\{C \Gamma[g2] F_1 F_\updownarrow \mathcal{G}^{[u/d] T}_{x,y} F_1 F_\updownarrow C \Gamma[g1] \mathcal{G}^{[u/d] }_{x,y}\right\}_{0}  ]*[f_\Gamma(g2,T) ]*/
  
    {
      Rcomplex contraction(1 , 0);
      contraction *= qdp_gcoeff(g2,true,false);
      contraction *= sink_phasefac(desired_mom_x,x_pos_vec,false);
    
      Rcomplex result_subdiag0(1.0);
      {
	SpinColorFlavorMatrix sdiag0_trset0_scfmat_prod(prop_ud_snk_x_src_y_trans_scfmat);
	sdiag0_trset0_scfmat_prod.pl(Fud);
	sdiag0_trset0_scfmat_prod.pl(F1);
	qdp_gl(sdiag0_trset0_scfmat_prod,g2);
	sdiag0_trset0_scfmat_prod.ccl(-1);
	sdiag0_trset0_scfmat_prod.pr(F1);
	sdiag0_trset0_scfmat_prod.pr(Fud);
	sdiag0_trset0_scfmat_prod.ccr(1);
	qdp_gr(sdiag0_trset0_scfmat_prod,g1);
	sdiag0_trset0_scfmat_prod *= prop_ud_snk_x_src_y_scfmat;
	Rcomplex sdiag0_trset0_cmplx;
	sdiag0_trset0_cmplx = sdiag0_trset0_scfmat_prod.Trace();
	result_subdiag0 *= sdiag0_trset0_cmplx;
      }
      contraction *= result_subdiag0;
    
    
      corrfunc(0,x_pos_vec[3]) += contraction;
    }
    /*Starting contraction 1*/
    /*[{\rm tr}_{scf,0}\left\{\Gamma[g2] C F_0 F_\updownarrow \mathcal{G}^{[u/d] T}_{x,y} F_1 F_\updownarrow C \Gamma[g1] \mathcal{G}^{[u/d] }_{x,y}\right\}_{0}  ]*/
  
    {
      Rcomplex contraction(1 , 0);
      contraction *= sink_phasefac(desired_mom_x,x_pos_vec,false);
    
      Rcomplex result_subdiag0(1.0);
      {
	SpinColorFlavorMatrix sdiag0_trset0_scfmat_prod(prop_ud_snk_x_src_y_trans_scfmat);
	sdiag0_trset0_scfmat_prod.pl(Fud);
	sdiag0_trset0_scfmat_prod.pl(F0);
	sdiag0_trset0_scfmat_prod.ccl(-1);
	qdp_gl(sdiag0_trset0_scfmat_prod,g2);
	sdiag0_trset0_scfmat_prod.pr(F1);
	sdiag0_trset0_scfmat_prod.pr(Fud);
	sdiag0_trset0_scfmat_prod.ccr(1);
	qdp_gr(sdiag0_trset0_scfmat_prod,g1);
	sdiag0_trset0_scfmat_prod *= prop_ud_snk_x_src_y_scfmat;
	Rcomplex sdiag0_trset0_cmplx;
	sdiag0_trset0_cmplx = sdiag0_trset0_scfmat_prod.Trace();
	result_subdiag0 *= sdiag0_trset0_cmplx;
      }
      contraction *= result_subdiag0;
    
    
      corrfunc(1,x_pos_vec[3]) += contraction;
    }
  }
  corrfunc.sumLattice();
}

void calc_pi_plus(CorrelationFunction &corrfunc, const char *prop, Lattice &lattice, const int *desired_mom_x){
  int g1(15), g2(15);
  //pi^+
  /*<<(\bar d,u)*(\bar u,d)>>*/
  /*Require a "CorrelationFunction &corrfunc"*/
  /*Require propagator "PropagatorContainer &prop_src_y_u_d_eitherflav_pcon corresponding to \mathcal{G}^{[u/d] }_{x,y} with source of either flavour (full prop matrix is generated using single flavour source). Source must be real.*/
  PropagatorContainer &prop_src_y_u_d_eitherflav_pcon = PropManager::getProp(prop);

  /*Fourier transform on sink index x*/
  /*Require a 3-component array 'desired_mom_x' representing the required momentum at this sink position*/
  {
    QuarkMomCombination momcomb;
    momcomb.add_prop(prop_src_y_u_d_eitherflav_pcon, true);
    momcomb.add_prop(prop_src_y_u_d_eitherflav_pcon, true);
    bool desired_mom_available(momcomb.contains(desired_mom_x));
    /*Create an appropriate error message if !desired_mom_available*/
    if(!desired_mom_available){ std::cout << "pion correlation functions with cosine sources test momentum error\n"; exit(-1); }
  }

  corrfunc.setNcontractions(2);
  for(int x=0;x<GJP.VolNodeSites();x++){
    int x_pos_vec[4];
    global_coord(x,x_pos_vec);
  
    /*Get all SpinColorFlavorMatrices needed*/
    SpinColorFlavorMatrix prop_ud_snk_x_src_y_hconj_scfmat(prop_src_y_u_d_eitherflav_pcon , lattice, x);
    prop_ud_snk_x_src_y_hconj_scfmat.hconj();
  
    SpinColorFlavorMatrix prop_ud_snk_x_src_y_cconj_scfmat(prop_src_y_u_d_eitherflav_pcon , lattice, x);
    prop_ud_snk_x_src_y_cconj_scfmat.cconj();
  
    /*Starting contraction 0*/
    /*[{\rm tr}_{scf,0}\left\{\gamma^5 C \Gamma[g1] \gamma^5 F_1 F_\updownarrow \mathcal{G}^{[u/d] *}_{x,y} F_1 F_\updownarrow \gamma^5 C \Gamma[g2] \gamma^5 \mathcal{G}^{[u/d] \dagger}_{x,y}\right\}_{0}  ]*[f_\Gamma(g1,T) ]*/
  
    {
      Rcomplex contraction(1 , 0);
      contraction *= qdp_gcoeff(g1,true,false);
      contraction *= sink_phasefac(desired_mom_x,x_pos_vec,false);
    
      Rcomplex result_subdiag0(1.0);
      {
	SpinColorFlavorMatrix sdiag0_trset0_scfmat_prod(prop_ud_snk_x_src_y_cconj_scfmat);
	sdiag0_trset0_scfmat_prod.pl(Fud);
	sdiag0_trset0_scfmat_prod.pl(F1);
	sdiag0_trset0_scfmat_prod.gl(-5);
	qdp_gl(sdiag0_trset0_scfmat_prod,g1);
	sdiag0_trset0_scfmat_prod.ccl(-1);
	sdiag0_trset0_scfmat_prod.gl(-5);
	sdiag0_trset0_scfmat_prod.pr(F1);
	sdiag0_trset0_scfmat_prod.pr(Fud);
	sdiag0_trset0_scfmat_prod.gr(-5);
	sdiag0_trset0_scfmat_prod.ccr(1);
	qdp_gr(sdiag0_trset0_scfmat_prod,g2);
	sdiag0_trset0_scfmat_prod.gr(-5);
	sdiag0_trset0_scfmat_prod *= prop_ud_snk_x_src_y_hconj_scfmat;
	Rcomplex sdiag0_trset0_cmplx;
	sdiag0_trset0_cmplx = sdiag0_trset0_scfmat_prod.Trace();
	result_subdiag0 *= sdiag0_trset0_cmplx;
      }
      contraction *= result_subdiag0;
    
    
      corrfunc(0,x_pos_vec[3]) += contraction;
    }
    /*Starting contraction 1*/
    /*[{\rm tr}_{scf,0}\left\{\gamma^5 C \Gamma[g1] \gamma^5 F_1 F_\updownarrow \mathcal{G}^{[u/d] *}_{x,y} F_0 F_\updownarrow \gamma^5 \Gamma[g2] \gamma^5 C \mathcal{G}^{[u/d] \dagger}_{x,y}\right\}_{0}  ]*[f_\Gamma(g1,T) ]*[f_\Gamma(g2,T) ]*/
  
    {
      Rcomplex contraction(1 , 0);
      contraction *= qdp_gcoeff(g1,true,false)*qdp_gcoeff(g2,true,false);
      contraction *= sink_phasefac(desired_mom_x,x_pos_vec,false);
    
      Rcomplex result_subdiag0(1.0);
      {
	SpinColorFlavorMatrix sdiag0_trset0_scfmat_prod(prop_ud_snk_x_src_y_cconj_scfmat);
	sdiag0_trset0_scfmat_prod.pl(Fud);
	sdiag0_trset0_scfmat_prod.pl(F1);
	sdiag0_trset0_scfmat_prod.gl(-5);
	qdp_gl(sdiag0_trset0_scfmat_prod,g1);
	sdiag0_trset0_scfmat_prod.ccl(-1);
	sdiag0_trset0_scfmat_prod.gl(-5);
	sdiag0_trset0_scfmat_prod.pr(F0);
	sdiag0_trset0_scfmat_prod.pr(Fud);
	sdiag0_trset0_scfmat_prod.gr(-5);
	qdp_gr(sdiag0_trset0_scfmat_prod,g2);
	sdiag0_trset0_scfmat_prod.gr(-5);
	sdiag0_trset0_scfmat_prod.ccr(1);
	sdiag0_trset0_scfmat_prod *= prop_ud_snk_x_src_y_hconj_scfmat;
	Rcomplex sdiag0_trset0_cmplx;
	sdiag0_trset0_cmplx = sdiag0_trset0_scfmat_prod.Trace();
	result_subdiag0 *= sdiag0_trset0_cmplx;
      }
      contraction *= result_subdiag0;
    
    
      corrfunc(1,x_pos_vec[3]) += contraction;
    }
  }
  corrfunc.sumLattice();
}


void calc_k0(CorrelationFunction &corrfunc, const char *prop_h, const char* prop_l, Lattice &lattice, const int *desired_mom_x){
  PropagatorContainer &prop_src_y_sprime_s_eitherflav_pcon = PropManager::getProp(prop_h);
  PropagatorContainer &prop_src_y_u_d_eitherflav_pcon = PropManager::getProp(prop_l);

  corrfunc.setNcontractions(1);
  for(int x=0;x<GJP.VolNodeSites();x++){
    int x_pos_vec[4];
    global_coord(x,x_pos_vec);
  
    /*Get all SpinColorFlavorMatrices needed*/
    SpinColorFlavorMatrix prop_sprimes_snk_x_src_y_hconj_scfmat(prop_src_y_sprime_s_eitherflav_pcon , lattice, x);
    prop_sprimes_snk_x_src_y_hconj_scfmat.hconj();
  
    SpinColorFlavorMatrix prop_ud_snk_x_src_y_scfmat(prop_src_y_u_d_eitherflav_pcon , lattice, x);
  
    {
      Rcomplex contraction(1 , 0);
      contraction *= sink_phasefac(desired_mom_x,x_pos_vec,false);
    
      Rcomplex result_subdiag0(1.0);
      {
	SpinColorFlavorMatrix sdiag0_trset0_scfmat_prod(prop_sprimes_snk_x_src_y_hconj_scfmat);
	sdiag0_trset0_scfmat_prod *= prop_ud_snk_x_src_y_scfmat;

	Rcomplex sdiag0_trset0_cmplx = sdiag0_trset0_scfmat_prod.Trace();
	result_subdiag0 *= sdiag0_trset0_cmplx;
      }
      contraction *= result_subdiag0;
      corrfunc(0,x_pos_vec[3]) += contraction;
    }
  }
  corrfunc.sumLattice();
}

void calc_a4k0(CorrelationFunction &corrfunc, const char *prop_h, const char* prop_l, Lattice &lattice, const int *desired_mom_x){
  PropagatorContainer &prop_src_y_sprime_s_eitherflav_pcon = PropManager::getProp(prop_h);
  PropagatorContainer &prop_src_y_u_d_eitherflav_pcon = PropManager::getProp(prop_l);

  corrfunc.setNcontractions(1);
  for(int x=0;x<GJP.VolNodeSites();x++){
    int x_pos_vec[4];
    global_coord(x,x_pos_vec);
  
    /*Get all SpinColorFlavorMatrices needed*/
    SpinColorFlavorMatrix prop_sprimes_snk_x_src_y_hconj_scfmat(prop_src_y_sprime_s_eitherflav_pcon , lattice, x);
    prop_sprimes_snk_x_src_y_hconj_scfmat.hconj();
  
    SpinColorFlavorMatrix prop_ud_snk_x_src_y_scfmat(prop_src_y_u_d_eitherflav_pcon , lattice, x);
  
    {
      Rcomplex contraction(1 , 0);
      contraction *= sink_phasefac(desired_mom_x,x_pos_vec,false);
    
      Rcomplex result_subdiag0(1.0);
      {
	SpinColorFlavorMatrix sdiag0_trset0_scfmat_prod(prop_sprimes_snk_x_src_y_hconj_scfmat);
	sdiag0_trset0_scfmat_prod.gr(3);
	sdiag0_trset0_scfmat_prod *= prop_ud_snk_x_src_y_scfmat;

	Rcomplex sdiag0_trset0_cmplx = sdiag0_trset0_scfmat_prod.Trace();
	result_subdiag0 *= sdiag0_trset0_cmplx;
      }
      contraction *= result_subdiag0;
      corrfunc(0,x_pos_vec[3]) += contraction;
    }
  }
  corrfunc.sumLattice();
}


void calc_etahat(CorrelationFunction &corrfunc, const char* prop_l, Lattice &lattice, const int *desired_mom_x, Rcomplex* txdep){
  PropagatorContainer &prop_src_y_u_d_eitherflav_pcon = PropManager::getProp(prop_l);

  for(int t=0;t<GJP.TnodeSites();t++){
    for(int x=0;x<GJP.XnodeSites();x++){
      txdep[t+GJP.TnodeSites()*x] = 0;
    }
  }

  corrfunc.setNcontractions(1);
  for(int x=0;x<GJP.VolNodeSites();x++){
    int x_pos_vec[4];
    global_coord(x,x_pos_vec);
  
    /*Get all SpinColorFlavorMatrices needed*/
    SpinColorFlavorMatrix prop_ud_snk_x_src_y_scfmat(prop_src_y_u_d_eitherflav_pcon , lattice, x);
    SpinColorFlavorMatrix prop_ud_snk_x_src_y_hconj_scfmat(prop_ud_snk_x_src_y_scfmat);
    prop_ud_snk_x_src_y_hconj_scfmat.hconj();
  
      {
      Rcomplex contraction(1 , 0);
      contraction *= sink_phasefac(desired_mom_x,x_pos_vec,false);
    
      Rcomplex result_subdiag0(1.0);
      {
	SpinColorFlavorMatrix sdiag0_trset0_scfmat_prod(prop_ud_snk_x_src_y_scfmat);
	sdiag0_trset0_scfmat_prod *= prop_ud_snk_x_src_y_hconj_scfmat;

	Rcomplex sdiag0_trset0_cmplx = sdiag0_trset0_scfmat_prod.Trace();
	result_subdiag0 *= sdiag0_trset0_cmplx;
      }
      contraction *= result_subdiag0;
      corrfunc(0,x_pos_vec[3]) += contraction;
     
      txdep[x_pos_vec[3]+GJP.TnodeSites()*x_pos_vec[0]] += result_subdiag0;
    }
  }
  corrfunc.sumLattice();
}



//CK: modified global sum routines used here such that QMP is not required
static void sum_double_array(Float* data, const int &size){
#ifdef USE_QMP
  QMP_sum_double_array((double *)data, size);  
#else
  slice_sum(data,size,99);
#endif
}


static void compute_coord(int x[4], const int size[4], int i)
{
    for(int j = 0; j < 4; ++j) {
        x[j] = i % size[j];
        i /= size[j];
    }
}
//Hantao's mres code
// void run_mres(const cps::QPropW &qp,
//               unsigned tsrc, const char fn[])

void run_mres(cps::QPropW &qp, std::vector<Rcomplex> &pion, std::vector<Rcomplex> &j5_q, const int &tsrc=0){
  if(qp.StoreMidprop() == 0) ERR.General("Hantao's code","mres","Requires midprop to be stored!\n");

  const int t_size_glb = GJP.TnodeSites() * GJP.Tnodes();
  
  const int lcl[4] = {GJP.XnodeSites(), GJP.YnodeSites(),
		      GJP.ZnodeSites(), GJP.TnodeSites(),};
  const int lcl_vol = lcl[0] * lcl[1] * lcl[2] * lcl[3];
  const int shift = GJP.TnodeSites() * GJP.TnodeCoor();
  
  pion.resize(t_size_glb, Rcomplex(0,0) );
  j5_q.resize(t_size_glb, Rcomplex(0,0) );

  //vector<Rcomplex> pion(t_size_glb, Rcomplex(0, 0));
  //vector<Rcomplex> j5_q(t_size_glb, Rcomplex(0, 0));

#pragma omp parallel
    {
        // threaded results
      std::vector<Rcomplex> pion_tmp(t_size_glb, Rcomplex(0, 0));
      std::vector<Rcomplex> j5_q_tmp(t_size_glb, Rcomplex(0, 0));

#pragma omp for
        for(int i = 0; i < lcl_vol; ++i) {
            int x[4];
            compute_coord(x, lcl, i);
            int t_glb = x[3] + shift;
            
            // J5 contraction (pion)
            WilsonMatrix p[2]  = {qp[i], qp[i]};
            p[1].hconj();
            // J5q contraction (midplane)
            WilsonMatrix q[2]  = {qp(i), qp(i)};
            q[1].hconj();
            
            pion_tmp[t_glb] += Trace(p[0], p[1]);
            j5_q_tmp[t_glb] += Trace(q[0], q[1]);
        } // sites
#pragma omp critical
        for(int t = 0; t < t_size_glb; ++t) {
            pion[t] += pion_tmp[(t+tsrc)%t_size_glb];
            j5_q[t] += j5_q_tmp[(t+tsrc)%t_size_glb];
        } // critical, for
    }//omp

    // FIXME
    if(GJP.Snodes() != 1) ERR.General("Hantaos code","mres","Requires snodes=1\n");
    sum_double_array((double *)pion.data(), 2 * t_size_glb);
    sum_double_array((double *)j5_q.data(), 2 * t_size_glb);

    // FILE *fp = Fopen(fn, "a+");
    // for(unsigned t = 0; t < t_size_glb; ++t) {
    //     Fprintf(fp, "%3u %3u %17.10e %17.10e %17.10e %17.10e\n", tsrc, t,
    //             pion[t].real(), pion[t].imag(),
    //             j5_q[t].real(), j5_q[t].imag());
    // }
    // Fclose(fp);
}




void nogparity_mrescalc_test(Lattice &lat){
  //Generate a simple propagator with midprop

  PropManager::clear();
  
  JobPropagatorArgs prop_args;
  SETUP_ARRAY(prop_args,props,PropagatorArg,1);
  
  PropagatorArg &parg = prop_args.props.props_val[0];
    
  parg.generics.tag = "prop";
  parg.generics.mass = 0.1;
  parg.generics.bc[0] = GJP.Xbc();
  parg.generics.bc[1] = GJP.Ybc();
  parg.generics.bc[2] = GJP.Zbc();
  parg.generics.bc[3] = GJP.Tbc();

  SETUP_ARRAY(parg,attributes,AttributeContainer,2);
    
  ELEM(parg,attributes,0).type = POINT_SOURCE_ATTR;
  PointSourceAttrArg &srcarg = ELEM(parg,attributes,0).AttributeContainer_u.point_source_attr;
  for(int i=0;i<4;i++) srcarg.pos[i] = 0;

  ELEM(parg,attributes,1).type = STORE_MIDPROP_ATTR;

  PropManager::setup(prop_args);
  PropManager::calcProps(lat);

  CommonArg c_arg;
  AlgGparityContract con(lat,c_arg);

  CorrelationFunction pion("pion",1,CorrelationFunction::THREADED);
  CorrelationFunction j5_q("j5q",1,CorrelationFunction::THREADED);
  ContractionTypeMres meas; meas.prop = "prop"; meas.file = "";

  {
    const char* first_prop = prop_args.props.props_val[0].generics.tag;
    PropagatorContainer &prop_pcon = PropManager::getProp(first_prop);
    QPropW & qp = prop_pcon.getProp(lat);
    printf("First prop %s site 0: ",first_prop);
    for(int ii=0;ii<18;ii++) printf("%f ",((Float*)&qp[0])[ii]);
    printf("\n");
  }


  //Test mres measurement against Hantao's version
  con.measure_mres(meas,pion,j5_q);

  std::vector<Rcomplex> pion_hantao, j5_q_hantao;
  run_mres(PropManager::getProp("prop").getProp(lat), pion_hantao,j5_q_hantao);
  
  bool fail(false);
  for(int t=0;t<GJP.Tnodes()*GJP.TnodeSites();++t){
    if( abs(pion(0,t).real() - pion_hantao[t].real()) > 1e-12 ){ printf("mres pion match fail real %.14e %.14e\n",pion(0,t).real(),pion_hantao[t].real()); fail = true; }
    else printf("mres pion match pass real %.14e %.14e\n",pion(0,t).real(),pion_hantao[t].real());

    if( abs(pion(0,t).imag() - pion_hantao[t].imag()) > 1e-12 ){ printf("mres pion match fail imag %.14e %.14e\n",pion(0,t).imag(),pion_hantao[t].imag()); fail = true; }
    else printf("mres pion match pass imag %.14e %.14e\n",pion(0,t).imag(),pion_hantao[t].imag());

    if( abs(j5_q(0,t).real() - j5_q_hantao[t].real()) > 1e-12 ){ printf("mres j5_q match fail real %.14e %.14e\n",j5_q(0,t).real(),j5_q_hantao[t].real()); fail = true; }
    else printf("mres j5_q match pass real %.14e %.14e\n",j5_q(0,t).real(),j5_q_hantao[t].real());

    if( abs(j5_q(0,t).imag() - j5_q_hantao[t].imag()) > 1e-12 ){ printf("mres j5_q match fail imag %.14e %.14e\n",j5_q(0,t).imag(),j5_q_hantao[t].imag()); fail = true; }
    else printf("mres j5_q match pass imag %.14e %.14e\n",j5_q(0,t).imag(),j5_q_hantao[t].imag());
  }
  if(fail){
    if(!UniqueID()) printf("mres comparison failed\n"); 
    exit(-1);
  }else if(!UniqueID()) printf("mres comparison passed\n"); 
}



#ifndef USE_BFM
void init_bfm(int *argc, char **argv[]){}
#else
void init_bfm(int *argc, char **argv[])
{
    const char* fname = "init_bfm()";

    Chroma::initialize(argc, argv);
    multi1d<int> nrow(Nd);

    for(int i = 0; i< Nd; ++i)
        nrow[i] = GJP.Sites(i);

    Layout::setLattSize(nrow);
    Layout::create();

    //This is regular DWF:
    // Fbfm::bfm_args[0].solver = DWF;
    // Fbfm::bfm_args[0].precon_5d = 1;

    //This is Mobius:
    Fbfm::bfm_args[0].solver = HmCayleyTanh;
    Fbfm::bfm_args[0].precon_5d = 0;
    
    //This is 4D even-odd preconditioned DWF:
    //It is identical to Mobius with Mobius parameter = 1
    // Fbfm::bfm_args[0].solver = HtCayleyTanh;
    // Fbfm::bfm_args[0].precon_5d = 0;


    Fbfm::bfm_args[0].Ls = GJP.SnodeSites();
    Fbfm::bfm_args[0].M5 = GJP.DwfHeight();
    Fbfm::bfm_args[0].mass = 0.1;
    Fbfm::bfm_args[0].residual = 1e-8;
    Fbfm::bfm_args[0].max_iter = 10000;
    Fbfm::bfm_args[0].Csw = 0.0;

    Fbfm::bfm_args[0].node_latt[0] = QDP::Layout::subgridLattSize()[0];
    Fbfm::bfm_args[0].node_latt[1] = QDP::Layout::subgridLattSize()[1];
    Fbfm::bfm_args[0].node_latt[2] = QDP::Layout::subgridLattSize()[2];
    Fbfm::bfm_args[0].node_latt[3] = QDP::Layout::subgridLattSize()[3];

    multi1d<int> procs = QDP::Layout::logicalSize();

    Fbfm::bfm_args[0].local_comm[0] = procs[0] > 1 ? 0 : 1;
    Fbfm::bfm_args[0].local_comm[1] = procs[1] > 1 ? 0 : 1;
    Fbfm::bfm_args[0].local_comm[2] = procs[2] > 1 ? 0 : 1;
    Fbfm::bfm_args[0].local_comm[3] = procs[3] > 1 ? 0 : 1;

    Fbfm::bfm_args[0].ncoor[0] = 0;
    Fbfm::bfm_args[0].ncoor[1] = 0;
    Fbfm::bfm_args[0].ncoor[2] = 0;
    Fbfm::bfm_args[0].ncoor[3] = 0;

    // mobius_scale = b + c in Andrew's notation
    bfmarg::mobius_scale = 4.;
#if TARGET == BGQ
    bfmarg::Threads(64);
#else
    bfmarg::Threads(1);
#endif

    bfmarg::Reproduce(0);
    bfmarg::ReproduceChecksum(0);
    bfmarg::ReproduceMasterCheck(0);
    bfmarg::Verbose(1);

    Fbfm::use_mixed_solver = true;

    VRB.Result("", fname, "init_bfm finished successfully\n");
}
#endif



int main(int argc,char *argv[])
{
  Start(&argc,&argv);
  CommandLine::is(argc,argv);

  bool gparity_X(false);
  bool gparity_Y(false);

  int arg0 = CommandLine::arg_as_int(0);
  printf("Arg0 is %d\n",arg0);
  if(arg0==0){
    gparity_X=true;
    printf("Doing G-parity test in X direction\n");
  }else if(arg0==1){
    printf("Doing G-parity test in X and Y directions\n");
    gparity_X = true;
    gparity_Y = true;
  }else if(arg0==2){
    printf("Doing standard periodic lattice\n");
  }

  bool dbl_latt_storemode(false);
  bool save_config(false);
  bool load_config(false);
  bool load_lrg(false);
  bool save_lrg(false);
  char *load_config_file;
  char *save_config_file;
  char *save_lrg_file;
  char *load_lrg_file;
  bool gauge_fix(false);
  bool verbose(false);

  int size[] = {2,2,2,2,2};
  int seed = 83209;

  bool do_mobius(false);

  int i=2;
  while(i<argc){
    char* cmd = argv[i];  
    if( strncmp(cmd,"-save_config",15) == 0){
      if(i==argc-1){ printf("-save_config requires an argument\n"); exit(-1); }
      save_config=true;
      save_config_file = argv[i+1];
      i+=2;
    }else if( strncmp(cmd,"-load_config",15) == 0){
      if(i==argc-1){ printf("-save_config requires an argument\n"); exit(-1); }
      load_config=true;
      load_config_file = argv[i+1];
      i+=2;
    }else if( strncmp(cmd,"-latt",10) == 0){
      if(i>argc-6){
	printf("Did not specify enough arguments for 'latt' (require 5 dimensions)\n"); exit(-1);
      }
      size[0] = CommandLine::arg_as_int(i); //CommandLine ignores zeroth input arg (i.e. executable name)
      size[1] = CommandLine::arg_as_int(i+1);
      size[2] = CommandLine::arg_as_int(i+2);
      size[3] = CommandLine::arg_as_int(i+3);
      size[4] = CommandLine::arg_as_int(i+4);
      i+=6;
    }else if( strncmp(cmd,"-save_double_latt",20) == 0){
      dbl_latt_storemode = true;
      i++;
    }else if( strncmp(cmd,"-mobius",20) == 0){
      do_mobius = true;
      i++;
    }else if( strncmp(cmd,"-load_lrg",15) == 0){
      if(i==argc-1){ printf("-load_lrg requires an argument\n"); exit(-1); }
      load_lrg=true;
      load_lrg_file = argv[i+1];
      i+=2;
    }else if( strncmp(cmd,"-save_lrg",15) == 0){
      if(i==argc-1){ printf("-save_lrg requires an argument\n"); exit(-1); }
      save_lrg=true;
      save_lrg_file = argv[i+1];
      i+=2;
    }else if( strncmp(cmd,"-gauge_fix",15) == 0){
      gauge_fix=true;
      i++;   
    }else if( strncmp(cmd,"-verbose",15) == 0){
      verbose=true;
      i++;
    }else if( strncmp(cmd,"-seed",15) == 0){
      std::stringstream os; os << argv[i+1]; os >> seed;
      printf("Seed set to %d\n",seed);
      i+=2;
    }else{
      if(UniqueID()==0) printf("Unrecognised argument: %s\n",cmd);
      exit(-1);
    }
  }
  
  printf("Lattice size is %d %d %d %d\n",size[0],size[1],size[2],size[3],size[4]);

  DoArg do_arg;
  do_arg.x_sites = size[0];
  do_arg.y_sites = size[1];
  do_arg.z_sites = size[2];
  do_arg.t_sites = size[3];
  do_arg.s_sites = size[4];
  do_arg.x_node_sites = 0;
  do_arg.y_node_sites = 0;
  do_arg.z_node_sites = 0;
  do_arg.t_node_sites = 0;
  do_arg.s_node_sites = 0;
  do_arg.x_nodes = 0;
  do_arg.y_nodes = 0;
  do_arg.z_nodes = 0;
  do_arg.t_nodes = 0;
  do_arg.s_nodes = 0;
  do_arg.updates = 0;
  do_arg.measurements = 0;
  do_arg.measurefreq = 0;
  do_arg.cg_reprod_freq = 10;
  do_arg.x_bc = BND_CND_PRD;
  do_arg.y_bc = BND_CND_PRD;
  do_arg.z_bc = BND_CND_PRD;
  do_arg.t_bc = BND_CND_APRD;
  do_arg.start_conf_kind = START_CONF_ORD;
  do_arg.start_conf_load_addr = 0x0;
  do_arg.start_seed_kind = START_SEED_INPUT;
  do_arg.start_seed_filename = "../rngs/ckpoint_rng.0";
  do_arg.start_conf_filename = "../configurations/ckpoint_lat.0";
  do_arg.start_conf_alloc_flag = 6;
  do_arg.wfm_alloc_flag = 2;
  do_arg.wfm_send_alloc_flag = 2;
  do_arg.start_seed_value = seed;
  do_arg.beta =   2.25;
  do_arg.c_1 =   -3.3100000000000002e-01;
  do_arg.u0 =   1.0000000000000000e+00;
  do_arg.dwf_height =   1.8000000000000000e+00;
  do_arg.dwf_a5_inv =   1.0000000000000000e+00;
  do_arg.power_plaq_cutoff =   0.0000000000000000e+00;
  do_arg.power_plaq_exponent = 0;
  do_arg.power_rect_cutoff =   0.0000000000000000e+00;
  do_arg.power_rect_exponent = 0;
  do_arg.verbose_level = -1202;//VERBOSE_DEBUG_LEVEL; //-1202;
  do_arg.checksum_level = 0;
  do_arg.exec_task_list = 0;
  do_arg.xi_bare =   1.0000000000000000e+00;
  do_arg.xi_dir = 3;
  do_arg.xi_v =   1.0000000000000000e+00;
  do_arg.xi_v_xi =   1.0000000000000000e+00;
  do_arg.clover_coeff =   0.0000000000000000e+00;
  do_arg.clover_coeff_xi =   0.0000000000000000e+00;
  do_arg.xi_gfix =   1.0000000000000000e+00;
  do_arg.gfix_chkb = 1;
  do_arg.asqtad_KS =   0.0000000000000000e+00;
  do_arg.asqtad_naik =   0.0000000000000000e+00;
  do_arg.asqtad_3staple =   0.0000000000000000e+00;
  do_arg.asqtad_5staple =   0.0000000000000000e+00;
  do_arg.asqtad_7staple =   0.0000000000000000e+00;
  do_arg.asqtad_lepage =   0.0000000000000000e+00;
  do_arg.p4_KS =   0.0000000000000000e+00;
  do_arg.p4_knight =   0.0000000000000000e+00;
  do_arg.p4_3staple =   0.0000000000000000e+00;
  do_arg.p4_5staple =   0.0000000000000000e+00;
  do_arg.p4_7staple =   0.0000000000000000e+00;
  do_arg.p4_lepage =   0.0000000000000000e+00;

  if(verbose) do_arg.verbose_level = VERBOSE_DEBUG_LEVEL;

  if(gparity_X) do_arg.x_bc = BND_CND_GPARITY;
  if(gparity_Y) do_arg.y_bc = BND_CND_GPARITY;

  GJP.Initialize(do_arg);

  SerialIO::dbl_latt_storemode = dbl_latt_storemode;
  
  LRG.Initialize();

  if(load_lrg){
    if(UniqueID()==0) printf("Loading RNG state from %s\n",load_lrg_file);
    LRG.Read(load_lrg_file,32);
  }
  if(save_lrg){
    if(UniqueID()==0) printf("Writing RNG state to %s\n",save_lrg_file);
    LRG.Write(save_lrg_file,32);
  }

  Lattice *lattice_p;
  if(do_mobius){
#ifndef USE_BFM
    ERR.General("","main","Cannot use Fbfm without bfm library!\n");
#endif
    init_bfm(&argc,&argv);
    lattice_p = new GnoneFbfm;
  }else{
    lattice_p = new GwilsonFdwf;
  }
  Lattice &lattice = *lattice_p;
  
  if(!load_config){
    printf("Creating gauge field\n"); fflush(stdout);
    lattice.SetGfieldDisOrd(); //unit gauge
    printf("Gauge checksum = %d\n", lattice.CheckSum());  fflush(stdout);
  }else{
    ReadLatticeParallel readLat;
    if(UniqueID()==0) printf("Reading: %s (NERSC-format)\n",load_config_file);
    readLat.read(lattice,load_config_file);
    if(UniqueID()==0) printf("Config read.\n");
  }
  if(do_mobius){
    GnoneFbfm* l =dynamic_cast<GnoneFbfm*>(lattice_p);
    l->BondCond(); //re-import gauge field to internal bfm object
  }

  if(save_config){
    if(UniqueID()==0) printf("Saving config to %s\n",save_config_file);

    QioArg wt_arg(save_config_file,0.001);
    
    wt_arg.ConcurIONumber=32;
    WriteLatticeParallel wl;
    wl.setHeader("disord_id","disord_label",0);
    wl.write(lattice,wt_arg);
    
    if(!wl.good()) ERR.General("main","()","Failed write lattice %s",save_config_file);

    if(UniqueID()==0) printf("Config written.\n");
  }

  if(gauge_fix){
    lattice.FixGaugeAllocate(FIX_GAUGE_COULOMB_T);
    lattice.FixGauge(1e-06,2000);
    if(!UniqueID()){ printf("Gauge fixing finished\n"); fflush(stdout); }
  }

  if(!gparity_X && !gparity_Y){
    if(!UniqueID()) printf("Running periodic lattice tests\n");

    nogparity_mrescalc_test(lattice);

    if(gauge_fix) lattice.FixGaugeFree();
    
    if(UniqueID()==0){
      printf("Main job complete\n"); 
      fflush(stdout);
    }
  
    return 0;
  }

  SpinColorFlavorMatrix one; _PropagatorBilinear_helper<SpinColorFlavorMatrix>::unit_matrix(one);

  bool fail(false);

  printf("Doing gcoeff test\n");
  //test gcoeff
  fail = false;
  for(int i=0;i<16;i++){
    SpinColorFlavorMatrix mat(one); 
    AlgGparityContract::qdp_gr(mat,i);  
    
    SpinColorFlavorMatrix mat_T(mat); 
    mat_T.transpose();
    mat_T *= AlgGparityContract::qdp_gcoeff(i,true,false);    //trans,conj

    if(!test_equals(mat_T,mat,1e-12)){ printf("Gamma^T err %d:\n",i); print(mat); print(mat_T); fail=true; }

    SpinColorFlavorMatrix mat_star(mat); 
    mat_star.cconj();
    mat_star *= AlgGparityContract::qdp_gcoeff(i,false,true); 
    
    if(!test_equals(mat_star,mat,1e-12)){ printf("Gamma^* err %d:\n",i); print(mat); print(mat_star); fail=true; }

    SpinColorFlavorMatrix mat_dag(mat); 
    mat_dag.hconj();
    mat_dag *= AlgGparityContract::qdp_gcoeff(i,true,true); 
    
    if(!test_equals(mat_dag,mat,1e-12)){ printf("Gamma^dag err %d:\n",i); print(mat); print(mat_dag); fail=true; }
  }
  if(fail){
    printf("gcoeff test fail\n");
    exit(0);
  }else{
    printf("gcoeff test passed\n");
  }

  printf("Doing pauli_coeff test\n");
  fail = false;
  const static FlavorMatrixType fmap[4] = {sigma0, sigma1, sigma2, sigma3};   

  for(int i=0;i<4;i++){
    SpinColorFlavorMatrix mat(one); 
    mat.pr(fmap[i]);
    
    SpinColorFlavorMatrix mat_T(mat); 
    mat_T.transpose();
    mat_T *= AlgGparityContract::pauli_coeff(i,true,false);    //trans,conj

    if(!test_equals(mat_T,mat,1e-12)){ printf("Sigma^T err %d:\n",i); print(mat); print(mat_T); fail=true; }

    SpinColorFlavorMatrix mat_star(mat); 
    mat_star.cconj();
    mat_star *= AlgGparityContract::pauli_coeff(i,false,true); 
    
    if(!test_equals(mat_star,mat,1e-12)){ printf("Sigma^* err %d:\n",i); print(mat); print(mat_star); fail=true; }

    SpinColorFlavorMatrix mat_dag(mat); 
    mat_dag.hconj();
    mat_dag *= AlgGparityContract::pauli_coeff(i,true,true); 
    
    if(!test_equals(mat_dag,mat,1e-12)){ printf("Sigma^dag err %d:\n",i); print(mat); print(mat_dag); fail=true; }
  }
  if(fail){
    printf("pauli_coeff test fail\n");
    exit(0);
  }else{
    printf("pauli_coeff test passed\n");
  }


  if(0){
    //test pion correlation functions with cosine sources is reproduced by new code using ContractionQuarkMomCombination
     printf("Starting pion correlation functions with cosine sources test\n");
    PropManager::clear();

    JobPropagatorArgs prop_args2;
    SETUP_ARRAY(prop_args2,props,PropagatorArg,2);
  
    char* names[2] = {"prop","prop_H"};
    BndCndType bndcnd[2] = {BND_CND_APRD,BND_CND_APRD};
    int psign[2] = {1,1};
    Float masses[2] = {0.1,0.5};

    for(int i=0;i<2;i++){
      PropagatorArg &parg = prop_args2.props.props_val[i];
    
      parg.generics.tag = names[i];
      parg.generics.mass = masses[i];
      parg.generics.bc[0] = GJP.Xbc();
      parg.generics.bc[1] = GJP.Ybc();
      parg.generics.bc[2] = GJP.Zbc();
      parg.generics.bc[3] = bndcnd[i];

      SETUP_ARRAY(parg,attributes,AttributeContainer,5);
    
      ELEM(parg,attributes,0).type = WALL_SOURCE_ATTR;
      WallSourceAttrArg &srcarg = ELEM(parg,attributes,0).AttributeContainer_u.wall_source_attr;
      srcarg.t = 0;

      ELEM(parg,attributes,1).type = GPARITY_FLAVOR_ATTR;
      GparityFlavorAttrArg &gparg = ELEM(parg,attributes,1).AttributeContainer_u.gparity_flavor_attr;
      gparg.flavor = 0;

      ELEM(parg,attributes,2).type = CG_ATTR;
      CGAttrArg &cgattr = ELEM(parg,attributes,2).AttributeContainer_u.cg_attr;
      cgattr.max_num_iter = 5000;
      cgattr.stop_rsd = 1e-08;
      cgattr.true_rsd = 1e-08;

      ELEM(parg,attributes,3).type = MOMENTUM_ATTR;
      MomentumAttrArg & momarg = ELEM(parg,attributes,3).AttributeContainer_u.momentum_attr;
      for(int ii=0;ii<3;ii++) 
	if(GJP.Bc(ii)==BND_CND_GPARITY) momarg.p[ii] = 1; //units of pi/2L
	else momarg.p[ii] = 0;

      ELEM(parg,attributes,4).type = MOM_COS_ATTR;
    }
    if(UniqueID()==0) printf("prop_args contains %d propagators\n", prop_args2.props.props_len);

    PropManager::setup(prop_args2);   
    PropManager::calcProps(lattice);

    int desired_mom_x[3] = {0,0,0};
    for(int i=0;i<3;i++) 
      if(GJP.Bc(i)==BND_CND_GPARITY) desired_mom_x[i] = 2;
      else desired_mom_x[i] = 0;

    CorrelationFunction pi_plus("pi^+");
    calc_pi_plus(pi_plus, names[0], lattice, desired_mom_x);
    
    const static Float Pi_const = 3.1415926535897932384626433832795;
    std::vector<Float> momvec(3); 
    for(int i=0;i<3;i++){
      momvec[i] = desired_mom_x[i]* Pi_const/(2.0*GJP.Nodes(i)*GJP.NodeSites(i));
    }
    
    ContractedBilinear<SpinColorFlavorMatrix> conbil;
    conbil.add_momentum(momvec);

    TimeStamp::set_file("timestamp.log");

    TimeStamp::stamp("Starting version 0");
    std::vector<Rcomplex> result_1_sigma3 = conbil.getBilinear(lattice,momvec,"prop",PropDFT::None,"prop",PropDFT::Dagger,
							       0,0,0,3);
    std::vector<Rcomplex> result_sigma3_sigma3 = conbil.getBilinear(lattice,momvec,"prop",PropDFT::None,"prop",PropDFT::Dagger,
								    0,3,0,3);
    TimeStamp::stamp("End of version 0");


    conbil.clear();
    conbil.add_momentum(momvec);

    TimeStamp::stamp("Starting version 1");
    std::vector<Rcomplex> resultv2_1_sigma3 = conbil.getBilinear(lattice,momvec,"prop",PropDFT::None,"prop",PropDFT::Dagger,
								 0,0,0,3,1);
    std::vector<Rcomplex> resultv2_sigma3_sigma3 = conbil.getBilinear(lattice,momvec,"prop",PropDFT::None,"prop",PropDFT::Dagger,
								      0,3,0,3,1);
    TimeStamp::stamp("End of version 1");

    TimeStamp::close_file();

    for(int t=0;t<GJP.TnodeSites();t++){
      Rcomplex pp = pi_plus(0,t) + pi_plus(1,t);
      Rcomplex pp2 = 0.5*result_1_sigma3[t] - 0.5*result_sigma3_sigma3[t];
      Rcomplex pp3 = 0.5*resultv2_1_sigma3[t] - 0.5*resultv2_sigma3_sigma3[t];
      printf("pi+:  orig (%f,%f) : v1 (%f,%f)  v2 (%f,%f)\n",pp.real(),pp.imag(),pp2.real(),pp2.imag(),pp3.real(),pp3.imag());      
    }
    ContractedBilinear<SpinColorFlavorMatrix> conbil_GlGh;
    conbil_GlGh.add_momentum(momvec);
    //conbil_GlGh.calculateBilinears(lattice,momvec,"prop",PropDFT::None,"prop_H",PropDFT::None);
    conbil_GlGh.calculateBilinears(lattice,"prop",PropDFT::None,"prop_H",PropDFT::None);
    
    ContractedBilinear<SpinColorFlavorMatrix> conbil_GhdagGldag;
    conbil_GhdagGldag.add_momentum(momvec);
    //conbil_GlGh.calculateBilinears(lattice,momvec,"prop",PropDFT::None,"prop_H",PropDFT::None);
    conbil_GhdagGldag.calculateBilinears(lattice,"prop_H",PropDFT::Dagger,"prop",PropDFT::Dagger);  

    ContractedBilinear<SpinColorFlavorMatrix> conbil_GldagGhdag;
    conbil_GldagGhdag.add_momentum(momvec);
    //conbil_GldagGhdag.calculateBilinears(lattice,momvec,"prop",PropDFT::Dagger,"prop_H",PropDFT::Dagger);  
    conbil_GldagGhdag.calculateBilinears(lattice,"prop",PropDFT::Dagger,"prop_H",PropDFT::Dagger);  

    ContractedBilinear<SpinColorFlavorMatrix> conbil_GltransGh;
    conbil_GltransGh.add_momentum(momvec);
    //conbil_GltransGh.calculateBilinears(lattice,momvec,"prop",PropDFT::Transpose,"prop_H",PropDFT::None);
    conbil_GltransGh.calculateBilinears(lattice,"prop",PropDFT::Transpose,"prop_H",PropDFT::None);
    
    ContractedBilinear<SpinColorFlavorMatrix> conbil_GlGhtrans;
    conbil_GlGhtrans.add_momentum(momvec);
    //conbil_GlGhtrans.calculateBilinears(lattice,momvec,"prop",PropDFT::None,"prop_H",PropDFT::Transpose);
    conbil_GlGhtrans.calculateBilinears(lattice,"prop",PropDFT::None,"prop_H",PropDFT::Transpose);

    ContractedBilinear<SpinColorFlavorMatrix> conbil_GhtransGl;
    conbil_GhtransGl.add_momentum(momvec);
    //conbil_GhtransGl.calculateBilinears(lattice,momvec,"prop_H",PropDFT::Transpose,"prop",PropDFT::None);
    conbil_GhtransGl.calculateBilinears(lattice,"prop_H",PropDFT::Transpose,"prop",PropDFT::None);

    ContractedBilinear<SpinColorFlavorMatrix> conbil_GlstarGhstar;
    conbil_GlstarGhstar.add_momentum(momvec);
    //conbil_GlstarGhstar.calculateBilinears(lattice,momvec,"prop",PropDFT::Conj,"prop_H",PropDFT::Conj);  
    conbil_GlstarGhstar.calculateBilinears(lattice,"prop",PropDFT::Conj,"prop_H",PropDFT::Conj);  



    //test trace identities that allow us to skip the calculation of certain combinations
    fail = false;
    for(int s1=0;s1<16;s1++){
      for(int f1=0;f1<4;f1++){
	for(int s2=0;s2<16;s2++){
	  for(int f2=0;f2<4;f2++){
	    //calculate tr(A Gl B Gh) by various means
	    std::vector<Rcomplex> GlGh = conbil_GlGh.getBilinear(lattice,momvec,"prop",PropDFT::None,"prop_H",PropDFT::None,
								 s1,f1,s2,f2);
	    std::vector<Rcomplex> from_GhdagGldag = conbil_GhdagGldag.getBilinear(lattice,momvec,"prop",PropDFT::None,"prop_H",PropDFT::None,
										  s1,f1,s2,f2);
	    std::vector<Rcomplex> from_GldagGhdag = conbil_GldagGhdag.getBilinear(lattice,momvec,"prop",PropDFT::None,"prop_H",PropDFT::None,
										  s1,f1,s2,f2);
	    std::vector<Rcomplex> from_GlstarGhstar = conbil_GlstarGhstar.getBilinear(lattice,momvec,"prop",PropDFT::None,"prop_H",PropDFT::None,
										      s1,f1,s2,f2);

	    std::vector<Rcomplex> GltransGh = conbil_GltransGh.getBilinear(lattice,momvec,"prop",PropDFT::Transpose,"prop_H",PropDFT::None,
									   s1,f1,s2,f2);
	    std::vector<Rcomplex> from_GlGhtrans = conbil_GlGhtrans.getBilinear(lattice,momvec,"prop",PropDFT::Transpose,"prop_H",PropDFT::None,
										s1,f1,s2,f2);

	    std::vector<Rcomplex> from_GhtransGl = conbil_GhtransGl.getBilinear(lattice,momvec,"prop",PropDFT::Transpose,"prop_H",PropDFT::None,
										s1,f1,s2,f2);

	    for(int t=0;t<GJP.TnodeSites();t++){
	      if( fabs( from_GhdagGldag[t].real()-GlGh[t].real() ) > 1e-12 || fabs( from_GhdagGldag[t].imag()-GlGh[t].imag() ) > 1e-12){
		printf("Gl Gh from  Gh^dag Gl^dag err %d %d %d %d: %d  -- (%f,%f) (%f,%f)\n",s1,f1,s2,f2,t,GlGh[t].real(),GlGh[t].imag(),
		       from_GhdagGldag[t].real(), from_GhdagGldag[t].imag());
		fail = true;
	      }
	      if( fabs( from_GldagGhdag[t].real()-GlGh[t].real() ) > 1e-12 || fabs( from_GldagGhdag[t].imag()-GlGh[t].imag() ) > 1e-12){
		printf("Gl Gh from  Gl^dag Gh^dag err %d %d %d %d: %d  -- (%f,%f) (%f,%f)\n",s1,f1,s2,f2,t,GlGh[t].real(),GlGh[t].imag(),
		       from_GldagGhdag[t].real(), from_GldagGhdag[t].imag());
		fail = true;
	      }
	      if( fabs( from_GlstarGhstar[t].real()-GlGh[t].real() ) > 1e-12 || fabs( from_GlstarGhstar[t].imag()-GlGh[t].imag() ) > 1e-12){
		printf("Gl Gh from  Gl^dag Gh^dag err %d %d %d %d: %d  -- (%f,%f) (%f,%f)\n",s1,f1,s2,f2,t,GlGh[t].real(),GlGh[t].imag(),
		       from_GlstarGhstar[t].real(), from_GlstarGhstar[t].imag());
		fail = true;
	      }
	      if( fabs( from_GlGhtrans[t].real()-GltransGh[t].real() ) > 1e-12 || fabs( from_GlGhtrans[t].imag()-GltransGh[t].imag() ) > 1e-12){
		printf("Gl^T Gh from  Gl Gh^T err %d %d %d %d: %d  -- (%f,%f) (%f,%f)\n",s1,f1,s2,f2,t,GltransGh[t].real(),GltransGh[t].imag(),
		       from_GlGhtrans[t].real(), from_GlGhtrans[t].imag());
		fail = true;
	      }	      
	      if( fabs( from_GhtransGl[t].real()-GltransGh[t].real() ) > 1e-12 || fabs( from_GhtransGl[t].imag()-GltransGh[t].imag() ) > 1e-12){
		printf("Gl^T Gh from  Gh^T Gl err %d %d %d %d: %d  -- (%f,%f) (%f,%f)\n",s1,f1,s2,f2,t,GltransGh[t].real(),GltransGh[t].imag(),
		       from_GhtransGl[t].real(), from_GhtransGl[t].imag());
		fail = true;
	      }
	    }
	  }
	}
      }
    }
    if(fail){
      printf("bilinear test fail\n");
      exit(0);
    }else{
      printf("bilinear test passed\n");
    }

  }



  {
    //test PropagatorBilinear
    PropManager::clear();

    JobPropagatorArgs prop_args2;
    SETUP_ARRAY(prop_args2,props,PropagatorArg,2);
  
    char* names[2] = {"prop","prop_H"};
    BndCndType bndcnd[2] = {BND_CND_APRD,BND_CND_APRD};
    int psign[2] = {1,1};
    Float masses[2] = {0.1,0.5};

    for(int i=0;i<2;i++){
      PropagatorArg &parg = prop_args2.props.props_val[i];
    
      parg.generics.tag = names[i];
      parg.generics.mass = masses[i];
      parg.generics.bc[0] = GJP.Xbc();
      parg.generics.bc[1] = GJP.Ybc();
      parg.generics.bc[2] = GJP.Zbc();
      parg.generics.bc[3] = bndcnd[i];

      SETUP_ARRAY(parg,attributes,AttributeContainer,5);
    
      ELEM(parg,attributes,0).type = WALL_SOURCE_ATTR;
      WallSourceAttrArg &srcarg = ELEM(parg,attributes,0).AttributeContainer_u.wall_source_attr;
      srcarg.t = 0;

      ELEM(parg,attributes,1).type = GPARITY_FLAVOR_ATTR;
      GparityFlavorAttrArg &gparg = ELEM(parg,attributes,1).AttributeContainer_u.gparity_flavor_attr;
      gparg.flavor = 0;

      ELEM(parg,attributes,2).type = CG_ATTR;
      CGAttrArg &cgattr = ELEM(parg,attributes,2).AttributeContainer_u.cg_attr;
      cgattr.max_num_iter = 5000;
      cgattr.stop_rsd = 1e-08;
      cgattr.true_rsd = 1e-08;

      ELEM(parg,attributes,3).type = MOMENTUM_ATTR;
      MomentumAttrArg & momarg = ELEM(parg,attributes,3).AttributeContainer_u.momentum_attr;
      for(int ii=0;ii<3;ii++) 
	if(GJP.Bc(ii)==BND_CND_GPARITY) momarg.p[ii] = 1; //units of pi/2L
	else momarg.p[ii] = 0;

      ELEM(parg,attributes,4).type = MOM_COS_ATTR;
    }
    if(UniqueID()==0) printf("prop_args contains %d propagators\n", prop_args2.props.props_len);

    PropManager::setup(prop_args2);   
    PropManager::calcProps(lattice);

    printf("Propagators calculated\n");
    
    int desired_mom_x[3] = {0,0,0};
    for(int i=0;i<3;i++) 
      if(GJP.Bc(i)==BND_CND_GPARITY) desired_mom_x[i] = 2;
      else desired_mom_x[i] = 0;

    const static Float Pi_const = 3.1415926535897932384626433832795;
    std::vector<Float> momvec(3); 
    for(int i=0;i<3;i++){
      momvec[i] = desired_mom_x[i]* Pi_const/(2.0*GJP.Nodes(i)*GJP.NodeSites(i));
    }
    
    PropagatorBilinear<SpinColorFlavorMatrix> pb;
    
    for(int Gamma=0;Gamma<16;Gamma++){
      for(int Sigma=0;Sigma<4;Sigma++){
	printf("Gamma = %d Sigma = %d\n",Gamma,Sigma);

	pb.clear();
	pb.add_momentum(momvec);
	printf("Calculating GlGh\n"); fflush(stdout);
	const std::vector<SpinColorFlavorMatrix> & GlGh = pb.getBilinear(lattice,momvec,
									 "prop",PropDFT::None,
									 "prop_H",PropDFT::None,
									 Gamma,Sigma);
	
	pb.clear();
	pb.add_momentum(momvec);
	printf("Calculating GltrGhtr\n"); fflush(stdout);
	const std::vector<SpinColorFlavorMatrix> & GhtrGltr = pb.getBilinear(lattice,momvec,
									 "prop_H",PropDFT::Transpose,
									 "prop",PropDFT::Transpose,
									 Gamma,Sigma);
	printf("Attempting to calculated GlGh from GhtrGltr using existing result\n");
	const std::vector<SpinColorFlavorMatrix> & GlGh_fromGhtrGltr = pb.getBilinear(lattice,momvec,
										      "prop",PropDFT::None,
										      "prop_H",PropDFT::None,
										      Gamma,Sigma);

      }
    }




  }








  // {
  //   //test pion correlation functions with zero-momentum wall sources
  //   printf("Starting pion correlation functions with zero-momentum wall sources test\n");
  //   PropManager::clear();

  //   JobPropagatorArgs prop_args2;
  //   SETUP_ARRAY(prop_args2,props,PropagatorArg,1);
  
  //   char* names[1] = {"prop"};
  //   BndCndType bndcnd[1] = {BND_CND_APRD};

  //   for(int i=0;i<1;i++){
  //     PropagatorArg &parg = prop_args2.props.props_val[i];
    
  //     parg.generics.tag = names[i];
  //     parg.generics.mass = 0.1;
  //     parg.generics.bc[0] = GJP.Xbc();
  //     parg.generics.bc[1] = GJP.Ybc();
  //     parg.generics.bc[2] = GJP.Zbc();
  //     parg.generics.bc[3] = bndcnd[i];

  //     SETUP_ARRAY(parg,attributes,AttributeContainer,3);
    
  //     ELEM(parg,attributes,0).type = WALL_SOURCE_ATTR;
  //     WallSourceAttrArg &srcarg = ELEM(parg,attributes,0).AttributeContainer_u.wall_source_attr;
  //     srcarg.t = 0;

  //     ELEM(parg,attributes,1).type = GPARITY_FLAVOR_ATTR;
  //     GparityFlavorAttrArg &gparg = ELEM(parg,attributes,1).AttributeContainer_u.gparity_flavor_attr;
  //     gparg.flavor = 0;

  //     ELEM(parg,attributes,2).type = CG_ATTR;
  //     CGAttrArg &cgattr = ELEM(parg,attributes,2).AttributeContainer_u.cg_attr;
  //     cgattr.max_num_iter = 5000;
  //     cgattr.stop_rsd = 1e-08;
  //     cgattr.true_rsd = 1e-08;
  //   }
  //   if(UniqueID()==0) printf("prop_args contains %d propagators\n", prop_args2.props.props_len);

  //   PropManager::setup(prop_args2);   
  //   PropManager::calcProps(lattice);

  //   {
  //     printf("0 momentum\n");

  //     int desired_mom_x[3] = {0,0,0};

  //     CorrelationFunction pi_minus("pi^-");
  //     CorrelationFunction pi_plus("pi^+");

  //     calc_pi_minus(pi_minus, names[0], lattice, desired_mom_x);
  //     calc_pi_plus(pi_plus, names[0], lattice, desired_mom_x);

  //     if(!UniqueID()){
  // 	for(int c=0;c<2;c++){
  // 	  printf("Contraction %d\n",c);
  // 	  for(int t=0;t<GJP.TnodeSites();t++){

  // 	    Rcomplex &pm = pi_minus(c,t);
  // 	    Rcomplex &pp = pi_plus(c,t);

  // 	    printf("pi+:pi-  (%f,%f) : (%f,%f)\n",pp.real(),pp.imag(),pm.real(),pm.imag());
  // 	  }
  // 	}
  //     }
  //   }
  // }

  // {
  //   //test real and imaginary parts of K^0 propagator
  //   printf("Starting test real and imaginary parts of K^0 and \hat eta propagator with point source\n");
  //   PropManager::clear();

  //   JobPropagatorArgs prop_args2;
  //   SETUP_ARRAY(prop_args2,props,PropagatorArg,2);
  
  //   char* names[2] = {"prop_h","prop_l"};
  //   BndCndType bndcnd[2] = {BND_CND_APRD,BND_CND_APRD};
  //   double mass[2] = {0.4,0.1};

  //   for(int i=0;i<2;i++){
  //     PropagatorArg &parg = prop_args2.props.props_val[i];
    
  //     parg.generics.tag = names[i];
  //     parg.generics.mass = mass[i];
  //     parg.generics.bc[0] = GJP.Xbc();
  //     parg.generics.bc[1] = GJP.Ybc();
  //     parg.generics.bc[2] = GJP.Zbc();
  //     parg.generics.bc[3] = bndcnd[i];

  //     SETUP_ARRAY(parg,attributes,AttributeContainer,3);
    
  //     ELEM(parg,attributes,0).type = POINT_SOURCE_ATTR;
  //     PointSourceAttrArg &srcarg = ELEM(parg,attributes,0).AttributeContainer_u.point_source_attr;
  //     srcarg.pos[0] = 0;
  //     srcarg.pos[1] = 0;
  //     srcarg.pos[2] = 0;
  //     srcarg.pos[3] = 0;

  //     ELEM(parg,attributes,1).type = GPARITY_FLAVOR_ATTR;
  //     GparityFlavorAttrArg &gparg = ELEM(parg,attributes,1).AttributeContainer_u.gparity_flavor_attr;
  //     gparg.flavor = 0;

  //     ELEM(parg,attributes,2).type = CG_ATTR;
  //     CGAttrArg &cgattr = ELEM(parg,attributes,2).AttributeContainer_u.cg_attr;
  //     cgattr.max_num_iter = 5000;
  //     cgattr.stop_rsd = 1e-08;
  //     cgattr.true_rsd = 1e-08;
  //   }
  //   if(UniqueID()==0) printf("prop_args contains %d propagators\n", prop_args2.props.props_len);

  //   PropManager::setup(prop_args2);   
  //   PropManager::calcProps(lattice);

  //   {
  //     printf("0 momentum\n");

  //     int desired_mom_x[3] = {0,0,0};
  //     CorrelationFunction kzero("K^0");
  //     calc_k0(kzero, names[0], names[1], lattice, desired_mom_x);

  //     if(!UniqueID()){
  // 	for(int t=0;t<GJP.TnodeSites();t++){
  // 	  Rcomplex &pp = kzero(0,t);
  // 	  printf("K^0[%d] (%f,%f)\n",t,pp.real(),pp.imag());
  // 	}
  //     }

  //     CorrelationFunction etahat("\\hat eta");
  //     Rcomplex etahat_txdep[GJP.TnodeSites()*GJP.XnodeSites()];
  //     calc_etahat(etahat,names[1],lattice,desired_mom_x, &(etahat_txdep[0]));

  //     if(!UniqueID()){
  // 	for(int t=0;t<GJP.TnodeSites();t++){
  // 	  Rcomplex &pp = etahat(0,t);
  // 	  printf("\\hat eta[%d] (%f,%f)\n",t,pp.real(),pp.imag());
  // 	}
  // 	printf("X dep:\n");
  // 	for(int t=0;t<GJP.TnodeSites();t++){
  // 	  printf("t=%d\n",t);
  // 	  for(int x=0;x<GJP.XnodeSites();x++){
  // 	    Rcomplex &pp = etahat_txdep[t+GJP.TnodeSites()*x];
  // 	    printf("\\hat eta[x=%d] (%f,%f)\n",x,pp.real(),pp.imag());
  // 	  }
  // 	}

  //     }
  //   }

  //   {
  //     printf("\\sqrt{n}\pi/L momentum\n");

  //     int desired_mom_x[3] = {0,0,0};
  //     if(GJP.Xbc()==BND_CND_GPARITY) desired_mom_x[0] = 2; //units of pi/2L
  //     if(GJP.Ybc()==BND_CND_GPARITY) desired_mom_x[1] = 2;

  //     CorrelationFunction kzero("K^0");
  //     calc_k0(kzero, names[0], names[1], lattice, desired_mom_x);

  //     if(!UniqueID()){
  // 	for(int t=0;t<GJP.TnodeSites();t++){
  // 	  Rcomplex &pp = kzero(0,t);
  // 	  printf("K^0[%d] (%f,%f)\n",t,pp.real(),pp.imag());
  // 	}
  //     }

  //     Rcomplex etahat_txdep[GJP.TnodeSites()*GJP.XnodeSites()];
  //     CorrelationFunction etahat("\\hat eta");
  //     calc_etahat(etahat,names[1],lattice,desired_mom_x,etahat_txdep);
  //     if(!UniqueID()){
  // 	for(int t=0;t<GJP.TnodeSites();t++){
  // 	  Rcomplex &pp = etahat(0,t);
  // 	  printf("\\hat eta[%d] (%f,%f)\n",t,pp.real(),pp.imag());
  // 	}
  // 	printf("X dep:\n");
  // 	for(int t=0;t<GJP.TnodeSites();t++){
  // 	  printf("t=%d\n",t);
  // 	  for(int x=0;x<GJP.XnodeSites();x++){
  // 	    Rcomplex &pp = etahat_txdep[t+GJP.TnodeSites()*x];
  // 	    printf("\\hat eta[x=%d] (%f,%f)\n",x,pp.real(),pp.imag());
  // 	  }
  // 	}
	
  //     }
  //   }
   
  // }


  // {
  //   //test real and imaginary parts of K^0 propagator
  //   printf("Starting test real and imaginary parts of K^0 propagator with cosine source\n");
  //   PropManager::clear();

  //   JobPropagatorArgs prop_args2;
  //   SETUP_ARRAY(prop_args2,props,PropagatorArg,2);
  
  //   char* names[2] = {"prop_h","prop_l"};
  //   BndCndType bndcnd[2] = {BND_CND_APRD,BND_CND_APRD};
  //   double mass[2] = {0.4,0.1};
    
  //   for(int i=0;i<2;i++){
  //     PropagatorArg &parg = prop_args2.props.props_val[i];
    
  //     parg.generics.tag = names[i];
  //     parg.generics.mass = mass[i];
  //     parg.generics.bc[0] = GJP.Xbc();
  //     parg.generics.bc[1] = GJP.Ybc();
  //     parg.generics.bc[2] = GJP.Zbc();
  //     parg.generics.bc[3] = bndcnd[i];

  //     SETUP_ARRAY(parg,attributes,AttributeContainer,5);
    
  //     ELEM(parg,attributes,0).type = WALL_SOURCE_ATTR;
  //     WallSourceAttrArg &srcarg = ELEM(parg,attributes,0).AttributeContainer_u.wall_source_attr;
  //     srcarg.t = 0;

  //     ELEM(parg,attributes,1).type = GPARITY_FLAVOR_ATTR;
  //     GparityFlavorAttrArg &gparg = ELEM(parg,attributes,1).AttributeContainer_u.gparity_flavor_attr;
  //     gparg.flavor = 0;

  //     ELEM(parg,attributes,2).type = CG_ATTR;
  //     CGAttrArg &cgattr = ELEM(parg,attributes,2).AttributeContainer_u.cg_attr;
  //     cgattr.max_num_iter = 5000;
  //     cgattr.stop_rsd = 1e-08;
  //     cgattr.true_rsd = 1e-08;

  //     ELEM(parg,attributes,3).type = MOMENTUM_ATTR;
  //     MomentumAttrArg & momarg = ELEM(parg,attributes,3).AttributeContainer_u.momentum_attr;
  //     for(int ii=0;ii<3;ii++) 
  // 	if(GJP.Bc(ii)==BND_CND_GPARITY) momarg.p[ii] = 1; //units of pi/2L
  // 	else momarg.p[ii] = 0;

  //     ELEM(parg,attributes,4).type = MOM_COS_ATTR;
  //   }
  //   if(UniqueID()==0) printf("prop_args contains %d propagators\n", prop_args2.props.props_len);

  //   PropManager::setup(prop_args2);   
  //   PropManager::calcProps(lattice);

  //   {
  //     printf("0 momentum\n");

  //     int desired_mom_x[3] = {0,0,0};
  //     CorrelationFunction kzero("K^0");
  //     calc_k0(kzero, names[0], names[1], lattice, desired_mom_x);

  //     if(!UniqueID()){
  // 	for(int t=0;t<GJP.TnodeSites();t++){
  // 	  Rcomplex &pp = kzero(0,t);
  // 	  printf("K^0[%d] (%f,%f)\n",t,pp.real(),pp.imag());
  // 	}
  //     }
  //   }

  // }


  // {
  //   //test real and imaginary parts of <A_4 K^0> Green's function
  //   printf("Starting test real and imaginary parts of <A_4 K^0> Green's function with cosine source\n");
  //   PropManager::clear();

  //   JobPropagatorArgs prop_args2;
  //   SETUP_ARRAY(prop_args2,props,PropagatorArg,2);
  
  //   char* names[2] = {"prop_h","prop_l"};
  //   BndCndType bndcnd[2] = {BND_CND_APRD,BND_CND_APRD};
  //   double mass[2] = {0.4,0.1};
    
  //   for(int i=0;i<2;i++){
  //     PropagatorArg &parg = prop_args2.props.props_val[i];
    
  //     parg.generics.tag = names[i];
  //     parg.generics.mass = mass[i];
  //     parg.generics.bc[0] = GJP.Xbc();
  //     parg.generics.bc[1] = GJP.Ybc();
  //     parg.generics.bc[2] = GJP.Zbc();
  //     parg.generics.bc[3] = bndcnd[i];

  //     SETUP_ARRAY(parg,attributes,AttributeContainer,5);
    
  //     ELEM(parg,attributes,0).type = WALL_SOURCE_ATTR;
  //     WallSourceAttrArg &srcarg = ELEM(parg,attributes,0).AttributeContainer_u.wall_source_attr;
  //     srcarg.t = 0;

  //     ELEM(parg,attributes,1).type = GPARITY_FLAVOR_ATTR;
  //     GparityFlavorAttrArg &gparg = ELEM(parg,attributes,1).AttributeContainer_u.gparity_flavor_attr;
  //     gparg.flavor = 0;

  //     ELEM(parg,attributes,2).type = CG_ATTR;
  //     CGAttrArg &cgattr = ELEM(parg,attributes,2).AttributeContainer_u.cg_attr;
  //     cgattr.max_num_iter = 5000;
  //     cgattr.stop_rsd = 1e-08;
  //     cgattr.true_rsd = 1e-08;

  //     ELEM(parg,attributes,3).type = MOMENTUM_ATTR;
  //     MomentumAttrArg & momarg = ELEM(parg,attributes,3).AttributeContainer_u.momentum_attr;
  //     for(int ii=0;ii<3;ii++) 
  // 	if(GJP.Bc(ii)==BND_CND_GPARITY) momarg.p[ii] = 1; //units of pi/2L
  // 	else momarg.p[ii] = 0;

  //     ELEM(parg,attributes,4).type = MOM_COS_ATTR;
  //   }
  //   if(UniqueID()==0) printf("prop_args contains %d propagators\n", prop_args2.props.props_len);

  //   PropManager::setup(prop_args2);   
  //   PropManager::calcProps(lattice);

  //   {
  //     printf("0 momentum\n");

  //     int desired_mom_x[3] = {0,0,0};
  //     CorrelationFunction a4kzero("K^0");
  //     calc_a4k0(a4kzero, names[0], names[1], lattice, desired_mom_x);

  //     if(!UniqueID()){
  // 	for(int t=0;t<GJP.TnodeSites();t++){
  // 	  Rcomplex &pp = a4kzero(0,t);
  // 	  printf("A_4 K^0[%d] (%f,%f)\n",t,pp.real(),pp.imag());
  // 	}
  //     }
  //   }

  // }


  if(gauge_fix) lattice.FixGaugeFree();

  if(UniqueID()==0){
    printf("Main job complete\n"); 
    fflush(stdout);
  }
  
  return 0;
}



