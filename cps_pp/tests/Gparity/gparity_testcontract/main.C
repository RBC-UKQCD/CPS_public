//In this code we demonstrate the relationship 
// \prop^(1-j,1-k) (y,z) = (-1)^|j-k| \gamma^5 C [ \prop^(j,k) (y,z) ]* C^\dagger \gamma^5

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

#include <util/spincolorflavormatrix.h>


#ifdef USE_OMP
#include <omp.h>
#endif

using namespace std;
USING_NAMESPACE_CPS

void print(const WilsonMatrix &w){
  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
      Complex c = w(i,0,j,0);
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
	  Complex ca = a(i,aa,j,bb);
	  Complex cb = b(i,aa,j,bb);
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
  for(int i=0;i<2;i++)
    for(int j=0;j<2;j++){
      printf("flavour indices %d %d\n",i,j);
      print(w(i,j));
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

void test_props(const JobPropagatorArgs &prop_args, Lattice &latt){
  //assume 2 props, one point source at 0 on flavour 0 and one point source at 0 on flavour 1
  PropagatorContainer &q_f0_pc = PropManager::getProp(prop_args.props.props_val[0].generics.tag);
  PropagatorContainer &q_f1_pc = PropManager::getProp(prop_args.props.props_val[1].generics.tag);
  
  QPropW &q_f0_qpw = q_f0_pc.getProp(latt);
  QPropW &q_f1_qpw = q_f1_pc.getProp(latt);
  
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

  for(int i=0;i<GJP.VolNodeSites();i++){
    int pos[4];
    pos[3] = i/spatial_vol + shift_t;
    pos[2] = (i%spatial_vol)/size_xy + shift_z;
    pos[1] = (i%size_xy)/size_x + shift_y;
    pos[0] = i%size_x + shift_x;

    printf("(%d,%d,%d,%d):\n",pos[0],pos[1],pos[2],pos[3]);

    WilsonMatrix g_0_0 = q_f0_qpw.SiteMatrix(i,0);
    WilsonMatrix g_1_0 = q_f0_qpw.SiteMatrix(i,1);
    WilsonMatrix g_0_1 = q_f1_qpw.SiteMatrix(i,0);
    WilsonMatrix g_1_1 = q_f1_qpw.SiteMatrix(i,1);

    //note, due to CPS oddity, A.ccl(1) = C^-1 A,  A.ccl(-1) = C A
    //                whereas, A.ccr(1) = A C,          A.ccr(-1) = A C^-1

    // \prop^(1-j,1-k) (y,z) = (-1)^|j-k| \gamma^5 C [ \prop^(j,k) (y,z) ]* C^\dagger \gamma^5
    {
      printf("Testing 1,1 = f(0,0):");
      WilsonMatrix f_0_0 = g_0_0;
      f_0_0.cconj();
      f_0_0.ccl(-1).gl(-5).ccr(-1).gr(-5);
      bool result = test_equals(g_1_1,f_0_0,1e-08);
      if(!result){
	printf(" false\n");
	print(g_1_1);
	print(f_0_0);
	exit(-1);
      }else printf(" true\n");
    }
    {
      printf("Testing 1,0 = f(0,1):");
      WilsonMatrix f_0_1 = g_0_1;
      f_0_1.cconj();
      f_0_1.ccl(-1).gl(-5).ccr(-1).gr(-5)*=Complex(-1.0,0.0);
      bool result = test_equals(g_1_0,f_0_1,1e-08);
      if(!result){
	printf(" false\n");
	print(g_1_0);
	print(f_0_1);
	exit(-1);
      }else printf(" true\n");
    }
    {
      printf("Testing 0,1 = f(1,0):");
      WilsonMatrix f_1_0 = g_1_0;
      f_1_0.cconj();
      f_1_0.ccl(-1).gl(-5).ccr(-1).gr(-5)*=Complex(-1.0,0.0);
      bool result = test_equals(g_0_1,f_1_0,1e-08);
      if(!result){
	printf(" false\n");
	print(g_0_1);
	print(f_1_0);
	exit(-1);
      }else printf(" true\n");
    }
    {
      printf("Testing 0,0 = f(1,1):");
      WilsonMatrix f_1_1 = g_1_1;
      f_1_1.cconj();
      f_1_1.ccl(-1).gl(-5).ccr(-1).gr(-5);
      bool result = test_equals(g_0_0,f_1_1,1e-08);
      if(!result){
	printf(" false\n");
	print(g_0_0);
	print(f_1_1);
	exit(-1);
      }else printf(" true\n");
    }
  }
    
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
    printf("Doing G-parity HMC test in X direction\n");
  }else if(arg0==1){
    printf("Doing G-parity HMC test in X and Y directions\n");
    gparity_X = true;
    gparity_Y = true;
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


  GwilsonFdwf lattice;
					       
  if(!load_config){
    printf("Creating gauge field\n");
    lattice.SetGfieldDisOrd(); //unit gauge
  }else{
    ReadLatticeParallel readLat;
    if(UniqueID()==0) printf("Reading: %s (NERSC-format)\n",load_config_file);
    readLat.read(lattice,load_config_file);
    if(UniqueID()==0) printf("Config read.\n");
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

#define SETUP_ARRAY(OBJ,ARRAYNAME,TYPE,SIZE)	\
  OBJ . ARRAYNAME . ARRAYNAME##_len = SIZE; \
  OBJ . ARRAYNAME . ARRAYNAME##_val = new TYPE [SIZE]

#define ELEM(OBJ,ARRAYNAME,IDX) OBJ . ARRAYNAME . ARRAYNAME##_val[IDX]

  //test with a point source on both flavours

  JobPropagatorArgs prop_args;
  SETUP_ARRAY(prop_args,props,PropagatorArg,2);
  
  char* names[2] = {"prop_f0","prop_f1"};
  
  for(int i=0;i<2;i++){
    PropagatorArg &parg = prop_args.props.props_val[i];
    
    parg.generics.tag = names[i];
    parg.generics.mass = 0.1;
    parg.generics.bc[0] = GJP.Xbc();
    parg.generics.bc[1] = GJP.Ybc();

    SETUP_ARRAY(parg,attributes,AttributeContainer,3);
    
    ELEM(parg,attributes,0).type = POINT_SOURCE_ATTR;
    PointSourceAttrArg &srcarg = ELEM(parg,attributes,0).AttributeContainer_u.point_source_attr;
    for(int j=0;j<4;j++) srcarg.pos[j] = 0;

    ELEM(parg,attributes,1).type = GPARITY_FLAVOR_ATTR;
    GparityFlavorAttrArg &gparg = ELEM(parg,attributes,1).AttributeContainer_u.gparity_flavor_attr;
    gparg.flavor = i;

    ELEM(parg,attributes,2).type = CG_ATTR;
    CGAttrArg &cgattr = ELEM(parg,attributes,2).AttributeContainer_u.cg_attr;
    cgattr.max_num_iter = 5000;
    cgattr.stop_rsd = 1e-08;
    cgattr.true_rsd = 1e-08;
  }
  if(UniqueID()==0) printf("prop_args contains %d propagators\n", prop_args.props.props_len);

  PropManager::setup(prop_args);

  PropagatorContainer &q_f0_pc = PropManager::getProp(prop_args.props.props_val[0].generics.tag);
  PropagatorContainer &q_f1_pc = PropManager::getProp(prop_args.props.props_val[1].generics.tag);

  QPropW & q_f1_qpw = q_f1_pc.getProp(lattice);
  QPropW & q_f0_qpw = q_f0_pc.getProp(lattice);

  for(int i=0;i<GJP.VolNodeSites();i++){
    printf("Starting site %d\n",i);
    SpinColorFlavorMatrix f(q_f0_pc,lattice,i);

    //first test generation
    
    //test f[0][1]
    WilsonMatrix &f01 = f(0,1);
    WilsonMatrix &q_0_1 = q_f1_qpw.SiteMatrix(i,0);
    
    printf("Testing f[0][1]\n");
    if(!test_equals(f01, q_0_1, 1e-08)){
      printf("Comparison of f[0][1] and q_0_1 failed on site %d:\n",i);
      print(f01);
      print(q_0_1);
      exit(-1);
    }
    printf("Passed\n");

    //test f[1][1]
    WilsonMatrix &f11 = f(1,1);
    WilsonMatrix &q_1_1 = q_f1_qpw.SiteMatrix(i,1);
    
    printf("Testing f[1][1]\n");
    if(!test_equals(f11, q_1_1, 1e-08)){
      printf("Comparison of f[1][1] and q_1_1 failed on site %d:\n",i);
      print(f11);
      print(q_1_1);
      exit(-1);
    }
    printf("Passed\n");

    //prop reln test
    {
      SpinColorFlavorMatrix g(f);
      g.cconj();
      g.pl(sigma3).ccl(-1).gl(-5).pl(Fud) .pr(sigma3).ccr(-1).gr(-5).pr(Fud);
      
      printf("Testing prop relation\n");
      if(!test_equals(g, f, 1e-08)){
	printf("Prop relation test failed on site %d:\n",i);
	print(g);
	print(f);
	exit(-1);
      }
      printf("Passed\n");
    }

    //prop transpose reln
    if(i==0){ //assuming point source at 0!
      SpinColorFlavorMatrix fT(f);
      fT.transpose();
      SpinColorFlavorMatrix g(f);
      g.pl(Fud).ccl(-1).pl(sigma3) .pr(Fud).ccr(-1).pr(sigma3);
      
      printf("Testing prop transpose relation\n");
      if(!test_equals(fT, g, 1e-08)){
	printf("Prop relation test failed on site %d:\n",i);
	print(fT);
	print(g);
	exit(-1);
      }
      printf("Passed\n");
    }

    {
      //Compare first contraction of pi^+ propagator
    
      //{\rm tr}_{scf,0}\left\{\gamma^5 C F_\updownarrow F_1 \mathcal{G}^{[u/d] T}_{x,z} F_\updownarrow F_1 \gamma^5 C \mathcal{G}^{[u/d] }_{x,z}\right\}_{0}

      //vs.

      //{\rm tr}\gamma^5 C G^{(1,0)}_{x,z}\gamma^5 C G^{(0,1) T} 
    
      SpinColorFlavorMatrix c1(f);
      c1.transpose();
      c1.pl(F1).pl(Fud).ccl(-1).gl(-5) . pr(Fud).pr(F1).gr(-5).ccr(1);
      c1 *= f;   
      Rcomplex con1 = c1.Trace();

      WilsonMatrix w1(q_f0_qpw.SiteMatrix(i,1));
      w1.ccl(-1).gl(-5) . gr(-5).ccr(1);
      WilsonMatrix w2(q_f1_qpw.SiteMatrix(i,0));
      w2.transpose();
      w1*=w2;
      Rcomplex con2 = w1.Trace();


      if( fabs(con1.real()-con2.real())>1e-08 || fabs(con1.imag()-con2.imag())>1e-08 ){
	printf("pi^+ correlator test failed on site %d:\n",i);
	printf("(%f,%f) : (%f,%f)\n",con1.real(),con1.imag(),con2.real(),con2.imag());
	exit(-1);
      }
    }
   
    {
      printf("Testing left and right projection\n");

      //Test left and right projection operators
      WilsonMatrix ww(q_f1_qpw.SiteMatrix(i,0));
      WilsonMatrix _1d2(0.0);
      for(int i=0;i<4;i++) for(int j=0;j<3;j++) _1d2(i,j,i,j) = Complex(1.0/2.0,0.0);

      WilsonMatrix _g5d2(_1d2); _g5d2.gr(-5);
      
      WilsonMatrix PL = _1d2-_g5d2;
      WilsonMatrix PR = _1d2+_g5d2;      

      printf("PL:\n");
      print(PL);
      printf("PR:\n");
      print(PR);

      WilsonMatrix PLww = PL*ww;
      WilsonMatrix PLww_f = ww; PLww_f.glPL();

      bool stop(false);

      if(!test_equals(PLww,PLww_f,1e-08)){
	printf("PLww test failed\n");
	print(PLww);
	print(PLww_f);
	stop=true;
      }


      WilsonMatrix wwPL = ww*PL;
      WilsonMatrix wwPL_f = ww; wwPL_f.grPL();
      if(!test_equals(wwPL,wwPL_f,1e-08)){
	printf("wwPL test failed\n");
	print(wwPL);
	print(wwPL_f);
	stop=true;
      }


      WilsonMatrix PRww = PR*ww;
      WilsonMatrix PRww_f = ww; PRww_f.glPR();
      if(!test_equals(PRww,PRww_f,1e-08)){
	printf("PRww test failed\n");
	print(PRww);
	print(PRww_f);
	stop=true;
      }


      WilsonMatrix wwPR = ww*PR;
      WilsonMatrix wwPR_f = ww; wwPR_f.grPR();
      if(!test_equals(wwPR,wwPR_f,1e-08)){
	printf("wwPR test failed\n");
	print(wwPR);
	print(wwPR_f);
	stop=true;
      }
      if(stop) exit(-1);


      printf("Passed\n");
	     
    }
 

  }

  {
    printf("Doing threaded correlation function test\n");

    //test threaded and unthreaded correlation function
    CorrelationFunction corrfunc_unthreaded("unthreaded");
    corrfunc_unthreaded.setNcontractions(1);

    for(int x=0;x<GJP.VolNodeSites();x++){
      int x_pos_vec[4];
      global_coord(x,x_pos_vec);

      WilsonMatrix A(q_f0_qpw.SiteMatrix(x,0));
      A.hconj();
      WilsonMatrix B(q_f0_qpw.SiteMatrix(x,0));
    
      Rcomplex result = Trace(A,B);
      corrfunc_unthreaded(0,x_pos_vec[3]) += result;
    }

    omp_set_num_threads(2);
    CorrelationFunction corrfunc_threaded("threaded",CorrelationFunction::THREADED);
    corrfunc_threaded.setNcontractions(1);

    printf("OMP num threads %d\n",omp_get_max_threads());
  
#pragma omp parallel for default(shared)
    for(int x=0;x<GJP.VolNodeSites();x++){
      printf("Thread %d, x %d start\n",omp_get_thread_num(),x); fflush(stdout);

      int x_pos_vec[4];
      global_coord(x,x_pos_vec);

      WilsonMatrix A(q_f0_qpw.SiteMatrix(x,0));
      A.hconj();
      WilsonMatrix B(q_f0_qpw.SiteMatrix(x,0));
    
      Rcomplex result = Trace(A,B);
      corrfunc_threaded(omp_get_thread_num(),0,x_pos_vec[3]) += result;
      printf("Thread %d, x %d end\n",omp_get_thread_num(),x); fflush(stdout);
    }
    //write does sumlattice which also sums over threads
    corrfunc_unthreaded.write("pooh1.dat");
    corrfunc_threaded.write("pooh2.dat");

    bool fail(false);
    
    for(int t=0;t<GJP.TnodeSites();t++){
      if( fabs( corrfunc_unthreaded(0,t).real() - corrfunc_threaded(0,t).real() )>1e-08 ||
	  fabs( corrfunc_unthreaded(0,t).imag() - corrfunc_threaded(0,t).imag() )>1e-08 ){
	printf("Fail at t=%d: (%f,%f) -- (%f,%f)\n",t,
	       corrfunc_unthreaded(0,t).real(),corrfunc_unthreaded(0,t).imag(),
	       corrfunc_threaded(0,t).real(),corrfunc_threaded(0,t).imag());
	fail = true;
      }
    }
    if(fail){
      printf("Failed\n");
      exit(-1);
    }

    printf("Passed\n");
  }


  //test P+A and P-A props
  {
    printf("Starting P+-A test\n");
    PropManager::clear();

    JobPropagatorArgs prop_args2;
    SETUP_ARRAY(prop_args2,props,PropagatorArg,4);
  
    char* names[4] = {"prop_P","prop_A","prop_F","prop_B"};
    BndCndType bndcnd[2] = {BND_CND_PRD,BND_CND_APRD};

    for(int i=0;i<2;i++){
      PropagatorArg &parg = prop_args2.props.props_val[i];
    
      parg.generics.tag = names[i];
      parg.generics.mass = 0.1;
      parg.generics.bc[0] = GJP.Xbc();
      parg.generics.bc[1] = GJP.Ybc();
      parg.generics.bc[2] = GJP.Zbc();
      parg.generics.bc[3] = bndcnd[i];

      SETUP_ARRAY(parg,attributes,AttributeContainer,3);
    
      ELEM(parg,attributes,0).type = POINT_SOURCE_ATTR;
      PointSourceAttrArg &srcarg = ELEM(parg,attributes,0).AttributeContainer_u.point_source_attr;
      for(int j=0;j<4;j++) srcarg.pos[j] = 0;

      ELEM(parg,attributes,1).type = GPARITY_FLAVOR_ATTR;
      GparityFlavorAttrArg &gparg = ELEM(parg,attributes,1).AttributeContainer_u.gparity_flavor_attr;
      gparg.flavor = 0;

      ELEM(parg,attributes,2).type = CG_ATTR;
      CGAttrArg &cgattr = ELEM(parg,attributes,2).AttributeContainer_u.cg_attr;
      cgattr.max_num_iter = 5000;
      cgattr.stop_rsd = 1e-08;
      cgattr.true_rsd = 1e-08;
    }
    PropCombination comb[2] = {A_PLUS_B,A_MINUS_B};

    for(int i=2;i<4;i++){
      PropagatorArg &parg = prop_args2.props.props_val[i];
    
      parg.generics.tag = names[i];
      parg.generics.mass = 0.1;
      parg.generics.bc[0] = GJP.Xbc();
      parg.generics.bc[1] = GJP.Ybc();

      SETUP_ARRAY(parg,attributes,AttributeContainer,1);

      ELEM(parg,attributes,0).type = PROP_COMBINATION_ATTR;
      PropCombinationAttrArg &carg = ELEM(parg,attributes,0).AttributeContainer_u.prop_combination_attr;
      carg.prop_A = names[0];
      carg.prop_B = names[1];
      carg.combination = comb[i-2];
    }

    if(UniqueID()==0) printf("prop_args contains %d propagators\n", prop_args2.props.props_len);

    PropManager::setup(prop_args2);
    
    PropagatorContainer &q_p = PropManager::getProp(names[0]);
    PropagatorContainer &q_a = PropManager::getProp(names[1]);
    PropagatorContainer &q_f = PropManager::getProp(names[2]);
    PropagatorContainer &q_b = PropManager::getProp(names[3]);
    
    QPropW & q_p_qpw = q_p.getProp(lattice);
    QPropW & q_a_qpw = q_a.getProp(lattice);
    QPropW & q_f_qpw = q_f.getProp(lattice);
    QPropW & q_b_qpw = q_b.getProp(lattice);

    q_f.printAttribs();
    q_b.printAttribs();

    bool fail(false);

    for(int x=0;x<GJP.VolNodeSites();x++){
      int x_pos_vec[4];
      global_coord(x,x_pos_vec);

      WilsonMatrix P(q_p_qpw.SiteMatrix(x,0));
      WilsonMatrix A(q_a_qpw.SiteMatrix(x,0));
      
      WilsonMatrix PpA = 0.5*(P+A);
      WilsonMatrix PmA = 0.5*(P-A);

      WilsonMatrix F(q_f_qpw.SiteMatrix(x,0));
      WilsonMatrix B(q_b_qpw.SiteMatrix(x,0));

      if(!test_equals(PpA,F,1e-08)){
	printf("P+A fail at x=%d: ",x);
	print(PpA);
	print(F);
	fail=true;
      }
      if(!test_equals(PmA,B,1e-08)){
	printf("P-A fail at x=%d: ",x);
	print(PmA);
	print(B);
	fail=true;
      }
    }
    if(fail){
      printf("P+-A test failed\n"); exit(-1);
    }else printf("P+-A Test Passed\n");

  }


  {
    //test pion correlation functions with cosine sources
     printf("Starting pion correlation functions with cosine sources test\n");
    PropManager::clear();

    JobPropagatorArgs prop_args2;
    SETUP_ARRAY(prop_args2,props,PropagatorArg,2);
  
    char* names[2] = {"prop+","prop-"};
    BndCndType bndcnd[2] = {BND_CND_APRD,BND_CND_APRD};
    int psign[2] = {1,-1};

    for(int i=0;i<2;i++){
      PropagatorArg &parg = prop_args2.props.props_val[i];
    
      parg.generics.tag = names[i];
      parg.generics.mass = 0.1;
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

    {
      printf("+ve momentum\n");

      int desired_mom_x[3] = {0,0,0};
      for(int i=0;i<3;i++) 
	if(GJP.Bc(i)==BND_CND_GPARITY) desired_mom_x[i] = 2;
	else desired_mom_x[i] = 0;

      CorrelationFunction pi_minus("pi^-");
      CorrelationFunction pi_plus("pi^+");

      calc_pi_minus(pi_minus, names[0], lattice, desired_mom_x);
      calc_pi_plus(pi_plus, names[0], lattice, desired_mom_x);

      if(!UniqueID()){
	for(int c=0;c<2;c++){
	  printf("Contraction %d\n",c);
	  for(int t=0;t<GJP.TnodeSites();t++){

	    Rcomplex &pm = pi_minus(c,t);
	    Rcomplex &pp = pi_plus(c,t);

	    printf("pi+:pi-  (%f,%f) : (%f,%f)\n",pp.real(),pp.imag(),pm.real(),pm.imag());
	  }
	}
      }
    }
    {
      printf("-ve momentum\n");

      int desired_mom_x[3] = {0,0,0};
      for(int i=0;i<3;i++) 
	if(GJP.Bc(i)==BND_CND_GPARITY) desired_mom_x[i] = -2;
	else desired_mom_x[i] = 0;

      CorrelationFunction pi_minus("pi^-");
      CorrelationFunction pi_plus("pi^+");

      calc_pi_minus(pi_minus, names[0], lattice, desired_mom_x);
      calc_pi_plus(pi_plus, names[0], lattice, desired_mom_x);

      if(!UniqueID()){
	for(int c=0;c<2;c++){
	  printf("Contraction %d\n",c);
	  for(int t=0;t<GJP.TnodeSites();t++){

	    Rcomplex &pm = pi_minus(c,t);
	    Rcomplex &pp = pi_plus(c,t);

	    printf("pi+:pi-  (%f,%f) : (%f,%f)\n",pp.real(),pp.imag(),pm.real(),pm.imag());
	  }
	}
      }
    }

  }

  {
    //test pion correlation functions with zero-momentum wall sources
    printf("Starting pion correlation functions with zero-momentum wall sources test\n");
    PropManager::clear();

    JobPropagatorArgs prop_args2;
    SETUP_ARRAY(prop_args2,props,PropagatorArg,1);
  
    char* names[1] = {"prop"};
    BndCndType bndcnd[1] = {BND_CND_APRD};

    for(int i=0;i<1;i++){
      PropagatorArg &parg = prop_args2.props.props_val[i];
    
      parg.generics.tag = names[i];
      parg.generics.mass = 0.1;
      parg.generics.bc[0] = GJP.Xbc();
      parg.generics.bc[1] = GJP.Ybc();
      parg.generics.bc[2] = GJP.Zbc();
      parg.generics.bc[3] = bndcnd[i];

      SETUP_ARRAY(parg,attributes,AttributeContainer,3);
    
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
    }
    if(UniqueID()==0) printf("prop_args contains %d propagators\n", prop_args2.props.props_len);

    PropManager::setup(prop_args2);   
    PropManager::calcProps(lattice);

    {
      printf("0 momentum\n");

      int desired_mom_x[3] = {0,0,0};

      CorrelationFunction pi_minus("pi^-");
      CorrelationFunction pi_plus("pi^+");

      calc_pi_minus(pi_minus, names[0], lattice, desired_mom_x);
      calc_pi_plus(pi_plus, names[0], lattice, desired_mom_x);

      if(!UniqueID()){
	for(int c=0;c<2;c++){
	  printf("Contraction %d\n",c);
	  for(int t=0;t<GJP.TnodeSites();t++){

	    Rcomplex &pm = pi_minus(c,t);
	    Rcomplex &pp = pi_plus(c,t);

	    printf("pi+:pi-  (%f,%f) : (%f,%f)\n",pp.real(),pp.imag(),pm.real(),pm.imag());
	  }
	}
      }
    }
  }

  {
    //test real and imaginary parts of K^0 propagator
    printf("Starting test real and imaginary parts of K^0 and \hat eta propagator with point source\n");
    PropManager::clear();

    JobPropagatorArgs prop_args2;
    SETUP_ARRAY(prop_args2,props,PropagatorArg,2);
  
    char* names[2] = {"prop_h","prop_l"};
    BndCndType bndcnd[2] = {BND_CND_APRD,BND_CND_APRD};
    double mass[2] = {0.4,0.1};

    for(int i=0;i<2;i++){
      PropagatorArg &parg = prop_args2.props.props_val[i];
    
      parg.generics.tag = names[i];
      parg.generics.mass = mass[i];
      parg.generics.bc[0] = GJP.Xbc();
      parg.generics.bc[1] = GJP.Ybc();
      parg.generics.bc[2] = GJP.Zbc();
      parg.generics.bc[3] = bndcnd[i];

      SETUP_ARRAY(parg,attributes,AttributeContainer,3);
    
      ELEM(parg,attributes,0).type = POINT_SOURCE_ATTR;
      PointSourceAttrArg &srcarg = ELEM(parg,attributes,0).AttributeContainer_u.point_source_attr;
      srcarg.pos[0] = 0;
      srcarg.pos[1] = 0;
      srcarg.pos[2] = 0;
      srcarg.pos[3] = 0;

      ELEM(parg,attributes,1).type = GPARITY_FLAVOR_ATTR;
      GparityFlavorAttrArg &gparg = ELEM(parg,attributes,1).AttributeContainer_u.gparity_flavor_attr;
      gparg.flavor = 0;

      ELEM(parg,attributes,2).type = CG_ATTR;
      CGAttrArg &cgattr = ELEM(parg,attributes,2).AttributeContainer_u.cg_attr;
      cgattr.max_num_iter = 5000;
      cgattr.stop_rsd = 1e-08;
      cgattr.true_rsd = 1e-08;
    }
    if(UniqueID()==0) printf("prop_args contains %d propagators\n", prop_args2.props.props_len);

    PropManager::setup(prop_args2);   
    PropManager::calcProps(lattice);

    {
      printf("0 momentum\n");

      int desired_mom_x[3] = {0,0,0};
      CorrelationFunction kzero("K^0");
      calc_k0(kzero, names[0], names[1], lattice, desired_mom_x);

      if(!UniqueID()){
	for(int t=0;t<GJP.TnodeSites();t++){
	  Rcomplex &pp = kzero(0,t);
	  printf("K^0[%d] (%f,%f)\n",t,pp.real(),pp.imag());
	}
      }

      CorrelationFunction etahat("\\hat eta");
      Rcomplex etahat_txdep[GJP.TnodeSites()*GJP.XnodeSites()];
      calc_etahat(etahat,names[1],lattice,desired_mom_x, &(etahat_txdep[0]));

      if(!UniqueID()){
	for(int t=0;t<GJP.TnodeSites();t++){
	  Rcomplex &pp = etahat(0,t);
	  printf("\\hat eta[%d] (%f,%f)\n",t,pp.real(),pp.imag());
	}
	printf("X dep:\n");
	for(int t=0;t<GJP.TnodeSites();t++){
	  printf("t=%d\n",t);
	  for(int x=0;x<GJP.XnodeSites();x++){
	    Rcomplex &pp = etahat_txdep[t+GJP.TnodeSites()*x];
	    printf("\\hat eta[x=%d] (%f,%f)\n",x,pp.real(),pp.imag());
	  }
	}

      }
    }

    {
      printf("\\sqrt{n}\pi/L momentum\n");

      int desired_mom_x[3] = {0,0,0};
      if(GJP.Xbc()==BND_CND_GPARITY) desired_mom_x[0] = 2; //units of pi/2L
      if(GJP.Ybc()==BND_CND_GPARITY) desired_mom_x[1] = 2;

      CorrelationFunction kzero("K^0");
      calc_k0(kzero, names[0], names[1], lattice, desired_mom_x);

      if(!UniqueID()){
	for(int t=0;t<GJP.TnodeSites();t++){
	  Rcomplex &pp = kzero(0,t);
	  printf("K^0[%d] (%f,%f)\n",t,pp.real(),pp.imag());
	}
      }

      Rcomplex etahat_txdep[GJP.TnodeSites()*GJP.XnodeSites()];
      CorrelationFunction etahat("\\hat eta");
      calc_etahat(etahat,names[1],lattice,desired_mom_x,etahat_txdep);
      if(!UniqueID()){
	for(int t=0;t<GJP.TnodeSites();t++){
	  Rcomplex &pp = etahat(0,t);
	  printf("\\hat eta[%d] (%f,%f)\n",t,pp.real(),pp.imag());
	}
	printf("X dep:\n");
	for(int t=0;t<GJP.TnodeSites();t++){
	  printf("t=%d\n",t);
	  for(int x=0;x<GJP.XnodeSites();x++){
	    Rcomplex &pp = etahat_txdep[t+GJP.TnodeSites()*x];
	    printf("\\hat eta[x=%d] (%f,%f)\n",x,pp.real(),pp.imag());
	  }
	}
	
      }
    }
   
  }


  {
    //test real and imaginary parts of K^0 propagator
    printf("Starting test real and imaginary parts of K^0 propagator with cosine source\n");
    PropManager::clear();

    JobPropagatorArgs prop_args2;
    SETUP_ARRAY(prop_args2,props,PropagatorArg,2);
  
    char* names[2] = {"prop_h","prop_l"};
    BndCndType bndcnd[2] = {BND_CND_APRD,BND_CND_APRD};
    double mass[2] = {0.4,0.1};
    
    for(int i=0;i<2;i++){
      PropagatorArg &parg = prop_args2.props.props_val[i];
    
      parg.generics.tag = names[i];
      parg.generics.mass = mass[i];
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

    {
      printf("0 momentum\n");

      int desired_mom_x[3] = {0,0,0};
      CorrelationFunction kzero("K^0");
      calc_k0(kzero, names[0], names[1], lattice, desired_mom_x);

      if(!UniqueID()){
	for(int t=0;t<GJP.TnodeSites();t++){
	  Rcomplex &pp = kzero(0,t);
	  printf("K^0[%d] (%f,%f)\n",t,pp.real(),pp.imag());
	}
      }
    }

  }


  {
    //test real and imaginary parts of <A_4 K^0> Green's function
    printf("Starting test real and imaginary parts of <A_4 K^0> Green's function with cosine source\n");
    PropManager::clear();

    JobPropagatorArgs prop_args2;
    SETUP_ARRAY(prop_args2,props,PropagatorArg,2);
  
    char* names[2] = {"prop_h","prop_l"};
    BndCndType bndcnd[2] = {BND_CND_APRD,BND_CND_APRD};
    double mass[2] = {0.4,0.1};
    
    for(int i=0;i<2;i++){
      PropagatorArg &parg = prop_args2.props.props_val[i];
    
      parg.generics.tag = names[i];
      parg.generics.mass = mass[i];
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

    {
      printf("0 momentum\n");

      int desired_mom_x[3] = {0,0,0};
      CorrelationFunction a4kzero("K^0");
      calc_a4k0(a4kzero, names[0], names[1], lattice, desired_mom_x);

      if(!UniqueID()){
	for(int t=0;t<GJP.TnodeSites();t++){
	  Rcomplex &pp = a4kzero(0,t);
	  printf("A_4 K^0[%d] (%f,%f)\n",t,pp.real(),pp.imag());
	}
      }
    }

  }


  if(gauge_fix) lattice.FixGaugeFree();

  if(UniqueID()==0){
    printf("Main job complete\n"); 
    fflush(stdout);
  }
  
  return 0;
}



