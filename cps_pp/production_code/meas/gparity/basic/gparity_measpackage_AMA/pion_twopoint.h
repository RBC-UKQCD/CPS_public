#ifndef _PION_TWOPOINT_H
#define _PION_TWOPOINT_H

#include "twopoint_function_generic.h"

CPS_START_NAMESPACE

//cf Eq 160 of GP paper. prop1 is the one that is daggered
//tsrc is used as the row index of the fmatrix. col index is (tsnk - tsrc + Lt) % Lt
//p1 is the source momentum of the first propagator (the one that will be daggered) and p2 is that of the second
//Option to use the wrong projection sign and wrong sink momentum (used for discussion in paper)
enum Pion2PtSinkOp { AX, AY, AZ, AT, P }; //Axial and pseudoscalar operator

inline SpinMatrixType snkOpMap(const Pion2PtSinkOp sink_op){
  SpinMatrixType snk_spn = spin_unit;
  switch(sink_op){
  case AX:
    snk_spn = gamma1; break;
  case AY:
    snk_spn = gamma2; break;
  case AZ:
    snk_spn = gamma3; break;
  case AT:
    snk_spn = gamma4; break;
  default:
    break;
  }
  return snk_spn;
}

//Pseudoscalar source and general sink
void pionTwoPointLWStandard(fMatrix<Rcomplex> &into, const int tsrc, const Pion2PtSinkOp sink_op, const ThreeMomentum &p_psibar, const ThreeMomentum &p_psi,
			    const PropSiteMatrixGetter &prop_dag, const PropSiteMatrixGetter &prop_undag){

  BasicOp src_op(spin_unit);
  BasicOp snk_op(snkOpMap(sink_op));

  Complex coeff(1.0);
  twoPointFunctionGeneric(into,tsrc,coeff,snk_op,src_op,p_psibar,p_psi,prop_dag,prop_undag);
}
void pionTwoPointLWGparity(fMatrix<Rcomplex> &into, const int tsrc, const Pion2PtSinkOp sink_op, const ThreeMomentum &p_psibar, const ThreeMomentum &p_psi,
			   const PropSiteMatrixGetter &prop_dag, const PropSiteMatrixGetter &prop_undag,
			   const PropSplane splane = SPLANE_BOUNDARY,
			   const bool use_wrong_proj_sign = false, const bool use_wrong_sink_mom = false){

  ThreeMomentum p_psi_src = p_psi;
  if(use_wrong_proj_sign) p_psi_src = -p_psi_src;

  GparityOpWithFlavorProject src_op(spin_unit,sigma3,p_psi_src);
  BasicGparityOp snk_op(snkOpMap(sink_op),sigma3);

  Complex coeff(0.5);
  twoPointFunctionGeneric(into,tsrc,coeff,snk_op,src_op,p_psibar,p_psi,prop_dag,prop_undag,splane,use_wrong_sink_mom);
}

//Time-component axial source and sink
void pionTwoPointA4A4LWGparity(fMatrix<Rcomplex> &into, const int tsrc, const ThreeMomentum &p_psibar, const ThreeMomentum &p_psi,
			   const PropSiteMatrixGetter &prop_dag, const PropSiteMatrixGetter &prop_undag){

  GparityOpWithFlavorProject src_op(gamma4,sigma3,p_psi);
  BasicGparityOp snk_op(gamma4,sigma3);

  Complex coeff(0.5);
  twoPointFunctionGeneric(into,tsrc,coeff,snk_op,src_op,p_psibar,p_psi,prop_dag,prop_undag);
}
void pionTwoPointA4A4LWStandard(fMatrix<Rcomplex> &into, const int tsrc, const ThreeMomentum &p_psibar, const ThreeMomentum &p_psi,
			   const PropSiteMatrixGetter &prop_dag, const PropSiteMatrixGetter &prop_undag){
  BasicOp src_op(gamma4);
  BasicOp snk_op(gamma4);

  Complex coeff(1.0);
  twoPointFunctionGeneric(into,tsrc,coeff,snk_op,src_op,p_psibar,p_psi,prop_dag,prop_undag);
}
 

//WW 2pt function
void pionTwoPointPPWWStandard(fMatrix<Rcomplex> &into, const int tsrc,
			     const WallSinkPropSiteMatrixGetter<WilsonMatrix> &prop_dag_W, const WallSinkPropSiteMatrixGetter<WilsonMatrix> &prop_undag_W){
  BasicOp src_op(spin_unit);
  BasicOp snk_op(spin_unit);

  Complex coeff(1.0);
  twoPointFunctionWallSinkGeneric(into, tsrc, coeff, snk_op, src_op, prop_dag_W, prop_undag_W);
}
//cf Eq 162 of GP paper. prop1 is the one that is daggered
void pionTwoPointPPWWGparity(fMatrix<Rcomplex> &into, const int tsrc, const ThreeMomentum &p_psi_snk, const ThreeMomentum &p_psi_src,
			     const WallSinkPropSiteMatrixGetter<SpinColorFlavorMatrix> &prop_dag_W, const WallSinkPropSiteMatrixGetter<SpinColorFlavorMatrix> &prop_undag_W){
  GparityOpWithFlavorProject src_op(spin_unit,sigma3,p_psi_src);
  GparityOpWithFlavorProject snk_op(spin_unit,sigma3,p_psi_snk);

  if(!UniqueID()) std::cout << "PPWW with source proj " << src_op.printProj() << " and sink proj " << snk_op.printProj() << '\n';

  Complex coeff(0.5);
  twoPointFunctionWallSinkGeneric(into, tsrc, coeff, snk_op, src_op, prop_dag_W, prop_undag_W);
}





//Pseudoscalar flavor singlet \bar\psi \gamma^5 \psi
//tsrc is used as the row index of the fmatrix. col index is (tsnk - tsrc + Lt) % Lt
void lightFlavorSingletLWGparity(fMatrix<Rcomplex> &into, const int tsrc, const ThreeMomentum &p_psibar, const ThreeMomentum &p_psi,
				 const PropSiteMatrixGetter &prop_dag, const PropSiteMatrixGetter &prop_undag){
  GparityOpWithFlavorProject src_op(spin_unit,sigma0,p_psi);
  BasicGparityOp snk_op(spin_unit,sigma0);

  Complex coeff(0.5,0);
  twoPointFunctionGeneric(into,tsrc,coeff,snk_op,src_op,p_psibar,p_psi,prop_dag,prop_undag);
}

//J5 or J5q
void J5Gparity(fMatrix<Rcomplex> &into, const int tsrc, const ThreeMomentum &p_psibar, const ThreeMomentum &p_psi,
	       const PropSiteMatrixGetter &prop_dag, const PropSiteMatrixGetter &prop_undag, const PropSplane splane = SPLANE_BOUNDARY, bool do_source_project = true){
  BasicGparityOp snk_op(spin_unit,sigma3);
  Complex coeff(0.5);

  if(do_source_project){
    GparityOpWithFlavorProject src_op(spin_unit,sigma3,p_psi);
    twoPointFunctionGeneric(into,tsrc,coeff,snk_op,src_op,p_psibar,p_psi,prop_dag,prop_undag,splane);
  }else{
    BasicGparityOp src_op(spin_unit,sigma3);
    twoPointFunctionGeneric(into,tsrc,coeff,snk_op,src_op,p_psibar,p_psi,prop_dag,prop_undag,splane);
  }
}

void J5Standard(fMatrix<Rcomplex> &into, const int tsrc, const ThreeMomentum &p_psibar, const ThreeMomentum &p_psi,
		const PropSiteMatrixGetter &prop_dag, const PropSiteMatrixGetter &prop_undag, const PropSplane splane = SPLANE_BOUNDARY){
  BasicOp src_op(spin_unit);
  BasicOp snk_op(spin_unit);

  Complex coeff(1.0);
  twoPointFunctionGeneric(into,tsrc,coeff,snk_op,src_op,p_psibar,p_psi,prop_dag,prop_undag,splane);
}






CPS_END_NAMESPACE
#endif
