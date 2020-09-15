#ifndef _KAON_TWOPOINT_H
#define _KAON_TWOPOINT_H

#include "twopoint_function_generic.h"

CPS_START_NAMESPACE


//Kaon two-point. prop_dag_h is the heavy quark propagator (the one to which g5-hermiticity is applied), and prop_undag_l the light-quark prop
//tsrc is used as the row index of the fmatrix. col index is (tsnk - tsrc + Lt) % Lt
void kaonTwoPointPPLWGparity(fMatrix<Rcomplex> &into, const int tsrc, const ThreeMomentum &p_psibar_l, const ThreeMomentum &p_psi_h,
			      const PropSiteMatrixGetter &prop_dag_h, const PropSiteMatrixGetter &prop_undag_l){
  GparityOpWithFlavorProject src_op(spin_unit,sigma0,p_psi_h);
  BasicGparityOp snk_op(spin_unit,sigma0);

  Complex coeff(0.5,0); //note positive sign because we define the creation/annihilation operator with a factor of i
  twoPointFunctionGeneric(into,tsrc,coeff,snk_op,src_op,p_psibar_l,p_psi_h,prop_dag_h,prop_undag_l);
}
//Kaon two-point with no Gparity
void kaonTwoPointPPLWStandard(fMatrix<Rcomplex> &into, const int tsrc, const ThreeMomentum &p_psibar_l, const ThreeMomentum &p_psi_h,
			      const PropSiteMatrixGetter &prop_dag_h, const PropSiteMatrixGetter &prop_undag_l){
  BasicOp src_op(spin_unit);
  BasicOp snk_op(spin_unit);

  Complex coeff(1.0,0); //note positive sign because we define the creation/annihilation operator with a factor of i
  twoPointFunctionGeneric(into,tsrc,coeff,snk_op,src_op,p_psibar_l,p_psi_h,prop_dag_h,prop_undag_l);
}


//*Physical* time-component axial operator sink   -i F0 g4 g5    (g5 is removed by g5-hermiticity)
void kaonTwoPointA4PhysPLWGparity(fMatrix<Rcomplex> &into, const int tsrc, const ThreeMomentum &p_psibar_l, const ThreeMomentum &p_psi_h,
			      const PropSiteMatrixGetter &prop_dag_h, const PropSiteMatrixGetter &prop_undag_l){
  GparityOpWithFlavorProject src_op(spin_unit,sigma0,p_psi_h);
  BasicGparityOp snk_op(gamma4,F0);

  Complex coeff(1.0,0);
  twoPointFunctionGeneric(into,tsrc,coeff,snk_op,src_op,p_psibar_l,p_psi_h,prop_dag_h,prop_undag_l);
}
//Physical time-component axial operator sink with no Gparity
void kaonTwoPointA4PhysPLWStandard(fMatrix<Rcomplex> &into, const int tsrc, const ThreeMomentum &p_psibar_l, const ThreeMomentum &p_psi_h,
			      const PropSiteMatrixGetter &prop_dag_h, const PropSiteMatrixGetter &prop_undag_l){
  BasicOp src_op(spin_unit);
  BasicOp snk_op(gamma4);

  Complex coeff(1.0,0);
  twoPointFunctionGeneric(into,tsrc,coeff,snk_op,src_op,p_psibar_l,p_psi_h,prop_dag_h,prop_undag_l);
}


//*Unphysical* time-component axial operator sink   -i F1 g4 g5, connects to unphysical kaon component    (g5 is removed by g5-hermiticity)
void kaonTwoPointA4UnphysPLWGparity(fMatrix<Rcomplex> &into, const int tsrc, const ThreeMomentum &p_psibar_l, const ThreeMomentum &p_psi_h,
			      const PropSiteMatrixGetter &prop_dag_h, const PropSiteMatrixGetter &prop_undag_l){
  GparityOpWithFlavorProject src_op(spin_unit,sigma0,p_psi_h);
  BasicGparityOp snk_op(gamma4,F1);

  Complex coeff(1.0,0);
  twoPointFunctionGeneric(into,tsrc,coeff,snk_op,src_op,p_psibar_l,p_psi_h,prop_dag_h,prop_undag_l);
}
//Time-component axial source and sink that connects to both the physical and unphysical components
void kaonTwoPointA4combA4combLWGparity(fMatrix<Rcomplex> &into, const int tsrc, const ThreeMomentum &p_psibar_l, const ThreeMomentum &p_psi_h,
			      const PropSiteMatrixGetter &prop_dag_h, const PropSiteMatrixGetter &prop_undag_l){
  GparityOpWithFlavorProject src_op(gamma4,sigma0,p_psi_h);
  BasicGparityOp snk_op(gamma4,sigma0);

  Complex coeff(0.5,0);
  twoPointFunctionGeneric(into,tsrc,coeff,snk_op,src_op,p_psibar_l,p_psi_h,prop_dag_h,prop_undag_l);
}

//Time-component axial source and sink with no Gparity
void kaonTwoPointA4PhysA4PhysLWStandard(fMatrix<Rcomplex> &into, const int tsrc, const ThreeMomentum &p_psibar_l, const ThreeMomentum &p_psi_h,
					const PropSiteMatrixGetter &prop_dag_h, const PropSiteMatrixGetter &prop_undag_l){
  BasicOp src_op(gamma4);
  BasicOp snk_op(gamma4);

  Complex coeff(1.0,0);
  twoPointFunctionGeneric(into,tsrc,coeff,snk_op,src_op,p_psibar_l,p_psi_h,prop_dag_h,prop_undag_l);
}


void kaonTwoPointPPWWGparity(fMatrix<Rcomplex> &into, const int tsrc, const ThreeMomentum &p_psi_l_snk, const ThreeMomentum &p_psi_h_src,
			     const WallSinkPropSiteMatrixGetter<SpinColorFlavorMatrix> &prop_dag_h_W, const WallSinkPropSiteMatrixGetter<SpinColorFlavorMatrix> &prop_undag_l_W){
  GparityOpWithFlavorProject src_op(spin_unit,sigma0,p_psi_h_src);
  GparityOpWithFlavorProject snk_op(spin_unit,sigma0,p_psi_l_snk);

  if(!UniqueID()) std::cout << "Kaon PPWW with source proj " << src_op.printProj() << " and sink proj " << snk_op.printProj() << '\n';

  Complex coeff(0.5,0);
  twoPointFunctionWallSinkGeneric(into, tsrc, coeff, snk_op, src_op, prop_dag_h_W, prop_undag_l_W);
}

void kaonTwoPointPPWWStandard(fMatrix<Rcomplex> &into, const int tsrc,
			     const WallSinkPropSiteMatrixGetter<WilsonMatrix> &prop_dag_h_W, const WallSinkPropSiteMatrixGetter<WilsonMatrix> &prop_undag_l_W){
  BasicOp src_op(spin_unit);
  BasicOp snk_op(spin_unit);

  Complex coeff(1.0,0);
  twoPointFunctionWallSinkGeneric(into, tsrc, coeff, snk_op, src_op, prop_dag_h_W, prop_undag_l_W);
}


CPS_END_NAMESPACE
#endif
