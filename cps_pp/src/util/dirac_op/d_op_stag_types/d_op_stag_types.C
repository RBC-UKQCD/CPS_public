#include <config.h>
CPS_START_NAMESPACE
/*! \file
  \brief  Definition of DiracOpStagTypes class constructor and destructor.

  $Id: d_op_stag_types.C,v 1.6 2004/09/02 16:59:01 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004/09/02 16:59:01 $
//  $Header: /space/cvs/cps/cps++/src/util/dirac_op/d_op_stag_types/d_op_stag_types.C,v 1.6 2004/09/02 16:59:01 zs Exp $
//  $Id: d_op_stag_types.C,v 1.6 2004/09/02 16:59:01 zs Exp $
//  $Name: v5_0_8_eigCG_Qi $
//  $Locker:  $
//  $RCSfile: d_op_stag_types.C,v $
//  $Revision: 1.6 $
//  $Source: /space/cvs/cps/cps++/src/util/dirac_op/d_op_stag_types/d_op_stag_types.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// d_op_stag_types.C
//
// Is derived from DiracOp and is relevant to
// all DiracOp classes with Staggered type fermions 
// These classes are derived from DiracOpStagTypes
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/dirac_op.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/error.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
/*!
  Only one instance of this class is allowed to be in existence at
  any time.
  \param latt The lattice on which this Dirac operator is defined
  \param f_field_out A (pointer to) a spin-colour field (optionally). 
  \param f_field_in A (pointer to) a spin-colour field (optionally).
  \param arg Parameters for the solver.
  \param convert Whether the lattice fields should be converted to
  to a new storage order appropriate for the type of fermion action.
  If this is ::CNV_FRM_NO, then just the gauge field is converted.
  If this is ::CNV_FRM_YES, then the fields \a f_field_out and \a f_field_in
  are also converted: This assumes they are initially in the same order as
  the gauge field.
 */
//------------------------------------------------------------------
DiracOpStagTypes::DiracOpStagTypes(Lattice & latt,
				   Vector *f_field_out,
				   Vector *f_field_in,
				   CgArg *arg,
				   CnvFrmType convert) :
                                   DiracOp(latt, 
				           f_field_out,
				           f_field_in, 
				           arg,
				           convert)
{
  cname = "DiracOpStagTypes";
  char *fname = "DiracOpStagTypes(L&,V*,V*,CgArg*,CnvFrmType)";
  VRB.Func(cname,fname);
}


//------------------------------------------------------------------
/*!
  If the storage order of any fields was changed by the constructor
  then they are changed to the canonical order by the destructor.
*/
//------------------------------------------------------------------
DiracOpStagTypes::~DiracOpStagTypes() {
  char *fname = "~DiracOpStagTypes()";
  VRB.Func(cname,fname);
}

#if 0 
\\ already in eigen_stag.C
//------------------------------------------------------------------
// RitzLatSize returns the size of a fermion on a node
// It uses the RitzMatType flag to determine the operator
// to use and relevant checkerboard sizes.
//------------------------------------------------------------------
int DiracOpStagTypes::RitzLatSize() {
  char *fname = "RitzLatSize()";
  VRB.Func(cname,fname);

  size_t f_size = GJP.VolNodeSites() * lat.FsiteSize();

  switch(dirac_arg->RitzMatOper)
  {
  case MAT_HERM:
  case MATDAG_MAT:
  case NEG_MATDAG_MAT:
    break;

  case MATPC_HERM:
  case MATPCDAG_MATPC:
  case MATPCDAG_MATPC_SHIFT:
  case NEG_MATPCDAG_MATPC:
    f_size >>= 1;
    break;

  default:
    ERR.General(cname,fname,"RitzMatOper %d not implemented\n",
                dirac_arg->RitzMatOper);
  }

  //printf("DiracOpStagTypes::RitzLatSize(): f_size = %d\n", f_size);

  return f_size;
}
#endif


//TIZB  PolynomialAccerelation
//
//  Q = [ -2 Ddag D + (alpha + beta) ] / [ alpha - beta ]
//
//  Output:  out =  T_n(Q) in
//
//  T_0 = 1,    T_1 = Q
//   T_{n+1}(Q) =  2 Q T_n(Q)  - T_{n-1}(Q)
//  
void DiracOpStagTypes::RitzMat(Vector *out, Vector *in,
		     MatrixPolynomialArg* cheby_arg) {
  char *fname = "RitzMat(V*,V*,I,F,F,F,V*,V*)";
  VRB.Func(cname,fname);

  const int Npol = cheby_arg-> Npol;
  const int size = RitzLatSize();
  //
  // Q = 2 / (alpha-beta)  Ddag D  -   (alpha+beta)/(alpha-beta)
  // 2 Q =   c1  (  c0 Ddag D  -   1 )
  //  c1 = 2 (alpha+beta)/(alpha-beta),   c0 =  2 / (alpha+beta)
  //
  const Float alpha =cheby_arg-> params.params_val[0];
  const Float beta =cheby_arg-> params.params_val[1];


  const Float c1 =   2.0*(alpha+beta)/(alpha-beta);
  const Float c0 =   2.0/(alpha+beta);

  Vector *tmp  = (Vector*)cheby_arg->tmp1;
  Vector *tmp2 = (Vector*)cheby_arg->tmp2;
  
  //  tmp2 =  T_0 v = v = in
  tmp2 -> CopyVec(in, size);
  //  tmp =  T_1 v = Q v = Q in
  //  QV = 0.5* (2Q)V = 0.5 c1 ( c0 Ddag D - 1)
  RitzMat(tmp, in);
  tmp->VecTimesEquFloat(c0, size);
  tmp->VecMinusEquVec(in,size);
  tmp->VecTimesEquFloat(0.5*c1, size);

  // debug
  out->CopyVec(tmp,size);
  //printf("cheby %f %f\n", alpha,beta);
  
  // loop over
  for(int i=2; i<=Npol; ++i){
    // out = 2 Q tmp
    RitzMat(out, tmp);
    out->VecTimesEquFloat(c0, size);
    out->VecMinusEquVec(tmp,size);
    out->VecTimesEquFloat(c1, size);

    // out = out - tmp2
    out->VecMinusEquVec(tmp2, size);
    if( i!=Npol) {
      // tmp2 = tmp
      tmp2->CopyVec(tmp, size);
      // tmp = out
      tmp->CopyVec(out, size);
    }
  }
}

//------------------------------------------------------------------
// RitzMat(Vector *out, Vector *in) :
// RitzMat is the base operator used in in Ritz.
// RitzMat works on the full or half lattice.
// The in, out fields are defined on the full or half lattice.
//------------------------------------------------------------------
void DiracOpStagTypes::RitzMat(Vector *out, Vector *in) {
  char *fname = "RitzMat(V*,V*)";
  VRB.Func(cname,fname);
  Float *dot=0;

  Float mass = dirac_arg->mass;
  Float c = 1.0/(64.0 + 4.0*mass*mass);

  switch(dirac_arg->RitzMatOper)
    {
    case MATDAG_MAT:
      MatPcDagMatPc(out, in, dot);
      break;

    case MATPCDAG_MATPC:
      MatPcDagMatPc(out, in, dot);
      //dbug: RitzLatSize();
      break;

    case MATPCDAG_MATPC_SHIFT:
      MatPcDagMatPc(out, in, dot);
      break;

    case NEG_MATPCDAG_MATPC:
      MatPcDagMatPc(out, in, dot);
      out->VecNegative(out, RitzLatSize());
      break;

    case NEG_MATDAG_MAT:
      MatPcDagMatPc(out, in, dot);
      out->VecNegative(out, RitzLatSize());
      break;

    case MATDAG_MAT_NORM:
      MatPcDagMatPc(out, in, dot);
      out->VecTimesEquFloat(c,RitzLatSize());
      break;

    case NEG_MATDAG_MAT_NORM:
      MatPcDagMatPc(out, in, dot);
      out->VecTimesEquFloat(-c,RitzLatSize());
      break;

    default:
      ERR.General(cname,fname,"RitzMatOper %d not implemented",
                  dirac_arg->RitzMatOper);
    }
}
CPS_END_NAMESPACE
