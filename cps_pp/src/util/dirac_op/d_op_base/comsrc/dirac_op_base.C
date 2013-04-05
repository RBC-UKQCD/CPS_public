#include <config.h>

CPS_START_NAMESPACE
/*! \file
  \brief  Definition of DiracOp class methods.
  
  $Id: dirac_op_base.C,v 1.13 2013-04-05 17:46:30 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2013-04-05 17:46:30 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_base/comsrc/dirac_op_base.C,v 1.13 2013-04-05 17:46:30 chulwoo Exp $
//  $Id: dirac_op_base.C,v 1.13 2013-04-05 17:46:30 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.13 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_base/comsrc/dirac_op_base.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------


CPS_END_NAMESPACE
#include <util/dirac_op.h>
#include "asq_data_types.h"
#include <util/pt.h>
#include <util/lattice.h>
#include <util/verbose.h>
#include <util/error.h>
#include <util/gjp.h>
#include <comms/glb.h>
#include <comms/cbuf.h>

CPS_START_NAMESPACE

//------------------------------------------------------------------
// Initialize static variables
//------------------------------------------------------------------
int DiracOp::scope_lock = 0;
IntFlopCounter DiracOp::CGflops = 0;
unsigned int DiracOp::CGcount = 0;

//------------------------------------------------------------------
// static variables
//------------------------------------------------------------------
static int nx[4];
static int bc[4] = {0,0,0,0};	// boundary condition on this node


static void BondCond(Lattice& lat, Matrix *u_base)
{
      int x[4];
  for(int u = 0; u < 4; ++u) {
    if(bc[u]) {
      for(x[0] = 0; x[0] < nx[0]; ++x[0]) {
        for(x[1] = 0; x[1] < nx[1]; ++x[1]) {
          for(x[2] = 0; x[2] < nx[2]; ++x[2]) {
            for(x[3] = 0; x[3] < nx[3]; ++x[3]) {
	      if(x[u]==nx[u]-1) {
	        Matrix *m = u_base+lat.GsiteOffset(x)+u;
	        *m *= -1;
	      }
	    }
          }
	}
      }
    }
  }
}


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
  If this is ::CNV_FRM_YES, then the fields \a f_field_out and f_field_in are
  also converted: This assumes they are initially in the same order as
  the gauge field.
 */
//------------------------------------------------------------------
DiracOp::DiracOp(Lattice & latt,           // Lattice object
		 Vector *f_field_out,      // Output fermion field ptr.
		 Vector *f_field_in,       // Input fermion field ptr.
		 CgArg *arg,               // Argument structure
		 CnvFrmType convert) : // Fermion conversion flag
		 lat(latt), 
		 f_out(f_field_out),
		 f_in(f_field_in), 
		 dirac_arg(arg),
		 cnv_frm(convert)
{
  cname = "DiracOp";
  char *fname = "DiracOp(L&,V*,V*,CgArg,CnvFrmType)";
  VRB.Clock(cname,fname,"Just entered\n");
  VRB.Func(cname,fname);

  //----------------------------------------------------------------
  // Check and set the scope_lock
  //----------------------------------------------------------------
  if(scope_lock != 0){
    ERR.General(cname,fname,
		"Only one DiracOp object is allowed to be on scope\n");
  }
  scope_lock = 1;

  //----------------------------------------------------------------
  // Set the gauge field pointer
  //----------------------------------------------------------------
  gauge_field = lat.GaugeField();

  //----------------------------------------------------------------
  // Save the circular buffer control register values
  //----------------------------------------------------------------
  saveCbufCntrlReg();

  //----------------------------------------------------------------
  // Initialize boundary condition on this node
  //----------------------------------------------------------------
  bc[0] = 0;
  bc[1] = 0;
  bc[2] = 0;
  bc[3] = 0;
  if(GJP.Xbc() == BND_CND_APRD)
  	bc[0] = GJP.XnodeCoor() == (GJP.Xnodes()-1) ? 1 : 0;
  if(GJP.Ybc() == BND_CND_APRD)
  	bc[1] = GJP.YnodeCoor() == (GJP.Ynodes()-1) ? 1 : 0;
  if(GJP.Zbc() == BND_CND_APRD)
  	bc[2] = GJP.ZnodeCoor() == (GJP.Znodes()-1) ? 1 : 0;
  if(GJP.Tbc() == BND_CND_APRD)
  	bc[3] = GJP.TnodeCoor() == (GJP.Tnodes()-1) ? 1 : 0;

  nx[0] = GJP.XnodeSites();
  nx[1] = GJP.YnodeSites();
  nx[2] = GJP.ZnodeSites();
  nx[3] = GJP.TnodeSites();


  //----------------------------------------------------------------
  // turn on the boundary condition
  //----------------------------------------------------------------
  if (ParTrans::scope_lock == 0)
  BondCond(latt, gauge_field);

  //???

  // Added in by Ping for anisotropic lattices
  //------------------------------------------------------------------
  // No data manipulations if the scaling factor is 1.0
  // [built into lat.MltFloat].
  {
    Float factor = GJP.XiVXi() / GJP.XiV();
    lat.MltFloat(factor, GJP.XiDir());
  }

}


//------------------------------------------------------------------
/*!
  If the storage order of any fields was changed by the constructor
  then they are changed to the canonical order by the destructor.
*/
//------------------------------------------------------------------
DiracOp::~DiracOp() {
  char *fname = "~DiracOp()";

  VRB.Func(cname,fname);

  //----------------------------------------------------------------
  // turn off the boundary condition
  //----------------------------------------------------------------
  if (ParTrans::scope_lock == 0)
  BondCond(lat, gauge_field);

  VRB.Clock(cname,fname,"Exiting\n");

  //----------------------------------------------------------------
  // Restore the circular buffer control register values
  //----------------------------------------------------------------
  restoreCbufCntrlReg();


  //----------------------------------------------------------------
  // Release the scope_lock
  //----------------------------------------------------------------
  scope_lock = 0;

  //???

  // Added in by Ping for anisotropic lattices
  //------------------------------------------------------------------
  // No data manipulations if the scaling factor is 1.0
  // [built into lat.MltFloat].
  {
    Float factor = GJP.XiV() / GJP.XiVXi();
    lat.MltFloat(factor, GJP.XiDir());
  }
}

//------------------------------------------------------------------
// DiracOpGlbSum(Float *): 
// The global sum used by InvCg. At the base class
// level it is just the usual global sum. It is overloaded
// by any definitions by the derived classes. This is needed
// for example by DiracOpDwf where the global sum has different
// meaning for s_nodes = 1 and s_nodes > 1.
//------------------------------------------------------------------

/*!
  \param float_p Pointer to the floating point number.
  \post \a float_p points to the global sum.
*/
void DiracOp::DiracOpGlbSum(Float *float_p) {
  glb_sum(float_p);
}


#if 1
void DiracOp::RitzMat(Vector *out, Vector *in,
		     MatrixPolynomialArg* cheby_arg) {
  ERR.NotImplemented(cname,"RitzMat(V*,V*,MP*)");
}
#else 
//TIZB  PolynomialAccerelation
//
//  Q = [ -2 Ddag D + (alpha + beta) ] / [ alpha - beta ]
//
//  Output:  out =  T_n(Q) in
//
//  T_0 = 1,    T_1 = Q
//   T_{n+1}(Q) =  2 Q T_n(Q)  - T_{n-1}(Q)
//  
// Calling virtual RitzMat(V*,V*)
//
void DiracOp::RitzMat(Vector *out, Vector *in,
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
  const Float alpha =cheby_arg-> params[0];
  const Float beta =cheby_arg-> params[1];


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
#endif

CPS_END_NAMESPACE
