#include <config.h>
CPS_START_NAMESPACE
/*! \file
  \brief Definition of DiracOpWilsonTypes class constructor and destructor.
  
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson_types/d_op_wilson_types.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// d_op_wilson_types.C
//
// Is derived from DiracOp and is relevant to
// all DiracOp classes with Wilson type fermions 
// (e.g DiracOpWilson, DiracOpClover, DiracOpDwf, ...). 
//  These classes are derived from DiracOpWilsonTypes
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
DiracOpWilsonTypes::DiracOpWilsonTypes(Lattice & latt,
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
  cname = "DiracOpWilsonTypes";
  char *fname = "DiracOpWilsonTypes(L&,V*,V*,CgArg*,CnvFrmType)";
  VRB.Func(cname,fname);
}


//------------------------------------------------------------------
/*!
  If the storage order of any fields was changed by the constructor
  then they are changed to the canonical order by the destructor.
*/
//------------------------------------------------------------------
DiracOpWilsonTypes::~DiracOpWilsonTypes() {
  char *fname = "~DiracOpWilsonTypes()";
  VRB.Func(cname,fname);
}


//------------------------------------------------------------------
/*!
  A vector is left multiplied by  gamma_5
  \param out The resulting vector
  \param in The vector to be multiplied.
  \param gamma_num This must be 15.
  \param nodevol The size, in terms of floating point numbers, of the vectors.
*/
//------------------------------------------------------------------
void DiracOpWilsonTypes::MultGamma(Vector *out,
				   const Vector *in, 
				   int gamma_num, 
				   int nodevol){
  char *fname = "MultGamma(V*,V*,i,i)";
  VRB.Func(cname,fname);

  // For the moment, only support gamma_5 * in
  if (gamma_num != 15)
    ERR.General(cname,fname, "only gamma_5 supported");

  IFloat *p = (IFloat *)out;
  IFloat *q = (IFloat *)in;
  for(int n = 0; n < nodevol; ++n)
    {
    int i;
    for(i = 0; i < 12; ++i)
      *p++ = *q++;

    for(i = 0; i < 12; ++i)
      *p++ = - *q++;
    }
}

//------------------------------------------------------------------
// RitzLatSize returns the size of a fermion on a node
// It uses the RitzMatType flag to determine the operator
// to use and relevant checkerboard sizes.
//------------------------------------------------------------------
int DiracOpWilsonTypes::RitzLatSize() {
  char *fname = "RitzLatSize()";
  VRB.Func(cname,fname);

  int f_size = GJP.VolNodeSites() * lat.FsiteSize();

  switch(dirac_arg->RitzMatOper)
  {
  case MAT_HERM:
  case MATDAG_MAT:
  case NEG_MATDAG_MAT:
    break;

  case MATPC_HERM:
  case MATPCDAG_MATPC:
  case NEG_MATPCDAG_MATPC:
    f_size >>= 1;
    break;

  default:
    ERR.General(cname,fname,"RitzMatOper %d not implemented\n",
		dirac_arg->RitzMatOper);
  }

  return f_size;
}

//------------------------------------------------------------------
// RitzMat(Vector *out, Vector *in) :
// RitzMat is the base operator used in in Ritz.
// RitzMat works on the full or half lattice.
// The in, out fields are defined on the full or half lattice.
//------------------------------------------------------------------
void DiracOpWilsonTypes::RitzMat(Vector *out, Vector *in) {
  char *fname = "RitzMat(V*,V*)";
  VRB.Func(cname,fname);

  switch(dirac_arg->RitzMatOper)
    {
    case MAT_HERM:
    case MATDAG_MAT:
      MatDagMat(out, in);
      break;
      
    case MATPCDAG_MATPC:
      MatPcDagMatPc(out, in);
      break;

    case NEG_MATPCDAG_MATPC:
      MatPcDagMatPc(out, in);
      out -> VecNegative(out, RitzLatSize());
      break;    
      
    case NEG_MATDAG_MAT:
      MatDagMat(out, in);
      out->VecNegative(out, RitzLatSize());
      break;
      
    default:
      ERR.General(cname,fname,"RitzMatOper %d not implemented\n",
		  dirac_arg->RitzMatOper);
    }
}

//------------------------------------------------------------------
// RitzEigMat(Vector *out, Vector *in) :
// RitzEigMat is the base operator used in in RitzEig.
// RitzEigMat works on the full or half lattice.
// The in, out fields are defined on the full or half lattice.
//------------------------------------------------------------------
void DiracOpWilsonTypes::RitzEigMat(Vector *out, Vector *in) {
  char *fname = "RitzEigMat(V*,V*)";
  VRB.Func(cname,fname);

  switch(dirac_arg->RitzMatOper)
  {
  case MAT_HERM:
    MatHerm(out, in);
    break;

  case MATDAG_MAT:
    MatDagMat(out, in);
    break;

  case MATPCDAG_MATPC:
    MatPcDagMatPc(out, in);
    break;

  case NEG_MATPCDAG_MATPC:
    MatPcDagMatPc(out, in);
    out->VecNegative(out, RitzLatSize());
    break;

  case NEG_MATDAG_MAT:
    MatDagMat(out, in);
    out->VecNegative(out, RitzLatSize());
    break;

  default:
    ERR.General(cname,fname,"RitzMatOper %d not implemented",
		dirac_arg->RitzMatOper);
  }
}


//------------------------------------------------------------------
/*!
  Multiplies a vector by \f$ M^\dagger M \f$ where \e M is the fermion
  matrix with no preconditioning, so the vectors are defined on the whole
  lattice.

  \param out The resulting vector.
  \param in The vector to be multiplied.
*/
//-----------------------------------------------------------------
void DiracOpWilsonTypes::MatDagMat(Vector *out, Vector *in) {
  char *fname = "MatDagMat(V*,V*)";
  VRB.Func(cname,fname);

  int temp_size = RitzLatSize();
  Vector *temp = (Vector *) smalloc(temp_size * sizeof(Float));
  if (temp == 0) 
    ERR.Pointer(cname, fname, "temp");
  VRB.Smalloc(cname,fname, "temp", temp, temp_size);

  Mat(temp, in);
  MatDag(out, temp);
  
  VRB.Sfree(cname, fname, "temp", temp);
  sfree(temp);
}


CPS_END_NAMESPACE
