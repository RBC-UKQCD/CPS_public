#include <config.h>
#include <stdio.h>
#include <util/dirac_op.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/error.h>
#include <util/wilson.h>
#include <util/time_cps.h>
#include <util/dwf.h>
//#include <mem/p2v.h>
#include <comms/glb.h>

CPS_START_NAMESPACE
int DiracOpDwf::eig_MatInv(Vector **V, const int vec_len, Float *M, const int nev, const int m, float **U, Rcomplex *invH, const int def_len, const Float *restart, const int restart_len,Vector *out, Vector *in, Float *true_res, PreserveType prs_in) 
{
  char *fname = "eig_MatInv(V*,V*,F*)";
  VRB.Func(cname,fname);

  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fdwf
  // and DiracOpDwf constructors but is done again in case the
  // user has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpDwf object.
  //----------------------------------------------------------------
  Dwf *dwf_arg = (Dwf *) dwf_lib_arg;
  dwf_arg->ls = GJP.SnodeSites();
  dwf_arg->dwf_kappa = 
    1.0 / ( 2 * (4 + GJP.DwfA5Inv() - GJP.DwfHeight()) );

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  Vector *temp2;
  int temp_size = GJP.VolNodeSites() * lat.FsiteSize() / 2;

//  printf("temp_size:%d\n",temp_size);
//  printf("MatInv : %e %e\n",in->NormSqNode(temp_size),out->NormSqNode(temp_size));



  // check out if converted
  //for (int ii = 0; ii < 2 * temp_size; ii++) {
  //  VRB.Result(cname, fname, "in[%d] = %e\n", ii, 
  //  *((IFloat *)in + ii));
  //  VRB.Result(cname, fname, "out[%d] = %e\n", ii, 
  //  *((IFloat *)out + ii));
  //}

  Vector *temp = (Vector *) smalloc(temp_size * sizeof(Float));
  if (temp == 0) ERR.Pointer(cname, fname, "temp");
  VRB.Smalloc(cname,fname, "temp", temp, temp_size * sizeof(Float));

  if(prs_in == PRESERVE_YES){
    temp2 = (Vector *) smalloc(temp_size * sizeof(Float));
    if (temp2 == 0) ERR.Pointer(cname, fname, "temp2");
    VRB.Smalloc(cname,fname, "temp2", temp2, temp_size * sizeof(Float));
  }

  // points to the even part of fermion source 
  Vector *even_in = (Vector *) ( (IFloat *) in + temp_size );

  // points to the even part of fermion solution
  Vector *even_out = (Vector *) ( (IFloat *) out + temp_size );

//  printf("MatInv : even : %e %e\n",even_in->NormSqNode(temp_size),even_out->NormSqNode(temp_size));

	
  Dslash(temp, even_in, CHKB_EVEN, DAG_NO);
//  printf("MatInv : even : Dslash : temp:%e even:%e\n",temp->NormSqNode(temp_size),even_in->NormSqNode(temp_size));

  fTimesV1PlusV2((IFloat *)temp, (IFloat) dwf_arg->dwf_kappa, (IFloat *)temp,
    (IFloat *)in, temp_size);

//  printf("MatInv : even : Dslash : temp:%e \n",temp->NormSqNode(temp_size));

  // save source
  if(prs_in == PRESERVE_YES){
    moveMem((IFloat *)temp2, (IFloat *)in, 
		temp_size * sizeof(IFloat) / sizeof(char));
  }


  int iter;
  switch (dirac_arg->Inverter) {
  case EIGCG:
	  MatPcDag(in,temp);
	  iter = InvEigCg(out,in,true_res,nev,m,V,vec_len,M,U,invH,def_len,restart,restart_len);
	  break;
//  case CG:
//    MatPcDag(in, temp);
//    iter = InvCg(out,in,true_res);
//    break;
//  case BICGSTAB:
//    iter = BiCGstab(out,temp,0.0,dirac_arg->bicgstab_n,true_res);
//    break;
  default:
    ERR.General(cname,fname,"InverterType %d not implemented in EigCG Code\n",
		dirac_arg->Inverter);
  }

  // restore source
  if(prs_in == PRESERVE_YES){
    moveMem((IFloat *)in, (IFloat *)temp2, 
		temp_size * sizeof(IFloat) / sizeof(char));
  }

  Dslash(temp, out, CHKB_ODD, DAG_NO);

  fTimesV1PlusV2((IFloat *)even_out, (IFloat) dwf_arg->dwf_kappa, (IFloat *)temp,
    (IFloat *) even_in, temp_size);

  VRB.Sfree(cname, fname, "temp", temp);
  sfree(temp);

  if(prs_in == PRESERVE_YES){
    VRB.Sfree(cname, fname, "temp2", temp2);
    sfree(temp2);
  }

  return iter;
}
CPS_END_NAMESPACE
