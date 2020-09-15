#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_mobius/noarch/mobius_dslash.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// mobius_dslash.C
//
// mobius_dslash is the derivative part of the fermion matrix. 
// The in, out fields are defined on the full lattice
// dag = 0/1 <--> Dslash/Dslash^dagger is calculated.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include<util/dwf.h>
#include<util/dirac_op.h>
#include<util/mobius.h>
#include<util/gjp.h>
#include<util/vector.h>
#include<util/verbose.h>
#include<util/error.h>


CPS_START_NAMESPACE
static void zTimesV1PluszTimesV2(Complex *a, Complex b, const Complex *c,
                    Complex d, const Complex *e, int len)
{
//#pragma omp parallel for
    for(int i = 0; i < len/2; ++i) {
            a[i] = b * c[i] + d*e[i];
                }
}

void mobius_unprec(Vector *out, 
		   Matrix *gauge_field, 
		   Vector *in, 
		   Float mass,
//		   int cb, 
		   int dag, 
		   Zmobus *mobius_lib_arg)
{
//------------------------------------------------------------------
// Apply 4-dimensional Dslash
//-b----------------------------------------------------------------- 

//  Vector  *frm_tmp2 = (Vector *) mobius_lib_arg->frm_tmp2;

  VRB.Debug("","mobius_unprec","mass=%e\n",mass);
  const unsigned long  f_size = (24/(2*6)) * mobius_lib_arg->vol_4d * mobius_lib_arg->ls;
  Vector  *frm_tmp2 = (Vector*)smalloc("","mobius_unpre()","frm_tmp2",2*f_size*(2*6)*sizeof(Float));
  memset(out,0,f_size*(2*6)*sizeof(Float));
  mobius_dslash_4(out+f_size, gauge_field, in, 0, dag, mobius_lib_arg, mass);
  mobius_dslash_4(out, gauge_field, in+f_size, 1, dag, mobius_lib_arg, mass);


//------------------------------------------------------------------
// Apply 5th-direction Dslash
//------------------------------------------------------------------
//  Complex* kappa_c = mobius_lib_arg->zmobius_kappa_c;
  Complex kappa_b = mobius_lib_arg->mobius_kappa_b;
//  Complex *one = new Complex[mobius_lib_arg->ls];
//  for(int i=0;i<mobius_lib_arg->ls;i++) one[i] = -1./kappa_c[i];
  Float kappa_c_inv = -1.0/mobius_lib_arg->mobius_kappa_c;

  memset(frm_tmp2,0,f_size*(2*6)*sizeof(Float));
  mobius_kappa_dslash_5_plus(frm_tmp2, in, mass, dag, mobius_lib_arg, kappa_c_inv);
  mobius_kappa_dslash_5_plus(frm_tmp2+f_size, in+f_size, mass, dag, mobius_lib_arg, kappa_c_inv);
  vecAddEquVec((IFloat*)out, (IFloat*) frm_tmp2, f_size*(2*6));
#if 1 
  // Multiply 2*kappa
  // do even / odd 
{
//  int local_ls = GJP.SnodeSites();
  int local_ls = mobius_lib_arg->ls;
  const int s_node_coor = GJP.SnodeCoor();
  const unsigned long  ls_stride = 24 * GJP.VolNodeSites()/2;
  const Complex half=-0.5;
//  unsigned long  size = GJP.VolNodeSites() * local_ls * 2 * Colors() * SpinComponents();
  unsigned long  size = GJP.VolNodeSites() * local_ls * 2 * 3 * 4;
  for(int ieo=0;ieo<2;++ieo){
//#pragma omp parallel for
    for(int s=0; s<local_ls;++s){
      int glb_s = s + local_ls*s_node_coor;
//      const Complex kappa_b =
//	2.0 / ( 2 * (GJP.ZMobius_b()[glb_s]
//		     *(4 - GJP.DwfHeight()) + GJP.DwfA5Inv()) );
      const Complex kb = 0.5 / kappa_b;
      int idx = s*ls_stride/2;// "/2" is for complex
      zTimesV1PluszTimesV2( (Complex*)out+idx+ieo*size/4,
			 kb, (Complex*)in+idx+ieo*size/4,
			 half, (Complex*)out+idx+ieo*size/4, ls_stride);
    }
  }
}
#endif

  sfree("","mobius_unpre()","frm_tmp2",frm_tmp2);


}

CPS_END_NAMESPACE
