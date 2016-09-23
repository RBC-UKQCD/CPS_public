//4d precond. mobius Dirac op of symmetric version 2nd kind but mult kappa_b in dslash_4:
//  (1 -  kappa_b M4eo M_5^-1 kappa_b M4oe M_5^-1)^dag
//= 1 - M_5^dag-1 M4oe^dag kappa_b^dag M_5^dag-1 M4oe^dag kappa_b^dag  
void  zmobius_mdag_sym3(Vector *out, 
		   Matrix *gauge_field, 
		   Vector *in, 
		   Float mass, 
		   Zmobus *mobius_lib_arg)
{
  const int dag=1;
  //------------------------------------------------------------------
  // Initializations
  //------------------------------------------------------------------
  const int f_size = 24 * mobius_lib_arg->vol_4d * mobius_lib_arg->ls / 2;
  const int ls=mobius_lib_arg->ls;
  const int vol_4d_cb = mobius_lib_arg->vol_4d / 2;
  const int ls_stride = 24 * vol_4d_cb;
  const int local_ls = GJP.SnodeSites();
  const int s_nodes = GJP.Snodes();
  const int global_ls = local_ls * s_nodes;
  const int s_node_coor = GJP.SnodeCoor();

  

  Vector  *frm_tmp2 = (Vector *) mobius_lib_arg->frm_tmp2;
  //Vector *temp = (Vector *) smalloc(f_size * sizeof(Float));
  Float norm;

  
  //  out = [ 1 -  M5inv^dag Moo^dag M5inv^dag Moe^dag kappa_b^dag ] in
  // (dslash_5 uses (1+-g5), not P_R,L, i.e. no factor of 1/2 which is here out front)
  //
  //    1. out = ( kappa_b Meo M5inv kappa_b Moe M5inv)^dag in
  //    2. -out +=  in
  

  //--------------------------------------------------------------
  //    1. ftmp2 = (kappa_b Meo M5inv kappa_b Moe M5inv)^dag in
  //--------------------------------------------------------------
  
  time_elapse();

  // out<- in
  DEBUG_MOBIUS_DSLASH("mdag start\n", time_elapse());
  moveFloat((IFloat*)out, (IFloat*)in, f_size);
  DEBUG_MOBIUS_DSLASH("out<-in %e\n", time_elapse());


  Complex *b_coeff = GJP.ZMobius_b();
  Complex *c_coeff = GJP.ZMobius_c();
  const Complex *kappa_b = mobius_lib_arg->zmobius_kappa_b;

  Complex Imag(0,1);
  DEBUG_MOBIUS_DSLASH("before new b,c coeff %e\n", time_elapse());
  for(int s=0;s<ls;++s){
    b_coeff[s] *= kappa_b[s];
    b_coeff[s] *= Imag;
    c_coeff[s] *= kappa_b[s];
    c_coeff[s] *= Imag;
  }
  DEBUG_MOBIUS_DSLASH("after new b,c coeff %e\n", time_elapse());
 
  // Apply Dslash^dag O <- E
  //------------------------------------------------------------------
  zmobius_dslash_4(frm_tmp2, gauge_field, out, 0, dag, mobius_lib_arg, mass);
  DEBUG_MOBIUS_DSLASH("mobius_dslash_4 dag 1 %e\n", time_elapse());

  //------------------------------------------------------------------
  // Apply [M_5^-1]^dag (hopping in 5th dir + diagonal)
  //------------------------------------------------------------------
  zmobius_m5inv(frm_tmp2, mass, dag, mobius_lib_arg,mobius_lib_arg->zmobius_kappa_ratio);
  DEBUG_MOBIUS_DSLASH("mobius_m5inv %e\n", time_elapse());
 
  //------------------------------------------------------------------
  // Apply Dslash E <- O dag
  //------------------------------------------------------------------
  zmobius_dslash_4(out, gauge_field, frm_tmp2, 1, dag, mobius_lib_arg, mass);
  DEBUG_MOBIUS_DSLASH("mobius_dslash_4  dag 1 %e\n", time_elapse());
  
  //------------------------------------------------------------------
  // Apply [M_5^-1]^dag (hopping in 5th dir + diagonal)
  //------------------------------------------------------------------
  zmobius_m5inv(out, mass, dag, mobius_lib_arg,mobius_lib_arg->zmobius_kappa_ratio);
  DEBUG_MOBIUS_DSLASH("mobius_m5inv %e\n", time_elapse());


  
  //------------------------------------------------------------------
  //    2. 
  //------------------------------------------------------------------
//#ifndef USE_BLAS
#if 1
  vecAddEquVec((IFloat*)out, (IFloat*)in, f_size);
  DEBUG_MOBIUS_DSLASH("out+=in %e\n", time_elapse());
#else
 !  SENTINEL !  USE_BLAS is not supported YET for complexified mobius
  cblas_dcopy(f_size, (IFloat*)in, 1, (IFloat*)out, 1);
#endif

  // Restore b,c coeffs
  for(int s=0;s<ls;++s){
    b_coeff[s] *= -Imag;
    b_coeff[s] /= kappa_b[s];
    c_coeff[s] *= -Imag;
    c_coeff[s] /= kappa_b[s];
  }
  DEBUG_MOBIUS_DSLASH("restore b,c coeffs %e\n", time_elapse());
  

}

