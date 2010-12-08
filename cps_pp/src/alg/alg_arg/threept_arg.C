/* from Matthew Lightman. generated via rpcvmlgen, but not automatically
should be integrated into normal rpcgenvml chain in the future
1/13/2010 CJ
*/
#include <alg/threept_arg.h>
CPS_START_NAMESPACE

bool ThreePtArg::Encode(char *filename,char *instance){
  VML vmls;
  if ( !vmls.Create(filename,VML_ENCODE)) return false;
  if ( !Vml(&vmls,instance) ) return false;
  vmls.Destroy(); return true;
}

bool ThreePtArg::Decode(char *filename,char *instance){
  VML vmls;
  if ( !vmls.Create(filename,VML_DECODE)) return false;
  if ( !Vml(&vmls,instance)) return false;
  vmls.Destroy(); return true;
}
bool ThreePtArg::Vml(VML *vmls,char *instance){
  if(!vml_ThreePtArg(vmls,instance,this)) return false;
  return true;
}


bool_t
vml_ThreePtArg (VML *vmls, char *name,ThreePtArg *objp)
{
  register int32_t *buf;

  vml_class_begin(vmls,"ThreePtArg",name);
  if (!vml_string (vmls, "results", &objp->results, ~0))
    return FALSE;
  if (!vml_string (vmls, "results_mres_ZA", &objp->results_mres_ZA, ~0))
    return FALSE;
  if (!vml_string (vmls, "results_pipi", &objp->results_pipi, ~0))
    return FALSE;
  if (!vml_CgArg (vmls, "cg", &objp->cg))
    return FALSE;
  if (!vml_RandomType (vmls, "rng", &objp->rng))
    return FALSE;
  if (!vml_int (vmls, "gauge_fix", &objp->gauge_fix))
    return FALSE;
  if (!vml_int (vmls, "t_src", &objp->t_src))
    return FALSE;
  if (!vml_int (vmls, "t_snk", &objp->t_snk))
    return FALSE;
  if (!vml_int (vmls, "t_shift", &objp->t_shift))
    return FALSE;
  if (!vml_int (vmls, "box_len", &objp->box_len))
    return FALSE;
  if (!vml_int (vmls, "t_op", &objp->t_op))
    return FALSE;
  if (!vml_int (vmls, "width", &objp->width))
    return FALSE;
  if (!vml_int (vmls, "num_hits", &objp->num_hits))
    return FALSE;
  if (!vml_int (vmls, "do_susy", &objp->do_susy))
    return FALSE;
  if (!vml_int (vmls, "num_light", &objp->num_light))
    return FALSE;
  if (!vml_vector (vmls, "l_mass", (char *)objp->l_mass, 20,
		   sizeof (Float), (vmlproc_t) vml_Float))
    return FALSE;
  if (!vml_int (vmls, "num_strange", &objp->num_strange))
    return FALSE;
  if (!vml_vector (vmls, "s_mass", (char *)objp->s_mass, 20,
		   sizeof (Float), (vmlproc_t) vml_Float))
    return FALSE;
  if (!vml_int (vmls, "num_heavy", &objp->num_heavy))
    return FALSE;
  if (!vml_vector (vmls, "h_mass", (char *)objp->h_mass, 20,
		   sizeof (Float), (vmlproc_t) vml_Float))
    return FALSE;
  if (!vml_int (vmls, "num_tK", &objp->num_tK))
    return FALSE;
  if (!vml_vector (vmls, "tK", (char *)objp->tK, 20,
		   sizeof (int), (vmlproc_t) vml_int))
    return FALSE;
  if (!vml_int (vmls, "do_zero_mom", &objp->do_zero_mom))
    return FALSE;
  if (!vml_int (vmls, "do_first_mom", &objp->do_first_mom))
    return FALSE;
  if (!vml_int (vmls, "do_second_mom", &objp->do_second_mom))
    return FALSE;
  if (!vml_int (vmls, "do_third_mom", &objp->do_third_mom))
    return FALSE;
  if (!vml_int (vmls, "do_pipi_non_zero_tot_mom", &objp->do_pipi_non_zero_tot_mom))
    return FALSE;
  if (!vml_int (vmls, "do_p_plus_a_kaon", &objp->do_p_plus_a_kaon))
    return FALSE;
  if (!vml_int (vmls, "do_kaon_at_walls", &objp->do_kaon_at_walls))
    return FALSE;
  if (!vml_int (vmls, "do_kaons_tK", &objp->do_kaons_tK))
    return FALSE;
  if (!vml_int (vmls, "chkpoints", &objp->chkpoints))
    return FALSE;
  if (!vml_string (vmls, "ensemble_label", &objp->ensemble_label, ~0))
    return FALSE;
  if (!vml_string (vmls, "ensemble_id", &objp->ensemble_id, ~0))
    return FALSE;
  if (!vml_int (vmls, "seqNum", &objp->seqNum))
    return FALSE;
  vml_class_end(vmls,"ThreePtArg",name);
  return TRUE;
}
CPS_END_NAMESPACE
