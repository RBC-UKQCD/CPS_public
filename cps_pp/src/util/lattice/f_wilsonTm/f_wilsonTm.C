#include<config.h>

CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of FwilsonTm class.

  $Id: f_wilsonTm.C,v 1.4 2013-03-18 19:33:14 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/f_wilsonTm/f_wilsonTm.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// f_wilsonTm.C
//
// FwilsonTm is derived from Fwilson and is relevant to
// twisted-mass wilson fermions
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/dirac_op.h>
#include <util/wilson.h>
#include <util/verbose.h>
#include <util/gjp.h>
#include <util/error.h>
#include <comms/scu.h>
#include <comms/glb.h>
#include <util/lattice/fforce_wilson_type.h>

CPS_START_NAMESPACE

//------------------------------------------------------------------
// This constructor does nothing.
// All initialization done by Fwilson constructor.
//------------------------------------------------------------------
FwilsonTm::FwilsonTm()
: Fwilson()
{
  cname = "FwilsonTm";
  char *fname = "FwilsonTm()";
  VRB.Func(cname,fname);
}

//------------------------------------------------------------------
// This destructor does nothing.
// All termination done by Fwilson destructor.
//------------------------------------------------------------------
FwilsonTm::~FwilsonTm()
{
  char *fname = "~FwilsonTm()";
  VRB.Func(cname,fname);
}

//------------------------------------------------------------------
// FclassType Fclass():
// It returns the type of fermion class.
//------------------------------------------------------------------
FclassType FwilsonTm::Fclass() const{
  return F_CLASS_WILSON_TM;
}

//------------------------------------------------------------------
// int FsiteSize() and int FchkbEvl() not changed for twisted-mass Wilson fermions
//------------------------------------------------------------------

//------------------------------------------------------------------
//~~ modified in f_wilsonTm to create wilsonTm fermions 
//~~ see full notes in f_wilson
//------------------------------------------------------------------
int FwilsonTm::FmatEvlInv(Vector *f_out, Vector *f_in, 
                          CgArg *cg_arg, 
                          Float *true_res,
                          CnvFrmType cnv_frm)
{
    char *fname = "FmatEvlInv(CgArg*,V*,V*,F*,CnvFrmType)";
    VRB.Func(cname,fname);

    DiracOpWilsonTm wilson(*this, f_out, f_in, cg_arg, cnv_frm);
    int iter = wilson.InvCg(&(cg_arg->true_rsd));
    if (true_res) *true_res = cg_arg ->true_rsd;

    return iter;
}

//------------------------------------------------------------------
// int FmatEvlMInv not changed for twisted-mass Wilson fermions
//------------------------------------------------------------------

//------------------------------------------------------------------
// void FminResExt not changed for twisted-mass Wilson fermions
//------------------------------------------------------------------

//------------------------------------------------------------------
// int FmatInv not changed for twisted-mass Wilson fermions
// NOTE: this call MatInv
//------------------------------------------------------------------

//------------------------------------------------------------------
// int FeigSolv not changed for twisted-mass Wilson fermions
//------------------------------------------------------------------

//------------------------------------------------------------------
// SetPhi(Vector *phi, Vector *frm1, Vector *frm2, Float mass,
//        Float epsilon, DagType dag):
// It sets the pseudofermion field phi from frm1, frm2.
// Note that frm2 is not used.
// Modified - now returns the (trivial) value of the action
// Now sets epsilon in cg_arg from new input parameter
//------------------------------------------------------------------
Float FwilsonTm::SetPhi(Vector *phi, Vector *frm1, Vector *frm2,
		      Float mass, Float epsilon, DagType dag){
  char *fname = "SetPhi(V*,V*,V*,F,F)";
  VRB.Func(cname,fname);

  CgArg cg_arg;
  cg_arg.mass = mass;
  cg_arg.epsilon = epsilon;

  if (phi == 0)
    ERR.Pointer(cname,fname,"phi") ;

  if (frm1 == 0)
    ERR.Pointer(cname,fname,"frm1") ;

  DiracOpWilsonTm wilson(*this, frm1, frm2, &cg_arg, CNV_FRM_NO) ;
  
  if (dag == DAG_YES) wilson.MatPcDag(phi, frm1) ;
  else wilson.MatPc(phi, frm1) ;

  return FhamiltonNode(frm1, frm1);
}

//------------------------------------------------------------------
// SetPhi(Vector *phi, Vector *frm1, Vector *frm2, Float mass, DagType dag):
// should never be called by wilsonTm fermions
//------------------------------------------------------------------
Float FwilsonTm::SetPhi(Vector *phi, Vector *frm1, Vector *frm2,
		      Float mass, DagType dag){
  char *fname = "SetPhi(V*,V*,V*,F,DagType)";
  VRB.Func(cname,fname);

  ERR.General(cname,fname,"f_wilsonTm: SetPhi(V*,V*,V*,F,DagType) not implemented here\n");

  return Float(0.0);
}

//------------------------------------------------------------------
// ForceArg RHMC_EvolveMomFforce not changed for twisted-mass Wilson fermions
//------------------------------------------------------------------

//------------------------------------------------------------------
// Float BhamiltonNode(Vector *boson, Float mass, Float epsilon):
// The boson Hamiltonian of the node sublattice.
// Now sets epsilon in cg_arg from new input parameter
//------------------------------------------------------------------
Float FwilsonTm::BhamiltonNode(Vector *boson, Float mass, Float epsilon){
  char *fname = "BhamiltonNode(V*,F,F)";
  VRB.Func(cname,fname);
  
  CgArg cg_arg;
  cg_arg.mass = mass;
  cg_arg.epsilon = epsilon;

  if (boson == 0)
    ERR.Pointer(cname,fname,"boson");

  int f_size = (GJP.VolNodeSites() * FsiteSize()) >> 1 ;

  Vector *bsn_tmp = (Vector *)
    smalloc(f_size*sizeof(Float));

  char *str_tmp = "bsn_tmp" ;

  if (bsn_tmp == 0)
    ERR.Pointer(cname,fname,str_tmp) ;

  VRB.Smalloc(cname,fname,str_tmp,bsn_tmp,f_size*sizeof(Float));

  DiracOpWilsonTm wilson(*this, boson, bsn_tmp, &cg_arg, CNV_FRM_NO) ;

  wilson.MatPc(bsn_tmp,boson);

  Float ret_val = bsn_tmp->NormSqNode(f_size) ;

  VRB.Sfree(cname,fname,str_tmp,bsn_tmp);

  sfree(bsn_tmp) ;

  return ret_val;
}

//------------------------------------------------------------------
// Float BhamiltonNode(Vector *boson, Float mass, Float epsilon):
// should never be called by wilsonTm fermions
//------------------------------------------------------------------
Float FwilsonTm::BhamiltonNode(Vector *boson, Float mass){
  char *fname = "BhamiltonNode(V*,F)";
  VRB.Func(cname,fname);
  
  ERR.General(cname,fname,"f_wilsonTm: BhamiltonNode(V*,F) not implemented here\n");

  return Float(0.0);
}

//------------------------------------------------------------------
// EvolveMomFforce(Matrix *mom, Vector *frm, Float mass, 
//                 Float epsilon, Float dt):
// It evolves the canonical momentum mom by dt
// using the fermion force.
// Now sets epsilon in cg_arg from new input parameter
// chi <- (M_f^\dag M_f)^{-1} M_f^\dag (RGV)
//------------------------------------------------------------------
ForceArg FwilsonTm::EvolveMomFforce(Matrix *mom, Vector *chi, 
                                    Float mass, Float epsilon, Float dt)
{
    const char *fname = "EvolveMomFforce(M*,V*,F,F,F)";
    if (mom == 0) ERR.Pointer(cname,fname,"mom");
    if (chi == 0) ERR.Pointer(cname,fname,"chi");

    const int f_size = FsiteSize() * GJP.VolNodeSites();
    Vector *v1 = (Vector *)smalloc(cname, fname, "v1", f_size*sizeof(Float));
    Vector *v2 = (Vector *)smalloc(cname, fname, "v2", f_size*sizeof(Float));

    {
        CgArg cg_arg ;
        cg_arg.mass = mass ;
        cg_arg.epsilon = epsilon;
      
        DiracOpWilsonTm wilson(*this, v1, v2, &cg_arg, CNV_FRM_YES) ;
        // chi <- (M_f^\dag M_f)^{-1} M_f^\dag (RGV)
        wilson.CalcHmdForceVecs(chi);
    }

    this->BondCond();
    FforceWilsonType cal_force(mom, this->GaugeField(),
                               (Float *)v2, (Float *)v1, 1, 2.0 * dt);
    ForceArg ret = cal_force.run();
    this->BondCond();

    sfree(cname, fname, "v1", v1);
    sfree(cname, fname, "v2", v2);
    return ret;
}

//------------------------------------------------------------------
// EvolveMomFforce(Matrix *mom, Vector *frm, Float mass, Float dt):
// should never be called by wilsonTm fermions
//------------------------------------------------------------------
ForceArg FwilsonTm::EvolveMomFforce(Matrix *mom, Vector *chi, 
			      Float mass, Float dt)
{
  char *fname = "EvolveMomFforce(M*,V*,F,F)";
  VRB.Func(cname,fname);

  ERR.General(cname,fname,"f_wilsonTm: EvolveMomFforce(M*,V*,F,F) not implemented here\n");

  return ForceArg(0.0,0.0,0.0);
}


//------------------------------------------------------------------
// EvolveMomFforce(Matrix *mom, Vector *frm, Vector *frm,
//                 Float mass, Float epsilon, Float dt):
// It evolves the canonical momentum mom by dt
// using the boson-component of the alg_quotient force.
// Now sets epsilon in cg_arg from new input parameter
// chi <- (M_f^\dag M_f)^{-1} M_f^\dag (RGV)
// phi <- M_b (M_b^\dag M_b)^{-1} M_f^\dag (RGV)
//------------------------------------------------------------------
ForceArg FwilsonTm::EvolveMomFforce(Matrix *mom, Vector *chi, Vector *eta,
                                    Float mass, Float epsilon, Float dt)
{
    const char *fname = "EvolveMomFforce(M*,V*,V*,F,F,F)";
    if (mom == 0)            { ERR.Pointer(cname,fname,"mom") ; }
    if (eta == 0)            { ERR.Pointer(cname,fname,"eta") ; }

    const int f_size = FsiteSize() * GJP.VolNodeSites();
    Vector *v1 = (Vector *)smalloc(cname, fname, "v1", f_size*sizeof(Float)) ;
    Vector *v2 = (Vector *)smalloc(cname, fname, "v2", f_size*sizeof(Float)) ;

    {
        CgArg cg_arg ;
        cg_arg.mass = mass ;
        cg_arg.epsilon = epsilon;

        DiracOpWilsonTm wilson(*this, v1, v2, &cg_arg, CNV_FRM_YES) ;

        //~~
        //~~ fermion version:  	wilson.CalcHmdForceVecs(chi)
        //~~ boson version:  	wilson.CalcBsnForceVecs(chi, eta)
        //~~
        // chi <- (M_f^\dag M_f)^{-1} M_f^\dag (RGV)
        // phi <- M_b (M_b^\dag M_b)^{-1} M_f^\dag (RGV)
        wilson.CalcBsnForceVecs(chi, eta) ;
    }

    this->BondCond();
    FforceWilsonType cal_force(mom, this->GaugeField(),
                               (Float *)v2, (Float *)v1, 1, -2.0 * dt);
    ForceArg ret = cal_force.run();
    this->BondCond();

    sfree(cname, fname, "v1", v1);
    sfree(cname, fname, "v2", v2);
    return ret;
}

//------------------------------------------------------------------
// ForceArg EvolveMomFforce(Matrix *mom, Vector *phi, Vector *eta,
//		      Float mass, Float epsilon, Float dt)
// should never be called by wilsonTm fermions
//------------------------------------------------------------------
ForceArg FwilsonTm::EvolveMomFforce(Matrix *mom, Vector *phi, Vector *eta,
		      Float mass, Float dt) {
  char *fname = "EvolveMomFforce(M*,V*,V*,F,F)";
  ERR.General(cname,fname,"f_wilsonTm: EvolveMomFForce(M*,V*,V*,F,F) not implemented here\n");

  return ForceArg(0.0,0.0,0.0);
}

CPS_END_NAMESPACE
