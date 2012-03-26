#include<config.h>

CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of FwilsonTm class.

  $Id: f_wilsonTm.C,v 1.3 2012-03-26 13:50:12 chulwoo Exp $
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

#define BENCHMARK
#ifdef BENCHMARK
#include <util/qcdio.h>
#include <sys/time.h>
unsigned long WfmFlopsTm;
#ifndef timersub
#define timersub(a, b, result)                                                \
  do {                                                                        \
    (result)->tv_sec = (a)->tv_sec - (b)->tv_sec;                             \
    (result)->tv_usec = (a)->tv_usec - (b)->tv_usec;                          \
    if ((result)->tv_usec < 0) {                                              \
      --(result)->tv_sec;                                                     \
      (result)->tv_usec += 1000000;                                           \
    }                                                                         \
  } while (0)
#endif
#endif

CPS_START_NAMESPACE

//------------------------------------------------------------------
// Initialize static variables.
//------------------------------------------------------------------

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
  int iter;
  char *fname = "FmatEvlInv(CgArg*,V*,V*,F*,CnvFrmType)";
  VRB.Func(cname,fname);

#if 1
{
  DiracOpWilsonTm wilson(*this, f_out, f_in, cg_arg, cnv_frm);

  WfmFlopsTm = 0;
  struct timeval t_start, t_stop;
  gettimeofday(&t_start,NULL);
  
  iter = wilson.InvCg(&(cg_arg->true_rsd));
  if (true_res) *true_res = cg_arg ->true_rsd;

  gettimeofday(&t_stop,NULL);
  timersub(&t_stop,&t_start,&t_start);
  double flops= (double)WfmFlopsTm;
  double secs = t_start.tv_sec + 1.E-6 *t_start.tv_usec;
//  printf("Wilson solve: %d iteratations %d flops %f Mflops per node\n",
//	 iter,WfmFlopsTm,flops/(secs*1000000) );
 }
#endif

  
  // Return the number of iterations
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
  char *fname = "EvolveMomFforce(M*,V*,F,F,F)";
  VRB.Func(cname,fname);

  Matrix *gauge = GaugeField() ;

  if (Colors() != 3) ERR.General(cname,fname,"Wrong nbr of colors.") ;
  if (SpinComponents() != 4) ERR.General(cname,fname,"Wrong nbr of spin comp.") ;
  if (mom == 0) ERR.Pointer(cname,fname,"mom") ;
  if (chi == 0) ERR.Pointer(cname,fname,"chi") ;

//------------------------------------------------------------------
// allocate space for two CANONICAL fermion fields.
//------------------------------------------------------------------
  int f_size = FsiteSize() * GJP.VolNodeSites() ;

  char *str_v1 = "v1" ;
  Vector *v1 = (Vector *)smalloc(cname, fname, str_v1, f_size*sizeof(Float)) ;

  char *str_v2 = "v2" ;
  Vector *v2 = (Vector *)smalloc(cname, fname, str_v2, f_size*sizeof(Float)) ;

//------------------------------------------------------------------
// allocate space for two CANONICAL fermion field on a site.
//------------------------------------------------------------------

  char *str_site_v1 = "site_v1";
  Float *site_v1 = (Float *)smalloc(cname, fname, str_site_v1, FsiteSize()*sizeof(Float));

  char *str_site_v2 = "site_v2";
  Float *site_v2 = (Float *)smalloc(cname, fname, str_site_v2, FsiteSize()*sizeof(Float));

  Float L1 = 0.0;
  Float L2 = 0.0;
  Float Linf = 0.0;


  {
    CgArg cg_arg ;
    cg_arg.mass = mass ;
    cg_arg.epsilon = epsilon;

    DiracOpWilsonTm wilson(*this, v1, v2, &cg_arg, CNV_FRM_YES) ;
    // chi <- (M_f^\dag M_f)^{-1} M_f^\dag (RGV)
    wilson.CalcHmdForceVecs(chi) ;
  }
#if 0
  VRB.Result(cname,fname,"Being skipped for debugging!");
#else

  int x, y, z, t, lx, ly, lz, lt ;

  lx = GJP.XnodeSites() ;
  ly = GJP.YnodeSites() ;
  lz = GJP.ZnodeSites() ;
  lt = GJP.TnodeSites() ;

//------------------------------------------------------------------
// start by summing first over direction (mu) and then over site
// to allow SCU transfers to happen face-by-face in the outermost
// loop.
//------------------------------------------------------------------

  int mu ;

  Matrix tmp, f ;

  for (mu=0; mu<4; mu++) {
    for (t=0; t<lt; t++)
    for (z=0; z<lz; z++)
    for (y=0; y<ly; y++)
    for (x=0; x<lx; x++) {
      int gauge_offset = x+lx*(y+ly*(z+lz*t)) ;
      int vec_offset = FsiteSize()*gauge_offset ;
      gauge_offset = mu+4*gauge_offset ;

      Float *v1_plus_mu ;
      Float *v2_plus_mu ;
      int vec_plus_mu_offset = FsiteSize() ;

      Float coeff = -2.0 * dt ;

      switch (mu) {
        case 0 :
          vec_plus_mu_offset *= (x+1)%lx+lx*(y+ly*(z+lz*t)) ;
          if ((x+1) == lx) {
            getPlusData( (IFloat *)site_v1,
                         (IFloat *)v1+vec_plus_mu_offset, FsiteSize(), mu) ;
            getPlusData( (IFloat *)site_v2,
                         (IFloat *)v2+vec_plus_mu_offset, FsiteSize(), mu) ;
            v1_plus_mu = site_v1 ;                        
            v2_plus_mu = site_v2 ;                        
            if (GJP.XnodeBc()==BND_CND_APRD) coeff = -coeff ;
          } else {
            v1_plus_mu = (Float *)v1+vec_plus_mu_offset ;
            v2_plus_mu = (Float *)v2+vec_plus_mu_offset ;
          }
          break ;
        case 1 :
          vec_plus_mu_offset *= x+lx*((y+1)%ly+ly*(z+lz*t)) ;
          if ((y+1) == ly) {
            getPlusData( (IFloat *)site_v1,
                         (IFloat *)v1+vec_plus_mu_offset, FsiteSize(), mu) ;
            getPlusData( (IFloat *)site_v2,
                         (IFloat *)v2+vec_plus_mu_offset, FsiteSize(), mu) ;
            v1_plus_mu = site_v1 ;                        
            v2_plus_mu = site_v2 ;                        
            if (GJP.YnodeBc()==BND_CND_APRD) coeff = -coeff ;
          } else {
            v1_plus_mu = (Float *)v1+vec_plus_mu_offset ;
            v2_plus_mu = (Float *)v2+vec_plus_mu_offset ;
          }
          break ;
        case 2 :
          vec_plus_mu_offset *= x+lx*(y+ly*((z+1)%lz+lz*t)) ;
          if ((z+1) == lz) {
            getPlusData( (IFloat *)site_v1,
                         (IFloat *)v1+vec_plus_mu_offset, FsiteSize(), mu) ;
            getPlusData( (IFloat *)site_v2,
                         (IFloat *)v2+vec_plus_mu_offset, FsiteSize(), mu) ;
            v1_plus_mu = site_v1 ;
            v2_plus_mu = site_v2 ;
            if (GJP.ZnodeBc()==BND_CND_APRD) coeff = -coeff ;
          } else {
            v1_plus_mu = (Float *)v1+vec_plus_mu_offset ;
            v2_plus_mu = (Float *)v2+vec_plus_mu_offset ;
          }
          break ;
        case 3 :
          vec_plus_mu_offset *= x+lx*(y+ly*(z+lz*((t+1)%lt))) ;
          if ((t+1) == lt) {
            getPlusData( (IFloat *)site_v1,
                         (IFloat *)v1+vec_plus_mu_offset, FsiteSize(), mu) ;
            getPlusData( (IFloat *)site_v2,
                         (IFloat *)v2+vec_plus_mu_offset, FsiteSize(), mu) ;
            v1_plus_mu = site_v1 ;
            v2_plus_mu = site_v2 ;
            if (GJP.TnodeBc()==BND_CND_APRD) coeff = -coeff ;
          } else {
            v1_plus_mu = (Float *)v1+vec_plus_mu_offset ;
            v2_plus_mu = (Float *)v2+vec_plus_mu_offset ;
          }
      } // end switch mu

      sproj_tr[mu](   (IFloat *)&tmp,
                      (IFloat *)v1_plus_mu,
                      (IFloat *)v2+vec_offset, 1, 0, 0);

      sproj_tr[mu+4]( (IFloat *)&f,
                      (IFloat *)v2_plus_mu,
                      (IFloat *)v1+vec_offset, 1, 0, 0);

      tmp += f ;

      f.DotMEqual(*(gauge+gauge_offset), tmp) ;

      tmp.Dagger(f) ;

      f.TrLessAntiHermMatrix(tmp) ;

      f *= coeff ;

      *(mom+gauge_offset) += f ;
      Float norm = f.norm();
      Float tmp = sqrt(norm);
      L1 += tmp;
      L2 += norm;
      Linf = (tmp>Linf ? tmp : Linf);
    }
  }
#endif

//------------------------------------------------------------------
// deallocate space for two CANONICAL fermion fields on a site.
//------------------------------------------------------------------
  VRB.Sfree(cname, fname, str_site_v2, site_v2) ;
  sfree(site_v2) ;

  VRB.Sfree(cname, fname, str_site_v1, site_v1) ;
  sfree(site_v1) ;

//------------------------------------------------------------------
// deallocate space for two CANONICAL fermion fields.
//------------------------------------------------------------------
  VRB.Sfree(cname, fname, str_v2, v2) ;
  sfree(v2) ;

  VRB.Sfree(cname, fname, str_v1, v1) ;
  sfree(v1) ;

  glb_sum(&L1);
  glb_sum(&L2);
  glb_max(&Linf);

  L1 /= 4.0*GJP.VolSites();
  L2 /= 4.0*GJP.VolSites();
  VRB.FuncEnd(cname,fname);

  return ForceArg(L1, sqrt(L2), Linf);
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
		      Float mass, Float epsilon, Float dt) {
  char *fname = "EvolveMomFforce(M*,V*,V*,F,F,F)";
  VRB.Func(cname,fname);

  Matrix *gauge = GaugeField() ;

  if (Colors() != 3)       { ERR.General(cname,fname,"Wrong nbr of colors.") ; }
  if (SpinComponents()!=4) { ERR.General(cname,fname,"Wrong nbr of spin comp.") ;}
  if (mom == 0)            { ERR.Pointer(cname,fname,"mom") ; }
  if (eta == 0)            { ERR.Pointer(cname,fname,"eta") ; }

//------------------------------------------------------------------
// allocate space for two CANONICAL fermion fields.
// these are all full fermion vector sizes ( i.e. *not* preconditioned )
//------------------------------------------------------------------

  const int f_size        ( FsiteSize() * GJP.VolNodeSites() );

  char *str_v1 = "v1" ;
  Vector *v1 = (Vector *)smalloc(f_size*sizeof(Float)) ;
  if (v1 == 0) ERR.Pointer(cname, fname, str_v1) ;
  VRB.Smalloc(cname, fname, str_v1, v1, f_size*sizeof(Float)) ;

  char *str_v2 = "v2" ;
  Vector *v2 = (Vector *)smalloc(f_size*sizeof(Float)) ;
  if (v2 == 0) ERR.Pointer(cname, fname, str_v2) ;
  VRB.Smalloc(cname, fname, str_v2, v2, f_size*sizeof(Float)) ;

//------------------------------------------------------------------
// two fermion vectors at a single position
//    - these will be used to store off-node
//      field components
//------------------------------------------------------------------

  char *str_site_v1 = "site_v1";
  Float *site_v1 = (Float *)smalloc(FsiteSize()*sizeof(Float));
  if (site_v1 == 0) ERR.Pointer(cname, fname, str_site_v1) ;
  VRB.Smalloc(cname, fname, str_site_v1, site_v1,
    FsiteSize()*sizeof(Float)) ;

  char *str_site_v2 = "site_v2";
  Float *site_v2 = (Float *)smalloc(FsiteSize()*sizeof(Float));
  if (site_v2 == 0) ERR.Pointer(cname, fname, str_site_v2) ;
  VRB.Smalloc(cname, fname, str_site_v2, site_v2,
    FsiteSize()*sizeof(Float)) ;
  
  Float L1 = 0.0;
  Float L2 = 0.0;
  Float Linf = 0.0;


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

#if 0
  VRB.Result(cname,fname,"Being skipped for debugging!");
#else

  // evolve the momenta by the fermion force
  int mu, x, y, z, t;
 
  const int lx(GJP.XnodeSites());
  const int ly(GJP.YnodeSites());
  const int lz(GJP.ZnodeSites());
  const int lt(GJP.TnodeSites());

//------------------------------------------------------------------
// start by summing first over direction (mu) and then over site
// to allow SCU transfers to happen face-by-face in the outermost
// loop.
//------------------------------------------------------------------

// (1-gamma_\mu) Tr_s[v1(x+\mu) v2^{\dagger}(x)] +          
// 		(1+gamma_\mu) Tr_s [v2(x+\mu) v1^{\dagger}(x)]
  for (mu=0; mu<4; mu++) {
    for (t=0; t<lt; t++)
      for (z=0; z<lz; z++)
        for (y=0; y<ly; y++)
          for (x=0; x<lx; x++) {
            // position offset
            int gauge_offset = x+lx*(y+ly*(z+lz*t)) ;
            
            // offset for vector field at this point
            int vec_offset = FsiteSize()*gauge_offset ;

            // offset for link in mu direction from this point
            gauge_offset = mu+4*gauge_offset ;

            Float *v1_plus_mu ;
            Float *v2_plus_mu ;
            int vec_plus_mu_offset = FsiteSize() ;

            // sign of coeff (look at momenta update)
            Float coeff = -2.0 * dt ;

            switch (mu) {
              case 0 :
                // next position in mu direction
                vec_plus_mu_offset *= (x+1)%lx+lx*(y+ly*(z+lz*t)) ;
                if ((x+1) == lx) {
                   // off-node
                   // fill site_v1 and site_v2 with v1 and v2 data
                   // from x=0 on next node, need loop because
                   // data is not contiguous in memory 
                   getPlusData( (IFloat *)site_v1,
                              (IFloat *)v1+vec_plus_mu_offset, FsiteSize(), mu) ;
                   getPlusData( (IFloat *)site_v2,
                              (IFloat *)v2+vec_plus_mu_offset, FsiteSize(), mu) ;

                   v1_plus_mu = site_v1 ;                        
                   v2_plus_mu = site_v2 ;                        

                   // GJP.XnodeBc() gives the forward boundary
                   // condition only (so this should work).
                   if (GJP.XnodeBc()==BND_CND_APRD) coeff = -coeff ;
                } else {
                   //
                   //  on - node: just add offset to v1 and v2
                   // (they are now 1 forward in the mu direction )
                   //
                   v1_plus_mu = (Float *)v1+vec_plus_mu_offset ;
                   v2_plus_mu = (Float *)v2+vec_plus_mu_offset ;
                }
                break ;

        // Repeat for the other directions
        case 1 :
          vec_plus_mu_offset *= x+lx*((y+1)%ly+ly*(z+lz*t)) ;
          if ((y+1) == ly) {
            getPlusData( (IFloat *)site_v1,
                         (IFloat *)v1+vec_plus_mu_offset, FsiteSize(), mu) ;
            getPlusData( (IFloat *)site_v2,
                         (IFloat *)v2+vec_plus_mu_offset, FsiteSize(), mu) ;
            v1_plus_mu = site_v1 ;                        
            v2_plus_mu = site_v2 ;                        
            if (GJP.YnodeBc()==BND_CND_APRD) coeff = -coeff ;
          } else {
            v1_plus_mu = (Float *)v1+vec_plus_mu_offset ;
            v2_plus_mu = (Float *)v2+vec_plus_mu_offset ;
          }
          break ;
        case 2 :
          vec_plus_mu_offset *= x+lx*(y+ly*((z+1)%lz+lz*t)) ;
          if ((z+1) == lz) {
            getPlusData( (IFloat *)site_v1,
                         (IFloat *)v1+vec_plus_mu_offset, FsiteSize(), mu) ;
            getPlusData( (IFloat *)site_v2,
                         (IFloat *)v2+vec_plus_mu_offset, FsiteSize(), mu) ;
            v1_plus_mu = site_v1 ;
            v2_plus_mu = site_v2 ;
            if (GJP.ZnodeBc()==BND_CND_APRD) coeff = -coeff ;
          } else {
            v1_plus_mu = (Float *)v1+vec_plus_mu_offset ;
            v2_plus_mu = (Float *)v2+vec_plus_mu_offset ;
          }
          break ;
        case 3 :
          vec_plus_mu_offset *= x+lx*(y+ly*(z+lz*((t+1)%lt))) ;
          if ((t+1) == lt) {
            getPlusData( (IFloat *)site_v1,
                         (IFloat *)v1+vec_plus_mu_offset, FsiteSize(), mu) ;
            getPlusData( (IFloat *)site_v2,
                         (IFloat *)v2+vec_plus_mu_offset, FsiteSize(), mu) ;
            v1_plus_mu = site_v1 ;
            v2_plus_mu = site_v2 ;
            if (GJP.TnodeBc()==BND_CND_APRD) coeff = -coeff ;
          } else {
            v1_plus_mu = (Float *)v1+vec_plus_mu_offset ;
            v2_plus_mu = (Float *)v2+vec_plus_mu_offset ;
          }
      } // end (the evil) mu switch 

      Matrix tmp_mat1, tmp_mat2;  

// ( 1 - gamma_\mu ) Tr_s [ v1(x+\mu) v2^{\dagger}(x) ]           
      sproj_tr[mu](   (IFloat *)&tmp_mat1,   	// output color matrix
                      (IFloat *)v1_plus_mu,		// row vector, NOT conjugated
                      (IFloat *)v2+vec_offset, 	// col vector, IS conjugated
                      1, 0, 0);				// 1 block, 0 strides

// (1 + gamma_\mu)  Tr_s [ v2(x+\mu) v1^{\dagger}(x) ]
      sproj_tr[mu+4]( (IFloat *)&tmp_mat2,		// output color matrix
                      (IFloat *)v2_plus_mu,		// row vector, NOT conjugated
                      (IFloat *)v1+vec_offset, 	// col vector, IS conjugated
                      1, 0, 0);				// 1 block, 0 strides

      // exactly what this sounds like
      tmp_mat1 += tmp_mat2 ;
            
      // multiply sum by the link in the \mu direction
      tmp_mat2.DotMEqual(*(gauge+gauge_offset), tmp_mat1) ;
      
      // take tracless antihermitian piece
      // TrLessAntiHermMatrix need to be passed
      // the dagger of the matrix in question
      tmp_mat1.Dagger(tmp_mat2) ;
      tmp_mat2.TrLessAntiHermMatrix(tmp_mat1) ;

      tmp_mat2 *= coeff ;
            
//~~
//~~ fermion version:  	(mom+gauge_offset) += f
//~~ boson version:  	(mom+gauge_offset) -= f
//~~
      *(mom+gauge_offset) -= tmp_mat2 ;

	 Float norm = tmp_mat2.norm();
	 Float tmp = sqrt(norm);
	 L1 += tmp;
	 L2 += norm;
	 Linf = (tmp>Linf ? tmp : Linf);
	 
    }
  }

//------------------------------------------------------------------
// deallocate space for two CANONICAL fermion fields on a site.
//------------------------------------------------------------------
  VRB.Sfree(cname, fname, str_site_v2, site_v2) ;
  sfree(site_v2) ;

  VRB.Sfree(cname, fname, str_site_v1, site_v1) ;
  sfree(site_v1) ;

//------------------------------------------------------------------
// deallocate space for two CANONICAL fermion fields.
//------------------------------------------------------------------
  VRB.Sfree(cname, fname, str_v2, v2) ;
  sfree(v2) ;

  VRB.Sfree(cname, fname, str_v1, v1) ;
  sfree(v1) ;

  glb_sum(&L1);
  glb_sum(&L2);
  glb_max(&Linf);

  L1 /= 4.0*GJP.VolSites();
  L2 /= 4.0*GJP.VolSites();
#endif

  VRB.FuncEnd(cname,fname);
  return ForceArg(L1, sqrt(L2), Linf);
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
