#include<config.h>

CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of FwilsonTm class.

*/
//------------------------------------------------------------------
//
// f_wilsonTm.C
//
// FwilsonTm is derived from Fwilson and is relevant to
// twisted-mass wilson fermions
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#ifdef USE_BFM_TM
#include <util/lattice/bfm_evo.h>
#define Printf if ( !UniqueID() ) printf
#endif
#include <util/lattice.h>
#include <util/qcdio.h>
#include <util/dirac_op.h>
#include <util/wilson.h>
#include <util/verbose.h>
#include <util/gjp.h>
#include <util/error.h>
#include <util/enum_func.h> //Added by CK for access to NumChkb
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
//CK: Well it should be! Added below
//------------------------------------------------------------------
int FwilsonTm::FmatEvlMInv(Vector **f_out, Vector *f_in, Float *shift, 
			 int Nshift, int isz, CgArg **cg_arg,
			 CnvFrmType cnv_frm, MultiShiftSolveType type, 
			 Float *alpha, Vector **f_out_d)
{
  char *fname = "FmatMInv(V*, V*, .....)";
  VRB.Func(cname,fname);

  size_t f_size = GJP.VolNodeSites() * FsiteSize() / (FchkbEvl()+1);
  if(GJP.Gparity()) f_size *= 2;
  Float dot = f_in -> NormSqGlbSum4D(f_size);

  Float *RsdCG = new Float[Nshift];
  for (int s=0; s<Nshift; s++) RsdCG[s] = cg_arg[s]->stop_rsd;

  //Fake the constructor
  DiracOpWilsonTm wilson(*this, f_out[0], f_in, cg_arg[0], cnv_frm);

  int return_value = wilson.MInvCG(f_out,f_in,dot,shift,Nshift,isz,RsdCG,type,alpha);  

  for (int s=0; s<Nshift; s++) cg_arg[s]->true_rsd = RsdCG[s];
  delete[] RsdCG;
  return return_value;
}


//------------------------------------------------------------------
// void FminResExt not changed for twisted-mass Wilson fermions
//------------------------------------------------------------------

//------------------------------------------------------------------
// int FmatInv not changed for twisted-mass Wilson fermions
// NOTE: this call MatInv
//------------------------------------------------------------------

//------------------------------------------------------------------
// int FeigSolv not changed for twisted-mass Wilson fermions
//CK: This is incorrect as the RitzEig function calls up to MatPcMatPcDag etc which differ between wilson and wilsonTm. We 
//    must therefore use the correct Dirac operator
//    I have also added a BFM version used when the compile switch USE_BFM_TM is active
//------------------------------------------------------------------

#ifndef BFM_GPARITY
int FwilsonTm::FeigSolv(Vector **f_eigenv, Float *lambda,
			Float *chirality, int *valid_eig,
			Float **hsum,
			EigArg *eig_arg, 
			CnvFrmType cnv_frm)
{
  char *fname = "FeigSolv(EigArg*,V*,F*,CnvFrmType)";
  VRB.Func(cname,fname);

  //Greg: I don't feel like added epsilon to EigArg
  VRB.Result(cname, fname, "Warning!! FwilsonTm::FeigSolv assumes epsilon = 0!!\n");

#if 1
  //CK: The regular CPS version

  int iter;

  CgArg cg_arg;
  cg_arg.mass = eig_arg->mass;
  cg_arg.epsilon = 0; //eig_arg->epsilon; //CK: passes down epsilon parameter for twisted mass Wilson fermions. Irrelevant here but this same function is used in FwilsonTm
  cg_arg.RitzMatOper = eig_arg->RitzMatOper;
  int N_eig = eig_arg->N_eig;
  int i;

  //=========================
  // convert fermion field
  //=========================

  if(cnv_frm == CNV_FRM_YES) //Fixed by CK to allow for single checkerboard input vectors as used in AlgActionRational. Previously it would always convert from CANONICAL to WILSON
  for(i=0; i < N_eig; ++i)  Fconvert(f_eigenv[i], WILSON, CANONICAL);

  //------------------------------------------------------------------
  //  we want both the eigenvalues of D_{hermitian} and
  //  D^{+}D.  To not change the arguments passed to RitzEig,
  //  we pass a float pointer which points to 2 * N_eig values
  //  and return both lambda and lambda^2 from RitzEig
  //------------------------------------------------------------------

  Float * lambda2 = (Float * ) smalloc (N_eig*2*sizeof(Float));
  if ( lambda2 == 0 ) ERR.Pointer(cname,fname, "lambda2");
  
  {
    DiracOpWilsonTm wilson(*this, (Vector*) 0 , (Vector*) 0, &cg_arg, CNV_FRM_NO);
    iter = wilson.RitzEig(f_eigenv, lambda2, valid_eig, eig_arg);
  }

  if(cnv_frm == CNV_FRM_YES) 
    for(i=0; i < N_eig; ++i) Fconvert(f_eigenv[i], CANONICAL, WILSON);


  /*
    the call to RitzEig returns a negative number if either the KS or CG maxes
    out, we wish to cope with this in alg_eig, so "pass it up". Clean up the
    storage order first in case we still want to use the eigenvectors as a
    guess.
  */
  if ( iter < 0 ) { return iter ; }


  // Compute chirality
  int Ncb = NumChkb(cg_arg.RitzMatOper);
  size_t f_size = (GJP.VolNodeSites() * FsiteSize()) * Ncb / 2; //CK: fixed

  Vector* v1 = (Vector *)smalloc(f_size*sizeof(Float));
  if (v1 == 0)
    ERR.Pointer(cname, fname, "v1");
  VRB.Smalloc(cname, fname, "v1", v1, f_size*sizeof(Float));

  int nspinvect = GJP.VolNodeSites() * Ncb/2;

  for(i=0; i < N_eig; ++i)
  {
    Gamma5(v1, f_eigenv[i], nspinvect);
    chirality[i] = f_eigenv[i]->ReDotProductGlbSum4D(v1, f_size);
  }

  VRB.Sfree(cname, fname, "v1", v1);
  sfree(v1);


  // rescale wilson eigenvalues to the convention  m + Dslash(U)
  Float factor = 4.0 + eig_arg->mass;
    
  
  FILE* fp=Fopen(eig_arg->fname,"a");
  for(i=0; i<N_eig; ++i)
    {
      lambda2[i] *= factor;	 		 //rescale eigenvalue
      lambda2[N_eig + i] *= ( factor * factor ); //rescale squared evalue
      lambda[i]=lambda2[i];                      //copy back
      
      //print out eigenvalue, eigenvalue^2, chirality 
      Fprintf(fp,"%d %g %g %g %d\n",i,
              (float)lambda2[i],
              (float)lambda2[N_eig + i],
	      (float)chirality[i],valid_eig[i]);
    }
  Fclose(fp);
  sfree(lambda2); 


  // Slice-sum the eigenvector density to make a 1D vector
  if (eig_arg->print_hsum){
    for(i=0; i < N_eig; ++i){
      //CK: The vector needs to be in canonical ordering. Thus if CNV_FRM_NO we need to convert to CANONICAL. If Ncb==1 we have
      //    only the odd part, so we will need to fill the even part with zeroes prior to conversion
      Vector* tosum = f_eigenv[i];
      if(cnv_frm == CNV_FRM_NO){
	//Create a temp copy of the eigenvector
	int alloc_size = f_size; if(Ncb==1) alloc_size *= 2;
	Float* full = (Float *)smalloc(alloc_size*sizeof(Float));
	for(int j=0;j<f_size;j++) full[j] = ((Float*)f_eigenv[i])[j];

	//Fill in even part with zero for Ncb==1
	if(Ncb==1) for(int j=f_size;j<alloc_size;j++) full[j] = 0; //zero even part

	//Convert
	if(cnv_frm == CNV_FRM_NO) Fconvert((Vector*)full, CANONICAL, WILSON);
	tosum = (Vector*)full;
      }
      tosum->NormSqArraySliceSum(hsum[i], FsiteSize(), eig_arg->hsum_dir);

      if(cnv_frm == CNV_FRM_NO) sfree(tosum);

    }
    
  }

  // The remaining part in QCDSP version are all about "downloading
  // eigenvectors", supposedly not applicable here.

  // Return the number of iterations
  return iter;
#endif
}

#else
int FwilsonTm::FeigSolv(Vector **f_eigenv, Float *lambda,
			Float *chirality, int *valid_eig,
			Float **hsum,
			EigArg *eig_arg, 
			CnvFrmType cnv_frm)
{
  char *fname = "FeigSolv(EigArg*,V*,F*,CnvFrmType)";
  VRB.Func(cname,fname);

#ifdef USE_BFM_TM
  //CK: Added a BFM implementation using Hantao's bfm_evo class

  // only 1 eigenvalue can be computed now.
  if(eig_arg->N_eig != 1) {
    ERR.General(cname,fname,"BFM FeigSolv can only calculate a single eigenvalue\n");
  }
  if(eig_arg->RitzMatOper != MATPCDAG_MATPC &&
     eig_arg->RitzMatOper != NEG_MATPCDAG_MATPC) {
    ERR.NotImplemented(cname, fname);
  }
    
  if(cnv_frm == CNV_FRM_YES) ERR.General(cname,fname,"BFM FeigSolv not implemented for non-Wilson ordered fermions");

  bfmarg bfm_arg;

  int nthread = GJP.SetNthreads();

  bfm_arg.solver = WilsonTM;
  for(int i=0;i<4;i++) bfm_arg.node_latt[i] = GJP.NodeSites(i);
  bfm_arg.verbose=1;
  bfm_arg.reproduce=0;

  bfmarg::Threads(nthread);
  bfmarg::Reproduce(0);
  bfmarg::ReproduceChecksum(0);
  bfmarg::ReproduceMasterCheck(0);
  bfmarg::Verbose(1);
  bfmarg::onepluskappanorm = 0;

  multi1d<int> ncoor = QDP::Layout::nodeCoord();
  multi1d<int> procs = QDP::Layout::logicalSize();

  if(GJP.Gparity()){
    bfm_arg.gparity = 1;
    Printf("G-parity directions: ");
    for(int d=0;d<3;d++)
      if(GJP.Bc(d) == BND_CND_GPARITY){ bfm_arg.gparity_dir[d] = 1; Printf("%d ",d); }
      else bfm_arg.gparity_dir[d] = 0;
    for(int d=0;d<4;d++){
      bfm_arg.nodes[d] = procs[d];
      bfm_arg.ncoor[d] = ncoor[d];
    }
    Printf("\n");
  }


  Printf("%d dim machine\n\t", procs.size());
  for(int mu=0;mu<4;mu++){
    Printf("%d ", procs[mu]);
    if ( procs[mu]>1 ) bfm_arg.local_comm[mu] = 0;
    else bfm_arg.local_comm[mu] = 1;
  }
  Printf("\nLocal comm = ");
  for(int mu=0;mu<4;mu++){
    Printf("%d ", bfm_arg.local_comm[mu]);
  }
  Printf("\n");

  double mq= eig_arg->mass;
  double epsilon = eig_arg->epsilon;
  Printf("mq=%g epsilon=%g\n",mq,epsilon);

  bfm_arg.precon_5d = 0;
  bfm_arg.Ls   = 1;
  bfm_arg.M5   = 0.0;
  bfm_arg.mass = toDouble(mq);
  bfm_arg.twistedmass = toDouble(epsilon);
  bfm_arg.Csw  = 0.0;
  bfm_arg.max_iter = 10000;
  bfm_arg.residual = eig_arg->Rsdlam;
  bfm_arg.max_iter = eig_arg->MaxCG;
  Printf("Initialising bfm operator\n");

  bfm_evo<double> bd;
  bd.init(bfm_arg);

  VRB.Result(cname, fname, "residual = %17.10e max_iter = %d mass = %17.10e\n",
	     bd.residual, bd.max_iter, bd.mass);
  BondCond();
  bd.cps_importGauge((Float*)GaugeField());

  Fermion_t in = bd.allocFermion();
  bd.cps_impexcbFermion((Float *)f_eigenv[0], in, 1, 1);

#pragma omp parallel
  {
    lambda[0] = bd.ritz(in, eig_arg->RitzMatOper == MATPCDAG_MATPC);
  }

  bd.cps_impexcbFermion((Float *)f_eigenv[0], in, 0, 1);

  // correct the eigenvalue for a dumb convention problem.
  if(eig_arg->RitzMatOper == NEG_MATPCDAG_MATPC) lambda[0] = -lambda[0];

  valid_eig[0] = 1;
  bd.freeFermion(in);

  BondCond();

  //Compute chirality and whatnot
  size_t f_size = (GJP.VolNodeSites() * FsiteSize())/ 2; //single checkerboard field
  if(GJP.Gparity()) f_size *= 2;
  Vector* v1 = (Vector *)pmalloc(f_size*sizeof(Float));

  int nspinvect = GJP.VolNodeSites()/2;
  if(GJP.Gparity()) nspinvect *= 2;

  Gamma5(v1, f_eigenv[0], nspinvect);
  chirality[0] = f_eigenv[0]->ReDotProductGlbSum4D(v1, f_size);
  pfree(v1);

  // Regular CPS code rescales wilson eigenvalues to the convention  m + Dslash(U), i.e. lambda *= (4+m)
  // There is also a normalization difference between CPS and BFM preconditioned matrices:
  //MdagM_BFM = 0.25/kappa^2 * MdagM CPS
  //Hence should multiply BFM eigenvalues by 4*kappa^2 = 1/(4+m)^2  as well as the above factor 
  
  //For twisted mass
  Float kappa = 1.0/2.0/sqrt( (eig_arg->mass + 4.0)*(eig_arg->mass + 4.0) + eig_arg->epsilon * eig_arg->epsilon );
  Float factor = 4*kappa*kappa*(4+eig_arg->mass);

//1.0/(4.0 + eig_arg->mass);
   
  FILE* fp=Fopen(eig_arg->fname,"a");
  lambda[0] *= factor;  //rescale eigenvalue
  lambda[1] =  lambda[0]*lambda[0]; //squared evalue
  
  //print out eigenvalue, eigenvalue^2, chirality 
  Fprintf(fp,"%d %g %g %g %d\n",0,
	  (float)lambda[0],
	  (float)lambda[1],
	  (float)chirality[0],valid_eig[0]);

  Fclose(fp);

  if (eig_arg->print_hsum) ERR.General(cname,fname,"Hsum print not yet implemented (too lazy to copy/paste code)");

  return 0;


#else
  //CK: The regular CPS version

  int iter;

  CgArg cg_arg;
  cg_arg.mass = eig_arg->mass;
  cg_arg.epsilon = eig_arg->epsilon; //CK: passes down epsilon parameter for twisted mass Wilson fermions. Irrelevant here but this same function is used in FwilsonTm
  cg_arg.RitzMatOper = eig_arg->RitzMatOper;
  int N_eig = eig_arg->N_eig;
  int i;

  //=========================
  // convert fermion field
  //=========================

  if(cnv_frm == CNV_FRM_YES) //Fixed by CK to allow for single checkerboard input vectors as used in AlgActionRational. Previously it would always convert from CANONICAL to WILSON
    for(i=0; i < N_eig; ++i)  Fconvert(f_eigenv[i], WILSON, CANONICAL);

  //------------------------------------------------------------------
  //  we want both the eigenvalues of D_{hermitian} and
  //  D^{+}D.  To not change the arguments passed to RitzEig,
  //  we pass a float pointer which points to 2 * N_eig values
  //  and return both lambda and lambda^2 from RitzEig
  //------------------------------------------------------------------

  Float * lambda2 = (Float * ) smalloc (cname,fname, "lambda2",N_eig*2*sizeof(Float));
  {
    DiracOpWilsonTm wilson(*this, (Vector*) 0 , (Vector*) 0, &cg_arg, CNV_FRM_NO);
    iter = wilson.RitzEig(f_eigenv, lambda2, valid_eig, eig_arg);
  }

  if(cnv_frm == CNV_FRM_YES) 
    for(i=0; i < N_eig; ++i) Fconvert(f_eigenv[i], CANONICAL, WILSON);


  /*
    the call to RitzEig returns a negative number if either the KS or CG maxes
    out, we wish to cope with this in alg_eig, so "pass it up". Clean up the
    storage order first in case we still want to use the eigenvectors as a
    guess.
  */
  if ( iter < 0 ) { return iter ; }


  // Compute chirality
  int Ncb = NumChkb(cg_arg.RitzMatOper);
  size_t f_size = (GJP.VolNodeSites() * FsiteSize()) * Ncb / 2; //CK: fixed
  if(GJP.Gparity()) f_size *= 2;

  Vector* v1 = (Vector *)smalloc(f_size*sizeof(Float));
  if (v1 == 0)
    ERR.Pointer(cname, fname, "v1");
  VRB.Smalloc(cname, fname, "v1", v1, f_size*sizeof(Float));

  int nspinvect = GJP.VolNodeSites() * Ncb/2;
  if(GJP.Gparity()) nspinvect *= 2;

  for(i=0; i < N_eig; ++i)
  {
    Gamma5(v1, f_eigenv[i], nspinvect);
    chirality[i] = f_eigenv[i]->ReDotProductGlbSum4D(v1, f_size);
  }

  VRB.Sfree(cname, fname, "v1", v1);
  sfree(v1);


  // rescale wilson eigenvalues to the convention  m + Dslash(U)
  Float factor = 4.0 + eig_arg->mass;
    
  
  FILE* fp=Fopen(eig_arg->fname,"a");
  for(i=0; i<N_eig; ++i)
    {
      lambda2[i] *= factor;	 		 //rescale eigenvalue
      lambda2[N_eig + i] *= ( factor * factor ); //rescale squared evalue
      lambda[i]=lambda2[i];                      //copy back
      
      //print out eigenvalue, eigenvalue^2, chirality 
      Fprintf(fp,"%d %g %g %g %d\n",i,
              (float)lambda2[i],
              (float)lambda2[N_eig + i],
	      (float)chirality[i],valid_eig[i]);
    }
  Fclose(fp);
  sfree(lambda2); 


  // Slice-sum the eigenvector density to make a 1D vector
  if (eig_arg->print_hsum){
    for(i=0; i < N_eig; ++i){
      //CK: The vector needs to be in canonical ordering. Thus if CNV_FRM_NO we need to convert to CANONICAL. If Ncb==1 we have
      //    only the odd part, so we will need to fill the even part with zeroes prior to conversion
      Vector* tosum = f_eigenv[i];
      if(cnv_frm == CNV_FRM_NO){
	//Create a temp copy of the eigenvector
	int alloc_size = f_size; if(Ncb==1) alloc_size *= 2;
	Float* full = (Float *)smalloc(alloc_size*sizeof(Float));
	for(int j=0;j<f_size;j++) full[j] = ((Float*)f_eigenv[i])[j];

	//Fill in even part with zero for Ncb==1
	if(Ncb==1) for(int j=f_size;j<alloc_size;j++) full[j] = 0; //zero even part

	//Convert
	if(cnv_frm == CNV_FRM_NO) Fconvert((Vector*)full, CANONICAL, WILSON);
	tosum = (Vector*)full;
      }
      tosum->NormSqArraySliceSum(hsum[i], FsiteSize(), eig_arg->hsum_dir);

      if(cnv_frm == CNV_FRM_NO) sfree(tosum);

    }
    
  }

  // The remaining part in QCDSP version are all about "downloading
  // eigenvectors", supposedly not applicable here.

  // Return the number of iterations
  return iter;
#endif
}
#endif

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
//CK: Of course it will not work, because calls to EvolveMomFforce will cause an error as they are not implemented unless the epsilon parameter is passed
//------------------------------------------------------------------

//CK: added one that passes down the epsilon parameter
ForceArg FwilsonTm::RHMC_EvolveMomFforce(Matrix *mom, Vector **sol, int degree,
					 int isz, Float *alpha, Float mass, Float epsilon,
					 Float dt, Vector **sol_d, 
					 ForceMeasure force_measure) {
  const char *fname = "RHMC_EvolveMomFforce";
  char *force_label;

  ForceArg Fdt;
  Float L1 = 0.0;
  Float L2 = 0.0;
  Float Linf = 0.0;

  int g_size = GJP.VolNodeSites() * GsiteSize();
  if(GJP.Gparity()) g_size *= 2;

  Matrix *mom_tmp;

  if (force_measure == FORCE_MEASURE_YES) {
    mom_tmp = (Matrix*)smalloc(g_size*sizeof(Float),cname, fname, "mom_tmp");
    ((Vector*)mom_tmp) -> VecZero(g_size);
    force_label = new char[100];
  } else {
    mom_tmp = mom;
  }

  for (int i=0; i<degree; i++) {
    ForceArg Fdt = EvolveMomFforce(mom_tmp,sol[i],mass,epsilon,dt*alpha[i]);
    if (force_measure == FORCE_MEASURE_YES) {
      sprintf(force_label, "Rational, mass = %e, pole = %d:", mass, i+isz);
      Fdt.print(dt, force_label);
    }
  }

  // If measuring the force, need to measure and then sum to mom
  if (force_measure == FORCE_MEASURE_YES) {
    for (int i=0; i<g_size/18; i++) {
      Float norm = (mom_tmp+i)->norm();
      Float tmp = sqrt(norm);
      L1 += tmp;
      L2 += norm;
      Linf = (tmp>Linf ? tmp : Linf);
    }
    glb_sum(&L1);
    glb_sum(&L2);
    glb_max(&Linf);

    L1 /= 4.0*GJP.VolSites();
    L2 /= 4.0*GJP.VolSites();

    fTimesV1PlusV2((IFloat*)mom, 1.0, (IFloat*)mom_tmp, (IFloat*)mom, g_size);

    delete[] force_label;
    sfree(mom_tmp, cname, fname, "mom_tmp");
  }

  return ForceArg(L1, sqrt(L2), Linf);
}


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

  size_t f_size = (GJP.VolNodeSites() * FsiteSize()) >> 1 ;
  if(GJP.Gparity()) f_size*=2;

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
#if 0
    const char *fname = "EvolveMomFforce(M*,V*,F,F,F)";
    if (mom == 0) ERR.Pointer(cname,fname,"mom");
    if (chi == 0) ERR.Pointer(cname,fname,"chi");

    const size_t f_size = FsiteSize() * GJP.VolNodeSites();
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
#else
  if(GJP.Gparity()) return EvolveMomFforceGparity(mom,chi,mass,epsilon,dt);

  char *fname = "EvolveMomFforce(M*,V*,F,F,F)";
  VRB.Func(cname,fname);

  Matrix *gauge = GaugeField() ;

  if (Colors() != 3) ERR.General(cname,fname,"Wrong nbr of colors.") ;
  if (SpinComponents() != 4) ERR.General(cname,fname,"Wrong nbr of spin comp.") ;
  if (mom == 0) ERR.Pointer(cname,fname,"mom") ;
  if (chi == 0) ERR.Pointer(cname,fname,"chi") ;
  if(GJP.Gparity()) ERR.General(cname,fname,"Use EvolveMomFforceGparity for G-parity BCs");

//------------------------------------------------------------------
// allocate space for two CANONICAL fermion fields.
//------------------------------------------------------------------
  size_t f_size = FsiteSize() * GJP.VolNodeSites() ;

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

  Matrix *gparity_1f_mombuf;
  if(GJP.Gparity1fX()){
    gparity_1f_mombuf = (Matrix *)fmalloc(cname,fname,"gparity_1f_mombuf",4 * GJP.VolNodeSites() * sizeof(Matrix) ) ;
    for(int i=0;i<4*GJP.VolNodeSites();i++) gparity_1f_mombuf[i].ZeroMatrix();
  }

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

  int lattsz[4] = {lx,ly,lz,lt};

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
      
      /* CK: Replaced nasty 300000 line switch statement with the following:*/
      {
	int pos_p_mu[] = {x,y,z,t};
	pos_p_mu[mu] = (pos_p_mu[mu]+1) % lattsz[mu];
	vec_plus_mu_offset *= pos_p_mu[0] + lx*(pos_p_mu[1] + ly*(pos_p_mu[2] + lz*pos_p_mu[3]));
      }      
      int pos[4] = {x,y,z,t};

      if ((pos[mu]+1) == lattsz[mu]) {
	getPlusData( (IFloat *)site_v1,
		     (IFloat *)v1+vec_plus_mu_offset, FsiteSize(), mu) ;
	getPlusData( (IFloat *)site_v2,
		     (IFloat *)v2+vec_plus_mu_offset, FsiteSize(), mu) ;
	v1_plus_mu = site_v1 ;                        
	v2_plus_mu = site_v2 ;                        
	if (GJP.NodeBc(mu)==BND_CND_APRD) coeff = -coeff ;
      } else {
	v1_plus_mu = (Float *)v1+vec_plus_mu_offset ;
	v2_plus_mu = (Float *)v2+vec_plus_mu_offset ;
      }

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

      if(GJP.Gparity1fX()) *(gparity_1f_mombuf+gauge_offset) = f;

      Float norm = f.norm();
      Float tmp = sqrt(norm);
      L1 += tmp;
      L2 += norm;
      Linf = (tmp>Linf ? tmp : Linf);
    }
  }

  //G-parity 1f add \delta p from other 'flavour'
  //Providing you set up the 1f prop source correctly (minus sign on UR quadrant) this code also works for
  //1f G-parity in both X and Y directions
  if(GJP.Gparity1fX()){
    int momsz = GsiteSize() * GJP.VolNodeSites();
    Matrix *buf2 = (Matrix *)fmalloc(cname,fname,"buf2",momsz * sizeof(Float) ) ;

    //Communicate \delta p from first half onto second half and vice versa
    Matrix *data_buf = gparity_1f_mombuf;
    Matrix *send_buf = data_buf;
    Matrix *recv_buf = buf2;

    int vol = GJP.VolNodeSites();
    int size[4] = {GJP.XnodeSites(),GJP.YnodeSites(),GJP.ZnodeSites(),GJP.TnodeSites()};

    if(GJP.Xnodes()>1){
      //pass between nodes
      for(int i=0;i<GJP.Xnodes()/2;i++){
	getMinusData((Float *)recv_buf, (Float *)send_buf, momsz , 0);
	data_buf = recv_buf;
	recv_buf = send_buf;
	send_buf = data_buf;
      }
    }else{
      //shift field by xsites/2
#ifdef USE_OMP
#pragma omp parallel for default(shared) private(mu)
      for(long i=0;i<vol*4;i++){
#else
      for(long i=0;i<vol*4;i++){
#endif
	  //i = mu + 4*(x + Lx*(y+Ly*(z+Lz*t) ) )
	  int mu = i%4;
	  int x = (i/4) % GJP.XnodeSites();
	  int pos_rem = i/4/GJP.XnodeSites(); //(y+Ly*(z+Lz*t)

	  int x_from = (x + GJP.XnodeSites()/2) % GJP.XnodeSites();
	  int i_from = mu + 4*(x_from + GJP.XnodeSites()*pos_rem);
	  buf2[i] = gparity_1f_mombuf[i_from];
	}
	data_buf = buf2;
      }

#ifdef USE_OMP
#pragma omp parallel for default(shared) private(mu)
      for (long i=0;i<vol*4;i++){
  	int pos[4];
	long rest=i;
	mu = rest%4; rest = rest/4;
	for(int j =0; j<4;j++){
	  pos[j]= rest%size[j]; rest = rest/size[j];
	}
#else
	int pos[4];
	for (mu=0; mu<4; mu++)
	  for (pos[3]=0; pos[3]<size[3]; pos[3]++)
	    for (pos[2]=0; pos[2]<size[2]; pos[2]++)
	      for (pos[1]=0; pos[1]<size[1]; pos[1]++)
		for (pos[0]=0; pos[0]<size[0]; pos[0]++){
#endif
		  int gauge_offset = mu + 4*(pos[0]+size[0]*(pos[1]+size[1]*(pos[2]+size[2]*pos[3])));
      
		  //complex conjugate the \delta p from the other flavour
		  Float *m = (Float*) &data_buf[gauge_offset];
		  for(int c=1;c<18;c+=2) m[c] *= -1;
      
		  //add it to the momentum at this site
		  *(mom+gauge_offset) += data_buf[gauge_offset];
  }
   ffree(buf2);
 }

    
 if(GJP.Gparity1fX()) ffree(gparity_1f_mombuf,cname,fname,"gparity_1f_mombuf");

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
// EvolveMomFforce(Matrix *mom, Vector *frm, Float mass, 
//                 Float epsilon, Float dt):
// It evolves the canonical momentum mom by dt
// using the fermion force.
// Now sets epsilon in cg_arg from new input parameter
// chi <- (M_f^\dag M_f)^{-1} M_f^\dag (RGV)
//------------------------------------------------------------------
ForceArg FwilsonTm::EvolveMomFforceGparity(Matrix *mom, Vector *chi, 
			      Float mass, Float epsilon, Float dt)
{
  char *fname = "EvolveMomFforce(M*,V*,F,F,F)";
  VRB.Func(cname,fname);

  Matrix *gauge = GaugeField() ;

  if (Colors() != 3) ERR.General(cname,fname,"Wrong nbr of colors.") ;
  if (SpinComponents() != 4) ERR.General(cname,fname,"Wrong nbr of spin comp.") ;
  if (mom == 0) ERR.Pointer(cname,fname,"mom") ;
  if (chi == 0) ERR.Pointer(cname,fname,"chi") ;
  if(!GJP.Gparity()) ERR.General(cname,fname,"G-parity boundary conditions are not active!");

//------------------------------------------------------------------
// allocate space for two CANONICAL fermion fields.
//------------------------------------------------------------------
  size_t f_size = FsiteSize() * GJP.VolNodeSites()*2;
  int f1_vectoff = f_size/2; //offset within fermion vector to f1 field

  char *str_v1 = "v1" ;
  Vector *v1 = (Vector *)smalloc(cname, fname, str_v1, f_size*sizeof(Float)) ;

  char *str_v2 = "v2" ;
  Vector *v2 = (Vector *)smalloc(cname, fname, str_v2, f_size*sizeof(Float)) ;

//------------------------------------------------------------------
// allocate space for two CANONICAL fermion field on a site.
// Each has 2 flavours
//------------------------------------------------------------------

  char *str_site_v1 = "site_v1";
  Float *site_v1 = (Float *)smalloc(cname, fname, str_site_v1, 2*FsiteSize()*sizeof(Float));

  char *str_site_v2 = "site_v2";
  Float *site_v2 = (Float *)smalloc(cname, fname, str_site_v2, 2*FsiteSize()*sizeof(Float));

  int f1_bufoff = FsiteSize(); //offset within buffer to f1 field

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

  Matrix tmp, tmp2, f ;

  for (mu=0; mu<4; mu++) {
    for (t=0; t<lt; t++)
    for (z=0; z<lz; z++)
    for (y=0; y<ly; y++)
    for (x=0; x<lx; x++) {
      int gauge_offset = x+lx*(y+ly*(z+lz*t)) ;
      int vec_offset = FsiteSize()*gauge_offset ;
      gauge_offset = mu+4*gauge_offset ;

      int pos[] = {x,y,z,t};
      int lattsz[] = {lx,ly,lz,lt}; 

      //Pointers for both flavours
      Float *v1_plus_mu[2] ;
      Float *v2_plus_mu[2] ;
      int vec_plus_mu_offset = FsiteSize() ;

      Float coeff[2] = {-2.0 * dt, -2.0 * dt}; 
      {
	int pos_p_mu[] = {x,y,z,t};
	pos_p_mu[mu] = (pos_p_mu[mu]+1) % lattsz[mu];
	vec_plus_mu_offset *= pos_p_mu[0] + lx*(pos_p_mu[1] + ly*(pos_p_mu[2] + lz*pos_p_mu[3]));
      }         

      if ((pos[mu]+1) == lattsz[mu]) {
	int buf_off[2] = {0,f1_bufoff};

	//If on global lattice boundary we must perform the G-parity flavor swap and multiply by the appropriate coefficient
	if(GJP.Bc(mu) == BND_CND_GPARITY && GJP.NodeCoor(mu)==GJP.Nodes(mu)-1){
	  //at global boundary
	  // d <- +CubarT,  CubarT <- -d
	  coeff[1] = -coeff[1];
	  buf_off[0] = f1_bufoff;
	  buf_off[1] = 0;
	}else if(GJP.NodeBc(mu)==BND_CND_APRD){
	  // d <- -d,   CubarT <- -CubarT
	  coeff[0] = -coeff[0];
	  coeff[1] = -coeff[1];
	}
	//Flavour 0
	getPlusData( (IFloat *)site_v1 + buf_off[0],
		     (IFloat *)v1+vec_plus_mu_offset, FsiteSize(), mu) ;
	getPlusData( (IFloat *)site_v2 + buf_off[0],
		     (IFloat *)v2+vec_plus_mu_offset, FsiteSize(), mu) ;
	v1_plus_mu[0] = site_v1 ;                        
	v2_plus_mu[0] = site_v2 ;                        

	//Flavour 1
	getPlusData( (IFloat *)site_v1 + buf_off[1],
		     (IFloat *)v1+vec_plus_mu_offset+f1_vectoff, FsiteSize(), mu) ;
	getPlusData( (IFloat *)site_v2 + buf_off[1],
		     (IFloat *)v2+vec_plus_mu_offset+f1_vectoff, FsiteSize(), mu) ;
	v1_plus_mu[1] = site_v1 + f1_bufoff;
	v2_plus_mu[1] = site_v2 + f1_bufoff;              

      } else {
	v1_plus_mu[0] = (Float *)v1+vec_plus_mu_offset ;
	v2_plus_mu[0] = (Float *)v2+vec_plus_mu_offset ;

	v1_plus_mu[1] = v1_plus_mu[0] + f1_vectoff;
	v2_plus_mu[1] = v2_plus_mu[0] + f1_vectoff;
      }

      //Contribution of f0 field is as standard
      sproj_tr[mu](   (IFloat *)&tmp,
                      (IFloat *)v1_plus_mu[0],
                      (IFloat *)v2+vec_offset, 1, 0, 0);

      sproj_tr[mu+4]( (IFloat *)&f,
                      (IFloat *)v2_plus_mu[0],
                      (IFloat *)v1+vec_offset, 1, 0, 0);

      tmp += f ;
      tmp *= coeff[0]; //different coefficients

      //For flavour 1 field we must flip the projection operator used for the force term
      //and swap the coordinates of the force vectors
      sproj_tr[mu+4](   (IFloat *)&tmp2,
                      (IFloat *)v1+vec_offset+f1_vectoff,
                      (IFloat *)v2_plus_mu[1], 1, 0, 0);

      sproj_tr[mu]( (IFloat *)&f,
		    (IFloat *)v2+vec_offset+f1_vectoff,
		    (IFloat *)v1_plus_mu[1], 1, 0, 0);      
      tmp2 += f ;
      f.Trans(tmp2);//must transpose outer product on colour index
      f *= coeff[1]; 

      tmp += f; //combined derivative part for both flavours

      //Multiply by gauge link and add to conj momentum
      f.DotMEqual(*(gauge+gauge_offset), tmp) ;
      tmp.Dagger(f) ;
      f.TrLessAntiHermMatrix(tmp) ;
      *(mom+gauge_offset) += f ;
      
      //Norm of the combined f0 and f1 force
      Float norm = f.norm();
      Float tmpf = sqrt(norm);
      L1 += tmpf;
      L2 += norm;
      Linf = (tmpf > Linf ? tmpf : Linf);

      //set force for U* link
      Matrix* mom_u = mom+gauge_offset;
      Matrix* mom_ustar = mom+gauge_offset+GJP.VolNodeSites()*4;
      mom_ustar->Conj((IFloat*)mom_u);
    }
  }
#endif
#endif

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
// chi = (M_f^\dag M_f)^{-1} M_f^\dag (RGV)
// phi = M_b (M_b^\dag M_b)^{-1} M_f^\dag (RGV)
//------------------------------------------------------------------
ForceArg FwilsonTm::EvolveMomFforce(Matrix *mom, Vector *chi, Vector *eta,
                                    Float mass, Float epsilon, Float dt)
{
#if 0
    const char *fname = "EvolveMomFforce(M*,V*,V*,F,F,F)";
    if (mom == 0)            { ERR.Pointer(cname,fname,"mom") ; }
    if (eta == 0)            { ERR.Pointer(cname,fname,"eta") ; }

    const size_t f_size = FsiteSize() * GJP.VolNodeSites();
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
#else
  if(GJP.Gparity()) return EvolveMomFforceGparity(mom,chi,eta,mass,epsilon,dt);

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

  size_t f_size        ( FsiteSize() * GJP.VolNodeSites() );

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
  
  Matrix *gparity_1f_mombuf;
  if(GJP.Gparity1fX()){
    gparity_1f_mombuf = (Matrix *)fmalloc(cname,fname,"gparity_1f_mombuf",4 * GJP.VolNodeSites() * sizeof(Matrix) ) ;
    for(int i=0;i<4*GJP.VolNodeSites();i++) gparity_1f_mombuf[i].ZeroMatrix();
  }

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
    // phi <- M_b (M_b^\dag M_b)^{-1} M_f^\dag (RGV) //CK: These appear to be wrong


    //Called from AlgActionQuotient:   chi = M_b (M_b^\dag M_b)^{-1} M_f^\dag (RGV)    eta = (M_f^\dag M_f)^{-1} M_f^\dag (RGV)
    wilson.CalcBsnForceVecs(chi, eta) ;

    //v1 = (chi, g5theta(ctheta,-stheta)D_eo chi)
    //v2 = ( -kappa^2 g5theta(ctheta,stheta)eta, -kappa^2 g5theta(ctheta,stheta) Deo^dag g5theta(ctheta,stheta)eta )
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

  int lattsz[4] = {lx,ly,lz,lt};
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

	    /* CK: Replaced nasty 300000 line switch statement with the following:*/
	    {
	      int pos_p_mu[] = {x,y,z,t};
	      pos_p_mu[mu] = (pos_p_mu[mu]+1) % lattsz[mu];
	      vec_plus_mu_offset *= pos_p_mu[0] + lx*(pos_p_mu[1] + ly*(pos_p_mu[2] + lz*pos_p_mu[3]));
	    }      
	    int pos[4] = {x,y,z,t};

	    if ((pos[mu]+1) == lattsz[mu]) {
	      getPlusData( (IFloat *)site_v1,
			   (IFloat *)v1+vec_plus_mu_offset, FsiteSize(), mu) ;
	      getPlusData( (IFloat *)site_v2,
			   (IFloat *)v2+vec_plus_mu_offset, FsiteSize(), mu) ;
	      v1_plus_mu = site_v1 ;                        
	      v2_plus_mu = site_v2 ;                        
	      if (GJP.NodeBc(mu)==BND_CND_APRD) coeff = -coeff ;
	    } else {
	      v1_plus_mu = (Float *)v1+vec_plus_mu_offset ;
	      v2_plus_mu = (Float *)v2+vec_plus_mu_offset ;
	    }

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

	    if(GJP.Gparity1fX()) *(gparity_1f_mombuf+gauge_offset) = tmp_mat2;

	    Float norm = tmp_mat2.norm();
	    Float tmp = sqrt(norm);
	    L1 += tmp;
	    L2 += norm;
	    Linf = (tmp>Linf ? tmp : Linf);
	 
	  }
  }


  //G-parity 1f add \delta p from other 'flavour'
  //Providing you set up the 1f prop source correctly (minus sign on UR quadrant) this code also works for
  //1f G-parity in both X and Y directions
  if(GJP.Gparity1fX()){
    int momsz = GsiteSize() * GJP.VolNodeSites();
    Matrix *buf2 = (Matrix *)fmalloc(cname,fname,"buf2",momsz * sizeof(Float) ) ;

    //Communicate \delta p from first half onto second half and vice versa
    Matrix *data_buf = gparity_1f_mombuf;
    Matrix *send_buf = data_buf;
    Matrix *recv_buf = buf2;

    int vol = GJP.VolNodeSites();
    int size[4] = {GJP.XnodeSites(),GJP.YnodeSites(),GJP.ZnodeSites(),GJP.TnodeSites()};

    if(GJP.Xnodes()>1){
      //pass between nodes
      for(int i=0;i<GJP.Xnodes()/2;i++){
	getMinusData((Float *)recv_buf, (Float *)send_buf, momsz , 0);
	data_buf = recv_buf;
	recv_buf = send_buf;
	send_buf = data_buf;
      }
    }else{
      //shift field by xsites/2
#ifdef USE_OMP
#pragma omp parallel for default(shared) private(mu)
      for(long i=0;i<vol*4;i++){
#else
      for(long i=0;i<vol*4;i++){
#endif
	  //i = mu + 4*(x + Lx*(y+Ly*(z+Lz*t) ) )
	  int mu = i%4;
	  int x = (i/4) % GJP.XnodeSites();
	  int pos_rem = i/4/GJP.XnodeSites(); //(y+Ly*(z+Lz*t)

	  int x_from = (x + GJP.XnodeSites()/2) % GJP.XnodeSites();
	  int i_from = mu + 4*(x_from + GJP.XnodeSites()*pos_rem);
	  buf2[i] = gparity_1f_mombuf[i_from];
	}
	data_buf = buf2;
      }

#ifdef USE_OMP
#pragma omp parallel for default(shared) private(mu)
      for (long i=0;i<vol*4;i++){
  	int pos[4];
	long rest=i;
	mu = rest%4; rest = rest/4;
	for(int j =0; j<4;j++){
	  pos[j]= rest%size[j]; rest = rest/size[j];
	}
#else
	int pos[4];
	for (mu=0; mu<4; mu++)
	  for (pos[3]=0; pos[3]<size[3]; pos[3]++)
	    for (pos[2]=0; pos[2]<size[2]; pos[2]++)
	      for (pos[1]=0; pos[1]<size[1]; pos[1]++)
		for (pos[0]=0; pos[0]<size[0]; pos[0]++){
#endif
		  int gauge_offset = mu + 4*(pos[0]+size[0]*(pos[1]+size[1]*(pos[2]+size[2]*pos[3])));
      
		  //complex conjugate the \delta p from the other flavour
		  Float *m = (Float*) &data_buf[gauge_offset];
		  for(int c=1;c<18;c+=2) m[c] *= -1;
      
		  //add it to the momentum at this site
		  //minus sign because it is a boson force
		  *(mom+gauge_offset) -= data_buf[gauge_offset];
		}
   ffree(buf2);
 }

    
 if(GJP.Gparity1fX()) ffree(gparity_1f_mombuf,cname,fname,"gparity_1f_mombuf");



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
#endif
  return ForceArg(L1, sqrt(L2), Linf);
}


ForceArg FwilsonTm::EvolveMomFforceGparity(Matrix *mom, Vector *chi, Vector *eta,
		      Float mass, Float epsilon, Float dt) {
  char *fname = "EvolveMomFforce(M*,V*,V*,F,F,F)";
  VRB.Func(cname,fname);

  Matrix *gauge = GaugeField() ;

  if (Colors() != 3)       { ERR.General(cname,fname,"Wrong nbr of colors.") ; }
  if (SpinComponents()!=4) { ERR.General(cname,fname,"Wrong nbr of spin comp.") ;}
  if (mom == 0)            { ERR.Pointer(cname,fname,"mom") ; }
  if (eta == 0)            { ERR.Pointer(cname,fname,"eta") ; }
  if (!GJP.Gparity())       { ERR.General(cname,fname,"G-parity boundary conditions not active!") ; }

//------------------------------------------------------------------
// allocate space for two CANONICAL fermion fields.
// these are all full fermion vector sizes ( i.e. *not* preconditioned )
//------------------------------------------------------------------

  size_t f_size( FsiteSize() * GJP.VolNodeSites()*2);
  int f1_vectoff = f_size/2; //offset within fermion vector to f1 field

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
  Float *site_v1 = (Float *)smalloc(2*FsiteSize()*sizeof(Float));
  if (site_v1 == 0) ERR.Pointer(cname, fname, str_site_v1) ;
  VRB.Smalloc(cname, fname, str_site_v1, site_v1,
    FsiteSize()*sizeof(Float)) ;

  char *str_site_v2 = "site_v2";
  Float *site_v2 = (Float *)smalloc(2*FsiteSize()*sizeof(Float));
  if (site_v2 == 0) ERR.Pointer(cname, fname, str_site_v2) ;
  VRB.Smalloc(cname, fname, str_site_v2, site_v2,
    FsiteSize()*sizeof(Float)) ;
  
  int f1_bufoff = FsiteSize(); //offset within buffer to f1 field

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

  int lattsz[4] = {lx,ly,lz,lt};

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

	    //Pointers for both flavours
	    Float *v1_plus_mu[2] ;
	    Float *v2_plus_mu[2] ;
	    int vec_plus_mu_offset = FsiteSize() ;
	    
	    Float coeff[2] = {-2.0 * dt, -2.0 * dt}; 

	    /* CK: Replaced nasty 300000 line switch statement with the following:*/
	    {
	      int pos_p_mu[] = {x,y,z,t};
	      pos_p_mu[mu] = (pos_p_mu[mu]+1) % lattsz[mu];
	      vec_plus_mu_offset *= pos_p_mu[0] + lx*(pos_p_mu[1] + ly*(pos_p_mu[2] + lz*pos_p_mu[3]));
	    }      
	    int pos[4] = {x,y,z,t};

	    if ((pos[mu]+1) == lattsz[mu]) {
	      int buf_off[2] = {0,f1_bufoff};

	      //If on global lattice boundary we must perform the G-parity flavor swap and multiply by the appropriate coefficient
	      if(GJP.Bc(mu) == BND_CND_GPARITY && GJP.NodeCoor(mu)==GJP.Nodes(mu)-1){
		//at global boundary
		// d <- +CubarT,  CubarT <- -d
		coeff[1] = -coeff[1];
		buf_off[0] = f1_bufoff;
		buf_off[1] = 0;
	      }else if(GJP.NodeBc(mu)==BND_CND_APRD){
		// d <- -d,   CubarT <- -CubarT
		coeff[0] = -coeff[0];
		coeff[1] = -coeff[1];
	      }

	      //Flavour 0
	      getPlusData( (IFloat *)site_v1 + buf_off[0],
			   (IFloat *)v1+vec_plus_mu_offset, FsiteSize(), mu) ;
	      getPlusData( (IFloat *)site_v2 + buf_off[0],
			   (IFloat *)v2+vec_plus_mu_offset, FsiteSize(), mu) ;
	      v1_plus_mu[0] = site_v1 ;                        
	      v2_plus_mu[0] = site_v2 ;                        

	      //Flavour 1
	      getPlusData( (IFloat *)site_v1 + buf_off[1],
			   (IFloat *)v1+vec_plus_mu_offset+f1_vectoff, FsiteSize(), mu) ;
	      getPlusData( (IFloat *)site_v2 + buf_off[1],
			   (IFloat *)v2+vec_plus_mu_offset+f1_vectoff, FsiteSize(), mu) ;
	      v1_plus_mu[1] = site_v1 + f1_bufoff;
	      v2_plus_mu[1] = site_v2 + f1_bufoff;       

	    } else {
	      v1_plus_mu[0] = (Float *)v1+vec_plus_mu_offset ;
	      v2_plus_mu[0] = (Float *)v2+vec_plus_mu_offset ;

	      v1_plus_mu[1] = v1_plus_mu[0] + f1_vectoff;
	      v2_plus_mu[1] = v2_plus_mu[0] + f1_vectoff;
	    }

	    Matrix tmp_mat1, tmp_mat2, tmp_mat3;  

	    // ( 1 - gamma_\mu ) Tr_s [ v1(x+\mu) v2^{\dagger}(x) ]           
	    sproj_tr[mu](   (IFloat *)&tmp_mat1,   	// output color matrix
			    (IFloat *)v1_plus_mu[0],		// row vector, NOT conjugated
			    (IFloat *)v2+vec_offset, 	// col vector, IS conjugated
			    1, 0, 0);				// 1 block, 0 strides

	    // (1 + gamma_\mu)  Tr_s [ v2(x+\mu) v1^{\dagger}(x) ]
	    sproj_tr[mu+4]( (IFloat *)&tmp_mat2,		// output color matrix
			    (IFloat *)v2_plus_mu[0],		// row vector, NOT conjugated
			    (IFloat *)v1+vec_offset, 	// col vector, IS conjugated
			    1, 0, 0);				// 1 block, 0 strides

	    tmp_mat1 += tmp_mat2 ;
	    tmp_mat1 *= coeff[0]; //different coefficients

	    //For flavour 1 field we must flip the projection operator used for the force term
	    //and swap the coordinates of the force vectors
	    sproj_tr[mu+4](   (IFloat *)&tmp_mat2,
			      (IFloat *)v1+vec_offset+f1_vectoff,
			      (IFloat *)v2_plus_mu[1], 1, 0, 0);

	    sproj_tr[mu]( (IFloat *)&tmp_mat3,
			  (IFloat *)v2+vec_offset+f1_vectoff,
			  (IFloat *)v1_plus_mu[1], 1, 0, 0);      
	    tmp_mat2 += tmp_mat3 ;
	    tmp_mat3.Trans(tmp_mat2);//must transpose outer product on colour index
	    tmp_mat3 *= coeff[1]; 

	    tmp_mat1 += tmp_mat3; //combined derivative part for both flavours

	    // multiply sum by the link in the \mu direction
	    tmp_mat2.DotMEqual(*(gauge+gauge_offset), tmp_mat1) ;
      
	    // take traceless antihermitian piece
	    // TrLessAntiHermMatrix need to be passed
	    // the dagger of the matrix in question
	    tmp_mat1.Dagger(tmp_mat2) ;
	    tmp_mat2.TrLessAntiHermMatrix(tmp_mat1) ;
            
	    //~~
	    //~~ fermion version:  	(mom+gauge_offset) += f
	    //~~ boson version:  	(mom+gauge_offset) -= f
	    //~~
	    *(mom+gauge_offset) -= tmp_mat2 ;

	    //Norm of the combined f0 and f1 force
	    Float norm = tmp_mat2.norm();
	    Float tmp = sqrt(norm);
	    L1 += tmp;
	    L2 += norm;
	    Linf = (tmp>Linf ? tmp : Linf);
	 
	    //set force for U* link
	    Matrix* mom_u = mom+gauge_offset;
	    Matrix* mom_ustar = mom+gauge_offset+GJP.VolNodeSites()*4;
	    mom_ustar->Conj((IFloat*)mom_u);
	  }
  }
#endif

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
