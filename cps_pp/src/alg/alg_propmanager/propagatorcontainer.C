#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <util/qcdio.h>
#ifdef PARALLEL
#include <comms/sysfunc_cps.h>
#endif
#include <comms/scu.h>
#include <comms/glb.h>

#include <util/lattice.h>
#include <util/time_cps.h>
#include <alg/do_arg.h>
#include <alg/no_arg.h>
#include <alg/common_arg.h>
#include <util/smalloc.h>

#include <util/ReadLatticePar.h>
#include <util/WriteLatticePar.h>

#include<unistd.h>
#include<config.h>

#include <alg/qpropw.h>
#include <alg/qpropw_arg.h>

#include<alg/propagatorcontainer.h>
#include<alg/propmanager.h>

#ifdef USE_BFM 

//CK: these are redefined by BFM (to the same values)
#undef ND
#undef SPINOR_SIZE
#undef HALF_SPINOR_SIZE
#undef GAUGE_SIZE
#endif

#ifdef USE_BFM
#include <util/lattice/bfm_evo.h>
#include <util/lattice/bfm_eigcg.h>
#include <util/lattice/fbfm.h>
#include <alg/eigen/Krylov_5d.h>
#endif

CPS_START_NAMESPACE

bool PropagatorContainer::fbfm_assume_bc_applied(true);

PropagatorContainer::PropagatorContainer(){ for(int i=0;i<50;i++) attributes[i]=NULL;}


void PropagatorContainer::setup(PropagatorArg &arg){
  for(int i=0;i<50;i++) if(attributes[i]) delete attributes[i];
  
  //special case, generic attributes are always present and maintained outside of the array
  AttributeContainer* p = new AttributeContainer;
  p->type = GENERIC_PROP_ATTR;
  p->AttributeContainer_u.generic_prop_attr.deep_copy(arg.generics);
  attributes[ (int)p->type ] = p;

  if(UniqueID()==0) printf("Got %d additional attributes\n",arg.attributes.attributes_len);

  for(int i=0;i<arg.attributes.attributes_len;i++){
    //no duplicates, later entries of same type overwrite earlier
    int idx = (int)arg.attributes.attributes_val[i].type;
    if(attributes[idx]!=NULL) delete attributes[idx];
    attributes[idx] = new AttributeContainer;
    attributes[idx]->deep_copy(arg.attributes.attributes_val[i]); //make a copy
  }
}

void PropagatorContainer::add(const AttributeContainer &p){
  int idx = (int)p.type;
  if(attributes[idx]!=NULL) delete attributes[idx]; //no duplicates
  attributes[idx] = new AttributeContainer;
  attributes[idx]->deep_copy(p);
}

AttributeContainer*  PropagatorContainer::findAttr(const AttrType &type) const{
  return attributes[(int)type];
}
PropagatorContainer::~PropagatorContainer(){
  for(int i=0;i<50;i++) if(attributes[i]) delete attributes[i];
}

bool PropagatorContainer::tagEquals(const char* what){
  GenericPropAttrArg *generics;
  if(!getAttr(generics)) ERR.General("PropagatorContainer","tagEquals(const char* what)","Propagator attribute list does not contain a GenericPropAttr\n");
  if(strcmp(generics->tag,what)==0) return true;
  return false;
}
char const* PropagatorContainer::tag() const{
  GenericPropAttrArg *generics;
  if(!getAttr(generics)) ERR.General("PropagatorContainer","tag()","Propagator attribute list does not contain a GenericPropAttr\n");
  return generics->tag;
}

void PropagatorContainer::printAttribs() const{
  printf("Propagator Attributes:\n");
  for(int i=0;i<50;i++){
    if(attributes[i]!=NULL){
      printf("%s:",AttrType_map[i].name);
      attributes[i]->print();
    }
  }
}
PropagatorType PropagatorContainer::type() const{
  if(hasAttr<GenericPropAttrArg>()) return getAttr<GenericPropAttrArg>()->type;
  else ERR.General("PropagatorContainer","type()","Propagator attribute list does not contain a GenericPropAttr\n");
}

PropagatorContainer* PropagatorContainer::create(const PropagatorType &ctype){
  if(ctype == QPROPW_TYPE){
    return new QPropWcontainer;
  }else if(ctype == A2A_PROP_TYPE){
    //return new A2ApropContainer;
  }else{
    ERR.General("PropagatorContainer","create(...)","Unknown type\n");
  }
}

void QPropWcontainer::readProp(Lattice &latt){
  if(prop!=NULL) return; //don't load if already inverted
  PropIOAttrArg *io;
  if(!getAttr(io)) return;
  if(!io->prop_on_disk) return;

  if(UniqueID()==0) printf("Prop is on disk, loading from file %s\n",io->qio_filename);
  
  //load prop with info from io
  CommonArg c_arg("label","filename");//find out what this does!
  
  prop = new QPropW(latt,&c_arg);
  prop->ReLoad(io->qio_filename);
}

void QPropWcontainer::calcProp(Lattice &latt){
  //function acts as factory for QPropW objects depending on the attribute objects
  if(prop!=NULL) return; //don't calculate twice

  const char *cname = "QPropWcontainer";
  const char *fname = "calcProp()";

  CommonArg c_arg("label","filename");//find out what this does!

  GenericPropAttrArg *generics;
  if(!getAttr(generics)) ERR.General(cname,fname,"Propagator attribute list does not contain a GenericPropAttr\n");

  //if propagator is a combination of other propagators, do this separately
  PropCombinationAttrArg *propcomb;
  if(getAttr(propcomb)){
    if(UniqueID()==0) printf("Propagator %s combining props %s and %s\n",generics->tag,propcomb->prop_A,propcomb->prop_B);

    PropagatorContainer &A = PropManager::getProp(propcomb->prop_A);
    PropagatorContainer &B = PropManager::getProp(propcomb->prop_B);

    if(A.type()!=QPROPW_TYPE || B.type()!=QPROPW_TYPE) ERR.General("QPropWcontainer","calcProp(Lattice &latt)","When combining propagators, either/both \"%s\" and \"%s\" are not a QPropWcontainer\n",propcomb->prop_A,propcomb->prop_B);

    //copy attributes from A. Does not check that attributes of A and B match
    propCombSetupAttrib();

    //calculate A and B if they have not yet been calculated
    QPropW &A_qpw = A.convert<QPropWcontainer>().getProp(latt);
    QPropW &B_qpw = B.convert<QPropWcontainer>().getProp(latt);

    prop = new QPropW(A_qpw);
    //perform the combination
    if(UniqueID()==0) printf("Propagator %s starting combination: prop %p, A_qpw %p B_qpw %p\n",generics->tag,prop,&A_qpw,&B_qpw);
    
    if(propcomb->combination == A_PLUS_B){
      prop->Average(B_qpw);
    }else if(propcomb->combination == A_MINUS_B){
      prop->LinComb(B_qpw,0.5,-0.5);
    }else{
      ERR.General(cname,fname,"Unknown PropagatorCombination");
    }
    if(UniqueID()==0) printf("Propagator %s finished combination\n",generics->tag);
    return;
  }

  //calculate the propagator
  if(UniqueID()==0) printf("Calculating propagator %s\n",generics->tag);

  CgArg cg;
  cg.mass =  generics->mass;
  cg.max_num_iter = 5000;
  cg.stop_rsd =   1.0000000000000000e-08;
  cg.true_rsd =   1.0000000000000000e-08;
  cg.RitzMatOper = NONE;
  cg.Inverter = CG;
  cg.bicgstab_n = 0;

  CGAttrArg *cgattr;
  if(getAttr(cgattr)){
    cg.max_num_iter = cgattr->max_num_iter;
    cg.stop_rsd = cgattr->stop_rsd;
    cg.true_rsd = cgattr->true_rsd; 
  }

  QPropWArg qpropw_arg;
  qpropw_arg.cg = cg;
  qpropw_arg.x = 0;
  qpropw_arg.y = 0;
  qpropw_arg.z = 0;
  qpropw_arg.t = 0;
  qpropw_arg.flavor = 0; //default on d field  
  qpropw_arg.ensemble_label = "ens";
  qpropw_arg.ensemble_id = "ens_id";
  qpropw_arg.StartSrcSpin = 0;
  qpropw_arg.EndSrcSpin = 4;
  qpropw_arg.StartSrcColor = 0;
  qpropw_arg.EndSrcColor = 3;
  qpropw_arg.gauge_fix_src = 0;
  qpropw_arg.gauge_fix_snk = 0;


  //Set the boundary conditions
  //For regular lattice classes it is sufficient to modify the boundary condition in GJP, as the boundary condition
  //is applied when the DiracOp class is instantiated within the inverter, and the lattice is restored afterwards.
  //However for Fbfm we must manually apply the BC and re-import the gauge field into the internal bfm objects
  bool is_fbfm = ( latt.Fclass() == F_CLASS_BFM );
  if(is_fbfm && PropagatorContainer::fbfm_assume_bc_applied)
    latt.BondCond(); //un-applies the existing BCs, reverting to periodic BCs and imports the gauge field into the bfm instances

  BndCndType init_bc[4];
  TwistedBcAttrArg *tbcarg;

  for(int i=0;i<4;i++){
    if(i<3 && generics->bc[i] != GJP.Bc(i) && !(generics->bc[i] == BND_CND_TWISTED || generics->bc[i] == BND_CND_GPARITY_TWISTED) )
      ERR.General(cname,fname,"Propagator %s: valence and sea spatial boundary conditions do not match (partially-twisted BCs are allowed)\n",generics->tag);

    if( GJP.Bc(i) != BND_CND_GPARITY && generics->bc[i] == BND_CND_GPARITY_TWISTED ) ERR.General(cname,fname,"Propagator %s: Cannot use twisted G-parity valence BCs in a non-Gparity direction");

    if(generics->bc[i] == BND_CND_TWISTED || generics->bc[i] == BND_CND_GPARITY_TWISTED){
      if(getAttr(tbcarg)) for(int j=0;j<3;j++) GJP.TwistAngle(j, tbcarg->theta[j]);
      else for(int j=0;j<3;j++) GJP.TwistAngle(j,0.0); //default twist angle is zero
    }
    init_bc[i] = GJP.Bc(i);
    GJP.Bc(i,generics->bc[i]);
  }
  if(is_fbfm) latt.BondCond(); //applies the new BCs and imports the gauge field into the bfm instances

  //fill out qpropw_arg arguments
  GparityFlavorAttrArg *flav;
  if(getAttr(flav)) qpropw_arg.flavor = flav->flavor;

  //mid-point correlator?
  if(hasAttr<StoreMidpropAttrArg>()) qpropw_arg.store_midprop = 1;

  PropIOAttrArg *io;
  if(getAttr(io)){ 
    qpropw_arg.save_prop = io->save_to_disk;
    #ifndef USE_QIO
    if(io->save_to_disk) ERR.General(cname,fname,"Cannot save propagators without QIO\n");
    #endif

    qpropw_arg.file = io->qio_filename;
  }

  GaugeFixAttrArg *gfix;
  if(getAttr(gfix)){
    if(latt.FixGaugeKind() == FIX_GAUGE_NONE && (gfix->gauge_fix_src==1 || gfix->gauge_fix_snk==1))
      ERR.General(cname,fname,"Gauge fixed source or sink requested but gauge fixing matrices have not been calculated\n");
    
    qpropw_arg.gauge_fix_src = gfix->gauge_fix_src; 
    qpropw_arg.gauge_fix_snk = gfix->gauge_fix_snk;
  }

  //Deal with deflated CG
  DeflatedCGAttrArg *defl_cg;
  if(getAttr(defl_cg)){
    if(!is_fbfm) ERR.General(cname,fname,"Currently, deflated CG only available via Fbfm\n");
#ifdef USE_BFM
    LanczosContainer &lanczos = PropManager::getLanczos(defl_cg->lanczos_tag);
    const LanczosContainerArg &lanc_args = lanczos.getArgs();
    int N_use = lanc_args.lanc_arg.N_true_get;

    if(Fbfm::use_mixed_solver){
      //Pass single-precision eigenvectors to Fbfm 
      //(note, this will fail if you haven't manually changed the precision of the eigenvectors in the LanczosContainer. I didn't want to do this
      // automatically as the change of precision introduces numerical errors. Better that it fails if you didn't explicitly want this!)
      BFM_Krylov::Lanczos_5d<float> &eigs = lanczos.getEigSinglePrec(latt);
      dynamic_cast<Fbfm&>(latt).set_deflation(&eigs.bq,&eigs.bl,N_use);
    }else{
      //Use double-precision eigenvectors
      lanczos.setPrecision(2);
      BFM_Krylov::Lanczos_5d<double> &eigs = lanczos.getEig(latt);
      dynamic_cast<Fbfm&>(latt).set_deflation(&eigs.bq,&eigs.bl,N_use);
    }
#else
    ERR.General(cname,fname,"Currently, deflated CG only available with Bfm/Fbfm\n");
#endif
  }

  PointSourceAttrArg *pt;
  WallSourceAttrArg *wl;
  VolumeSourceAttrArg *vl;
  if(getAttr(pt) && getAttr(wl) || getAttr(pt) && getAttr(vl) || getAttr(vl) && getAttr(wl) )
    ERR.General(cname,fname,"Propagator %s: Must specify only one source type attribute\n",generics->tag);

  //NOTE: Momentum units are:
  //                         Periodic  2\pi/L
  //                     Antiperiodic   \pi/L
  //                         G-parity   \pi/2L
  //In G-parity directions, the propagators are antiperiodic in 2L

  if(getAttr(pt)){
    if(!UniqueID()){ printf("Doing a point source propagator\n"); }
    qpropw_arg.x = pt->pos[0];
    qpropw_arg.y = pt->pos[1];
    qpropw_arg.z = pt->pos[2];
    qpropw_arg.t = pt->pos[3];
    
    //assemble the arg objects
    prop = new QPropWPointSrc(latt,&qpropw_arg,&c_arg);

    MomentumAttrArg *mom;
    if(getAttr(mom)){
      if(!UniqueID()){ printf("Adding momentum to point source\n"); }
      //for a point source with momentum, apply the e^-ipx factor at the source location 
      //such that it can be treated in the same way as a momentum source

      int local_site[4];
      bool on_node = true;
      for(int i=0;i<4;i++){ 
	local_site[i] = pt->pos[i] - GJP.NodeSites(i)*GJP.NodeCoor(i);
	if(local_site[i] < 0 || local_site[i] >= GJP.NodeSites(i)) on_node=false;
      }
      if(on_node){
	ThreeMom tmom(mom->p);
	Site s(local_site[0],local_site[1],local_site[2],local_site[3]);
	Complex phase;

	MomCosAttrArg *cos;
	if(getAttr(cos)) //a cosine point source
	  phase = tmom.FactCos(s);
	else
	  phase = tmom.Fact(s);
	
	int wmat_shift = GJP.Gparity() ? qpropw_arg.flavor * GJP.VolNodeSites() : 0;
	(*prop)[s.Index()+wmat_shift] *= phase; 
	printf("Propagator %s has 3-momentum (%d,%d,%d) and position (%d,%d,%d,%d), source phase factor is (%e,%e)\n",
			     generics->tag,mom->p[0],mom->p[1],mom->p[2],pt->pos[0],pt->pos[1],pt->pos[2],pt->pos[3],phase.real(),phase.imag());
      }
    }
  }else if(getAttr(wl)){
    qpropw_arg.t = wl->t;
    
    if(qpropw_arg.gauge_fix_src == 1 && latt.FixGaugeKind() != FIX_GAUGE_COULOMB_T)
      ERR.General(cname,fname,"Gauge fixed wall/mom source requested, but gauge is not Coulomb-T\n");

    MomentumAttrArg *mom;
    if(getAttr(mom)){
      MomCosAttrArg *cos;
      if(getAttr(cos)){ //a cosine source
	prop = new QPropWMomCosSrc(latt,&qpropw_arg,mom->p,&c_arg); //knows about correct G-parity units of momentum (cf. ThreeMom::CalcLatMom)
      }else{ //a regular momentum source
	prop = new QPropWMomSrc(latt,&qpropw_arg,mom->p,&c_arg); //knows about correct G-parity units of momentum
      }
    }else{
      prop = new QPropWWallSrc(latt,&qpropw_arg,&c_arg);
    }
  }else if(getAttr(vl)){
    if(qpropw_arg.gauge_fix_src == 1 && latt.FixGaugeKind() == FIX_GAUGE_NONE)
      ERR.General(cname,fname,"Gauge fixed volume source requested, but no gauge fixing has been performed\n");
    
    MomentumAttrArg *mom;
    if(getAttr(mom)){
      MomCosAttrArg *cos;
      if(getAttr(cos)){
	ERR.General(cname,fname,"No volume cosine source implemented");
      }else{
	prop = new QPropWVolMomSrc(latt,&qpropw_arg,mom->p,&c_arg);
      }
    }else{
      prop = new QPropWVolSrc(latt,&qpropw_arg,&c_arg);
    }
  }else{
    ERR.General(cname,fname,"Propagator %s has no source type AttrArg\n",generics->tag);
  }
  
  if(getAttr(io) && io->save_to_disk){
    io->prop_on_disk = true; //QPropW saves the prop
  }
  
  if(is_fbfm)    
    latt.BondCond(); //un-apply the new BCs and reimport the gauge field into the bfm instances
  

  //restore boundary conditions to original
  for(int i=0;i<4;i++) GJP.Bc(i,init_bc[i]);
  for(int j=0;j<3;j++) GJP.TwistAngle(j,0);

  if(is_fbfm && PropagatorContainer::fbfm_assume_bc_applied) latt.BondCond(); //restore original BCs
}


int QPropWcontainer::flavor() const{
  GparityFlavorAttrArg *flav;
  if(getAttr(flav)) return flav->flavor;
  return 0;
}
void QPropWcontainer::momentum(int *into) const{
    MomentumAttrArg *mom;
    if(getAttr(mom)){
      for(int i=0;i<3;i++) into[i] = mom->p[i];
    }else{
      for(int i=0;i<3;i++) into[i] = 0;
    }
}

void QPropWcontainer::setProp(QPropW *to, bool _own){
  if(own && prop != NULL) delete prop;
  prop = to;
  own = _own;
}

QPropW & QPropWcontainer::getProp(Lattice &latt){
  if(prop==NULL){
    readProp(latt);
    calcProp(latt); //will calculate if prop was not read
  }
  return *prop;
}


QPropWcontainer & QPropWcontainer::verify_convert(PropagatorContainer &pc, const char* cname, const char* fname){
  if(pc.type()!=QPROPW_TYPE) ERR.General(cname,fname,"Expect propagator \"%s\" to be QPropW type\n",pc.tag());
  return pc.convert<QPropWcontainer>();
}
const QPropWcontainer & QPropWcontainer::verify_convert(const PropagatorContainer &pc, const char* cname, const char* fname){
  if(pc.type()!=QPROPW_TYPE) ERR.General(cname,fname,"Expect propagator \"%s\" to be QPropW type\n",pc.tag());
  return pc.convert<QPropWcontainer>();
}


void QPropWcontainer::deleteProp(){
  if(!own) ERR.General("QPropWcontainer","deleteProp","Cannot delete prop that is not owned by this container\n");
  if(prop!=NULL){ delete prop; prop=NULL; }
}



void QPropWcontainer::propCombSetupAttrib(){
  PropCombinationAttrArg *propcomb;
  if(getAttr(propcomb)){
    PropagatorContainer &A = PropManager::getProp(propcomb->prop_A);
    //copy attributes from A. Does not check that attributes of A and B match
    for(int i=0;i<50;i++) 
      if((AttrType)i != GENERIC_PROP_ATTR && attributes[i]==NULL && A.findAttr( (AttrType)i )!=NULL) 
  	add(*A.findAttr( (AttrType)i ) );
  }
}


std::vector<std::vector<int> > QPropWcontainer::get_allowed_momenta() const{
  std::vector<std::vector<int> > out;

  MomCosAttrArg *cos;
  MomentumAttrArg *mom;

  if(getAttr(mom)){
    std::vector<int> p(3); p[0] = mom->p[0]; p[1] = mom->p[1]; p[2] = mom->p[2];
    out.push_back(p);

    if(getAttr(cos)){
      std::vector<int> mp(p); for(int i=0;i<3;i++) mp[i]*=-1;
      out.push_back(mp);
    }
  }else{
    out.push_back(std::vector<int>(3,0));
  }
  return out;
}


int QPropWcontainer::getSourceTimeslice(){
  PointSourceAttrArg *pt;
  if(getAttr(pt)) return pt->pos[3];
  WallSourceAttrArg *wl;
  if(getAttr(wl)) return wl->t;

  ERR.General("QPropWcontainer","getSourceTimeslice with prop %s: Could not determine source timeslice",tag());
}

void QPropWcontainer::setSourceTimeslice(const int &t){
  PointSourceAttrArg *pt;
  if(getAttr(pt)){
    pt->pos[3] = t;
    return;
  }
  WallSourceAttrArg *wl;
  if(getAttr(wl)){
    wl->t = t;
    return;
  }
  ERR.General("QPropWcontainer","setSourceTimeslice with prop %s: Could not determine source timeslice",tag());
}







// A2ApropContainer::~A2ApropContainer(){ if(prop!=NULL) delete prop; }


// void A2ApropContainer::readProp(Lattice &latt){
//   //Is there a way to read back in the eigenvectors?

//   return;
// }

// void A2ApropContainer::calcProp(Lattice &latt){
//   if(prop!=NULL) return; //don't calculate twice

//   const char *cname = "A2ApropContainer";
//   const char *fname = "calcProp()";

//   GenericPropAttrArg *generics;
//   if(!getAttr(generics)) ERR.General(cname,fname,"Propagator attribute list does not contain a GenericPropAttr\n");

//   A2AAttrArg *a2a_arg;
//   if(!getAttr(a2a_arg)) ERR.General(cname,fname,"Propagator attribute list does not contain a A2AAttr\n");

//   //Get eigenvectors
//   LanczosContainer &lanczos = PropManager::getLanczos(a2a_arg->lanczos_tag);
//   Lanczos_5d<double> &eig = lanczos.getEig(latt);
//   if(a2a_arg->nl !=0 && eig.dop.mass != generics->mass) ERR.General(cname,fname,"Eigenvalues provided were calculated with a different mass than that specified in propagator arguments\n");

//   CommonArg c_arg("label","/dev/null");//find out what this does!

//   prop = new A2APropbfm(latt,*a2a_arg,c_arg,&eig);

//   bfm_evo<double> &dwf = eig.dop;

//   //If the Lanczos mass is different from the prop mass (only allowed if we are not using any low-modes from the Lanczos)
//   //then we must change the mass used in the Dirac operator
//   double restore_mass; bool do_restore(false);
//   if(dwf.mass != generics->mass){
//     restore_mass = dwf.mass;
//     dwf.mass = generics->mass;
//     dwf.GeneralisedFiveDimEnd();
//     dwf.GeneralisedFiveDimInit();
//     do_restore=true;
//   }

//   //Calculate the vectors and FFT vectors
//   prop->allocate_vw();
//   if(!UniqueID()) printf("%s::%s computing A2A low modes component for prop %s\n",cname,fname,generics->tag); 
//   prop->compute_vw_low(dwf);
//   if(!UniqueID()) printf("%s::%s computing A2A high modes component for prop %s\n",cname,fname,generics->tag); 
//   prop->compute_vw_high(dwf);
//   if(!UniqueID()) printf("%s::%s computing FFT of V and W for A2A prop %s\n",cname,fname,generics->tag); 
//   prop->fft_vw();
//   if(!UniqueID()) printf("%s::%s Finished computing A2A prop %s\n",cname,fname,generics->tag); 

//   if(do_restore){
//     dwf.mass = restore_mass;
//     dwf.GeneralisedFiveDimEnd();
//     dwf.GeneralisedFiveDimInit();
//   }
// }


// A2APropbfm & A2ApropContainer::getProp(Lattice &latt){
//   if(prop==NULL){
//     readProp(latt);
//     calcProp(latt); //will calculate if prop was not read
//   }
//   return *prop;
// }


// A2ApropContainer & A2ApropContainer::verify_convert(PropagatorContainer &pc, const char* cname, const char* fname){
//   if(pc.type()!=A2A_PROP_TYPE) ERR.General(cname,fname,"Expect propagator \"%s\" to be A2A type\n",pc.tag());
//   return pc.convert<A2ApropContainer>();
// }
// const A2ApropContainer & A2ApropContainer::verify_convert(const PropagatorContainer &pc, const char* cname, const char* fname){
//   if(pc.type()!=A2A_PROP_TYPE) ERR.General(cname,fname,"Expect propagator \"%s\" to be A2A type\n",pc.tag());
//   return pc.convert<A2ApropContainer>();
// }


// void A2ApropContainer::deleteProp(){
//   if(prop!=NULL){ delete prop; prop=NULL; }
// }



#ifdef USE_BFM

void LanczosContainer::deleteEig(){
  if(lanczos!=NULL){ delete lanczos; lanczos=NULL;}
  if(lanczos_f!=NULL){ delete lanczos_f; lanczos_f=NULL;}
}

void LanczosContainer::setupBfm(const int &prec){
  //Setup bfmarg
  bfmarg dwfa;
  
  for(int d=0;d<4;d++){
#ifdef BFM_GPARITY
    dwfa.nodes[d] = GJP.Nodes(d);
#endif
    dwfa.ncoor[d] = GJP.NodeCoor(d);
    dwfa.node_latt[d] = GJP.NodeSites(d);
  }

#ifdef BFM_GPARITY
  if(GJP.Gparity()){
    dwfa.gparity = 1;
    for(int d=0;d<3;d++) dwfa.gparity_dir[d] = (GJP.Bc(d) == BND_CND_GPARITY ? 1 : 0);
  }
#endif
  dwfa.verbose=1;
  dwfa.reproduce=0;

  for(int mu=0;mu<4;mu++)
    if ( GJP.Nodes(mu)>1 ) dwfa.local_comm[mu] = 0;
    else dwfa.local_comm[mu] = 1;

  dwfa.Ls   = GJP.SnodeSites();
  dwfa.M5   = toDouble(GJP.DwfHeight());

  //Now for user specified arguments
  dwfa.mass = args.lanc_arg.mass;
  dwfa.max_iter = args.cg_max_iter;
  dwfa.residual = args.cg_residual;
  dwfa.precon_5d = args.cg_precon_5d;
  dwfa.mobius_scale = args.mobius_scale;

  if(args.solver == BFM_DWF) dwfa.solver = DWF;
  else if(args.solver == BFM_HmCayleyTanh) dwfa.solver = HmCayleyTanh;
  else ERR.General("LanczosContainer","setupBfm","Enum for chosen solver not in factory, add it!\n");

  //Create and initialise bfm object
  if(prec == 2){
    dwf = new bfm_evo<double>();
    dwf->init(dwfa);
  }else if(prec == 1){
    dwf_f = new bfm_evo<float>();
    dwf_f->init(dwfa);
  }else ERR.General("LanczosContainer","setupBfm","Invalid precision, %d\n",prec);
}

//Calculate in double precision always
void LanczosContainer::calcEig(Lattice &latt){
  if(lanczos!=NULL) return; //no need to calculate twice
  
  if(lanczos_f != NULL) delete lanczos_f;

  if(dwf==NULL){
    setupBfm(2);
    reload_gauge=true;
  }

  //boundary conditions
  BndCndType init_tbc = GJP.Bc(3);

  if(reload_gauge){    //reimport gauge
    GJP.Bc(3,args.tbc);
    latt.BondCond(); //Don't forget to apply the boundary conditions!
    Float* gauge = (Float*) latt.GaugeField();
    dwf->cps_importGauge(gauge); 
  }  
  
  lanczos = new BFM_Krylov::Lanczos_5d<double>(*dwf,args.lanc_arg);

  lanczos->Run();

  if(reload_gauge){
    latt.BondCond(); //Don't forget to un-apply the boundary conditions!
    GJP.Bc(3,init_tbc);
  }
  precision = 2;
}

BFM_Krylov::Lanczos_5d<double> & LanczosContainer::getEig(Lattice &latt){
  if(precision != 2) ERR.General("LanczosContainer","getEig","Current precision is not double\n");
  if(lanczos == NULL) calcEig(latt);
  return *lanczos;
}
BFM_Krylov::Lanczos_5d<float> & LanczosContainer::getEigSinglePrec(Lattice &latt){
  if(precision != 1) ERR.General("LanczosContainer","getEigSinglePrec","Current precision is not single\n");
  if(lanczos_f == NULL){
    calcEig(latt);
    setPrecision(1);
  }
  return *lanczos_f;
}


void LanczosContainer::set_lanczos(BFM_Krylov::Lanczos_5d<double> *to){
  lanczos = to;
  if(to!=NULL) dwf = &lanczos->dop;
  else dwf = NULL;
  
  precision = 2;
  if(lanczos_f) delete lanczos_f;
}

template<typename FloatType>
void setupLanczosNoAlloc(BFM_Krylov::Lanczos_5d<FloatType> &lanc){ //just resize the vectors for the Lanczos output
  lanc.bq.resize(lanc.get);
  lanc.bl.resize(lanc.get);
}

//ASSUMED TO BE EXECUTED IN PARALLEL REGION!
template<typename FloatOut, typename FloatIn>
void moveAndChangePrecision(Fermion_t &out, Fermion_t &in, bfm_evo<FloatOut> &bfm_out,  bfm_evo<FloatIn> &bfm_in){
  //Allocate memory on receiving end
  out = bfm_out.threadedAllocCompactFermion();
  //Copy and change precision
  mixed_cg::threaded_convFermion_fast(out,in,bfm_out,bfm_in);
  //Deallocate input fermion
  bfm_in.threadedFreeFermion(in); in = NULL;
}


void LanczosContainer::setPrecision(const int &prec){ //convert from float to double or double to float
  if(precision == prec) return;
  const char *cname = "LanczosContainer";

  if(!UniqueID()){ printf("LanczosContainer with tag %s converting from precision %d to %d\n", tag(), precision, prec); fflush(stdout); }

  if(lanczos_f == NULL && lanczos == NULL){
    precision = prec; return;
  }

  if(dwf_f == NULL) setupBfm(1);
  if(dwf == NULL) setupBfm(2);
  
  if(precision == 2 && prec == 1){
  
    if(lanczos_f != NULL) ERR.General("LanczosContainer","setPrecision(1) [currently 2]","Expect lanczos_f pointer to be NULL!");
    lanczos_f = new BFM_Krylov::Lanczos_5d<float>(*dwf_f,args.lanc_arg);
    setupLanczosNoAlloc<float>(*lanczos_f);
    
    for(int i=0;i<lanczos_f->get;i++){
#pragma omp parallel
      {
	if(!lanczos->prec) moveAndChangePrecision<float,double>(lanczos_f->bq[i][0], lanczos->bq[i][0], *dwf_f, *dwf);
	moveAndChangePrecision<float,double>(lanczos_f->bq[i][1], lanczos->bq[i][1], *dwf_f, *dwf);
      }
      lanczos_f->bl[i] = lanczos->bl[i];
    }
    if(!UniqueID()){ printf("Deleting lanczos\n"); fflush(stdout); }

    delete lanczos; //because we nulled the pointers we deallocated, the lanczos destructor should not attempt to destroy them

  }else if(precision == 1 && prec == 2){

    if(lanczos != NULL) ERR.General("LanczosContainer","setPrecision(2) [currently 1]","Expect lanczos pointer to be NULL!");
    lanczos = new BFM_Krylov::Lanczos_5d<double>(*dwf,args.lanc_arg);
    setupLanczosNoAlloc<double>(*lanczos);

    for(int i=0;i<lanczos->get;i++){
#pragma omp parallel
      {
	if(!lanczos_f->prec) moveAndChangePrecision<double,float>(lanczos->bq[i][0], lanczos_f->bq[i][0], *dwf, *dwf_f);
	moveAndChangePrecision<double,float>(lanczos->bq[i][1], lanczos_f->bq[i][1], *dwf, *dwf_f);
      }
      lanczos->bl[i] = lanczos_f->bl[i];
    }
    delete lanczos_f;

  }else ERR.General("LanczosContainer","setPrecision","Invalid precision, %d",prec);

  precision = prec;
}

#endif

CPS_END_NAMESPACE


