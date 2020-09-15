#ifndef _COMPUTE_PROPS_H
#define _COMPUTE_PROPS_H

#define SETUP_ARRAY(OBJ,ARRAYNAME,TYPE,SIZE)	\
  OBJ . ARRAYNAME . ARRAYNAME##_len = SIZE; \
  OBJ . ARRAYNAME . ARRAYNAME##_val = new TYPE [SIZE]

#define ELEM(OBJ,ARRAYNAME,IDX) OBJ . ARRAYNAME . ARRAYNAME##_val[IDX]

void setupGparityPropagators(const std::string &tag, const int &tsrc, const int &mass, const double &resid, const bool &gauge_fix){
  //Setup  wall-momentum source propagators with flavor 0 and with both momentum p and -p, where p = pi/2L G  where G is the vector of G-parity directions. Together these can be used to form the 2x2 flavor matrix propagator
  //The +p propagator is given the tag specified, and the -p propagator the tag with _negmom appended
  
  std::string names[2] = {tag, tag+"_negmom"};
  for(int i=0;i<2;i++){
    PropagatorArg parg;

    GenericPropAttrArg & gen = parg.generics;
    gen.type = QPROPW_TYPE;
    gen.tag = strdup(names[i].c_str());
    gen.mass = mass;
    for(int d=0;d<4;d++)
      gen.bc[d] = GJP.Bc(d);

    SETUP_ARRAY(parg,attributes,AttributeContainer,5 + (gauge_fix ? 1:0) );
    
    ELEM(parg,attributes,0).type = WALL_SOURCE_ATTR;
    WallSourceAttrArg &srcarg = ELEM(parg,attributes,0).AttributeContainer_u.wall_source_attr;
    srcarg.t = tsrc;
    
    ELEM(parg,attributes,1).type = GPARITY_FLAVOR_ATTR;
    GparityFlavorAttrArg &gparg = ELEM(parg,attributes,1).AttributeContainer_u.gparity_flavor_attr;
    gparg.flavor = 0;

    ELEM(parg,attributes,2).type = MOMENTUM_ATTR;
    MomentumAttrArg &momarg = ELEM(parg,attributes,2).AttributeContainer_u.momentum_attr;
    for(int pp=0;pp<3;pp++)
      momarg.p[pp] = ( GJP.Bc(pp) == BND_CND_GPARITY ? 1 : 0 ); //units of pi/2L
    
    ELEM(parg,attributes,3).type = CG_ATTR;
    CGAttrArg &cgarg = ELEM(parg,attributes,3).AttributeContainer_u.cg_attr;
    cgarg.max_num_iter = 10000;
    cgarg.stop_rsd = resid;
    cgarg.true_rsd = resid;

    ELEM(parg,attributes,4).type = GPARITY_COMPLEX_CONJ_SOURCE_PARTNER_PROP_ATTR; //set this prop as partner to its complex conjugate
    CGAttrArg &oarg = ELEM(parg,attributes,4).AttributeContainer_u.gparity_complex_conj_source_partner_prop_attr;
    oarg.tag = strdup( names[(i+1) % 2].c_str() );

    if(gauge_fix){
      ELEM(parg,attributes,5).type = GAUGE_FIX_ATTR;
      GaugeFixAttrArg &gfarg = ELEM(parg,attributes,5).AttributeContainer_u.gauge_fix_attr;
      gfarg.gauge_fix_src = 1;
      gfarg.gauge_fix_snk = 0;
    }
    
    cout << "Adding propagator " << names[i] << endl;
    PropManager::addProp(prg);
  }
}

void setupStandardPropagator(const std::string &tag, const int &tsrc, const int &mass, const double &resid, const bool &gauge_fix){
  //Setup  wall source propagators with zero momentum
  
  PropagatorArg parg;

  GenericPropAttrArg & gen = parg.generics;
  gen.type = QPROPW_TYPE;
  gen.tag = strdup(tag.c_str());
  gen.mass = mass;
  for(int d=0;d<4;d++)
    gen.bc[d] = GJP.Bc(d);
  
  SETUP_ARRAY(parg,attributes,AttributeContainer,2 + (gauge_fix ? 1:0) );
    
  ELEM(parg,attributes,0).type = WALL_SOURCE_ATTR;
  WallSourceAttrArg &srcarg = ELEM(parg,attributes,0).AttributeContainer_u.wall_source_attr;
  srcarg.t = tsrc;
    
  ELEM(parg,attributes,1).type = CG_ATTR;
  CGAttrArg &cgarg = ELEM(parg,attributes,1).AttributeContainer_u.cg_attr;
  cgarg.max_num_iter = 10000;
  cgarg.stop_rsd = resid;
  cgarg.true_rsd = resid;

  if(gauge_fix){
    ELEM(parg,attributes,2).type = GAUGE_FIX_ATTR;
    GaugeFixAttrArg &gfarg = ELEM(parg,attributes,2).AttributeContainer_u.gauge_fix_attr;
    gfarg.gauge_fix_src = 1;
    gfarg.gauge_fix_snk = 0;
  }
    
  cout << "Adding propagator " << tag << endl;
  PropManager::addProp(prg);
}

void setupPropagator(const std::string &tag, const int &tsrc, const int &mass, const double &resid, const bool &gauge_fix){
  if(GJP.Gparity()) return setupGparityPropagators(tag,tsrc,mass,resid,gauge_fix);
  else setupStandardPropagator(tag,tsrc,mass,resid,gauge_fix);
}


#endif
