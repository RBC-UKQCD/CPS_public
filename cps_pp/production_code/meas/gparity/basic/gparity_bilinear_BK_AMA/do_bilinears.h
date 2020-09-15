#ifndef _DO_BILINEARS_H
#define _DO_BILINEARS_H

// struct ContractionTypeAllBilinears{
//   string prop_1<>;
//   string prop_2<>;
//   MomArg momenta<>;
//   string file<>;

//  rpccommand GENERATE_PRINT_METHOD;
//  rpccommand GENERATE_DEEPCOPY_METHOD;
// };

void setupBilinearArgsGparity(std::vector<ContractionTypeAllBilinears> &into, const int &tsrc, const bool &sloppy){
  //We want ll, lh and hh bilinears with momenta  pi/L G,  -pi/L G and 0, where G is the vector of G-parity directions (with integer entries)
  std::string base = sloppy ? "sloppy" : "exact";
  std::ostringstream tstr; tstr << "_t" << tsrc;

  MomArg mom[3];
  Float mul[3] = { 1.0, -1.0, 0.0 };  
  for(int m=0;m<3;m++){
    for(int d=0;d<3;d++){
      int L = GJP.NodeSites(d) * GJP.Nodes(d); 
      mom[m][d] = ( GJP.Bc(d) == BND_CND_GPARITY ? mul[m] * M_PI/Float(L) : 0.0 );
    }
  }

  std::string prop_names[] = { base + "_l" + tstr.str(), base + "_h" + tstr.str() };
  for(int p1=0;p1<2;p1++){
    for(int p2=0;p2<2;p2++){
      into.resize(into.size()+1);
      ContractionTypeAllBilinears &c = into.back();
      c.prop_1 = strdup(prop_names[p1].c_str());
      c.prop_2 = strdup(prop_names[p2].c_str());
      
      c.momenta.momenta_len = 3;
      c.momenta.momenta_val = new MomArg(3);
      for(int m=0;m<3;m++)
	c.momenta.momenta_val[m].deep_copy(mom[m]);
    }
  }
}
void setupBilinearArgsStandard(std::vector<ContractionTypeAllBilinears> &into, const int &tsrc, const bool &sloppy){
  //We want ll, lh and hh bilinears with momenta  0
  std::string base = sloppy ? "sloppy" : "exact";
  std::ostringstream tstr; tstr << "_t" << tsrc;

  MomArg mom;
  mom.p[0] = 0; mom.p[1] = 0; mom.p[2] = 0;
  
  std::string prop_names[] = { base + "_l" + tstr.str(), base + "_h" + tstr.str() };
  for(int p1=0;p1<2;p1++){
    for(int p2=0;p2<2;p2++){
      into.resize(into.size()+1);
      ContractionTypeAllBilinears &c = into.back();
      c.prop_1 = strdup(prop_names[p1].c_str());
      c.prop_2 = strdup(prop_names[p2].c_str());
      
      c.momenta.momenta_len = 1;
      c.momenta.momenta_val = new MomArg(1);
      c.momenta.momenta_val[0].deep_copy(mom);
    }
  }
}




template<typename MatrixType>
void do_bilinears(Lattice &latt, const int &conf_idx){
  int global_T = GJP.Tnodes()*GJP.TnodeSites();
  
  //Create storage for sums of sloppy and 'rest' (exact-sloppy) bilinear sets for each propagator pair
  int nbils = ama_arg.bilinear_args.bilinear_args_len;
  std::vector< ContractedBilinearSimple<MatrixType> > sloppy_sum(nbils);
  std::vector< ContractedBilinearSimple<MatrixType> > exact_sum(nbils);

  for(int tsrc = 0; tsrc < global_T; tsrc++){ //timeslice of source
    //Do all sloppy calculations
    std::vector< ContractedBilinearSimple<MatrixType> > bils_t(nbils); 
    
    for(int bil_meas = 0; bil_meas < nbils; bil_meas++){
      if(!UniqueID()){ printf("OUTPUT: Calculating sloppy propagators %d on timeslice %d\n",bil_meas,tsrc); fflush(stdout); }

      ContractionTypeAllBilinears &bil_arg = ama_arg.bilinear_args.bilinear_args_val[bil_meas];
      PropagatorContainer &prop_A = PropManager::getProp(bil_arg.prop_1);
      PropagatorContainer &prop_B = PropManager::getProp(bil_arg.prop_2);
    
      //Get the CGAttr (add if not already present) to set sloppy solve residuals
      CGAttrArg* cg_attr_A;
      CGAttrArg* cg_attr_B;
      if(!prop_A.getAttr(cg_attr_A)) cg_attr_A = addCGattr(prop_A);
      if(!prop_B.getAttr(cg_attr_B)) cg_attr_B = addCGattr(prop_B);

      //Set the sloppy residuals
      cg_attr_A->stop_rsd = ama_arg.sloppy_precision;
      cg_attr_A->true_rsd = ama_arg.sloppy_precision;
    
      cg_attr_B->stop_rsd = ama_arg.sloppy_precision;
      cg_attr_B->true_rsd = ama_arg.sloppy_precision;
    
      //Set source timeslices
      dynamic_cast<QPropWcontainer&>(prop_A).setSourceTimeslice(tsrc);
      dynamic_cast<QPropWcontainer&>(prop_B).setSourceTimeslice(tsrc);

      //Calculate props
      prop_A.calcProp(latt);	
      prop_B.calcProp(latt);

      if(!UniqueID()){ printf("OUTPUT: Calculating sloppy bilinears %d on timeslice %d\n",bil_meas,tsrc); fflush(stdout); }
      
      //Calculate sloppy bilinear
      ContractedBilinearSimple<MatrixType> &bil = bils_t[bil_meas];
      _multimom_helper<ContractedBilinearSimple<MatrixType> >::add_momenta(bil, bil_arg.momenta.momenta_val, bil_arg.momenta.momenta_len);
      bil.calculateBilinears(latt, bil_arg.prop_1, PropDFT::Dagger, bil_arg.prop_2, PropDFT::None); 

      //Temporal shift so prop A always lives at 0 (take advantage of time translation symmetry)
      if(tsrc>0) bil.Tshift(-tsrc);

      bil_sum[bil_meas] += bil;

      //Don't discard the propagators yet as they may need to be re-used for other bilinear measurements
    }
    //Discard all sloppy propagators
    for(int l=0;l<prop_arg.props.props_len;l++)
      PropManager::getProp(prop_arg.props.props_val[l].generics.tag).deleteProp();

    //Now check if we are on a timeslice for an exact solve. If so repeat the above with exact props and compute the 
    bool do_exact = false;
    for(int i=0; i< ama_arg.exact_solve_timeslices.exact_solve_timeslices_len; i++)
      if(tsrc == ama_arg.exact_solve_timeslices.exact_solve_timeslices_val[i]){ do_exact = true; break; }

    if(do_exact){
      for(int bil_meas = 0; bil_meas < nbils; bil_meas++){
	if(!UniqueID()) printf("OUTPUT: Calculating exact propagators %d on timeslice %d\n",bil_meas,tsrc);

	ContractionTypeAllBilinears &bil_arg = ama_arg.bilinear_args.bilinear_args_val[bil_meas];
	PropagatorContainer &prop_A = PropManager::getProp(bil_arg.prop_1);
	PropagatorContainer &prop_B = PropManager::getProp(bil_arg.prop_2);
	
	CGAttrArg* cg_attr_A; prop_A.getAttr(cg_attr_A);
	CGAttrArg* cg_attr_B; prop_B.getAttr(cg_attr_B);

	cg_attr_A->stop_rsd = ama_arg.exact_precision;
	cg_attr_A->true_rsd = ama_arg.exact_precision;
    
	cg_attr_B->stop_rsd = ama_arg.exact_precision;
	cg_attr_B->true_rsd = ama_arg.exact_precision;
    
	prop_A.calcProp(latt);	
	prop_B.calcProp(latt);

	if(!UniqueID()) printf("OUTPUT: Calculating exact bilinears %d on timeslice %d\n",bil_meas,tsrc);

	//Calculate exact bilinear
	ContractedBilinearSimple<MatrixType> ebil;
	_multimom_helper<ContractedBilinearSimple<MatrixType> >::add_momenta(ebil, bil_arg.momenta.momenta_val, bil_arg.momenta.momenta_len);
	ebil.calculateBilinears(latt, bil_arg.prop_1, PropDFT::Dagger, bil_arg.prop_2, PropDFT::None); 

	if(tsrc>0) ebil.Tshift(-tsrc);

	ebil -= bils_t[bil_meas]; //just the 'rest' part

	rest_sum[bil_meas] += ebil;
      }
      for(int l=0;l<prop_arg.props.props_len;l++)
	PropManager::getProp(prop_arg.props.props_val[l].generics.tag).deleteProp();
    }
  }

  //Take timeslice average of the sums
  for(int bil_meas = 0; bil_meas < nbils; bil_meas++){
    bil_sum[bil_meas]/=Float(global_T); //timeslice average
    rest_sum[bil_meas]/=Float(ama_arg.exact_solve_timeslices.exact_solve_timeslices_len);
  }
  //Write out results
  for(int bil_meas = 0; bil_meas < nbils; bil_meas++){
    std::string avg_file;
    { 
      std::ostringstream filestr; 
      filestr << ama_arg.bilinear_args.bilinear_args_val[bil_meas].file << "." << conf_idx; 
      avg_file = filestr.str();
    }
    std::string rest_file;
    { 
      std::ostringstream filestr; 
      filestr << ama_arg.bilinear_args.bilinear_args_val[bil_meas].file << "_rest." << conf_idx; 
      rest_file = filestr.str();
    }   
    bil_sum[bil_meas].write(avg_file);
    rest_sum[bil_meas].write(rest_file);
  }
}
    


#endif
