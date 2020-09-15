//This file is included from prop_dft.h, do not attempt to include it on its own

template<typename SiteObjectType>
void DFT<SiteObjectType>::transform(std::vector<SiteObjectType> &result) const{ //result will be a vector of size Lt
  //create vector for individual threads to accumulate into
  int global_T = GJP.TnodeSites()*GJP.Tnodes();
  int local_T = GJP.TnodeSites();
  int local_toff = GJP.TnodeCoor()*local_T;

  std::vector<std::vector<SiteObjectType> > thread_result(local_T);
  result.resize(local_T);

  for(int t=0;t<local_T;t++){
    setZero(result[t]);

    thread_result[t].resize(omp_get_max_threads());
    for(int thr=0;thr<omp_get_max_threads();thr++) setZero(thread_result[t][thr]);
  }

  //threaded site loop
#pragma omp parallel for default(shared)
  for(int x=0;x<GJP.VolNodeSites();x++){
    int x_pos_vec[4];
    {
      int rem = x;
      for(int i=0;i<4;i++){
	x_pos_vec[i] = rem % GJP.NodeSites(i) + GJP.NodeCoor(i)*GJP.NodeSites(i);
	rem /= GJP.NodeSites(i);
      }
    }

    int local_t = x_pos_vec[3] - local_toff;
    SiteObjectType tmp;
    evaluate(tmp, x, x_pos_vec, local_t);
    accumulate(thread_result[local_t][omp_get_thread_num()], tmp);
  }

  //accumulate over threads
  int work = threadAccumulateWork();
#pragma omp parallel for default(shared)
  for(int i=0;i<work;i++){
    threadAccumulate(result, thread_result, i);
  }

  //lattice sum, output vector is of size global_T
  std::vector<SiteObjectType> sum_globT(global_T);

  for(int global_t=0;global_t<global_T;global_t++){
    SiteObjectType &m = sum_globT[global_t];

    int local_t = global_t - local_toff;
    if(local_t >= 0 && local_t < local_T) m = result[local_t];
    else setZero(m);

    latticeSum(m);
  }
  result = sum_globT;
}



template<typename MatrixType>
void PropDFT::conjugateMomReorder(std::vector<std::vector<MatrixType> > &to, const std::vector<std::vector<MatrixType> > &from){
  //vector<vector<MatrixType> >[nmom][nt]
  for(mom_idx_map_type::iterator mom_it = mom_idx_map.begin(); mom_it != mom_idx_map.end(); mom_it++){
    std::vector<Float> mom = mom_it->first;
    int pos = mom_it->second;
    int conjpos = pos;
    if(mom[0]!=0.0 || mom[1]!=0.0 || mom[2]!=0.0){
      mom[0]*=-1; mom[1]*=-1; mom[2]*=-1;
      mom_idx_map_type::iterator conj_it = mom_idx_map.find(mom);
      if(conj_it == mom_idx_map.end()) ERR.General("PropDFT","conjugateMomReorder","Minus mom not present in map, why?");
      conjpos = conj_it->second;
    }
    to[conjpos] = from[pos];
  }
}

template<typename MatrixType>
struct doFlipMomentum{};

template<>
struct doFlipMomentum<WilsonMatrix>{
  static void flip(WilsonMatrix &mat);
};
template<>
struct doFlipMomentum<SpinColorFlavorMatrix>{
  static void flip(SpinColorFlavorMatrix &mat);
};


template<typename MatrixType>
void PropDFT::do_superscript(MatrixType &mat, const PropSuperscript &ss){
  if(ss == OpNone) return;
  else if(ss == OpTranspose) mat.transpose();
  else if(ss == OpConj) mat.cconj();
  else if(ss == OpDagger) mat.hconj();
  else if(ss == OpTransposeFlipMomentum){ doFlipMomentum<MatrixType>::flip(mat); mat.transpose(); }
  else if(ss == OpConjFlipMomentum){ doFlipMomentum<MatrixType>::flip(mat); mat.cconj(); }
  else if(ss == OpDaggerFlipMomentum){ doFlipMomentum<MatrixType>::flip(mat); mat.hconj(); }
  else if(ss == OpFlipMomentum){ doFlipMomentum<MatrixType>::flip(mat); }
  else ERR.General("PropDFT","do_superscript","Unknown superscript\n");
}

template<typename MatrixType>
void FourierProp<MatrixType>::add_momentum(std::vector<Float> sink_mom){
  if(!cosine_sink) return PropDFT::add_momentum(sink_mom);
  else{
    this->mom_idx_map[ sink_mom ] = this->nmom++; //do not add minus mom for cosine sink
  }
}

template<typename MatrixType>
void FourierProp<MatrixType>::enableCosineSink(){
  if(calculation_started) ERR.General("PropDFT","enableCosineSink(..)","Cannot modify cosine sink status after calculations have begun");
  if(nmom!=0) ERR.General("PropDFT","enableCosineSink(..)","Cannot modify cosine sink status after momenta have been added");
  cosine_sink = true;
}

template<typename MatrixType>
void FourierProp<MatrixType>::gaugeFixSink(const bool &tf){ 
  if(calculation_started) ERR.General("PropDFT","gaugeFixSink(..)","Cannot modify gauge fixing status after calculations have begun");
  gauge_fix_sink = tf; 
}

template<typename MatrixType>
void FourierProp<MatrixType>::clear(){
  calculation_started = false;
  gauge_fix_sink = false;
  cosine_sink = false;
  props.clear();
  PropDFT::clear();
}

template<typename MatrixType>
void FourierProp<MatrixType>::calcProp(const std::string &tag, Lattice &lat){
  calculation_started = true;
  std::vector<std::vector<MatrixType> > &mats = props[tag]; //vector in the sink momentum index
  mats.resize(nmom);
        
  int global_T = GJP.TnodeSites()*GJP.Tnodes();
  int local_T = GJP.TnodeSites();
  int local_toff = GJP.TnodeCoor()*local_T;

  std::vector<std::vector<std::vector<MatrixType> > > thread_mats(nmom); //[p][t][thread]
  for(int p=0;p<nmom;p++){
    mats[p].resize(local_T); //size to local_T for now. When lattice sum is performed, replace the vector with one of size global_T
    thread_mats[p].resize(local_T);
    for(int t=0;t<local_T;t++){
      mats[p][t] = 0.0;
      thread_mats[p][t].resize(omp_get_max_threads());
      for(int thr=0;thr<omp_get_max_threads();thr++) thread_mats[p][t][thr] = 0.0;
    }
  }
  QPropWcontainer &prop = QPropWcontainer::verify_convert(PropManager::getProp(tag.c_str()),"FourierProp<MatrixType>","calcProp(const std::string &tag, Lattice &lat)");

#pragma omp parallel for default(shared)
  for(int x=0;x<GJP.VolNodeSites();x++){
    int x_pos_vec[4];
    global_coord(x,x_pos_vec);
    int local_t = x_pos_vec[3] - local_toff;

    /*Get all SpinColorFlavorMatrices needed*/
    MatrixType mat; _FourierProp_helper<MatrixType>::site_matrix(mat, prop, lat, x);
      
    //Gauge fix sink
    if(gauge_fix_sink && lat.FixGaugeKind() != FIX_GAUGE_NONE) _FourierProp_helper<MatrixType>::mult_gauge_fix_mat(mat,x,lat);

    if(GJP.Gparity1f2fComparisonCode()){ 
      //For G-parity in 2 directions comparison, just Fourier transform over the two lower quadrants: the upper quadrants are simply
      //copies of these but the upper-right quadrant has the opposite sign:
      //  | C\bar u^T |    -d     |
      //  |    d      | C\bar u^T |
      if(GJP.Gparity1fX() && GJP.Gparity1fY() && x_pos_vec[1] >= GJP.YnodeSites()*GJP.Ynodes()/2) continue;

      //testing of 1f versus 2f code, if in 1f mode, where the second half of the lattice represents the second flavour,
      //shift position vectors on second half of lattice back onto first half
      //to ensure correct coordinate used in Fourier transform
      if(GJP.Gparity1fX() && x_pos_vec[0] >= GJP.XnodeSites()*GJP.Xnodes()/2) x_pos_vec[0] -= GJP.XnodeSites()*GJP.Xnodes()/2;
      if(GJP.Gparity1fY() && x_pos_vec[1] >= GJP.YnodeSites()*GJP.Ynodes()/2) x_pos_vec[1] -= GJP.YnodeSites()*GJP.Ynodes()/2;
    }

    //accumulate Fourier transform into vector of mats
    for(mom_idx_map_type::iterator mom_it = mom_idx_map.begin(); mom_it != mom_idx_map.end(); ++mom_it){
      const std::vector<Float> & mom = mom_it->first;
      const int &vec_pos = mom_it->second;

      Rcomplex phase;

      if(cosine_sink){
	// cos(p1*x1)*cos(p2*x2)*cos(p3*x3)
	phase = 1.0;
	for(int i=0;i<3;i++) phase *= cos(mom[i]*x_pos_vec[i]);
      }else{
	// e^{i p.x}
	Float pdotx = 0.0;
	for(int i=0;i<3;i++) pdotx += mom[i]*x_pos_vec[i];
	phase = Rcomplex(cos(pdotx),sin(pdotx));
      }
      thread_mats[ vec_pos ][local_t][omp_get_thread_num()] += mat * phase;
    }
  }
  //accumulate thread results; thread this too
#pragma omp parallel for default(shared)
  for(int i=0;i<nmom*local_T;i++){
    int rem = i;
    int t = rem % local_T; rem/=local_T;
    int p = rem;

    for(int thr=0;thr<omp_get_max_threads();thr++){
      mats[p][t] += thread_mats[p][t][thr];
    }
  }
    
  //lattice sum, output vector is of size global_T
  for(int p=0;p<nmom;p++){
    std::vector<MatrixType> matsum_globT(global_T);
    for(int global_t=0;global_t<global_T;global_t++){
      MatrixType &m = matsum_globT[global_t];

      int local_t = global_t - local_toff;
      if(local_t >= 0 && local_t < local_T) m = mats[p][local_t];
      else m = 0.0;

      _FourierProp_helper<MatrixType>::lattice_sum(m);
    }
    mats[p] = matsum_globT;
  }
    
}

template<typename MatrixType>
const std::vector<MatrixType> & FourierProp<MatrixType>::getFTProp(Lattice &lat,
						       const std::vector<Float> &sink_momentum, 
						       char const* tag
						       ){
  if(!mom_idx_map.count( sink_momentum )) ERR.General("FourierProp","getFTProp","Desired sink momentum is not in the list of those momentum components generated");
  int momentum_idx = mom_idx_map[sink_momentum];

  typename map_info_type::iterator it = props.find(tag);
  if(it != props.end()){ 	//ok, we have found the matching matrix
    return it->second[momentum_idx];
  }
  calcProp(tag, lat);
  return props.at(tag)[momentum_idx];
}

template<typename MatrixType>
void FourierProp<MatrixType>::write(const std::string &tag, const char *file, Lattice &lat){
  FILE *fp;
  if ((fp = Fopen(file, "w")) == NULL) {
    ERR.FileW("FourierProp","write(const char *file)",file);
  }
  write(tag, fp, lat);
  Fclose(fp);
}

template<typename MatrixType>
void FourierProp<MatrixType>::write(const std::string &tag, FILE *fp, Lattice &lat){
  typename map_info_type::const_iterator propset = props.find(tag);
  if(propset == props.end()) calcProp(tag, lat);

  const std::vector<std::vector<MatrixType> > &mats = props[tag];

  //Find all p^2 and sort
  std::vector< std::pair<Float,int> > p2list;
  std::map<int,std::vector<Float> > p2map;
  find_p2sorted(p2list,p2map);
    
  for(int p=0;p<p2list.size();p++){
    int pidx = p2list[p].second;
    const std::vector<Float> &mom = p2map[pidx];
    const Float &p2 = p2list[p].first;

    for(int t=0;t<GJP.TnodeSites()*GJP.Tnodes();t++){
      _FourierProp_helper<MatrixType>::write(fp,mats[pidx][t], p2, mom, t);
    }
  }
    
}

template<typename MatrixType>
void PropagatorBilinear<MatrixType>::calcAllBilinears(const prop_info_pair &props, Lattice &lat){
  nmom_fixed = true;
  const int &nidx = _PropagatorBilinear_helper<MatrixType>::nidx;
  std::vector< std::vector<std::vector<MatrixType> >* > mats(nidx); //internal vector<vector<matrix> > are indexed by [p][t] where t is global time
  for(int idx=0;idx<nidx;idx++){
    if(all_mats[idx] == NULL){
      all_mats[idx] = new map_info_type(); 
    }
    mats[idx] = &(*all_mats[idx])[props];
  }
  QPropWcontainer &prop_A = QPropWcontainer::verify_convert(PropManager::getProp(props.first.first.c_str()),"PropagatorBilinear<MatrixType>","calcAllBilinears(const prop_info_pair &props, Lattice &lat)");
  QPropWcontainer &prop_B = QPropWcontainer::verify_convert(PropManager::getProp(props.second.first.c_str()),"PropagatorBilinear<MatrixType>","calcAllBilinears(const prop_info_pair &props, Lattice &lat)");

  int global_T = GJP.Tnodes()*GJP.TnodeSites();
  int local_T = GJP.TnodeSites();
  int local_toff = GJP.TnodeCoor()*local_T;

  std::vector<std::vector<std::vector<std::vector<MatrixType> > > > thread_mats(nidx); //[nidx][nmom][nt][nthread]
  for(int i=0;i<nidx;i++){
    thread_mats[i].resize(nmom);
    mats[i]->resize(nmom);
    for(int p=0;p<nmom;p++){
      mats[i]->at(p).resize(local_T); //initially of size local_T, resize to global_T on lattice sum
      thread_mats[i][p].resize(local_T);
      for(int t=0;t<local_T;t++){
	mats[i]->at(p)[t] = 0.0;
	thread_mats[i][p][t].resize(omp_get_max_threads());
	for(int thr=0;thr<omp_get_max_threads();thr++) thread_mats[i][p][t][thr] = 0.0;
      }
    }
  }
#pragma omp parallel for default(shared)
  for(int x=0;x<GJP.VolNodeSites();x++){
    int x_pos_vec[4];
    global_coord(x,x_pos_vec);
    int local_t = x_pos_vec[3] - local_toff;

    /*Get all SpinColorFlavorMatrices needed*/
    MatrixType mat_A; _PropagatorBilinear_helper<MatrixType>::site_matrix(mat_A, prop_A, lat, x);
    do_superscript<MatrixType>(mat_A, props.first.second);
  
    MatrixType mat_B; _PropagatorBilinear_helper<MatrixType>::site_matrix(mat_B, prop_B, lat, x);
    do_superscript<MatrixType>(mat_B, props.second.second);

    for(int i=0;i<nidx;i++){
      MatrixType idx_mat_A = mat_A;
      //right-multiply with spin and flavour matrices
      std::pair<int,int> spin_flav = _PropagatorBilinear_helper<MatrixType>::unmap(i);
      _PropagatorBilinear_helper<MatrixType>::rmult_matrix(mat_A, spin_flav);
      
      idx_mat_A *= mat_B;

      //accumulate Fourier transform into vector of mats
      for(mom_idx_map_type::iterator mom_it = mom_idx_map.begin(); mom_it != mom_idx_map.end(); ++mom_it){
	const std::vector<Float> & mom = mom_it->first;
	const int &vec_pos = mom_it->second;

	Float pdotx = 0.0;
	for(int j=0;j<3;j++) pdotx += mom[j]*x_pos_vec[j];
	Rcomplex phase(cos(pdotx),sin(pdotx));
	thread_mats[ i ][ vec_pos ][ local_t ][omp_get_thread_num()] += idx_mat_A * phase;
      }
    }
  }
  //accumulate thread results, this can also be threaded
  int work = nidx * nmom * local_T; //= idx + nidx*p + nidx*nmom*t
#pragma omp parallel for default(shared)
  for(int i=0;i<work;i++){
    int rem = i;
    int idx = rem % nidx; rem/=nidx;
    int p = rem %nmom; rem/=nmom;
    int t = rem;

    for(int thr=0;thr<omp_get_max_threads();thr++){
      mats[idx]->at(p)[t] += thread_mats[idx][p][t][thr];
    }
  }

  //lattice sum, output vector is of size global_T
  for(int i=0;i<nidx;i++){
    for(int p=0;p<nmom;p++){
      std::vector<MatrixType> matsum_globT(global_T);
      for(int global_t=0;global_t<global_T;global_t++){
	MatrixType &m = matsum_globT[global_t];
	
	int local_t = global_t - local_toff;
	if(local_t >= 0 && local_t < local_T) m = mats[i]->at(p)[local_t];
	else m = 0.0;

	_PropagatorBilinear_helper<MatrixType>::lattice_sum(m);
      }
      mats[i]->at(p) = matsum_globT;
    }
  }

}

template<typename MatrixType>
void PropagatorBilinear<MatrixType>::calcBilinear(const int &idx, const prop_info_pair &props, Lattice &lat){
  nmom_fixed = true;
  bool first(false);
  if(all_mats[idx] == NULL){ all_mats[idx] = new map_info_type(); first = true; }

  std::vector<std::vector<MatrixType> > &mats = (*all_mats[idx])[props]; //vector in the sink momentum index and t
  mats.resize(nmom);

  QPropWcontainer &prop_A = QPropWcontainer::verify_convert(PropManager::getProp(props.first.first.c_str()), "PropagatorBilinear<MatrixType>","calcBilinear(const int &idx, const prop_info_pair &props, Lattice &lat)");
  QPropWcontainer &prop_B = QPropWcontainer::verify_convert(PropManager::getProp(props.second.first.c_str()),"PropagatorBilinear<MatrixType>","calcBilinear(const int &idx, const prop_info_pair &props, Lattice &lat)");
  std::pair<int,int> spin_flav = _PropagatorBilinear_helper<MatrixType>::unmap(idx); //if MatrixType is WilsonMatrix, the second index can be ignored
        
  int global_T = GJP.Tnodes()*GJP.TnodeSites();
  int local_T = GJP.TnodeSites();
  int local_toff = GJP.TnodeCoor()*local_T;

  std::vector<std::vector<std::vector<MatrixType> > > thread_mats(nmom); //[nmom][nt][nthread]
  for(int p=0;p<nmom;p++){
    mats[p].resize(local_T); //initially of size local_T, resize to global_T on lattice sum
    thread_mats[p].resize(local_T);
    for(int t=0;t<local_T;t++){
      mats[p][t] = 0.0;
      thread_mats[p][t].resize(omp_get_max_threads());
      for(int thr=0;thr<omp_get_max_threads();thr++) thread_mats[p][t][thr] = 0.0;
    }
  }

#pragma omp parallel for default(shared)
  for(int x=0;x<GJP.VolNodeSites();x++){
    int x_pos_vec[4];
    global_coord(x,x_pos_vec);
    int local_t = x_pos_vec[3] - local_toff;

    /*Get all SpinColorFlavorMatrices needed*/
    MatrixType mat_A; _PropagatorBilinear_helper<MatrixType>::site_matrix(mat_A, prop_A, lat, x);
    do_superscript<MatrixType>(mat_A, props.first.second);
  
    MatrixType mat_B; _PropagatorBilinear_helper<MatrixType>::site_matrix(mat_B, prop_B, lat, x);
    do_superscript<MatrixType>(mat_B, props.second.second);
      
    //right-multiply with spin and flavour matrices
    _PropagatorBilinear_helper<MatrixType>::rmult_matrix(mat_A, spin_flav);
      
    mat_A *= mat_B;

    //accumulate Fourier transform into vector of mats
    for(mom_idx_map_type::iterator mom_it = mom_idx_map.begin(); mom_it != mom_idx_map.end(); ++mom_it){
      const std::vector<Float> & mom = mom_it->first;
      const int &vec_pos = mom_it->second;

      Float pdotx = 0.0;
      for(int i=0;i<3;i++) pdotx += mom[i]*x_pos_vec[i];
      Rcomplex phase(cos(pdotx),sin(pdotx));
      thread_mats[ vec_pos ][ local_t ][omp_get_thread_num()] += mat_A * phase;
    }
  }
  //accumulate thread results, this can also be threaded
  int work = nmom * local_T; //= p + nmom*t
#pragma omp parallel for default(shared)
  for(int i=0;i<work;i++){
    int rem = i; 
    int p = rem %nmom; rem/=nmom;
    int t = rem;

    for(int thr=0;thr<omp_get_max_threads();thr++){
      mats[p][t] += thread_mats[p][t][thr];
    }
  }

  //lattice sum, output vector is of size global_T
  for(int p=0;p<nmom;p++){
    std::vector<MatrixType> matsum_globT(global_T);
    for(int global_t=0;global_t<global_T;global_t++){
      MatrixType &m = matsum_globT[global_t];
	
      int local_t = global_t - local_toff;
      if(local_t >= 0 && local_t < local_T) m = mats[p][local_t];
      else m = 0.0;

      _PropagatorBilinear_helper<MatrixType>::lattice_sum(m);
    }
    mats[p] = matsum_globT;
  }

}

//Get a (locally) Fourier transformed bilinear as a function of time
//Sigma is ignored if MatrixType is WilsonMatrix
template<typename MatrixType>
const std::vector<MatrixType> & PropagatorBilinear<MatrixType>::getBilinear(Lattice &lat,
									    const std::vector<Float> &sink_momentum, 
									    char const* tag_A, const PropSuperscript &ss_A,  
									    char const* tag_B, const PropSuperscript &ss_B, 
									    const int &Gamma, const int &Sigma
									    ){
  if(!mom_idx_map.count( sink_momentum )) ERR.General("PropagatorBilinear","getBilinear","Desired sink momentum is not in the list of those momentum components generated");
  int momentum_idx = mom_idx_map[sink_momentum];

  prop_info_pair props( prop_info(tag_A,ss_A), prop_info(tag_B,ss_B) );
  int scf_idx = _PropagatorBilinear_helper<MatrixType>::scf_map(Gamma,Sigma);
    
  if(scf_idx >= _PropagatorBilinear_helper<MatrixType>::nidx) ERR.General("PropagatorBilinear","getBilinear","Invalid Spin/Flavour matrix indices");

  if(all_mats[scf_idx] != NULL){
    typename map_info_type::iterator it = all_mats[scf_idx]->find(props);
    if(it != all_mats[scf_idx]->end()){ 	//ok, we have found the matching matrix
      return it->second[momentum_idx];
    }
  }
  calcBilinear(scf_idx, props, lat);
  return all_mats[scf_idx]->at(props)[momentum_idx];
}

template<typename MatrixType>
void PropagatorBilinear<MatrixType>::calcAllBilinears(Lattice &lat,
						      char const* tag_A, const PropSuperscript &ss_A,  
						      char const* tag_B, const PropSuperscript &ss_B){
  prop_info_pair props( prop_info(tag_A,ss_A), prop_info(tag_B,ss_B) );
  calcAllBilinears(props,lat);
}
  
template<typename MatrixType>
void PropagatorBilinear<MatrixType>::clear(){
  nmom_fixed = false;
  for(int i=0;i<_PropagatorBilinear_helper<MatrixType>::nidx;i++){
    if(all_mats[i]!=NULL){
      delete all_mats[i];
      all_mats[i] = NULL;
    }
    PropDFT::clear();
  }
}

template<typename MatrixType>
void PropagatorBilinear<MatrixType>::add_momentum(std::vector<Float> sink_mom){
  if(nmom_fixed) ERR.General("PropagatorBilinear","add_momentum","Cannot add momentum after bilinears have begun being calculated"); 
  return PropDFT::add_momentum(sink_mom);
}

template<typename MatrixType>
void PropagatorBilinear<MatrixType>::write(char const* tag_A, const PropSuperscript &ss_A,  
					   char const* tag_B, const PropSuperscript &ss_B,
					   const int &Gamma, const int &Sigma,
					   const char *file, Lattice &lat){
  FILE *fp;
  if ((fp = Fopen(file, "w")) == NULL) {
    ERR.FileW("PropagatorBilinear","write(const char *file)",file);
  }
  write(tag_A,ss_A,tag_B,ss_B, Gamma, Sigma, fp, lat);
  Fclose(fp);
}

template<typename MatrixType>
void PropagatorBilinear<MatrixType>::write(char const* tag_A, const PropSuperscript &ss_A,  
					   char const* tag_B, const PropSuperscript &ss_B,
					   const int &Gamma, const int &Sigma,
					   FILE *fp, Lattice &lat){
  static const char* fname = "write(...)";
  prop_info_pair props( prop_info(tag_A,ss_A), prop_info(tag_B,ss_B) );
  int scf_idx = _PropagatorBilinear_helper<MatrixType>::scf_map(Gamma,Sigma);

  if(all_mats[scf_idx] == NULL || all_mats[scf_idx]->count(props) ==0 ) calcBilinear(scf_idx, props, lat);
      
  const std::vector<std::vector<MatrixType> > &data = (*all_mats[scf_idx])[props];

  //Find all p^2 and sort
  std::vector< std::pair<Float,int> > p2list;
  std::map<int,std::vector<Float> > p2map;
  find_p2sorted(p2list,p2map);

  for(int p=0;p<p2list.size();p++){
    int pidx = p2list[p].second;
    const std::vector<Float> &mom = p2map[pidx];
    const Float &p2 = p2list[p].first;

    for(int t=0;t<GJP.TnodeSites()*GJP.Tnodes();t++){
      _PropagatorBilinear_helper<MatrixType>::write(fp,data[pidx][t], scf_idx, p2, mom, t);
    }
  }
}
    
template<typename MatrixType>
void PropagatorBilinear<MatrixType>::write(char const* tag_A, const PropSuperscript &ss_A,  
					   char const* tag_B, const PropSuperscript &ss_B,
					   const char *file, Lattice &lat){
  FILE *fp;
  if ((fp = Fopen(file, "w")) == NULL) {
    ERR.FileW("PropagatorBilinear","write(const char *file)",file);
  }
  write(tag_A,ss_A,tag_B,ss_B, fp, lat);
  Fclose(fp);
}

template<typename MatrixType>
void PropagatorBilinear<MatrixType>::write(char const* tag_A, const PropSuperscript &ss_A,  
					   char const* tag_B, const PropSuperscript &ss_B,
					   FILE *fp, Lattice &lat){
  static const char* fname = "write(...)";
  prop_info_pair props( prop_info(tag_A,ss_A), prop_info(tag_B,ss_B) );
    
  bool calc_needed(false);
  for(int i=0;i<_PropagatorBilinear_helper<MatrixType>::nidx;i++){
    if(all_mats[i] == NULL || all_mats[i]->count(props) ==0 ){ calc_needed = true; break; }
  }
  if(calc_needed) calcAllBilinears(props,lat);

  //Find all p^2 and sort
  std::vector< std::pair<Float,int> > p2list;
  std::map<int,std::vector<Float> > p2map;
  find_p2sorted(p2list,p2map);

  for(int scf_idx = 0; scf_idx < _PropagatorBilinear_helper<MatrixType>::nidx; scf_idx++){
    const std::vector<std::vector<MatrixType> > &data = (*all_mats[scf_idx])[props];

    for(int p=0;p<p2list.size();p++){
      int pidx = p2list[p].second;
      const std::vector<Float> &mom = p2map[pidx];
      const Float &p2 = p2list[p].first;
	
      for(int t=0;t<GJP.TnodeSites()*GJP.Tnodes();t++){
	_PropagatorBilinear_helper<MatrixType>::write(fp,data[pidx][t], scf_idx, p2, mom, t);
      }
    }
  }
}





template<typename MatrixType>
int ContractedBilinear<MatrixType>::idx_map(const int &mat1, const int &mat2, const int & mom_idx, const int &t) const{    
  return mat1 + nmat * (mat2 + nmat*(mom_idx + nmom*t));
}

template<typename MatrixType>
void ContractedBilinear<MatrixType>::idx_unmap(const int &idx, int &mat1, int &mat2, int & mom_idx, int &t) const{
  int rem = idx;
  mat1 = rem % nmat;  rem /= nmat;
  mat2 = rem % nmat;  rem /= nmat;
  mom_idx = rem % nmom; rem /= nmom;
  t = rem;
}

//We calculate \sum_x tr( M1 A(x) M2 B(x) )  where M1 and M2 are spin(-flavor) matrices and A,B are propagors with the superscript (dagger,conj,etc) specified in the prop_info_pair
//In this version we first determine   \sum_x A(x) M2 B(x) in a threaded loop over the sites, then afterwards loop over M1 and form the trace.

template<typename MatrixType>
void ContractedBilinear<MatrixType>::calcAllContractedBilinears1(const prop_info_pair &props, Lattice &lat){
  Timer::reset(); Elapsed time;
  if(!UniqueID()){ printf("ContractedBilinear::calcAllContractedBilinears1 starting with %d threads\n",omp_get_max_threads()); }

  if(array_size==-1) array_size = nmat*nmat*nmom*GJP.Tnodes()*GJP.TnodeSites(); //also acts as a lock to prevent further momenta from being added

  if(!results.count(props)) results[props] = new Rcomplex[array_size];
  Rcomplex *into = results[props];
  for(int i=0;i<array_size;i++) into[i] = 0.0;

  int local_T = GJP.TnodeSites();
  int global_T = local_T*GJP.Tnodes();
  int local_toff = GJP.TnodeCoor() * local_T;

  //This version uses more memory but has a reduced number of flops/site, so should run faster
  int mat_per_thread = nmat*nmom*local_T;
  MatrixType* bil = new MatrixType[mat_per_thread]; //idx = mat2 + nmat*(p + nmom*t))
  MatrixType* thread_bil = new MatrixType[mat_per_thread*omp_get_max_threads()]; //idx = mat2 + nmat*(p + nmom*(t + local_T*thread))
    
  for(int i=0;i<mat_per_thread*omp_get_max_threads();i++){
    if(i<mat_per_thread) bil[i] = 0.0;
    thread_bil[i] = 0.0;
  }
  QPropWcontainer &prop_A = QPropWcontainer::verify_convert(PropManager::getProp(props.first.first.c_str()), "ContractedBilinear<MatrixType>","calcAllContractedBilinears1(const prop_info_pair &props, Lattice &lat)");
  QPropWcontainer &prop_B = QPropWcontainer::verify_convert(PropManager::getProp(props.second.first.c_str()),"ContractedBilinear<MatrixType>","calcAllContractedBilinears1(const prop_info_pair &props, Lattice &lat)");
  time = Timer::relative_time(); time.print("Starting site loop at time ");

#pragma omp parallel for default(shared)
  for(int x=0;x<GJP.VolNodeSites();x++){
    int me = omp_get_thread_num();

    int x_pos_vec[4];
    global_coord(x,x_pos_vec);
    int local_t = x_pos_vec[3] - local_toff;

    //Get phases
    Rcomplex phases[nmom];
    for(mom_idx_map_type::iterator mom_it = mom_idx_map.begin(); mom_it != mom_idx_map.end(); ++mom_it){
      const std::vector<Float> & mom = mom_it->first;
      const int &vec_pos = mom_it->second;

      Float pdotx = 0.0;
      for(int i=0;i<3;i++) pdotx += mom[i]*x_pos_vec[i];
      phases[vec_pos].real(cos(pdotx));
      phases[vec_pos].imag(sin(pdotx));	
    }
    //Get propagators and act with superscripts
    MatrixType mat_A; _PropagatorBilinear_helper<MatrixType>::site_matrix(mat_A, prop_A, lat, x);
    do_superscript<MatrixType>(mat_A, props.first.second);
  
    MatrixType mat_B; _PropagatorBilinear_helper<MatrixType>::site_matrix(mat_B, prop_B, lat, x);
    do_superscript<MatrixType>(mat_B, props.second.second);
      
    //Loop over gamma matrix
    MatrixType tmp;
    for(int mat2=0;mat2<nmat;mat2++){  //spin + 16*flav
      tmp = mat_A;

      std::pair<int,int> spin_flav = _PropagatorBilinear_helper<MatrixType>::unmap(mat2);
      //right-multiply with spin and flavour matrices
      _PropagatorBilinear_helper<MatrixType>::rmult_matrix(tmp, spin_flav);

      tmp*=mat_B;
      for(int mom=0;mom<nmom;mom++){
	int off = mat2 + nmat*(mom + nmom*(local_t + local_T*omp_get_thread_num() ));
	thread_bil[off] += tmp * phases[mom];
      }
    }

    if(me == 0 && x==0){
      time = Timer::relative_time(); time.print("First site time ");
    }
  }//end of position loop
  
  time = Timer::relative_time(); time.print("Finished site loop and starting thread accumulation ");
  
  //thread the accumulate too, sum over data from all threads on each thread to prevent multiple writes to same location
#pragma omp parallel for default(shared)
  for(int i=0; i< mat_per_thread; i++){
    //idx_into = mat2 + nmat*(p + nmom*t))
    int mat2, p, t, rem(i);
    mat2 = rem%nmat; rem/=nmat;
    p = rem%nmom; rem/=nmom;
    t = rem;

    for(int thr=0;thr<omp_get_max_threads();thr++){
      //idx_from = mat2 + nmat*(p + nmom*(t + T*thread)) = idx_into + nmat*nmom*T*thread
      int idx_from = mat2 + nmat*(p + nmom*(t + local_T*thr));
      bil[i] += thread_bil[idx_from];
    }
  }
    
  delete[] thread_bil;
  
  time = Timer::relative_time(); time.print("Finished thread accumulation and starting outer matrix mult ");

  //now we have the bilinear for each inner mat, loop over outer (gamma1,sigma1) and inner (gamma2,sigma2) mats, p and t (threaded) and do trace
  int ntrace = nmat*nmat*nmom*GJP.TnodeSites(); //do traces local to this node and poke onto required array element

#pragma omp parallel for default(shared)
  for(int i=0;i<ntrace;i++){
    int mat1,mat2,p,t; //t is *local* here
    idx_unmap(i,mat1,mat2,p,t);

    int global_t = t + local_toff;

    std::pair<int,int> spin_flav = _PropagatorBilinear_helper<MatrixType>::unmap(mat1);
    MatrixType outer; _PropagatorBilinear_helper<MatrixType>::unit_matrix(outer);
    _PropagatorBilinear_helper<MatrixType>::rmult_matrix(outer, spin_flav);
      
    int tr_off = idx_map(mat1,mat2,p,global_t);
    int bil_off = mat2 + nmat*(p + nmom*t);
    into[tr_off] = Trace(outer, bil[bil_off]);
  }

  delete[] bil;

  time = Timer::relative_time(); time.print("Finished outer matrix mult and performing global sum ");

  //lattice sum
  slice_sum( (Float*)into, 2*array_size, 99); //2 for re/im, 99 is a *magic* number (we are abusing slice_sum here)

  time = Timer::relative_time(); time.print("Finished contraction");
}

//We calculate \sum_x tr( M1 A(x) M2 B(x) )  where M1 and M2 are spin(-flavor) matrices and A,B are propagors with the superscript (dagger,conj,etc) specified in the prop_info_pair
//In this version we calculate the full trace within the loop over sites. This version is slower than version 1

template<typename MatrixType>
void ContractedBilinear<MatrixType>::calcAllContractedBilinears2(const prop_info_pair &props, Lattice &lat){
  if(array_size==-1) array_size = nmat*nmat*nmom*GJP.Tnodes()*GJP.TnodeSites(); //also acts as a lock to prevent further momenta from being added

  if(!results.count(props)) results[props] = new Rcomplex[array_size];
  Rcomplex *into = results[props];

  //This version uses less memory but has more flops/site
  int thread_result_sz = omp_get_num_threads()*array_size; //idx = thread + nthread*(mat1 + nmat*(mat2+nmat*(p + nmom*t)))  t is *global*
  Rcomplex* thread_result = new Rcomplex[thread_result_sz];
    
  for(int i=0;i<thread_result_sz;i++){
    if(i<array_size) into[i] = 0.0;
    thread_result[i] =  0.0;
  }
  QPropWcontainer &prop_A = QPropWcontainer::verify_convert(PropManager::getProp(props.first.first.c_str()),"ContractedBilinear<MatrixType>","calcAllContractedBilinears2(const prop_info_pair &props, Lattice &lat)");
  QPropWcontainer &prop_B = QPropWcontainer::verify_convert(PropManager::getProp(props.second.first.c_str()),"ContractedBilinear<MatrixType>","calcAllContractedBilinears2(const prop_info_pair &props, Lattice &lat)");

#pragma omp parallel for default(shared)
  for(int x=0;x<GJP.VolNodeSites();x++){
    int x_pos_vec[4];
    global_coord(x,x_pos_vec);
  
    //Get phases
    Rcomplex phases[nmom];
    for(mom_idx_map_type::iterator mom_it = mom_idx_map.begin(); mom_it != mom_idx_map.end(); ++mom_it){
      const std::vector<Float> & mom = mom_it->first;
      const int &vec_pos = mom_it->second;

      Float pdotx = 0.0;
      for(int i=0;i<3;i++) pdotx += mom[i]*x_pos_vec[i];
      phases[vec_pos].real(cos(pdotx));
      phases[vec_pos].imag(sin(pdotx));
    }

    //Get propagators and act with superscripts
    MatrixType mat_A; _PropagatorBilinear_helper<MatrixType>::site_matrix(mat_A, prop_A, lat, x);
    do_superscript<MatrixType>(mat_A, props.first.second);
  
    MatrixType mat_B; _PropagatorBilinear_helper<MatrixType>::site_matrix(mat_B, prop_B, lat, x);
    do_superscript<MatrixType>(mat_B, props.second.second);

    MatrixType tmp;

    for(int mat1=0;mat1<nmat;mat1++){
      //use commutativity of trace to start with mat_B
      tmp = mat_B;
      std::pair<int,int> spin_flav = _PropagatorBilinear_helper<MatrixType>::unmap(mat1);
      _PropagatorBilinear_helper<MatrixType>::rmult_matrix(tmp, spin_flav);

      for(int mat2=0;mat2<nmat;mat2++){
	spin_flav = _PropagatorBilinear_helper<MatrixType>::unmap(mat2);
	MatrixType inner(mat_A);
	_PropagatorBilinear_helper<MatrixType>::rmult_matrix(inner, spin_flav); 
	Rcomplex tr = Trace(tmp,inner);

	for(int mom=0;mom<nmom;mom++){//idx = thread + nthread*(mat1 + nmat*(mat2+nmat*(p + nmom*t)))
	  int off = omp_get_thread_num() + omp_get_max_threads()*(mat1 + nmat*(mat2 + nmat*(mom + nmom*x_pos_vec[3])));
	  thread_result[off] += tr * phases[mom];
	}
      }
    }

  }//end of site loop

  //accumulate thread results, thread this!
#pragma omp parallel for default(shared)
  for(int i=0;i<array_size;i++){
    for(int thr=0;thr<omp_get_max_threads();thr++){
      int idx_from = thr + omp_get_max_threads()*i;
      into[i] += thread_result[idx_from];
    }
  }

  delete[] thread_result;

  //lattice sum
  slice_sum( (Float*)into, 2*array_size, 99); //2 for re/im, 99 is a *magic* number (we are abusing slice_sum here)
}

template<typename MatrixType>
void ContractedBilinear<MatrixType>::calculateBilinears(Lattice &lat,
							const prop_info_pair &props,
							const int &version
							){
  if(!results.count(props)){
    if(version==0) calcAllContractedBilinears1(props,lat);
    else calcAllContractedBilinears2(props,lat);
  }
}

template<typename MatrixType>
void ContractedBilinear<MatrixType>::calculateBilinears(Lattice &lat,
							char const* tag_A, const PropSuperscript &ss_A,  
							char const* tag_B, const PropSuperscript &ss_B, 
							const int &version
							){
  prop_info_pair props( prop_info(tag_A,ss_A), prop_info(tag_B,ss_B) );

  calculateBilinears(lat,props,version);
}


//Get a Fourier transformed bilinear correlation function as a function of time
//Sigma is ignored if MatrixType is WilsonMatrix
template<typename MatrixType>
std::vector<Rcomplex> ContractedBilinear<MatrixType>::getBilinear(Lattice &lat,
								  const std::vector<Float> &sink_momentum, 
								  char const* tag_A, const PropSuperscript &ss_A,  
								  char const* tag_B, const PropSuperscript &ss_B, 
								  const int &Gamma1, const int &Sigma1,
								  const int &Gamma2, const int &Sigma2,
								  const int &version
								  ){
  if(!mom_idx_map.count( sink_momentum )) ERR.General("PropagatorBilinear","getBilinear","Desired sink momentum is not in the list of those momentum components generated");
  int momentum_idx = mom_idx_map[sink_momentum];     
  prop_info_pair props( prop_info(tag_A,ss_A), prop_info(tag_B,ss_B) );
  calculateBilinears(lat,props,version);

  int scf_idx1 = _PropagatorBilinear_helper<MatrixType>::scf_map(Gamma1,Sigma1);
  int scf_idx2 = _PropagatorBilinear_helper<MatrixType>::scf_map(Gamma2,Sigma2);

  Rcomplex *con = results[props];
  std::vector<Rcomplex> out(GJP.TnodeSites()*GJP.Tnodes());
  for(int t=0;t<out.size();t++){
    out[t] = con[ idx_map(scf_idx1,scf_idx2,momentum_idx,t) ];
  }
  return out;
}

//for use with WilsonMatrix where sigma (flavour matrix idx) does not play a role
template<typename MatrixType>
std::vector<Rcomplex> ContractedBilinear<MatrixType>::getBilinear(Lattice &lat,
								  const std::vector<Float> &sink_momentum, 
								  char const* tag_A, const PropSuperscript &ss_A,  
								  char const* tag_B, const PropSuperscript &ss_B, 
								  const int &Gamma1, const int &Gamma2, const int &version
								  ){
  return getBilinear(lat,sink_momentum,tag_A,ss_A,tag_B,ss_B,Gamma1,0,Gamma2,0,version);
}
    
template<typename MatrixType>
void ContractedBilinear<MatrixType>::add_momentum(const std::vector<Float> &sink_mom){
  //If array_size is set we cannot add any more momenta without having to do lots of extra work to fill in gaps for existing bilinears
  if(array_size!=-1) ERR.General("ContractedBilinear","add_momentum","Cannot add momentum after bilinears have begun being calculated"); 
  return PropDFT::add_momentum(sink_mom);
}

template<typename MatrixType>
void ContractedBilinear<MatrixType>::clear(){
  for(typename map_info_type::iterator it = results.begin(); it!= results.end(); ++it){
    delete[] it->second;
  }
  results.clear();
  array_size = -1;
  PropDFT::clear();
}

template<typename MatrixType>
void ContractedBilinear<MatrixType>::write(char const* tag_A, const PropSuperscript &ss_A,  
					   char const* tag_B, const PropSuperscript &ss_B, 
					   const int &Gamma1, const int &Sigma1,
					   const int &Gamma2, const int &Sigma2,
					   const char *file, Lattice &lat){
#ifndef USE_OFSTREAM
  FILE *fp;
  if ((fp = Fopen(file, "w")) == NULL) {
    ERR.FileW("ContractedBilinear","write(...)",file);
  }
  write(tag_A,ss_A,tag_B,ss_B, Gamma1, Sigma1, Gamma2, Sigma2, fp, lat);
  Fclose(fp);
#else
  if(!UniqueID()) printf("ContractedBilinear File write using stream version\n");
  std::ostream *str = _DFTutils::open_ofstream(file);
  write(tag_A,ss_A,tag_B,ss_B, Gamma1, Sigma1, Gamma2, Sigma2, *str, lat);
  _DFTutils::close_ofstream(str);
#endif

}

template<typename MatrixType>
template<typename OutputType>
void ContractedBilinear<MatrixType>::write(char const* tag_A, const PropSuperscript &ss_A,  
					   char const* tag_B, const PropSuperscript &ss_B, 
					   const int &Gamma1, const int &Sigma1,
					   const int &Gamma2, const int &Sigma2,
					   OutputType fp, Lattice &lat){
  prop_info_pair props( prop_info(tag_A,ss_A), prop_info(tag_B,ss_B) );

  calculateBilinears(lat,props,0); //only calculates if not yet done

  int scf_idx1 = _PropagatorBilinear_helper<MatrixType>::scf_map(Gamma1,Sigma1);
  int scf_idx2 = _PropagatorBilinear_helper<MatrixType>::scf_map(Gamma2,Sigma2);

  Rcomplex *con = results[props];

  //Find all p^2 and sort
  std::vector< std::pair<Float,int> > p2list;
  std::map<int,std::vector<Float> > p2map;
  find_p2sorted(p2list,p2map);
    
  for(int p=0;p<p2list.size();p++){
    int pidx = p2list[p].second;
    const std::vector<Float> &mom = p2map[pidx];
    const Float &p2 = p2list[p].first;

    for(int t=0;t<GJP.TnodeSites()*GJP.Tnodes();t++){
      int off = idx_map(scf_idx1,scf_idx2,pidx,t);
      const Rcomplex &val = con[ off ];
      _ContractedBilinear_helper<MatrixType>::write(fp,val,scf_idx1,scf_idx2,p2, mom, t);
    }
  }
}
  
//write all combinations
template<typename MatrixType>
void ContractedBilinear<MatrixType>::write(char const* tag_A, const PropSuperscript &ss_A,  
					   char const* tag_B, const PropSuperscript &ss_B,
					   const std::string &file, Lattice &lat, const bool &binary){
  if(!UniqueID()) printf("ContractedBilinear writing to file \"%s\"\n",file.c_str());
  if(!UniqueID() && binary) printf("ContractedBilinear binary mode enabled\n");
#ifndef USE_OFSTREAM
  FILE *fp;
  if ((fp = Fopen(file.c_str(), binary ? "wb" : "w")) == NULL) {
    ERR.FileW("ContractedBilinear","write(...)",file.c_str());
  }

  if(binary && !UniqueID()){
    //To make sure that the file is read correctly later (correct endedness, etc), write a known number at the beginning which can be checked
    const static double test = 3.141592654;
    fwrite(&test, sizeof(double), 1, fp);
  }
  write(tag_A,ss_A,tag_B,ss_B,fp, lat, binary);
  Fclose(fp);
#else
  if(binary){
    ERR.General("ContractedBilinear","write(...)","Binary write not implemented for stream version");
  }
  if(!UniqueID()) printf("ContractedBilinear File write using stream version\n");
  std::ostream *str = _DFTutils::open_ofstream(file.c_str());
  write(tag_A,ss_A,tag_B,ss_B,*str, lat);
   _DFTutils::close_ofstream(str);
#endif
}

template<typename MatrixType>
template<typename OutputType>
void ContractedBilinear<MatrixType>::write(char const* tag_A, const PropSuperscript &ss_A,  
					   char const* tag_B, const PropSuperscript &ss_B, 
					   OutputType fp, Lattice &lat, const bool &binary){
  prop_info_pair props( prop_info(tag_A,ss_A), prop_info(tag_B,ss_B) );

  calculateBilinears(lat,props,0); //only calculates if not yet done

  Timer::reset(); Elapsed time;
  time = Timer::relative_time(); time.print("Starting file write setup");

  Rcomplex *con = results[props];

  //Find all p^2 and sort
  std::vector< std::pair<Float,int> > p2list;
  std::map<int,std::vector<Float> > p2map;
  find_p2sorted(p2list,p2map);

  std::vector< std::vector<Float> > mom_list;
  for(int p=0;p<p2list.size();p++){
    int pidx = p2list[p].second;
    mom_list.push_back(p2map[pidx]);
  }
  time = Timer::relative_time(); time.print("Starting write");
  for(int scf_idx1 = 0; scf_idx1 < nmat; scf_idx1++){
    for(int scf_idx2 = 0; scf_idx2 < nmat; scf_idx2++){
      for(int p=0;p<p2list.size();p++){
	int pidx = p2list[p].second;
	const std::vector<Float> &mom = mom_list[p];  //p2map[pidx];
	const Float &p2 = p2list[p].first;
	for(int t=0;t<GJP.TnodeSites()*GJP.Tnodes();t++){
	  int off = idx_map(scf_idx1,scf_idx2,pidx,t);
	  const Rcomplex &val = con[ off ];
	  binary ? 
	    _ContractedBilinear_helper<MatrixType>::binary_write(fp,val,scf_idx1,scf_idx2,p2, mom, t) : 
	    _ContractedBilinear_helper<MatrixType>::write(fp,val,scf_idx1,scf_idx2,p2, mom, t);
	}
      }
    }
  }
  time = Timer::relative_time(); time.print("Finished write");   
}




template<typename MatrixType>
void ContractedBilinearSimple<MatrixType>::add_momentum(const std::vector<Float> &sink_mom){
  //If array_size is set we cannot add any more momenta without having to do lots of extra work to fill in gaps for existing bilinears
  if(array_size!=-1) ERR.General("ContractedBilinearSimple","add_momentum","Cannot add momentum after bilinears have begun being calculated"); 
  return PropDFT::add_momentum(sink_mom);
}
template<typename MatrixType>
void ContractedBilinearSimple<MatrixType>::clear(){
  if(array_size != -1) delete[] results;
  array_size = -1;
  PropDFT::clear();
}
//We calculate \sum_x tr( M1 A(x) M2 B(x) )  where M1 and M2 are spin(-flavor) matrices and A,B are propagors with the superscript (dagger,conj,etc) specified in the prop_info_pair
//In this version we first determine   \sum_x A(x) M2 B(x) in a threaded loop over the sites, then afterwards loop over M1 and form the trace.

template<typename MatrixType>
void ContractedBilinearSimple<MatrixType>::calculateBilinears(Lattice &lat,
							      char const* tag_A, const PropSuperscript &ss_A,  
							      char const* tag_B, const PropSuperscript &ss_B){
  if(array_size==-1){
    array_size = nmat*nmat*nmom*GJP.Tnodes()*GJP.TnodeSites(); //also acts as a lock to prevent further momenta from being added
    results = new Rcomplex[array_size];
  }
  Rcomplex *into = results;
  for(int i=0;i<array_size;i++) into[i] = 0.0;

  int local_T = GJP.TnodeSites();
  int global_T = local_T*GJP.Tnodes();
  int local_toff = GJP.TnodeCoor() * local_T;

  //This version uses more memory but has a reduced number of flops/site, so should run faster
  int mat_per_thread = nmat*nmom*local_T;
  MatrixType* bil = new MatrixType[mat_per_thread]; //idx = mat2 + nmat*(p + nmom*t))
  MatrixType* thread_bil = new MatrixType[mat_per_thread*omp_get_max_threads()]; //idx = mat2 + nmat*(p + nmom*(t + local_T*thread))
    
  for(int i=0;i<mat_per_thread*omp_get_max_threads();i++){
    if(i<mat_per_thread) bil[i] = 0.0;
    thread_bil[i] = 0.0;
  }
  QPropWcontainer &prop_A = QPropWcontainer::verify_convert(PropManager::getProp(tag_A), "ContractedBilinearSimple<MatrixType>","calculateBilinears(..)");
  QPropWcontainer &prop_B = QPropWcontainer::verify_convert(PropManager::getProp(tag_B),"ContractedBilinearSimple<MatrixType>","calculateBilinears(..)");

#pragma omp parallel for default(shared)
  for(int x=0;x<GJP.VolNodeSites();x++){
    int x_pos_vec[4];
    global_coord(x,x_pos_vec);
    int local_t = x_pos_vec[3] - local_toff;

    //Get phases
    Rcomplex phases[nmom];
    for(mom_idx_map_type::iterator mom_it = mom_idx_map.begin(); mom_it != mom_idx_map.end(); ++mom_it){
      const std::vector<Float> & mom = mom_it->first;
      const int &vec_pos = mom_it->second;

      Float pdotx = 0.0;
      for(int i=0;i<3;i++) pdotx += mom[i]*x_pos_vec[i];
      phases[vec_pos].real(cos(pdotx));
      phases[vec_pos].imag(sin(pdotx));	
    }
    //Get propagators and act with superscripts
    MatrixType mat_A; _PropagatorBilinear_helper<MatrixType>::site_matrix(mat_A, prop_A, lat, x);
    do_superscript<MatrixType>(mat_A, ss_A);
  
    MatrixType mat_B; _PropagatorBilinear_helper<MatrixType>::site_matrix(mat_B, prop_B, lat, x);
    do_superscript<MatrixType>(mat_B, ss_B);
      
    //Loop over gamma matrix
    MatrixType tmp;
    for(int mat2=0;mat2<nmat;mat2++){  //spin + 16*flav
      tmp = mat_A;

      std::pair<int,int> spin_flav = _PropagatorBilinear_helper<MatrixType>::unmap(mat2);
      //right-multiply with spin and flavour matrices
      _PropagatorBilinear_helper<MatrixType>::rmult_matrix(tmp, spin_flav);

      tmp*=mat_B;
      for(int mom=0;mom<nmom;mom++){
	int off = mat2 + nmat*(mom + nmom*(local_t + local_T*omp_get_thread_num() ));
	thread_bil[off] += tmp * phases[mom];
      }
    }
  }//end of position loop
  

  //thread the accumulate too, sum over data from all threads on each thread to prevent multiple writes to same location
#pragma omp parallel for default(shared)
  for(int i=0; i< mat_per_thread; i++){
    //idx_into = mat2 + nmat*(p + nmom*t))
    int mat2, p, t, rem(i);
    mat2 = rem%nmat; rem/=nmat;
    p = rem%nmom; rem/=nmom;
    t = rem;

    for(int thr=0;thr<omp_get_max_threads();thr++){
      //idx_from = mat2 + nmat*(p + nmom*(t + T*thread)) = idx_into + nmat*nmom*T*thread
      int idx_from = mat2 + nmat*(p + nmom*(t + local_T*thr));
      bil[i] += thread_bil[idx_from];
    }
  }
    
  delete[] thread_bil;
  
  //now we have the bilinear for each inner mat, loop over outer (gamma1,sigma1) and inner (gamma2,sigma2) mats, p and t (threaded) and do trace
  int ntrace = nmat*nmat*nmom*GJP.TnodeSites(); //do traces local to this node and poke onto required array element

#pragma omp parallel for default(shared)
  for(int i=0;i<ntrace;i++){
    int mat1,mat2,p,t; //t is *local* here
    idx_unmap(i,mat1,mat2,p,t);

    int global_t = t + local_toff;

    std::pair<int,int> spin_flav = _PropagatorBilinear_helper<MatrixType>::unmap(mat1);
    MatrixType outer; _PropagatorBilinear_helper<MatrixType>::unit_matrix(outer);
    _PropagatorBilinear_helper<MatrixType>::rmult_matrix(outer, spin_flav);
      
    int tr_off = idx_map(mat1,mat2,p,global_t);
    int bil_off = mat2 + nmat*(p + nmom*t);
    into[tr_off] = Trace(outer, bil[bil_off]);
  }

  delete[] bil;

  //lattice sum
  slice_sum( (Float*)into, 2*array_size, 99); //2 for re/im, 99 is a *magic* number (we are abusing slice_sum here)
}

//Get a Fourier transformed bilinear correlation function as a function of time
//Sigma is ignored if MatrixType is WilsonMatrix
template<typename MatrixType>
std::vector<Rcomplex> ContractedBilinearSimple<MatrixType>::getBilinear(const std::vector<Float> &sink_momentum, 
									const int &Gamma1, const int &Sigma1,
									const int &Gamma2, const int &Sigma2){
  if(!mom_idx_map.count( sink_momentum )) ERR.General("ContractedBilinearSimple","getBilinear","Desired sink momentum is not in the list of those momentum components generated");
  if(array_size==-1) ERR.General("ContractedBilinearSimple","getBilinear","Must calculate the bilinears before trying to retrieve them");

  int momentum_idx = mom_idx_map[sink_momentum];     

  int scf_idx1 = _PropagatorBilinear_helper<MatrixType>::scf_map(Gamma1,Sigma1);
  int scf_idx2 = _PropagatorBilinear_helper<MatrixType>::scf_map(Gamma2,Sigma2);

  std::vector<Rcomplex> out(GJP.TnodeSites()*GJP.Tnodes());
  for(int t=0;t<out.size();t++){
    out[t] = results[ idx_map(scf_idx1,scf_idx2,momentum_idx,t) ];
  }
  return out;
}

//for use with WilsonMatrix where sigma (flavour matrix idx) does not play a role
template<typename MatrixType>
std::vector<Rcomplex> ContractedBilinearSimple<MatrixType>::getBilinear(const std::vector<Float> &sink_momentum, 
									const int &Gamma1, const int &Gamma2){
  return getBilinear(sink_momentum,Gamma1,0,Gamma2,0);
}

//write all combinations
template<typename MatrixType>
void ContractedBilinearSimple<MatrixType>::write(FILE *fp){
  if(array_size==-1) ERR.General("ContractedBilinearSimple","write","Must calculate the bilinears before trying to retrieve them");

  //Find all p^2 and sort
  std::vector< std::pair<Float,int> > p2list;
  std::map<int,std::vector<Float> > p2map;
  find_p2sorted(p2list,p2map);

  for(int scf_idx1 = 0; scf_idx1 < nmat; scf_idx1++){
    for(int scf_idx2 = 0; scf_idx2 < nmat; scf_idx2++){
      for(int p=0;p<p2list.size();p++){
	int pidx = p2list[p].second;
	const std::vector<Float> &mom = p2map[pidx];
	const Float &p2 = p2list[p].first;
	  
	for(int t=0;t<GJP.TnodeSites()*GJP.Tnodes();t++){
	  int off = idx_map(scf_idx1,scf_idx2,pidx,t);
	  const Rcomplex &val = results[ off ];
	  _ContractedBilinear_helper<MatrixType>::write(fp,val,scf_idx1,scf_idx2,p2, mom, t);
	}
      }
    }
  }
    
}

template<typename MatrixType>
void ContractedBilinearSimple<MatrixType>::write(const std::string &file){
  if(array_size==-1) ERR.General("ContractedBilinearSimple","write","Must calculate the bilinears before trying to retrieve them");
  if(!UniqueID()) printf("ContractedBilinearSimple writing to file \"%s\"\n",file.c_str());
  FILE *fp;
  if ((fp = Fopen(file.c_str(), "w")) == NULL) {
    ERR.FileW("ContractedBilinear","write(...)",file.c_str());
  }
  write(fp);
  Fclose(fp);
}

//Contents become the sum of the contents of this object and r
template<typename MatrixType>
ContractedBilinearSimple<MatrixType> & ContractedBilinearSimple<MatrixType>::operator+=(const ContractedBilinearSimple<MatrixType> &r){
  //If this object is empty, we simply copy the contents of r
  if(array_size == -1){
    copy_momenta(r);
    array_size = r.array_size;
    results = new Rcomplex[array_size];
    for(int i=0;i<array_size;i++) results[i] = r.results[i];
  }else{
    if(r.array_size != array_size) ERR.General("ContractedBilinearSimple","operator+=","Cannot combine two objects with different results sizes");
    for(int i=0;i<array_size;i++) results[i] += r.results[i];
  }
  return *this;
}
template<typename MatrixType>
ContractedBilinearSimple<MatrixType> & ContractedBilinearSimple<MatrixType>::operator-=(const ContractedBilinearSimple<MatrixType> &r){
  //If this object is empty, we simply copy the contents of r with a minus sign
  if(array_size == -1){
    copy_momenta(r);
    array_size = r.array_size;
    results = new Rcomplex[array_size];
    for(int i=0;i<array_size;i++) results[i] = -r.results[i];
  }else{
    if(r.array_size != array_size) ERR.General("ContractedBilinearSimple","operator-=","Cannot combine two objects with different results sizes");
    for(int i=0;i<array_size;i++) results[i] -= r.results[i];
  }
  return *this;
}

//Contents are divided by a float
template<typename MatrixType>
ContractedBilinearSimple<MatrixType> & ContractedBilinearSimple<MatrixType>::operator/=(const Float &r){
  if(array_size==-1) return *this;
  for(int i=0;i<array_size;i++) results[i] /= r;
  return *this;
}

//Shift the data:   data[t] -> data[t+dt]  (modulo lattice size)
template<typename MatrixType>
void ContractedBilinearSimple<MatrixType>::Tshift(const int &dt){
  if(array_size == -1) return;
  Rcomplex *tmp = new Rcomplex[array_size];
  int global_T = GJP.Tnodes()*GJP.TnodeSites();

#pragma omp parallel for
  for(int i=0;i<array_size;i++){
    int mat1,mat2,mom_idx,t;
    idx_unmap(i,mat1,mat2,mom_idx,t);
    int tnew = (t+dt+global_T) % global_T;
    tmp[ idx_map(mat1,mat2,mom_idx,tnew) ] = results[i];
  }
  results = tmp;
}



template<typename MatrixType>
void PropagatorQuadrilinear<MatrixType>::unmap(int pair_idx, int &idx1, int &idx2) const{
  idx1 = pair_idx % nmat; pair_idx/=nmat;
  idx2 = pair_idx;
}

template<typename MatrixType>
void PropagatorQuadrilinear<MatrixType>::calcQuadrilinear(const int &idx, const prop_info_quad &props, Lattice &lat){
  nmom_fixed = true;
  std::pair<prop_info_quad,int> key(props,idx);
    
  std::vector<std::vector<TensorType> > &result = all_mats[key]; //[p][t]
  std::vector<std::vector<std::vector<TensorType> > > thread_result; //[p][t][thread]

  int global_T = GJP.Tnodes()*GJP.TnodeSites();
  int local_T = GJP.TnodeSites();
  int local_toff = GJP.TnodeCoor()*local_T;

  result.resize(nmom);
  thread_result.resize(nmom);
  for(int p=0;p<nmom;p++){
    result[p].resize(local_T); //initially of size local_T, resize to global_T on lattice sum
    thread_result[p].resize(local_T);
    for(int t=0;t<local_T;t++){
      result[p][t] = 0.0;
      thread_result[p][t].resize(omp_get_max_threads());
      for(int thr=0;thr<omp_get_max_threads();thr++){
	thread_result[p][t][thr] = 0.0;
      }
    }
  }
  const char* cname = "PropagatorQuadrilinear<MatrixType>";
  const char* fname = "calcQuadrilinear(const int &idx, const prop_info_quad &props, Lattice &lat)";
  //prop_info_quad:  [bil][prop][first=tag, second=superscript]
  QPropWcontainer &prop_A = QPropWcontainer::verify_convert(PropManager::getProp(props.first.first.first.c_str()),cname,fname);
  QPropWcontainer &prop_B = QPropWcontainer::verify_convert(PropManager::getProp(props.first.second.first.c_str()),cname,fname);
  QPropWcontainer &prop_C = QPropWcontainer::verify_convert(PropManager::getProp(props.second.first.first.c_str()),cname,fname);
  QPropWcontainer &prop_D = QPropWcontainer::verify_convert(PropManager::getProp(props.second.second.first.c_str()),cname,fname);

  int idx1, idx2;
  unmap(idx,idx1,idx2);

  std::pair<int,int> spin_flav_1 = _PropagatorBilinear_helper<MatrixType>::unmap(idx1);
  std::pair<int,int> spin_flav_2 = _PropagatorBilinear_helper<MatrixType>::unmap(idx2);
        
#pragma omp parallel for default(shared)
  for(int x=0;x<GJP.VolNodeSites();x++){
    int x_pos_vec[4];
    global_coord(x,x_pos_vec);
    int local_t = x_pos_vec[3] - local_toff;

    /*Get all SpinColorFlavorMatrices needed*/
    MatrixType mat_A; _PropagatorBilinear_helper<MatrixType>::site_matrix(mat_A, prop_A, lat, x);
    do_superscript<MatrixType>(mat_A, props.first.first.second);
  
    MatrixType mat_B; _PropagatorBilinear_helper<MatrixType>::site_matrix(mat_B, prop_B, lat, x);
    do_superscript<MatrixType>(mat_B, props.first.second.second);
      
    MatrixType mat_C; _PropagatorBilinear_helper<MatrixType>::site_matrix(mat_C, prop_C, lat, x);
    do_superscript<MatrixType>(mat_C, props.second.first.second);
  
    MatrixType mat_D; _PropagatorBilinear_helper<MatrixType>::site_matrix(mat_D, prop_D, lat, x);
    do_superscript<MatrixType>(mat_D, props.second.second.second);

    //right-multiply with spin and flavour matrices
    _PropagatorBilinear_helper<MatrixType>::rmult_matrix(mat_A, spin_flav_1);
    _PropagatorBilinear_helper<MatrixType>::rmult_matrix(mat_C, spin_flav_2);
      
    mat_A *= mat_B;
    mat_C *= mat_D;

    if(GJP.Gparity1f2fComparisonCode()){ 
      //For G-parity in 2 directions comparison, just Fourier transform over the two lower quadrants: the upper quadrants are simply
      //copies of these but the upper-right quadrant has the opposite sign:
      //  | C\bar u^T |    -d     |
      //  |    d      | C\bar u^T |
      if(GJP.Gparity1fX() && GJP.Gparity1fY() && x_pos_vec[1] >= GJP.YnodeSites()*GJP.Ynodes()/2) continue;

      //testing of 1f versus 2f code, if in 1f mode, where the second half of the lattice represents the second flavour,
      //shift position vectors on second half of lattice back onto first half
      //to ensure correct coordinate used in Fourier transform
      if(GJP.Gparity1fX() && x_pos_vec[0] >= GJP.XnodeSites()*GJP.Xnodes()/2) x_pos_vec[0] -= GJP.XnodeSites()*GJP.Xnodes()/2;
      if(GJP.Gparity1fY() && x_pos_vec[1] >= GJP.YnodeSites()*GJP.Ynodes()/2) x_pos_vec[1] -= GJP.YnodeSites()*GJP.Ynodes()/2;
    }

    //accumulate Fourier transform into vector of mats
    for(mom_idx_map_type::iterator mom_it = mom_idx_map.begin(); mom_it != mom_idx_map.end(); ++mom_it){
      const std::vector<Float> & mom = mom_it->first;
      const int &vec_pos = mom_it->second;

      Float pdotx = 0.0;
      for(int i=0;i<3;i++) pdotx += mom[i]*x_pos_vec[i];
      Rcomplex phase(cos(pdotx),sin(pdotx));
      thread_result[ vec_pos ][ local_t ][omp_get_thread_num()].add(mat_A,mat_C,phase);
    }
  }
  //accumulate thread results, this can also be threaded
  int work = nmom * local_T; //= p + nmom*t
#pragma omp parallel for default(shared)
  for(int i=0;i<work;i++){
    int rem = i; 
    int p = rem %nmom; rem/=nmom;
    int t = rem;

    for(int thr=0;thr<omp_get_max_threads();thr++){
      result[p][t] += thread_result[p][t][thr];
    }
  }

  //lattice sum, output vector is of size global_T
  for(int p=0;p<nmom;p++){
    std::vector<TensorType> matsum_globT(global_T);
    for(int global_t=0;global_t<global_T;global_t++){
      TensorType &m = matsum_globT[global_t];
	
      int local_t = global_t - local_toff;
      if(local_t >= 0 && local_t < local_T) m = result[p][local_t];
      else{ m = 0.0; }

      _PropagatorQuadrilinear_helper<MatrixType>::lattice_sum(m);
    }
    result[p] = matsum_globT;
  }

}

template<typename MatrixType>
void PropagatorQuadrilinear<MatrixType>::calcAllQuadrilinears(const prop_info_quad &props, Lattice &lat){
  //NOTE:
  //RUNNING THIS FUNCTION IS A VERY BAD IDEA FOR G-PARITY BCS:
  //SCFTensor storage 24^4 * 2 * 8 = 5.0625 MB
  //64^2 of these is 20.25 GB    *TOO DAMN BIG*

  //SCTensor storage 12^4 * 2 * 8 = 0.3164 MB
  //16^2 of these is 81 MB  is fine

  //For G-parity I would run with a restricted set of matrix choices
  //(For many applications, the lattice discrete symmetries reduces the number of non-zero vertices after lattice sum)

  nmom_fixed = true;

  int global_T = GJP.Tnodes()*GJP.TnodeSites();
  int local_T = GJP.TnodeSites();
  int local_toff = GJP.TnodeCoor()*local_T;

  int nmatpairs = nmat*nmat;
  std::vector<std::vector<TensorType> >* into[nmatpairs];
  for(int i=0;i<nmatpairs;i++){
    std::pair<prop_info_quad,int> key(props,i);
    into[i] = &all_mats[key]; //[p][t]
    into[i]->resize(nmom);
    for(int p=0;p<nmom;p++){
      into[i]->at(p).resize(local_T); //initially of size local_T, resize to global_T on lattice sum
      for(int t=0;t<local_T;t++){
	into[i]->at(p)[t] = 0.0;
      }
    }
  }
  static const char* cname = "PropagatorQuadrilinear<MatrixType>";
  static const char* fname = "calcAllQuadrilinears(const prop_info_quad &props, Lattice &lat)";

  QPropWcontainer &prop_A = QPropWcontainer::verify_convert(PropManager::getProp(props.first.first.first.c_str()),cname,fname); //prop_info_quad:  [bil][prop][first=tag, second=superscript]
  QPropWcontainer &prop_B = QPropWcontainer::verify_convert(PropManager::getProp(props.first.second.first.c_str()),cname,fname);
  QPropWcontainer &prop_C = QPropWcontainer::verify_convert(PropManager::getProp(props.second.first.first.c_str()),cname,fname);
  QPropWcontainer &prop_D = QPropWcontainer::verify_convert(PropManager::getProp(props.second.second.first.c_str()),cname,fname);

#pragma omp parallel for default(shared)
  for(int idx=0; idx<nmatpairs; idx++){

    int idx1, idx2;
    unmap(idx,idx1,idx2);

    std::pair<int,int> spin_flav_1 = _PropagatorBilinear_helper<MatrixType>::unmap(idx1);
    std::pair<int,int> spin_flav_2 = _PropagatorBilinear_helper<MatrixType>::unmap(idx2);
        
    for(int x=0;x<GJP.VolNodeSites();x++){
      int x_pos_vec[4];
      global_coord(x,x_pos_vec);
      int local_t = x_pos_vec[3] - local_toff;

      /*Get all SpinColorFlavorMatrices needed*/
      MatrixType mat_A; _PropagatorBilinear_helper<MatrixType>::site_matrix(mat_A, prop_A, lat, x);
      do_superscript<MatrixType>(mat_A, props.first.first.second);
  
      MatrixType mat_B; _PropagatorBilinear_helper<MatrixType>::site_matrix(mat_B, prop_B, lat, x);
      do_superscript<MatrixType>(mat_B, props.first.second.second);
      
      MatrixType mat_C; _PropagatorBilinear_helper<MatrixType>::site_matrix(mat_C, prop_C, lat, x);
      do_superscript<MatrixType>(mat_C, props.second.first.second);
  
      MatrixType mat_D; _PropagatorBilinear_helper<MatrixType>::site_matrix(mat_D, prop_D, lat, x);
      do_superscript<MatrixType>(mat_D, props.second.second.second);

      //right-multiply with spin and flavour matrices
      _PropagatorBilinear_helper<MatrixType>::rmult_matrix(mat_A, spin_flav_1);
      _PropagatorBilinear_helper<MatrixType>::rmult_matrix(mat_C, spin_flav_2);
      
      mat_A *= mat_B;
      mat_C *= mat_D;

      if(GJP.Gparity1f2fComparisonCode()){ 
	//For G-parity in 2 directions comparison, just Fourier transform over the two lower quadrants: the upper quadrants are simply
	//copies of these but the upper-right quadrant has the opposite sign:
	//  | C\bar u^T |    -d     |
	//  |    d      | C\bar u^T |
	if(GJP.Gparity1fX() && GJP.Gparity1fY() && x_pos_vec[1] >= GJP.YnodeSites()*GJP.Ynodes()/2) continue;

	//testing of 1f versus 2f code, if in 1f mode, where the second half of the lattice represents the second flavour,
	//shift position vectors on second half of lattice back onto first half
	//to ensure correct coordinate used in Fourier transform
	if(GJP.Gparity1fX() && x_pos_vec[0] >= GJP.XnodeSites()*GJP.Xnodes()/2) x_pos_vec[0] -= GJP.XnodeSites()*GJP.Xnodes()/2;
	if(GJP.Gparity1fY() && x_pos_vec[1] >= GJP.YnodeSites()*GJP.Ynodes()/2) x_pos_vec[1] -= GJP.YnodeSites()*GJP.Ynodes()/2;
      }

      //accumulate Fourier transform into vector of mats
      for(mom_idx_map_type::iterator mom_it = mom_idx_map.begin(); mom_it != mom_idx_map.end(); ++mom_it){
	const std::vector<Float> & mom = mom_it->first;
	const int &vec_pos = mom_it->second;

	Float pdotx = 0.0;
	for(int i=0;i<3;i++) pdotx += mom[i]*x_pos_vec[i];
	Rcomplex phase(cos(pdotx),sin(pdotx));
	into[ idx ]->at( vec_pos )[ local_t ].add(mat_A,mat_C,phase);
      }
    }
  }

  //lattice sum, output vector is of size global_T
  for(int idx=0; idx<nmatpairs; idx++){
    for(int p=0;p<nmom;p++){
      std::vector<TensorType> matsum_globT(global_T);
      for(int global_t=0;global_t<global_T;global_t++){
	TensorType &m = matsum_globT[global_t];
	  
	int local_t = global_t - local_toff;
	if(local_t >= 0 && local_t < local_T) m = into[idx]->at(p)[local_t];
	else{ m = 0.0; }
	  
	_PropagatorQuadrilinear_helper<MatrixType>::lattice_sum(m);
      }
      into[idx]->at(p) = matsum_globT;
    }
  }

}



//Get a Fourier transformed quadrilinear correlation function as a function of time
//Sigma is ignored if MatrixType is WilsonMatrix
//quadrilinear form is    A Gamma1 Sigma1 B  \otimes  C Gamma2 Sigma2 D        where \otimes is an outer product and A,B,C,D are propagators
template<typename MatrixType>
const std::vector<typename _PropagatorQuadrilinear_helper<MatrixType>::TensorType> &PropagatorQuadrilinear<MatrixType>::getQuadrilinear(Lattice &lat,
										   const std::vector<Float> &sink_momentum, 
										   char const* tag_A, const PropSuperscript &ss_A,  
										   char const* tag_B, const PropSuperscript &ss_B, 
										   char const* tag_C, const PropSuperscript &ss_C,  
										   char const* tag_D, const PropSuperscript &ss_D, 
										   const int &Gamma1, const int &Sigma1,
										   const int &Gamma2, const int &Sigma2){				   
  if(!mom_idx_map.count( sink_momentum )) ERR.General("PropagatorQuadrilinear","getQuadrilinear","Desired sink momentum is not in the list of those momentum components generated");
  int momentum_idx = mom_idx_map[sink_momentum];     
  prop_info_pair props_1( prop_info(tag_A,ss_A), prop_info(tag_B,ss_B) );
  prop_info_pair props_2( prop_info(tag_C,ss_C), prop_info(tag_D,ss_D) );
  prop_info_quad props(props_1,props_2);

  int scf_idx1 = _PropagatorBilinear_helper<MatrixType>::scf_map(Gamma1,Sigma1);
  int scf_idx2 = _PropagatorBilinear_helper<MatrixType>::scf_map(Gamma2,Sigma2);
  int idx = scf_idx1 + nmat*scf_idx2;

  std::pair<prop_info_quad,int> key(props,idx);
  if(!all_mats.count(key)) calcQuadrilinear(idx,props,lat);
  return all_mats[key][momentum_idx];
}

//for use with WilsonMatrix where sigma (flavour matrix idx) does not play a role
template<typename MatrixType>
const std::vector<typename _PropagatorQuadrilinear_helper<MatrixType>::TensorType> & PropagatorQuadrilinear<MatrixType>::getQuadrilinear(Lattice &lat,
										    const std::vector<Float> &sink_momentum, 
										    char const* tag_A, const PropSuperscript &ss_A,  
										    char const* tag_B, const PropSuperscript &ss_B, 
										    char const* tag_C, const PropSuperscript &ss_C,  
										    char const* tag_D, const PropSuperscript &ss_D, 
										    const int &Gamma1, const int &Gamma2
										    ){
  return getQuadrilinear(lat,sink_momentum,tag_A,ss_A,tag_B,ss_B,tag_C,ss_C,tag_D,ss_D,Gamma1,0,Gamma2,0);
}
    
template<typename MatrixType>
void PropagatorQuadrilinear<MatrixType>::add_momentum(std::vector<Float> sink_mom){
  //If array_size is set we cannot add any more momenta without having to do lots of extra work to fill in gaps for existing bilinears
  if(nmom_fixed) ERR.General("PropagatorQuadrilinear","add_momentum","Cannot add momentum after vertices have begun being calculated"); 
  return PropDFT::add_momentum(sink_mom);
}

template<typename MatrixType>
void PropagatorQuadrilinear<MatrixType>::clear(){
  all_mats.clear();
  nmom_fixed = false;
  PropDFT::clear();
}

template<typename MatrixType>
void PropagatorQuadrilinear<MatrixType>::write(char const* tag_A, const PropSuperscript &ss_A,  
					       char const* tag_B, const PropSuperscript &ss_B,
					       char const* tag_C, const PropSuperscript &ss_C,  
					       char const* tag_D, const PropSuperscript &ss_D, 
					       const int &Gamma1, const int &Sigma1,
					       const int &Gamma2, const int &Sigma2,
					       const char *file, Lattice &lat){
  FILE *fp;
  if ((fp = Fopen(file, "w")) == NULL) {
    ERR.FileW("PropagatorQuadrlinear","write(...) %s",file);
  }
  write(tag_A,ss_A,tag_B,ss_B,tag_C,ss_C,tag_D,ss_D,Gamma1,Sigma1,Gamma2,Sigma2,fp,lat);
  Fclose(fp);
}

template<typename MatrixType>
void PropagatorQuadrilinear<MatrixType>::write(char const* tag_A, const PropSuperscript &ss_A,  
					       char const* tag_B, const PropSuperscript &ss_B,
					       char const* tag_C, const PropSuperscript &ss_C,  
					       char const* tag_D, const PropSuperscript &ss_D,
					       const int &Gamma1, const int &Sigma1,
					       const int &Gamma2, const int &Sigma2,
					       FILE *fp, Lattice &lat){
  static const char* fname = "write(....)";
  prop_info_pair props_1( prop_info(tag_A,ss_A), prop_info(tag_B,ss_B) );
  prop_info_pair props_2( prop_info(tag_C,ss_C), prop_info(tag_D,ss_D) );
  prop_info_quad props(props_1,props_2);

  int scf_idx1 = _PropagatorBilinear_helper<MatrixType>::scf_map(Gamma1,Sigma1);
  int scf_idx2 = _PropagatorBilinear_helper<MatrixType>::scf_map(Gamma2,Sigma2);
  int idx = scf_idx1 + nmat*scf_idx2;

  std::pair<prop_info_quad,int> key(props,idx);
  if(!all_mats.count(key)) calcQuadrilinear(idx,props,lat);

  //Find all p^2 and sort
  std::vector< std::pair<Float,int> > p2list;
  std::map<int,std::vector<Float> > p2map;
  find_p2sorted(p2list,p2map);

  const std::vector<std::vector<TensorType> > &data = all_mats[key];

  for(int p=0;p<p2list.size();p++){
    int pidx = p2list[p].second;
    const std::vector<Float> &mom = p2map[pidx];
    const Float &p2 = p2list[p].first;
      
    for(int t=0;t<GJP.TnodeSites()*GJP.Tnodes();t++){
      _PropagatorQuadrilinear_helper<MatrixType>::write(fp,data[pidx][t], scf_idx1, scf_idx2, p2, mom, t);
    }
  }
}

template<typename MatrixType>
void PropagatorQuadrilinear<MatrixType>::write(char const* tag_A, const PropSuperscript &ss_A,  
					       char const* tag_B, const PropSuperscript &ss_B,
					       char const* tag_C, const PropSuperscript &ss_C,  
					       char const* tag_D, const PropSuperscript &ss_D, 
					       const char *file, Lattice &lat){
  FILE *fp;
  if ((fp = Fopen(file, "w")) == NULL) {
    ERR.FileW("PropagatorQuadrilinear","write(...) %s",file);
  }
  write(tag_A,ss_A,tag_B,ss_B,tag_C,ss_C,tag_D,ss_D,fp,lat);
  Fclose(fp);
}

template<typename MatrixType>
void PropagatorQuadrilinear<MatrixType>::write(char const* tag_A, const PropSuperscript &ss_A,  
					       char const* tag_B, const PropSuperscript &ss_B,
					       char const* tag_C, const PropSuperscript &ss_C,  
					       char const* tag_D, const PropSuperscript &ss_D,
					       FILE *fp, Lattice &lat){
  static const char* fname = "write(....)";
  prop_info_pair props_1( prop_info(tag_A,ss_A), prop_info(tag_B,ss_B) );
  prop_info_pair props_2( prop_info(tag_C,ss_C), prop_info(tag_D,ss_D) );
  prop_info_quad props(props_1,props_2);

  bool do_calc(false);
  for(int i=0;i<nmat*nmat;i++){ 
    std::pair<prop_info_quad,int> key(props,i);
    if(!all_mats.count(key)){ do_calc = true; break; }
  }
  if(do_calc) calcAllQuadrilinears(props,lat);

  //Find all p^2 and sort
  std::vector< std::pair<Float,int> > p2list;
  std::map<int,std::vector<Float> > p2map;
  find_p2sorted(p2list,p2map);

  for(int i=0;i<nmat*nmat;i++){
    std::pair<prop_info_quad,int> key(props,i);
    const std::vector<std::vector<TensorType> > &data = all_mats[key];
      
    int scf_idx2 = i;
    int scf_idx1 = scf_idx2 % nmat; scf_idx2/=nmat;

    for(int p=0;p<p2list.size();p++){
      int pidx = p2list[p].second;
      const std::vector<Float> &mom = p2map[pidx];
      const Float &p2 = p2list[p].first;
      
      for(int t=0;t<GJP.TnodeSites()*GJP.Tnodes();t++){
	_PropagatorQuadrilinear_helper<MatrixType>::write(fp,data[pidx][t], scf_idx1, scf_idx2, p2, mom, t);
      }
    }
  }
}

template<typename MatrixType>
int ContractedWallSinkBilinearSpecMomentum<MatrixType>::idx_map(const int &mat1, const int &mat2, const int & mom_idx, const int &t) const{    
  return mat1 + nmat * (mat2 + nmat*(mom_idx + this->nmompairs*t));
}

template<typename MatrixType>
void ContractedWallSinkBilinearSpecMomentum<MatrixType>::idx_unmap(const int &idx, int &mat1, int &mat2, int & mom_idx, int &t) const{
  int rem = idx;
  mat1 = rem % nmat;  rem /= nmat;
  mat2 = rem % nmat;  rem /= nmat;
  mom_idx = rem % this->nmompairs; rem /= this->nmompairs;
  t = rem;
}

template<typename MatrixType>
void ContractedWallSinkBilinearSpecMomentum<MatrixType>::do_superscript(MatrixType &mat, const PropSuperscript &ss){
  return PropDFT::do_superscript<MatrixType>(mat,ss);
}


template<typename MatrixType>
void ContractedWallSinkBilinearSpecMomentum<MatrixType>::enableCosineSink(){
  if(array_size!=-1) ERR.General("ContractedWallSinkBilinearSpecMomentum","enableCosineSink","Cannot enable cosine sink after bilinears have begun being calculated"); 
  if(nmompairs!=0) ERR.General("ContractedWallSinkBilinearSpecMomentum","enableCosineSink","Must enable momentum sink prior to adding momentum"); 
  cosine_sink = true;
  fprop.enableCosineSink();
}

inline bool OpFlipsMomentumPhase(const PropSuperscript ss){
  switch(ss){
  case OpNone:
  case OpTranspose:
  case OpConjFlipMomentum:
  case OpDaggerFlipMomentum:
    return false;
  case OpConj:
  case OpFlipMomentum:
  case OpDagger:
  case OpTransposeFlipMomentum:
    return true;
  default:
    ERR.General("","OpFlipsMomentumPhase","Unknown operation\n");
  }
}

template<typename MatrixType>
void ContractedWallSinkBilinearSpecMomentum<MatrixType>::calcAllContractedBilinears(const prop_info_pair &props, Lattice &lat){
  Timer::reset(); Elapsed time;
  time = Timer::relative_time(); time.print("ContractedWallSinkBilinearSpecMomentum::calcAllContractedBilinears starting");

  int global_T = GJP.Tnodes()*GJP.TnodeSites();
  if(array_size==-1) array_size = nmat*nmat*nmompairs*global_T; //also acts as a lock to prevent further momenta from being added

  if(!results.count(props)) results[props] = new Rcomplex[array_size];
  Rcomplex *into = results[props];
  for(int i=0;i<array_size;i++) into[i] = 0.0;

  time = Timer::relative_time(); time.print("Starting outer mom loop");

  for(mom_pair_idx_map_type::iterator mom_it = mom_pair_idx_map.begin(); mom_it != mom_pair_idx_map.end(); ++mom_it){
    const std::vector<Float> & mom1 = mom_it->first.first;
    std::vector<Float> minus_mom1(mom1); for(int i=0;i<3;i++) minus_mom1[i]*=-1;
    const std::vector<Float> & mom2 = mom_it->first.second;
    std::vector<Float> minus_mom2(mom2); for(int i=0;i<3;i++) minus_mom2[i]*=-1;

    const int &p_vec_pos = mom_it->second;
    std::vector<MatrixType> prop_1;
    std::vector<MatrixType> prop_2;

    time = Timer::relative_time(); time.print("Starting computation of FT props");
    if(cosine_sink){ //fprop is put in cosine sink mode when enableCosineSink() is called
      prop_1 = fprop.getFTProp(lat,mom1,props.first.first.c_str());
      prop_2 = fprop.getFTProp(lat,mom2,props.second.first.c_str());
    }else{
      //We choose a convention that the momentum supplied for the propagators is applied after the PropSuperscript operation.
      //As the Fourier transform is applied prior to the operation being performed, and some of the operations flip the momentum, we need to choose the opposite phase for the FT
      if(OpFlipsMomentumPhase(props.first.second)) prop_1 = fprop.getFTProp(lat,minus_mom1,props.first.first.c_str());
      else prop_1 = fprop.getFTProp(lat,mom1,props.first.first.c_str());

      if(OpFlipsMomentumPhase(props.second.second)) prop_2 = fprop.getFTProp(lat,minus_mom2,props.second.first.c_str());
      else prop_2 = fprop.getFTProp(lat,mom2,props.second.first.c_str());
    }
	
    for(int t=0;t<global_T;t++){ //apply superscripts
      do_superscript(prop_1[t], props.first.second);
      do_superscript(prop_2[t], props.second.second);
    }
    time = Timer::relative_time(); time.print("Starting contraction loop");
	
    //loop over outer (gamma1,sigma1) and inner (gamma2,sigma2) mats, t and do trace
    int ntrace = nmat*nmat*global_T; 

#pragma omp parallel for default(shared)
    for(int i=0;i<ntrace;i++){
      int rem = i;
      int mat1 = rem % nmat; rem/=nmat;
      int mat2 = rem % nmat; rem/=nmat;
      int t = rem;

      std::pair<int,int> spin_flav_1 = _PropagatorBilinear_helper<MatrixType>::unmap(mat1);
      std::pair<int,int> spin_flav_2 = _PropagatorBilinear_helper<MatrixType>::unmap(mat2);

      //tr( A G_1 B G_2 ) = tr( G_1 B G_2 A )
      MatrixType tmp1 = prop_1[t]; _PropagatorBilinear_helper<MatrixType>::rmult_matrix(tmp1, spin_flav_2);
      MatrixType tmp2 = prop_2[t]; _PropagatorBilinear_helper<MatrixType>::rmult_matrix(tmp2, spin_flav_1);

      int tr_off = idx_map(mat1,mat2,p_vec_pos,t);
      into[tr_off] = Trace(tmp1,tmp2);
    }
    time = Timer::relative_time(); time.print("Finished mom comb");
  }
  time = Timer::relative_time(); time.print("Finished contraction");
}

template<typename MatrixType>
void ContractedWallSinkBilinearSpecMomentum<MatrixType>::calculateBilinears(Lattice &lat,
									    const prop_info_pair &props){
  if(!results.count(props)){
    calcAllContractedBilinears(props,lat);
  }
}

template<typename MatrixType>
bool ContractedWallSinkBilinearSpecMomentum<MatrixType>::mom_sort_pred(const std::pair<Float,int> &i, const std::pair<Float,int> &j){ 
  return (i.first<j.first);
}

template<typename MatrixType>
void ContractedWallSinkBilinearSpecMomentum<MatrixType>::find_p2sorted(std::vector< std::pair<Float,int> > &p2list, std::map<int,std::pair<std::vector<Float>,std::vector<Float> > > &p2map){
  //Find all p^2 and sort, output into p2list and p2map
  p2list.clear(); p2map.clear(); p2list.reserve(nmompairs);
  for(mom_pair_idx_map_type::iterator mom_it = mom_pair_idx_map.begin(); mom_it != mom_pair_idx_map.end(); ++mom_it){
    const std::vector<Float> & mom1 = mom_it->first.first;
    const std::vector<Float> & mom2 = mom_it->first.second;
    const int &vec_pos = mom_it->second;

    //total momentum
    std::vector<Float> mom(3); for(int i=0;i<3;i++) mom[i] = mom1[i]+mom2[i];

    std::pair<Float,int> p2idx( mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2], vec_pos );
    p2list.push_back(p2idx);
    p2map[vec_pos] = mom_it->first;
  }
  std::sort(p2list.begin(),p2list.end(),mom_sort_pred);
}

template<typename MatrixType>
void ContractedWallSinkBilinearSpecMomentum<MatrixType>::calculateBilinears(Lattice &lat,
									    char const* tag_A, const PropSuperscript &ss_A,  
									    char const* tag_B, const PropSuperscript &ss_B
									    ){
  prop_info_pair props( prop_info(tag_A,ss_A), prop_info(tag_B,ss_B) );

  calculateBilinears(lat,props);
}


//Get a Fourier transformed wall sink bilinear correlation function as a function of time
//Sigma is ignored if MatrixType is WilsonMatrix
template<typename MatrixType>
std::vector<Rcomplex> ContractedWallSinkBilinearSpecMomentum<MatrixType>::getBilinear(Lattice &lat,
										      const std::pair<std::vector<Float>,std::vector<Float> > &sink_momenta, 
										      char const* tag_A, const PropSuperscript &ss_A,  
										      char const* tag_B, const PropSuperscript &ss_B, 
										      const int &Gamma1, const int &Sigma1,
										      const int &Gamma2, const int &Sigma2
										      ){
  if(!this->mom_pair_idx_map.count( sink_momenta )) ERR.General("ContractedWallSinkBilinearSpecMomentum","getBilinear","Desired sink momentum combination is not in the list of those generated");
  int momentum_idx = this->mom_pair_idx_map[sink_momenta];     
  prop_info_pair props( prop_info(tag_A,ss_A), prop_info(tag_B,ss_B) );
  calculateBilinears(lat,props);

  int scf_idx1 = _PropagatorBilinear_helper<MatrixType>::scf_map(Gamma1,Sigma1);
  int scf_idx2 = _PropagatorBilinear_helper<MatrixType>::scf_map(Gamma2,Sigma2);

  Rcomplex *con = results[props];
  std::vector<Rcomplex> out(GJP.TnodeSites()*GJP.Tnodes());
  for(int t=0;t<out.size();t++){
    out[t] = con[ idx_map(scf_idx1,scf_idx2,momentum_idx,t) ];
  }
  return out;
}

//for use with WilsonMatrix where sigma (flavour matrix idx) does not play a role
template<typename MatrixType>
std::vector<Rcomplex> ContractedWallSinkBilinearSpecMomentum<MatrixType>::getBilinear(Lattice &lat,
										      const std::pair<std::vector<Float>,std::vector<Float> > &sink_momenta,
										      char const* tag_A, const PropSuperscript &ss_A,  
										      char const* tag_B, const PropSuperscript &ss_B, 
										      const int &Gamma1, const int &Gamma2
										      ){
  return getBilinear(lat,sink_momenta,tag_A,ss_A,tag_B,ss_B,Gamma1,0,Gamma2,0);
}
    
template<typename MatrixType>
void ContractedWallSinkBilinearSpecMomentum<MatrixType>::add_momentum(const std::pair<std::vector<Float>,std::vector<Float> > &sink_momenta){
  //If array_size is set we cannot add any more momenta without having to do lots of extra work to fill in gaps for existing bilinears
  if(array_size!=-1) ERR.General("ContractedWallSinkBilinearSpecMomentum","add_momentum","Cannot add momentum after bilinears have begun being calculated"); 

  mom_pair_idx_map[ sink_momenta ] = nmompairs++;
  fprop.add_momentum(sink_momenta.first);
  fprop.add_momentum(sink_momenta.second);
}

template<typename MatrixType>
void ContractedWallSinkBilinearSpecMomentum<MatrixType>::clear(){
  for(typename map_info_type::iterator it = results.begin(); it!= results.end(); ++it){
    delete[] it->second;
  }
  results.clear();
  array_size = -1;
  mom_pair_idx_map.clear();
  nmompairs = 0;
  cosine_sink = false;
}

template<typename MatrixType>
void ContractedWallSinkBilinearSpecMomentum<MatrixType>::write(char const* tag_A, const PropSuperscript &ss_A,  
							       char const* tag_B, const PropSuperscript &ss_B, 
							       const int &Gamma1, const int &Sigma1,
							       const int &Gamma2, const int &Sigma2,
							       const char *file, Lattice &lat){
#ifndef USE_OFSTREAM
  FILE *fp;
  if ((fp = Fopen(file, "w")) == NULL) {
    ERR.FileW("ContractedWallSinkBilinearSpecMomentum","write(...)",file);
  }
  write(tag_A,ss_A,tag_B,ss_B, Gamma1, Sigma1, Gamma2, Sigma2, fp, lat);
  Fclose(fp);
#else
  std::ostream *str = _DFTutils::open_ofstream(file);
  write(tag_A,ss_A,tag_B,ss_B, Gamma1, Sigma1, Gamma2, Sigma2, *str, lat);
  _DFTutils::close_ofstream(str);
#endif
}

template<typename MatrixType>
template<typename OutputType>
void ContractedWallSinkBilinearSpecMomentum<MatrixType>::_writeit(OutputType into,Rcomplex *con,const int &scf_idx1, const int &scf_idx2, 
								  const std::vector< std::pair<Float,int> > &p2list, const std::map<int,std::pair<std::vector<Float>,std::vector<Float> > > &p2map, const bool &binary){
  for(int p=0;p<p2list.size();p++){
    const Float &p2 = p2list[p].first;
    int pidx = p2list[p].second;

    std::map<int,std::pair<std::vector<Float>,std::vector<Float> > >::const_iterator map_loc = p2map.find(pidx);

    const std::vector<Float> &mom1 = map_loc->second.first;
    const std::vector<Float> &mom2 = map_loc->second.second;

    for(int t=0;t<GJP.TnodeSites()*GJP.Tnodes();t++){
      int off = idx_map(scf_idx1,scf_idx2,pidx,t);
      const Rcomplex &val = con[ off ];
      binary ?
	_ContractedWallSinkBilinearSpecMomentum_helper<MatrixType>::binary_write(into,val,scf_idx1,scf_idx2,p2, mom1, mom2, t) :
	_ContractedWallSinkBilinearSpecMomentum_helper<MatrixType>::write(into,val,scf_idx1,scf_idx2,p2, mom1, mom2, t);
    }
  }
}

template<typename MatrixType>
template<typename OutputType>
void ContractedWallSinkBilinearSpecMomentum<MatrixType>::write(char const* tag_A, const PropSuperscript &ss_A,  
							       char const* tag_B, const PropSuperscript &ss_B, 
							       const int &Gamma1, const int &Sigma1,
							       const int &Gamma2, const int &Sigma2,
							       OutputType into, Lattice &lat){
  prop_info_pair props( prop_info(tag_A,ss_A), prop_info(tag_B,ss_B) );

  calculateBilinears(lat,props); //only calculates if not yet done

  int scf_idx1 = _PropagatorBilinear_helper<MatrixType>::scf_map(Gamma1,Sigma1);
  int scf_idx2 = _PropagatorBilinear_helper<MatrixType>::scf_map(Gamma2,Sigma2);

  Rcomplex *con = results[props];

  //Find all total sink p^2 and sort by the quark momenta
  std::vector< std::pair<Float,int> > p2list;
  std::map<int,std::pair<std::vector<Float>,std::vector<Float> > > p2map;
  find_p2sorted(p2list,p2map);
  _writeit(into,con,scf_idx1,scf_idx2,p2list,p2map);
}
  
//write all combinations
template<typename MatrixType>
void ContractedWallSinkBilinearSpecMomentum<MatrixType>::write(char const* tag_A, const PropSuperscript &ss_A,  
							       char const* tag_B, const PropSuperscript &ss_B,
							       const std::string &file, Lattice &lat, const bool &binary){
  if(!UniqueID()) printf("ContractedWallSinkBilinearSpecMomentum writing to file \"%s\"\n",file.c_str());
  if(!UniqueID() && binary) printf("ContractedWallSinkBilinearSpecMomentum binary mode enabled\n");
#ifndef USE_OFSTREAM
  FILE *fp;
  if ((fp = Fopen(file.c_str(), binary ? "wb" : "w")) == NULL) {
    ERR.FileW("ContractedWallSinkBilinearSpecMomentum","write(...)",file.c_str());
  }
  if(binary && !UniqueID()){
    //To make sure that the file is read correctly later (correct endedness, etc), write a known number at the beginning which can be checked
    const static double test = 3.141592654;
    fwrite(&test, sizeof(double), 1, fp);
  }
  write(tag_A,ss_A,tag_B,ss_B,fp, lat, binary);
  Fclose(fp);
#else
  if(!UniqueID()) printf("ContractedWallSinkBilinearSpecMomentum File write using stream version\n");
  if(binary){
    ERR.General("ContractedWallSinkBilinearSpecMomentum","write(...)","Binary write not implemented for stream version");
  }
  std::ostream *str = _DFTutils::open_ofstream(file.c_str());
  write(tag_A,ss_A,tag_B,ss_B,*str, lat);
  _DFTutils::close_ofstream(str);
#endif

}

template<typename MatrixType>
template<typename OutputType>
void ContractedWallSinkBilinearSpecMomentum<MatrixType>::write(char const* tag_A, const PropSuperscript &ss_A,  
							       char const* tag_B, const PropSuperscript &ss_B, 
							       OutputType into, Lattice &lat, const bool &binary){
  prop_info_pair props( prop_info(tag_A,ss_A), prop_info(tag_B,ss_B) );

  calculateBilinears(lat,props); //only calculates if not yet done

  Timer::reset(); Elapsed time;
  time = Timer::relative_time(); time.print("Starting file write ");

  Rcomplex *con = results[props];

  //Find all p^2 and sort
  std::vector< std::pair<Float,int> > p2list;
  std::map<int,std::pair<std::vector<Float>,std::vector<Float> > > p2map;
  find_p2sorted(p2list,p2map);

  for(int scf_idx1 = 0; scf_idx1 < nmat; scf_idx1++){
    for(int scf_idx2 = 0; scf_idx2 < nmat; scf_idx2++){
      _writeit(into,con,scf_idx1,scf_idx2,p2list,p2map,binary);
    }
  }
  time = Timer::relative_time(); time.print("Finished file write ");
}
