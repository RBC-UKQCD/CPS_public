#ifndef CPS_FIELD_UTILS_H
#define CPS_FIELD_UTILS_H

CPS_START_NAMESPACE

inline void compareFermion(const CPSfermion5D<ComplexD> &A, const CPSfermion5D<ComplexD> &B, const std::string &descr = "Ferms", const double tol = 1e-9){
  double fail = 0.;
  for(int i=0;i<GJP.VolNodeSites()*GJP.SnodeSites();i++){
    int x[5]; int rem = i;
    for(int ii=0;ii<5;ii++){ x[ii] = rem % GJP.NodeSites(ii); rem /= GJP.NodeSites(ii); }
    
    for(int f=0;f<GJP.Gparity()+1;f++){
      for(int sc=0;sc<24;sc++){
	double vbfm = *((double*)A.site_ptr(i,f) + sc);
	double vgrid = *((double*)B.site_ptr(i,f) + sc);
	    
	double diff_rat = fabs( 2.0 * ( vbfm - vgrid )/( vbfm + vgrid ) );
	double rat_grid_bfm = vbfm/vgrid;
	if(vbfm == 0.0 && vgrid == 0.0){ diff_rat = 0.;	 rat_grid_bfm = 1.; }
	if( (vbfm == 0.0 && fabs(vgrid) < 1e-50) || (vgrid == 0.0 && fabs(vbfm) < 1e-50) ){ diff_rat = 0.;	 rat_grid_bfm = 1.; }

	if(diff_rat > tol){
	  printf("Fail: (%d,%d,%d,%d,%d; %d; %d) A %g B %g rat_A_B %g fracdiff %g\n",x[0],x[1],x[2],x[3],x[4],f,sc,vbfm,vgrid,rat_grid_bfm,diff_rat);
	  fail = 1.0;
	}//else printf("Pass: (%d,%d,%d,%d,%d; %d; %d) A %g B %g rat_A_B %g fracdiff %g\n",x[0],x[1],x[2],x[3],x[4],f,sc,vbfm,vgrid,rat_grid_bfm,diff_rat);
      }
    }
  }
  glb_max(&fail);
  
  if(fail!=0.0){
    if(!UniqueID()){ printf("Failed %s check\n", descr.c_str()); fflush(stdout); } 
    exit(-1);
  }else{
    if(!UniqueID()){ printf("Passed %s check\n", descr.c_str()); fflush(stdout); }
  }
}

template<typename FieldType, typename my_enable_if<_equal<typename ComplexClassify<typename FieldType::FieldSiteType>::type, complex_double_or_float_mark>::value,int>::type = 0>
inline void compareField(const FieldType &A, const FieldType &B, const std::string &descr = "Field", const double tol = 1e-9, bool print_all = false){
  typedef typename FieldType::FieldSiteType::value_type value_type;
  
  double fail = 0.;
  for(int xf=0;xf<A.nfsites();xf++){
    int f; int x[FieldType::FieldDimensionPolicy::EuclideanDimension];
    A.fsiteUnmap(xf, x,f);

    for(int i=0;i<FieldType::FieldSiteSize;i++){
      value_type const* av = (value_type const*)(A.fsite_ptr(xf)+i);
      value_type const* bv = (value_type const*)(B.fsite_ptr(xf)+i);
      for(int reim=0;reim<2;reim++){
	value_type diff_rat = (av[reim] == 0.0 && bv[reim] == 0.0) ? 0.0 : fabs( 2.*(av[reim]-bv[reim])/(av[reim]+bv[reim]) );
	if(diff_rat > tol || print_all){
	  if(!print_all) std::cout << "Fail: (";
	  else std::cout << "Pass: (";
	  
	  for(int xx=0;xx<FieldType::FieldDimensionPolicy::EuclideanDimension-1;xx++)
	    std::cout << x[xx] << ", ";
	  std::cout << x[FieldType::FieldDimensionPolicy::EuclideanDimension-1];

	  std::cout << ") f=" << f << " reim " << reim << " A " << av[reim] << " B " << bv[reim] << " fracdiff " << diff_rat << std::endl;
	  if(!print_all) fail = 1.;
	}
      }
    }
  }
  glb_max(&fail);
  
  if(fail!=0.0){
    if(!UniqueID()){ printf("Failed %s check\n", descr.c_str()); fflush(stdout); } 
    exit(-1);
  }else{
    if(!UniqueID()){ printf("Passed %s check\n", descr.c_str()); fflush(stdout); }
  }
}








#ifdef USE_BFM
inline void exportBFMcb(CPSfermion5D<ComplexD> &into, Fermion_t from, bfm_evo<double> &dwf, int cb, bool singleprec_evec = false){
  Fermion_t zero_a = dwf.allocFermion();
#pragma omp parallel
  {   
    dwf.set_zero(zero_a); 
  }
  Fermion_t etmp = dwf.allocFermion(); 
  Fermion_t tmp[2];
  tmp[!cb] = zero_a;
  if(singleprec_evec){
    const int len = 24 * dwf.node_cbvol * (1 + dwf.gparity) * dwf.cbLs;
#pragma omp parallel for
    for(int j = 0; j < len; j++) {
      ((double*)etmp)[j] = ((float*)(from))[j];
    }
    tmp[cb] = etmp;
  }else tmp[cb] = from;

  dwf.cps_impexFermion(into.ptr(),tmp,0);
  dwf.freeFermion(zero_a);
  dwf.freeFermion(etmp);
}
#endif

#ifdef USE_GRID
template<typename GridPolicies>
inline void exportGridcb(CPSfermion5D<ComplexD> &into, typename GridPolicies::GridFermionField &from, typename GridPolicies::FgridFclass &latg){
  Grid::GridCartesian *FGrid = latg.getFGrid();
  typename GridPolicies::GridFermionField tmp_g(FGrid);
  tmp_g = Grid::zero;

  setCheckerboard(tmp_g, from);
  latg.ImportFermion((Vector*)into.ptr(), tmp_g);
}
#endif

#ifdef USE_QMP

//Cyclic permutation of *4D* CPSfield with std::complex type and FourDpolicy dimension policy
//Conventions are direction of *data flow*: For shift n in direction +1   f'(x) = f(x-\hat i)  so data is sent in the +x direction.

#define CONDITION _equal<typename ComplexClassify<mf_Complex>::type, complex_double_or_float_mark>::value && (_equal<DimensionPolicy,FourDpolicy>::value || _equal<DimensionPolicy,SpatialPolicy>::value)

template< typename mf_Complex, int SiteSize, typename DimensionPolicy, typename FlavorPolicy, typename AllocPolicy>
void cyclicPermute(CPSfield<mf_Complex,SiteSize,DimensionPolicy,FlavorPolicy,AllocPolicy> &to, const CPSfield<mf_Complex,SiteSize,DimensionPolicy,FlavorPolicy,AllocPolicy> &from,
		   const int dir, const int pm, const int n,
		   typename my_enable_if<CONDITION , const int>::type dummy = 0){
  enum {Dimension = DimensionPolicy::EuclideanDimension};
  assert(dir < Dimension);
  assert(n < GJP.NodeSites(dir));
  assert(pm == 1 || pm == -1);
	   
  if(&to == &from){
    if(n==0) return;    
    CPSfield<mf_Complex,SiteSize,DimensionPolicy,FlavorPolicy,AllocPolicy> tmpfrom(from);
    return cyclicPermute(to,tmpfrom,dir,pm,n);
  }
  if(n == 0){
    to = from;
    return;
  }

  QMP_barrier();
  
  //Prepare face to send. If we send in the + direction we need to collect the slice starting {L-n ... L-1} (inclusive), and if we send in the - dir we collect the slice {0... n-1}
  int bsites = n; //sites on boundary
  int bsizes[Dimension]; bsizes[dir] = n;
  int boff[Dimension]; boff[dir] = (pm == 1 ? GJP.NodeSites(dir)-n : 0);
		 
  for(int i=0;i<Dimension;i++)
    if(i != dir){
      bsizes[i] = GJP.NodeSites(i);
      bsites *= bsizes[i];
      boff[i] = 0;
    }
  int flav_off = from.flav_offset();
  int nf = from.nflavors();
  
  int bufsz = bsites * SiteSize * nf;
  int halfbufsz = bufsz/2;

  QMP_mem_t *recv_mem = QMP_allocate_memory(bufsz * sizeof(mf_Complex));
  mf_Complex *recv_buf = (mf_Complex *)QMP_get_memory_pointer(recv_mem);

  QMP_mem_t *send_mem = QMP_allocate_memory(bufsz * sizeof(mf_Complex));
  mf_Complex *send_buf = (mf_Complex *)QMP_get_memory_pointer(send_mem);

#pragma omp parallel for
  for(int i=0;i<bsites;i++){
    int rem = i;
    int coor[Dimension];
    for(int d=0;d<Dimension;d++){ coor[d] = rem % bsizes[d] + boff[d]; rem/=bsizes[d]; }

    mf_Complex const* site_ptr = from.site_ptr(coor);
    mf_Complex* bp = send_buf + i*SiteSize;
    memcpy(bp,site_ptr,SiteSize*sizeof(mf_Complex));
    if(nf == 2){
      site_ptr += flav_off;
      bp += halfbufsz;
      memcpy(bp,site_ptr,SiteSize*sizeof(mf_Complex));
    }
  }
  QMP_barrier();
 
  //Copy remaining sites from on-node data with shift
  int rsizes[Dimension]; rsizes[dir] = GJP.NodeSites(dir) - n;
  int rsites = GJP.NodeSites(dir) - n;
  //if we sent in the + direction we need to shift the remaining L-n sites {0...L-n-1} forwards by n to make way for a new slice at the left side
  //if we sent in the - direction we need to shift the remaining L-n sites {n ... L-1} backwards by n to make way for a new slice at the right side
  
  int roff[Dimension]; roff[dir] = (pm == 1 ? 0 : n);  
  for(int i=0;i<Dimension;i++)
    if(i != dir){
      rsizes[i] = GJP.NodeSites(i);
      rsites *= rsizes[i];
      roff[i] = 0;
    }

#pragma omp parallel for
  for(int i=0;i<rsites;i++){
    int rem = i;
    int from_coor[Dimension];
    for(int d=0;d<Dimension;d++){ from_coor[d] = rem % rsizes[d] + roff[d]; rem/=rsizes[d]; }
    
    int to_coor[Dimension]; memcpy(to_coor,from_coor,Dimension*sizeof(int));
    to_coor[dir] = (pm == +1 ? from_coor[dir] + n : from_coor[dir] - n);
    
    mf_Complex const* from_ptr = from.site_ptr(from_coor);
    mf_Complex * to_ptr = to.site_ptr(to_coor);

    memcpy(to_ptr,from_ptr,SiteSize*sizeof(mf_Complex));
    if(nf == 2){
      from_ptr += flav_off;
      to_ptr += flav_off;
      memcpy(to_ptr,from_ptr,SiteSize*sizeof(mf_Complex));
    }
  }
  
  //Send/receive
  QMP_msgmem_t send_msg = QMP_declare_msgmem(send_buf,bufsz * sizeof(mf_Complex));
  QMP_msgmem_t recv_msg = QMP_declare_msgmem(recv_buf,bufsz * sizeof(mf_Complex));
  
  QMP_msghandle_t send = QMP_declare_send_relative(send_msg, dir, pm, 0);
  QMP_msghandle_t recv = QMP_declare_receive_relative(recv_msg, dir, -pm, 0);
  QMP_start(recv);
  QMP_start(send);
  
  QMP_status_t send_status = QMP_wait(send);
  if (send_status != QMP_SUCCESS) 
    QMP_error("Send failed in cyclicPermute: %s\n", QMP_error_string(send_status));
  QMP_status_t rcv_status = QMP_wait(recv);
  if (rcv_status != QMP_SUCCESS) 
    QMP_error("Receive failed in PassDataT: %s\n", QMP_error_string(rcv_status));

  //Copy received face into position. For + shift the origin we copy into is the left-face {0..n-1}, for a - shift its the right-face {L-n .. L-1}
  boff[dir] = (pm == 1 ? 0 : GJP.NodeSites(dir)-n);
#pragma omp parallel for
  for(int i=0;i<bsites;i++){
    int rem = i;
    int coor[Dimension];
    for(int d=0;d<Dimension;d++){ coor[d] = rem % bsizes[d] + boff[d]; rem/=bsizes[d]; }
    
    mf_Complex * site_ptr = to.site_ptr(coor);
    mf_Complex const* bp = recv_buf + i*SiteSize;
    memcpy(site_ptr,bp,SiteSize*sizeof(mf_Complex));
    if(nf == 2){
      site_ptr += flav_off;
      bp += halfbufsz;
      memcpy(site_ptr,bp,SiteSize*sizeof(mf_Complex));
    }
  }

  QMP_free_msghandle(send);
  QMP_free_msghandle(recv);
  QMP_free_msgmem(send_msg);
  QMP_free_msgmem(recv_msg);
  QMP_free_memory(send_mem);
  QMP_free_memory(recv_mem);
  QMP_barrier();
}
#undef CONDITION

# ifdef USE_GRID

#define CONDITION _equal<typename ComplexClassify<mf_Complex>::type, grid_vector_complex_mark>::value && (_equal<DimensionPolicy,FourDSIMDPolicy>::value || _equal<DimensionPolicy,ThreeDSIMDPolicy>::value)

//Version with SIMD vectorized data
template< typename mf_Complex, int SiteSize, typename DimensionPolicy, typename FlavorPolicy, typename AllocPolicy>
void cyclicPermute(CPSfield<mf_Complex,SiteSize,DimensionPolicy,FlavorPolicy,AllocPolicy> &to, const CPSfield<mf_Complex,SiteSize,DimensionPolicy,FlavorPolicy,AllocPolicy> &from,
		   const int dir, const int pm, const int n,
		   typename my_enable_if<CONDITION, const int>::type dummy = 0){
  enum {Dimension = DimensionPolicy::EuclideanDimension};
  assert(dir < Dimension);
  assert(n < GJP.NodeSites(dir));
  assert(pm == 1 || pm == -1);
  
  if(&to == &from){
    if(n==0) return;    
    CPSfield<mf_Complex,SiteSize,DimensionPolicy,FlavorPolicy,AllocPolicy> tmpfrom(from);
    return cyclicPermute(to,tmpfrom,dir,pm,n);
  }
  if(n == 0){
    to = from;
    return;
  }

  const int nsimd = mf_Complex::Nsimd();
  
  //Use notation c (combined index), o (outer index) i (inner index)
  
  int bcsites = n; //sites on boundary
  int bcsizes[Dimension]; bcsizes[dir] = n;
  int bcoff[Dimension]; bcoff[dir] = (pm == 1 ? GJP.NodeSites(dir)-n : 0);
  int bcoff_postcomms[Dimension]; bcoff_postcomms[dir] = (pm == 1 ? 0 : GJP.NodeSites(dir)-n);
  
  for(int i=0;i<Dimension;i++)
    if(i != dir){
      bcsizes[i] = GJP.NodeSites(i);
      bcsites *= bcsizes[i];
      bcoff[i] = 0;
      bcoff_postcomms[i] = 0;
    }

  //Build table of points on face (both outer and inner index)
  int nf = from.nflavors();
  int flav_off = from.flav_offset();

  typedef typename Grid::GridTypeMapper<mf_Complex>::scalar_type scalarType;
  
  int bufsz = bcsites * SiteSize * nf;

  QMP_mem_t *recv_mem = QMP_allocate_memory(bufsz * sizeof(scalarType));
  scalarType *recv_buf = (scalarType *)QMP_get_memory_pointer(recv_mem);

  QMP_mem_t *send_mem = QMP_allocate_memory(bufsz * sizeof(scalarType));
  scalarType *send_buf = (scalarType *)QMP_get_memory_pointer(send_mem);

  int osites = from.nsites();
  std::vector<int> to_oi_buf_map(nf * osites * nsimd); //map from outer and inner index of destination site to offset within buffer, used *after* comms.
  //map i + nsimd*(o + osites*f) as index
  
#pragma omp parallel for
  for(int c=0;c<bcsites;c++){
    int rem = c;
    int coor[Dimension];
    for(int d=0;d<Dimension;d++){ coor[d] = rem % bcsizes[d]; rem/=bcsizes[d]; }

    int coor_dest[Dimension];
    for(int d=0;d<Dimension;d++){
      coor_dest[d] = coor[d] + bcoff_postcomms[d];
      coor[d] += bcoff[d];
    }
    
    int i = from.SIMDmap(coor);
    int o = from.siteMap(coor);

    int i_dest = from.SIMDmap(coor_dest);
    int o_dest = from.siteMap(coor_dest);

    Grid::Vector<scalarType> ounpacked(nsimd);
    for(int f=0;f<nf;f++){
      mf_Complex const *osite_ptr = from.site_ptr(o,f);
      int send_buf_off = (c + bcsites*f)*SiteSize;
      scalarType* bp = send_buf + send_buf_off;
      to_oi_buf_map[ i_dest + nsimd*(o_dest+osites*f) ] = send_buf_off;
      
      for(int s=0;s<SiteSize;s++){
	vstore(*(osite_ptr++), ounpacked.data());
	*(bp++) = ounpacked[i];
      }      
    }
  }

  //Send/receive
  QMP_msgmem_t send_msg = QMP_declare_msgmem(send_buf,bufsz * sizeof(scalarType));
  QMP_msgmem_t recv_msg = QMP_declare_msgmem(recv_buf,bufsz * sizeof(scalarType));
  
  QMP_msghandle_t send = QMP_declare_send_relative(send_msg, dir, pm, 0);
  QMP_msghandle_t recv = QMP_declare_receive_relative(recv_msg, dir, -pm, 0);
  QMP_start(recv);
  QMP_start(send);
  
  QMP_status_t send_status = QMP_wait(send);
  if (send_status != QMP_SUCCESS) 
    QMP_error("Send failed in cyclicPermute: %s\n", QMP_error_string(send_status));
  QMP_status_t rcv_status = QMP_wait(recv);
  if (rcv_status != QMP_SUCCESS) 
    QMP_error("Receive failed in PassDataT: %s\n", QMP_error_string(rcv_status));


  
  //Copy remaining sites from on-node data with shift and pull in data from buffer simultaneously
  //if we sent in the + direction we need to shift the remaining L-n sites {0...L-n-1} forwards by n to make way for a new slice at the left side
  //if we sent in the - direction we need to shift the remaining L-n sites {n ... L-1} backwards by n to make way for a new slice at the right side
  //Problem is we don't want two threads writing to the same AVX register at the same time. Therefore we thread the loop over the destination SIMD vectors and work back
  std::vector< std::vector<int> > lane_offsets(nsimd,  std::vector<int>(Dimension) );
  for(int i=0;i<nsimd;i++) from.SIMDunmap(i, lane_offsets[i].data() );

#pragma omp parallel for
  for(int oto = 0;oto < osites; oto++){
    int oto_base_coor[Dimension]; to.siteUnmap(oto,oto_base_coor);

    //For each destination lane compute the source site index and lane
    int from_lane[nsimd];
    int from_osite_idx[nsimd]; //also use for recv_buf offsets for sites pulled over boundary
    for(int lane = 0; lane < nsimd; lane++){
      int offrom_coor[Dimension];
      for(int d=0;d<Dimension;d++) offrom_coor[d] = oto_base_coor[d] + lane_offsets[lane][d];
      offrom_coor[dir] += (pm == 1 ? -n : n);

      if(offrom_coor[dir] < 0 || offrom_coor[dir] >= GJP.NodeSites(dir)){
	from_lane[lane] = -1; //indicates data is in recv_buf	
	from_osite_idx[lane] = to_oi_buf_map[ lane + nsimd*oto ]; //here is for flavor 0 - remember to offset for second flav
      }else{
	from_lane[lane] = from.SIMDmap(offrom_coor);
	from_osite_idx[lane] = from.siteMap(offrom_coor);
      }
    }

    //Now loop over flavor and element within the site as well as SIMD lanes of the destination vector and gather what we need to poke - then poke it
    Grid::Vector<scalarType> towrite(nsimd);
    Grid::Vector<scalarType> unpack(nsimd);
    
    for(int f=0;f<nf;f++){
      for(int s=0;s<SiteSize;s++){
	for(int tolane=0;tolane<nsimd;tolane++){	  
	  if(from_lane[tolane] != -1){
	    mf_Complex const* from_osite_ptr = from.site_ptr(from_osite_idx[tolane], f) + s;
	    vstore(*from_osite_ptr,unpack.data());
	    towrite[tolane] = unpack[ from_lane[tolane] ];
	  }else{
	    //data is in buffer
	    towrite[tolane] = recv_buf[ from_osite_idx[tolane] + s + f*bcsites*SiteSize ];
	  }
	    
	}
	mf_Complex* to_osite_ptr = to.site_ptr(oto,f) + s;
	vset(*to_osite_ptr, towrite.data());	
      }
    }
  }

  QMP_free_msghandle(send);
  QMP_free_msghandle(recv);
  QMP_free_msgmem(send_msg);
  QMP_free_msgmem(recv_msg);
  QMP_free_memory(send_mem);
  QMP_free_memory(recv_mem);
  QMP_barrier();
}
#undef CONDITION

# endif //ifdef USE_GRID

#else //ifdef USE_QMP

#define CONDITION _equal<typename ComplexClassify<mf_Complex>::type, complex_double_or_float_mark>::value && (_equal<DimensionPolicy,FourDpolicy>::value || _equal<DimensionPolicy,SpatialPolicy>::value)

template< typename mf_Complex, int SiteSize, typename DimensionPolicy, typename FlavorPolicy, typename AllocPolicy>
void cyclicPermute(CPSfield<mf_Complex,SiteSize,DimensionPolicy,FlavorPolicy,AllocPolicy> &to, const CPSfield<mf_Complex,SiteSize,DimensionPolicy,FlavorPolicy,AllocPolicy> &from,
		   const int dir, const int pm, const int n,
		   typename my_enable_if<CONDITION , const int>::type dummy = 0){
  enum {Dimension = DimensionPolicy::EuclideanDimension};
  assert(dir < Dimension);
  assert(n < GJP.NodeSites(dir));
  assert(pm == 1 || pm == -1);
  
  if(&to == &from){
    if(n==0) return;    
    CPSfield<mf_Complex,SiteSize,DimensionPolicy,FlavorPolicy,AllocPolicy> tmpfrom(from);
    return cyclicPermute(to,tmpfrom,dir,pm,n);
  }
  if(n == 0){
    to = from;
    return;
  }
  const int nodes = GJP.Xnodes()*GJP.Ynodes()*GJP.Znodes()*GJP.Tnodes()*GJP.Snodes();
  if(nodes != 1) ERR.General("","cyclicPermute","Parallel implementation requires QMP\n");

#pragma omp parallel for
  for(int i=0;i<from.nfsites();i++){
    int f; int x[Dimension];
    from.fsiteUnmap(i,x,f);
    x[dir] = (x[dir] + pm * n + 5*GJP.NodeSites(dir) ) % GJP.NodeSites(dir);
    const mf_Complex* from_ptr = from.fsite_ptr(i);
    mf_Complex* to_ptr = to.site_ptr(x,f);
    memcpy(to_ptr,from_ptr,SiteSize*sizeof(mf_Complex));
  }
}
#undef CONDITION

# ifdef USE_GRID

#define CONDITION _equal<typename ComplexClassify<mf_Complex>::type, grid_vector_complex_mark>::value && (_equal<DimensionPolicy,FourDSIMDPolicy>::value || _equal<DimensionPolicy,ThreeDSIMDPolicy>::value)

//Version with SIMD vectorized data
template< typename mf_Complex, int SiteSize, typename DimensionPolicy, typename FlavorPolicy, typename AllocPolicy>
void cyclicPermute(CPSfield<mf_Complex,SiteSize,DimensionPolicy,FlavorPolicy,AllocPolicy> &to, const CPSfield<mf_Complex,SiteSize,DimensionPolicy,FlavorPolicy,AllocPolicy> &from,
		   const int dir, const int pm, const int n,
		   typename my_enable_if<CONDITION, const int>::type dummy = 0){
  enum {Dimension = DimensionPolicy::EuclideanDimension};
  assert(dir < Dimension);
  assert(n < GJP.NodeSites(dir));
  assert(pm == 1 || pm == -1);
  
  if(&to == &from){
    if(n==0) return;    
    CPSfield<mf_Complex,SiteSize,DimensionPolicy,FlavorPolicy,AllocPolicy> tmpfrom(from);
    return cyclicPermute(to,tmpfrom,dir,pm,n);
  }
  if(n == 0){
    to = from;
    return;
  }
  const int nodes = GJP.Xnodes()*GJP.Ynodes()*GJP.Znodes()*GJP.Tnodes()*GJP.Snodes();
  if(nodes != 1) ERR.General("","cyclicPermute","Parallel implementation requires QMP\n");
  
  const int nsimd = mf_Complex::Nsimd();

  typedef typename mf_Complex::scalar_type scalar_type;
  const int nthr = omp_get_max_threads();
  scalar_type* tmp_store_thr[nthr]; for(int i=0;i<nthr;i++) tmp_store_thr[i] = (scalar_type*)memalign(128,nsimd*sizeof(scalar_type));
  
#pragma omp parallel for
  for(int ofto=0;ofto<to.nfsites();ofto++){ //loop over outer site index
    const int me = omp_get_thread_num();
    int f; int oxto[Dimension];
    to.fsiteUnmap(ofto,oxto,f);

    mf_Complex* to_base_ptr = to.fsite_ptr(ofto);
    
    scalar_type* tmp_store = tmp_store_thr[me];

    //indexed by destination lane
    mf_Complex const* from_base_ptrs[nsimd];
    int from_lane_idx[nsimd];
      
    for(int tolane = 0; tolane < nsimd; tolane++){
      int ixto_off[Dimension];
      to.SIMDunmap(tolane,ixto_off); //get offset of inner site on tolane

      int xfrom[Dimension]; for(int d=0;d<Dimension;d++) xfrom[d] = oxto[d] + ixto_off[d]; //full coord corresponding to tolane + outer site
      xfrom[dir] = (xfrom[dir] - pm * n + 5*GJP.NodeSites(dir) ) % GJP.NodeSites(dir);

      from_base_ptrs[tolane] = from.site_ptr(xfrom,f);
      from_lane_idx[tolane] = from.SIMDmap(xfrom);
    }

    for(int s=0;s<SiteSize;s++){
      for(int tolane = 0; tolane < nsimd; tolane++)
	tmp_store[tolane] = *( (scalar_type*)(from_base_ptrs[tolane] + s) + from_lane_idx[tolane] ); //cast SIMD type to scalar type pointer
      vset(*(to_base_ptr + s), tmp_store);
    }
  }            
  for(int i=0;i<nthr;i++) free(tmp_store_thr[i]);  
}
#undef CONDITION
  
# endif //ifdef USE_GRID

#endif //ifdef USE_QMP


inline int getShiftSign(const int of){ return of > 0 ? +1 : -1; }

//Invoke multiple independent permutes to offset field by vector 'shift' assuming field is periodic
template<typename FieldType>
void shiftPeriodicField(FieldType &to, const FieldType &from, const std::vector<int> &shift){
  int nd = shift.size(); //assume ascending: x,y,z,t
  int nshift_dirs = 0;
  for(int i=0;i<nd;i++) if(shift[i]!=0) ++nshift_dirs;

  if(nshift_dirs == 0){
    if(&to != &from) to = from;
    return;
  }else if(nshift_dirs == 1){
    for(int d=0;d<nd;d++){
      if(shift[d] != 0){
	cyclicPermute(to,from,d,getShiftSign(shift[d]),abs(shift[d]) );
	return;
      }
    }    
  }else{
    FieldType tmp1 = from;
    FieldType tmp2 = from;
    FieldType * send = &tmp1;
    FieldType * recv = &tmp2;

    int shifts_done = 0;
    for(int d=0;d<nd;d++){
      if(shift[d] != 0){
	cyclicPermute(shifts_done < nshift_dirs-1 ? *recv : to,*send,d,getShiftSign(shift[d]),abs(shift[d]) );
	++shifts_done;
	if(shifts_done < nshift_dirs) std::swap(send,recv);
	else return;
      }
    }   
  }
}







template<typename CPSfieldType>
void fft(CPSfieldType &into, const CPSfieldType &from, const bool* do_dirs, const bool inverse_transform = false,
	 typename my_enable_if<_equal<typename ComplexClassify<typename CPSfieldType::FieldSiteType>::type, complex_double_or_float_mark>::value, const int>::type = 0
	 ){
  typedef typename LocalToGlobalInOneDirMap<typename CPSfieldType::FieldDimensionPolicy>::type DimPolGlobalInOneDir;
  typedef CPSfieldGlobalInOneDir<typename CPSfieldType::FieldSiteType, CPSfieldType::FieldSiteSize, DimPolGlobalInOneDir, typename CPSfieldType::FieldFlavorPolicy, typename CPSfieldType::FieldAllocPolicy> CPSfieldTypeGlobalInOneDir;

  int dcount = 0;
  
  for(int mu=0;mu<CPSfieldType::FieldDimensionPolicy::EuclideanDimension;mu++)
    if(do_dirs[mu]){
      CPSfieldTypeGlobalInOneDir tmp_dbl(mu);
      tmp_dbl.gather( dcount==0 ? from : into );
      tmp_dbl.fft(inverse_transform);
      tmp_dbl.scatter(into);
      dcount ++;
    }
}

#ifdef USE_GRID
template<typename CPSfieldType>
void fft(CPSfieldType &into, const CPSfieldType &from, const bool* do_dirs, const bool inverse_transform = false,
	 typename my_enable_if<_equal<typename ComplexClassify<typename CPSfieldType::FieldSiteType>::type, grid_vector_complex_mark>::value, const int>::type = 0
	 ){
  typedef typename Grid::GridTypeMapper<typename CPSfieldType::FieldSiteType>::scalar_type ScalarType;
  typedef typename CPSfieldType::FieldDimensionPolicy::EquivalentScalarPolicy ScalarDimPol;
  typedef CPSfield<ScalarType, CPSfieldType::FieldSiteSize, ScalarDimPol, typename CPSfieldType::FieldFlavorPolicy, StandardAllocPolicy> ScalarFieldType;

  NullObject null_obj;
  ScalarFieldType tmp_in(null_obj);
  ScalarFieldType tmp_out(null_obj);
  tmp_in.importField(from);
  fft(tmp_out, tmp_in, do_dirs, inverse_transform);
  tmp_out.exportField(into);
}
#endif
  
template<typename CPSfieldType>
void fft(CPSfieldType &fftme, const bool* do_dirs){
  fft(fftme,fftme,do_dirs);
}

template<typename CPSfieldType>
void fft_opt(CPSfieldType &into, const CPSfieldType &from, const bool* do_dirs, const bool inverse_transform = false,
	     typename my_enable_if<_equal<typename ComplexClassify<typename CPSfieldType::FieldSiteType>::type, complex_double_or_float_mark>::value, const int>::type = 0
	     ){
#ifndef USE_MPI
  fft(into,from,do_dirs,inverse_transform);
#else
  
  enum { Dimension = CPSfieldType::FieldDimensionPolicy::EuclideanDimension };
  int ndirs_fft = 0; for(int i=0;i<Dimension;i++) if(do_dirs[i]) ++ndirs_fft;
  if(! ndirs_fft ) return;

  //Need info on the MPI node mapping
  assert(GJP.Snodes() == 1);
  std::vector<int> node_map;
  getMPIrankMap(node_map);

  CPSfieldType tmp(from.getDimPolParams());

  //we want the last fft to end up in 'into'. Intermediate FFTs cycle between into and tmp as temp storage. Thus for odd ndirs_fft, the first fft should output to 'into', for even it should output to 'tmp'
  CPSfieldType *tmp1, *tmp2;
  if(ndirs_fft % 2 == 1){
    tmp1 = &into; tmp2 = &tmp;
  }else{
    tmp1 = &tmp; tmp2 = &into;
  }
  
  CPSfieldType* src = tmp2;
  CPSfieldType* out = tmp1;

  int fft_count = 0;
  for(int mu=0; mu<Dimension; mu++){
    if(do_dirs[mu]){
      CPSfieldType const *msrc = fft_count == 0 ? &from : src;
      fft_opt_mu(*out, *msrc, mu, node_map, inverse_transform);
      ++fft_count;
      std::swap(src,out);      
    }
  }
#endif
}

#ifdef USE_MPI
template<typename CPSfieldType>
void fft_opt_mu(CPSfieldType &into, const CPSfieldType &from, const int mu, const std::vector<int> &node_map, const bool inverse_transform,
	     typename my_enable_if<_equal<typename ComplexClassify<typename CPSfieldType::FieldSiteType>::type, complex_double_or_float_mark>::value, const int>::type = 0
	     ){
  enum {SiteSize = CPSfieldType::FieldSiteSize, Dimension = CPSfieldType::FieldDimensionPolicy::EuclideanDimension };
  typedef typename CPSfieldType::FieldSiteType ComplexType;
  typedef typename ComplexType::value_type FloatType;
  typedef typename FFTWwrapper<FloatType>::complexType FFTComplex;
  const int nf = from.nflavors();
  const int foff = from.flav_offset();
  const int nthread = omp_get_max_threads();
  
  //Eg for fft in X-direction, divide up Y,Z,T work over nodes in X-direction doing linear FFTs.
  const int munodesites = GJP.NodeSites(mu);
  const int munodes = GJP.Nodes(mu);
  const int mutotalsites = munodesites*munodes;
  const int munodecoor = GJP.NodeCoor(mu);
  const int n_orthdirs = Dimension - 1;
  FloatType Lmu(mutotalsites);
  
  int orthdirs[n_orthdirs]; //map of orthogonal directions to mu
  int total_work_munodes = 1; //sites orthogonal to FFT direction
  int o=0;
  for(int i=0;i< Dimension;i++)
    if(i!=mu){
      total_work_munodes *= GJP.NodeSites(i);
      orthdirs[o++] = i;
    }

  //Divvy up work over othogonal directions
  int munodes_work[munodes];
  int munodes_off[munodes];
  for(int i=0;i<munodes;i++)
    thread_work(munodes_work[i],munodes_off[i], total_work_munodes, i, munodes); //use for node work instead :)

  //Get MPI ranks of nodes in mu direction
  int my_node_coor[4];
  for(int i=0;i<4;i++) my_node_coor[i] = GJP.NodeCoor(i);
  
  int munodes_mpiranks[munodes];
  for(int i=0;i<munodes;i++){
    int munode_coor[4]; memcpy(munode_coor,my_node_coor,4*sizeof(int));
    munode_coor[mu] = i;

    const int munode_lex = node_lex( munode_coor, 4 );
    munodes_mpiranks[i] = node_map[munode_lex];
  }

  //Gather send data
  ComplexType* send_bufs[munodes];
  int send_buf_sizes[munodes];
  for(int i=0;i<munodes;i++){
    send_buf_sizes[i] = munodes_work[i] * munodesites * nf * SiteSize;
    send_bufs[i] = (ComplexType*)malloc( send_buf_sizes[i] * sizeof(ComplexType) );

    for(int w = 0; w < munodes_work[i]; w++){ //index of orthogonal site within workload for i'th node in mu direction
      const int orthsite = munodes_off[i] + w;
      int coor_base[Dimension] = {0};
	  
      //Unmap orthsite into a base coordinate
      int rem = orthsite;
      for(int a=0;a<n_orthdirs;a++){
	const int dir_a = orthdirs[a];
	coor_base[dir_a] = rem % GJP.NodeSites(dir_a); rem /= GJP.NodeSites(dir_a);
      }

      for(int f=0;f<nf;f++){
	for(int xmu=0;xmu<munodesites;xmu++){
	  ComplexType* to = send_bufs[i] + SiteSize * (w + munodes_work[i]*( f + nf*xmu ) );  //with musite changing slowest
	  coor_base[mu] = xmu;
	  ComplexType const* frm = from.site_ptr(coor_base,f);

	  memcpy(to,frm,SiteSize*sizeof(ComplexType));
	}
      }
    }
  }
  MPI_Request send_req[munodes];
  MPI_Request recv_req[munodes];
  MPI_Status status[munodes];

  //Prepare recv buf
  const int bufsz = munodes_work[munodecoor] * mutotalsites * nf * SiteSize; //complete line in mu for each orthogonal coordinate
  ComplexType* recv_buf = (ComplexType*)malloc(bufsz * sizeof(ComplexType) );

  //Setup send/receive    
  for(int i=0;i<munodes;i++){ //works fine to send to all nodes, even if this involves a send to self.
    int sret = MPI_Isend(send_bufs[i], send_buf_sizes[i]*sizeof(ComplexType), MPI_CHAR, munodes_mpiranks[i], 0, MPI_COMM_WORLD, &send_req[i]);
    assert(sret == MPI_SUCCESS);

    int rret = MPI_Irecv(recv_buf + i*munodes_work[munodecoor]*nf*SiteSize*munodesites, send_buf_sizes[i]*sizeof(ComplexType), MPI_CHAR, munodes_mpiranks[i], MPI_ANY_TAG, MPI_COMM_WORLD, &recv_req[i]);
    assert(rret == MPI_SUCCESS);
  }

      
  int wret = MPI_Waitall(munodes,recv_req,status);
  assert(wret == MPI_SUCCESS);
      
  //Do FFT
  const int howmany = munodes_work[munodecoor] * nf * SiteSize;
  const int howmany_per_thread_base = howmany / nthread;
  //Divide work orthogonal to mu, 'howmany', over threads. Note, this may not divide howmany equally. The difference is made up by adding 1 unit of work to threads in ascending order until total work matches. Thus we need 2 plans: 1 for the base amount and one for the base+1

  //if(!UniqueID()) printf("FFT work per site %d, divided over %d threads with %d work each. Remaining work %d allocated to ascending threads\n", howmany, nthread, howmany_per_thread_base, howmany - howmany_per_thread_base*nthread);

  int fft_phase = inverse_transform ? FFTW_BACKWARD : FFTW_FORWARD;
  
  static FFTplanContainer<FloatType> plan_f_base[Dimension]; //destructors deallocate plans
  static FFTplanContainer<FloatType> plan_f_base_p1[Dimension];
      
  static int plan_howmany[Dimension];
  static bool plan_init = false;
  static int plan_fft_phase;
  
  if(!plan_init || plan_howmany[mu] != howmany || fft_phase != plan_fft_phase){
    if(!plan_init) for(int i=0;i<Dimension;i++) plan_howmany[i] = -1;

    typename FFTWwrapper<FloatType>::complexType *tmp_f; //I don't think it actually does anything with this

    plan_fft_phase = fft_phase;
    const int fft_work_per_musite = howmany_per_thread_base;
    const int musite_stride = howmany; //stride between musites
    
    plan_f_base[mu].setPlan(1, &mutotalsites, fft_work_per_musite, 
			    tmp_f, NULL, musite_stride, 1,
			    tmp_f, NULL, musite_stride, 1,
			    plan_fft_phase, FFTW_ESTIMATE);
    plan_f_base_p1[mu].setPlan(1, &mutotalsites, fft_work_per_musite+1, 
			       tmp_f, NULL, musite_stride, 1,
			       tmp_f, NULL, musite_stride, 1,
			       plan_fft_phase, FFTW_ESTIMATE);	
    plan_init = true; //other mu's will still init later
  }
  FFTComplex*fftw_mem = (FFTComplex*)recv_buf;

#pragma omp parallel
  {
    assert(nthread == omp_get_num_threads()); //plans will be messed up if not true
    const int me = omp_get_thread_num();
    int thr_work, thr_off;
    thread_work(thr_work, thr_off, howmany, me, nthread);

    const FFTplanContainer<FloatType>* thr_plan_ptr;
    
    if(thr_work == howmany_per_thread_base) thr_plan_ptr = &plan_f_base[mu];
    else if(thr_work == howmany_per_thread_base + 1) thr_plan_ptr = &plan_f_base_p1[mu];
    else assert(0); //catch if logic for thr_work changes

    FFTWwrapper<FloatType>::execute_dft(thr_plan_ptr->getPlan(), fftw_mem + thr_off, fftw_mem + thr_off); 
  }

  wret = MPI_Waitall(munodes,send_req,status);
  assert(wret == MPI_SUCCESS);
      
  //Send back out. Reuse the old send buffers as receive buffers and vice versa
  for(int i=0;i<munodes;i++){ //works fine to send to all nodes, even if this involves a send to self
    int sret = MPI_Isend(recv_buf + i*munodes_work[munodecoor]*nf*SiteSize*munodesites, send_buf_sizes[i]*sizeof(ComplexType), MPI_CHAR, munodes_mpiranks[i], 0, MPI_COMM_WORLD, &send_req[i]);
    assert(sret == MPI_SUCCESS);

    int rret = MPI_Irecv(send_bufs[i], send_buf_sizes[i]*sizeof(ComplexType), MPI_CHAR, munodes_mpiranks[i], MPI_ANY_TAG, MPI_COMM_WORLD, &recv_req[i]);
    assert(rret == MPI_SUCCESS);
  }

  wret = MPI_Waitall(munodes,recv_req,status);
  assert(wret == MPI_SUCCESS);

  //Poke into output
  for(int i=0;i<munodes;i++){
#pragma omp parallel for
    for(int w = 0; w < munodes_work[i]; w++){ //index of orthogonal site within workload for i'th node in mu direction
      const int orthsite = munodes_off[i] + w;
      int coor_base[Dimension] = {0};
	  
      //Unmap orthsite into a base coordinate
      int rem = orthsite;
      for(int a=0;a<n_orthdirs;a++){
	int dir_a = orthdirs[a];
	coor_base[dir_a] = rem % GJP.NodeSites(dir_a); rem /= GJP.NodeSites(dir_a);
      }

      for(int f=0;f<nf;f++){
	for(int xmu=0;xmu<munodesites;xmu++){	      
	  coor_base[mu] = xmu;
	  ComplexType* to = into.site_ptr(coor_base,f);
	  ComplexType const* frm = send_bufs[i] + SiteSize * (w + munodes_work[i]*( f + nf*xmu ) );
	  if(!inverse_transform) memcpy(to,frm,SiteSize*sizeof(ComplexType));
	  else for(int s=0;s<SiteSize;s++) to[s] = frm[s]/Lmu;
	}
      }
    }
  }

  wret = MPI_Waitall(munodes,send_req,status);
  assert(wret == MPI_SUCCESS);
      
  free(recv_buf);
  for(int i=0;i<munodes;i++) free(send_bufs[i]);
}
#endif


#ifdef USE_GRID
template<typename CPSfieldType>
void fft_opt(CPSfieldType &into, const CPSfieldType &from, const bool* do_dirs, const bool inverse_transform = false,
	     typename my_enable_if<_equal<typename ComplexClassify<typename CPSfieldType::FieldSiteType>::type, grid_vector_complex_mark>::value, const int>::type = 0
	     ){ //we can avoid the copies below but with some effort - do at some point
# ifdef USE_MPI
  fft(into,from,do_dirs,inverse_transform);
# else
  typedef typename Grid::GridTypeMapper<typename CPSfieldType::FieldSiteType>::scalar_type ScalarType;
  typedef typename CPSfieldType::FieldDimensionPolicy::EquivalentScalarPolicy ScalarDimPol;
  typedef CPSfield<ScalarType, CPSfieldType::FieldSiteSize, ScalarDimPol, typename CPSfieldType::FieldFlavorPolicy, StandardAllocPolicy> ScalarFieldType;

  NullObject null_obj;
  ScalarFieldType tmp_in(null_obj);
  ScalarFieldType tmp_out(null_obj);
  tmp_in.importField(from);
  fft_opt(tmp_out, tmp_in, do_dirs, inverse_transform);
  tmp_out.exportField(into);
# endif
}
#endif

CPS_END_NAMESPACE
#endif
