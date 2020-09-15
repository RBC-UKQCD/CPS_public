#ifndef MESONFIELD_IO
#define MESONFIELD_IO

struct hostEndian{
  enum EndianType { BIG, LITTLE };
  inline static EndianType get(){ //copied from fpconv
    char end_check[4] = {1,0,0,0};
    uint32_t *lp = (uint32_t *)end_check;
    if ( *lp == 0x1 ) { 
      return LITTLE;
    } else {
      return BIG;
    }
  }
};

template<typename T>
struct FPformat{
  inline static FP_FORMAT get(){ //also taken from fpconv
    assert(sizeof(T) == 4 || sizeof(T) == 8);
    static const hostEndian::EndianType endian = hostEndian::get();
    
    if(sizeof(T) == 8){
      return endian == hostEndian::LITTLE ? FP_IEEE64LITTLE : FP_IEEE64BIG;
    }else {  // 32 bits
      union { 
	float pinum;
	char pichar[4];
      }cpspi;

      FP_FORMAT format;
      
      cpspi.pinum = FPConv_PI;
      if(endian == hostEndian::BIG) {
	format = FP_IEEE32BIG;
	for(int i=0;i<4;i++) {
	  if(cpspi.pichar[i] != FPConv_ieee32pi_big[i]) {
	    format = FP_TIDSP32;
	    break;
	  }
	}
      }
      else {
	format = FP_IEEE32LITTLE;
	for(int i=0;i<4;i++) {
	  if(cpspi.pichar[i] != FPConv_ieee32pi_big[3-i]) {
	    format = FP_TIDSP32;
	    break;
	  }
	}
      }
      return format;
    } // end of 32 bits
  }   
};

template<typename T>
struct FPformat<std::complex<T> >{
  inline static FP_FORMAT get(){ return FPformat<T>::get(); }
};

template<typename mf_Policies, template <typename> class A2AfieldL,  template <typename> class A2AfieldR>
void A2AmesonField<mf_Policies,A2AfieldL,A2AfieldR>::write(const std::string &filename, FP_FORMAT fileformat) const{
  if(!UniqueID()) printf("Writing meson field of size %d kB to file %s\n",byte_size()/1024,filename.c_str());
  std::ofstream *file = !UniqueID() ? new std::ofstream(filename.c_str(),std::ofstream::out) : NULL;

  write(file,fileformat);

  if(!UniqueID()){
    file->close();
    delete file;
  }
}  
  
//ostream pointer should only be open on node 0 - should be NULL otherwise
template<typename mf_Policies, template <typename> class A2AfieldL,  template <typename> class A2AfieldR>
void A2AmesonField<mf_Policies,A2AfieldL,A2AfieldR>::write(std::ostream *file_ptr, FP_FORMAT fileformat) const{
  if(!UniqueID()) assert(file_ptr != NULL);
  else assert(file_ptr == NULL);

  if(!UniqueID()){
    assert(!file_ptr->fail());
    file_ptr->exceptions ( std::ofstream::failbit | std::ofstream::badbit );
  }

#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  
  if(node_mpi_rank != -1){
#ifdef USE_MPI    
    int my_rank;
    int ret = MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    if(ret != MPI_SUCCESS) ERR.General("A2AmesonField","write","Comm_rank failed\n");
    
    if(!UniqueID() && node_mpi_rank != my_rank)
      ERR.General("A2AmesonField","write","Mesonfield must be present on head node\n");
#else
    ERR.General("A2AmesonField","write","Mesonfield must be present on head node (no mpi)\n");
#endif
    
  }

  if(!UniqueID()){
    FP_FORMAT dataformat = FPformat<ScalarComplexType>::get();
    FPConv conv;
    if(fileformat == FP_AUTOMATIC)
      fileformat = dataformat;
    else
      assert(conv.size(fileformat) == conv.size(dataformat));
    
    conv.setHostFormat(dataformat);
    conv.setFileFormat(fileformat);

    int dsize = conv.size(dataformat);
    unsigned int checksum = conv.checksum( (char*)mf, 2*fsize, dataformat);

    //Header
    std::ostream &file = *file_ptr;     

    file << "BEGIN_HEADER\n";
    file << "HDR_VERSION = 1\n";
    file << "DATA_FORMAT = " << conv.name(fileformat) << '\n';
    file << "CHECKSUM = " << checksum << "\n";
    file << "END_HEADER\n";
    
    //Parameters    
    char* la2aparams_buf = (char*)malloc(10000 * sizeof(char));
    char* ra2aparams_buf = (char*)malloc(10000 * sizeof(char));
    {
      VML vml;
      vml.Create(la2aparams_buf,10000,VML_ENCODE);
      A2AArg &a2a_args_l = const_cast<A2AArg &>(lindexdilution.getArgs());
      assert( a2a_args_l.Vml(&vml,"A2AARGS_L") );
      vml.Destroy();
    }
    {
      VML vml;
      vml.Create(ra2aparams_buf,10000,VML_ENCODE);
      A2AArg &a2a_args_r = const_cast<A2AArg &>(rindexdilution.getArgs());
      assert( a2a_args_r.Vml(&vml,"A2AARGS_R") );
      vml.Destroy();
    }		
      
    file << "BEGIN_PARAMS\n";
    //int nmodes_l, nmodes_r;  //derived from index dilutions
    file << "FSIZE = " << fsize << "\n";
    file << "TL = " << tl << "\n";
    file << "TR = " << tr << "\n";
    file << "LINDEXDILUTION = " << lindexdilution.name() << "\n";
    file << "RINDEXDILUTION = " << rindexdilution.name() << "\n";
    file << "STRLEN_A2AARGS_L = " << strlen(la2aparams_buf) << '\n';
    file << "STRLEN_A2AARGS_R = " << strlen(ra2aparams_buf) << '\n';
    file << la2aparams_buf << ra2aparams_buf;
    file << "END_PARAMS\n";      

    free(la2aparams_buf);
    free(ra2aparams_buf);

    //Data
    file << "BEGIN_DATA\n";
      
    static const int chunk = 32768; //32kb chunks
    assert(chunk % dsize == 0);
    int fdinchunk = chunk/dsize;
    char* wbuf = (char*)malloc(chunk * sizeof(char)); 
      
    char const* dptr = (char const*)mf;

    int off = 0;
    int nfd = 2*fsize;
    while(off < nfd){
      int grab = std::min(nfd-off, fdinchunk); //How many data elements to grab
      int grabchars = grab * dsize;
      conv.host2file(wbuf,dptr,grab);
      file.write(wbuf,grabchars);
      off += grab;
      dptr += grabchars;
    }
    file << "\nEND_DATA\n";
    free(wbuf);
  }

#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

template<typename mf_Policies, template <typename> class A2AfieldL,  template <typename> class A2AfieldR>
void A2AmesonField<mf_Policies,A2AfieldL,A2AfieldR>::read(const std::string &filename){
  if(!UniqueID()) printf("Reading meson field from file %s\n",filename.c_str());
  std::ifstream *file = !UniqueID() ?  new std::ifstream(filename.c_str()) : NULL;
    
  read(file);

  if(!UniqueID()){
    file->close();
    delete file;
  }
}

//istream pointer should only be open on node 0 - should be NULL otherwise
template<typename mf_Policies, template <typename> class A2AfieldL,  template <typename> class A2AfieldR>
void A2AmesonField<mf_Policies,A2AfieldL,A2AfieldR>::read(std::istream *file_ptr){
  if(!UniqueID()) assert(file_ptr != NULL);
  else assert(file_ptr == NULL);
  
  if(!UniqueID()){
    assert(!file_ptr->fail());
    file_ptr->exceptions ( std::ofstream::failbit | std::ofstream::badbit );
  }
#ifdef USE_MPI
  int my_mpi_rank = getMyMPIrank(); //Get this node's mpi rank  
  int head_mpi_rank = getHeadMPIrank(); //Broadcast to all nodes the mpi rank of the head node (UniqueID() == 0)
#else
  if(GJP.Xnodes()*GJP.Ynodes()*GJP.Znodes()*GJP.Tnodes()*GJP.Snodes() != 1) ERR.General("A2AmesonField","read","Parallel implementation requires MPI\n");
#endif
  
  int read_fsize;
  unsigned int checksum;

  int a2aparams_l_buflen, a2aparams_r_buflen;
  char *a2aparams_l_buf, *a2aparams_r_buf;
  char dformatbuf[256];
  
  if(UniqueID() == 0){
    std::istream &file = *file_ptr;
    assert(!file.fail());

    std::string str;
    
    //Header
    getline(file,str); assert(str == "BEGIN_HEADER");
    getline(file,str); assert(str == "HDR_VERSION = 1");
    getline(file,str); assert( sscanf(str.c_str(),"DATA_FORMAT = %s",dformatbuf) == 1 );
    getline(file,str); assert( sscanf(str.c_str(),"CHECKSUM = %u",&checksum) == 1 );
    getline(file,str); assert(str == "END_HEADER");

    //Params
    getline(file,str); assert(str == "BEGIN_PARAMS");    
    getline(file,str); assert( sscanf(str.c_str(),"FSIZE = %d",&read_fsize) == 1 );
    getline(file,str); assert( sscanf(str.c_str(),"TL = %d",&tl) == 1 );
    getline(file,str); assert( sscanf(str.c_str(),"TR = %d",&tr) == 1 );

    char nmbuf[256];
    getline(file,str); assert( sscanf(str.c_str(),"LINDEXDILUTION = %s",nmbuf) == 1 );
    assert( std::string(nmbuf) == lindexdilution.name() );

    getline(file,str); assert( sscanf(str.c_str(),"RINDEXDILUTION = %s",nmbuf) == 1 );
    assert( std::string(nmbuf) == rindexdilution.name() );
    
    getline(file,str); assert( sscanf(str.c_str(),"STRLEN_A2AARGS_L = %d",&a2aparams_l_buflen) == 1 );
    getline(file,str); assert( sscanf(str.c_str(),"STRLEN_A2AARGS_R = %d",&a2aparams_r_buflen) == 1 );
    ++a2aparams_l_buflen; //leave room for null character!
    ++a2aparams_r_buflen;
    
    //Read the VMLs for the left and right A2A params and squirt out for each node to decode
    a2aparams_l_buf = (char*)malloc(a2aparams_l_buflen * sizeof(char));
    a2aparams_r_buf = (char*)malloc(a2aparams_r_buflen * sizeof(char));
    file.get(a2aparams_l_buf,a2aparams_l_buflen,EOF);
    file.get(a2aparams_r_buf,a2aparams_r_buflen,EOF);

    getline(file,str); assert(str == "END_PARAMS");
  }
#ifdef USE_MPI  
  MPI_Barrier(MPI_COMM_WORLD);
  
  //Squirt A2Aparams and whatnot over to other nodes for data setup
  int ret = MPI_Bcast(&checksum, 1, MPI_UNSIGNED, head_mpi_rank, MPI_COMM_WORLD);
  if(ret != MPI_SUCCESS) ERR.General("A2AmesonField","read","Squirt 1 fail\n");
  
  ret = MPI_Bcast(&read_fsize, 1, MPI_INT, head_mpi_rank, MPI_COMM_WORLD);
  if(ret != MPI_SUCCESS) ERR.General("A2AmesonField","read","Squirt 2 fail\n");
  
  ret = MPI_Bcast(&tl, 1, MPI_INT, head_mpi_rank, MPI_COMM_WORLD);
  if(ret != MPI_SUCCESS) ERR.General("A2AmesonField","read","Squirt 3 fail\n");

  ret = MPI_Bcast(&tr, 1, MPI_INT, head_mpi_rank, MPI_COMM_WORLD);
  if(ret != MPI_SUCCESS) ERR.General("A2AmesonField","read","Squirt 4 fail\n");
  
  ret = MPI_Bcast(&a2aparams_l_buflen, 1, MPI_INT, head_mpi_rank, MPI_COMM_WORLD);
  if(ret != MPI_SUCCESS) ERR.General("A2AmesonField","read","Squirt 5 fail\n");

  ret = MPI_Bcast(&a2aparams_r_buflen, 1, MPI_INT, head_mpi_rank, MPI_COMM_WORLD);
  if(ret != MPI_SUCCESS) ERR.General("A2AmesonField","read","Squirt 6 fail\n");

  if(UniqueID() != 0){
    //Other nodes create bufs for a2a params
    a2aparams_l_buf = (char*)malloc(a2aparams_l_buflen * sizeof(char));
    a2aparams_r_buf = (char*)malloc(a2aparams_r_buflen * sizeof(char));    
  }

  ret = MPI_Bcast(a2aparams_l_buf, a2aparams_l_buflen, MPI_CHAR, head_mpi_rank, MPI_COMM_WORLD);
  if(ret != MPI_SUCCESS) ERR.General("A2AmesonField","read","Squirt 7 fail\n");

  ret = MPI_Bcast(a2aparams_r_buf, a2aparams_r_buflen, MPI_CHAR, head_mpi_rank, MPI_COMM_WORLD);
  if(ret != MPI_SUCCESS) ERR.General("A2AmesonField","read","Squirt 8 fail\n");
#endif
  
  //Every node parse the params buffers
  A2AArg read_a2a_args_l, read_a2a_args_r;
  {
    VML vml;
    vml.Create(a2aparams_l_buf,a2aparams_l_buflen,VML_DECODE);
    assert( read_a2a_args_l.Vml(&vml,"A2AARGS_L") );
    vml.Destroy();
  }
  {
    VML vml;
    vml.Create(a2aparams_r_buf,a2aparams_r_buflen,VML_DECODE);
    assert( read_a2a_args_r.Vml(&vml,"A2AARGS_R") );
    vml.Destroy();
  }
  free(a2aparams_l_buf); free(a2aparams_r_buf);

  
  //Setup the data buffer
  this->setup(read_a2a_args_l,read_a2a_args_r,tl,tr);

  assert(read_fsize == fsize);

  if(!UniqueID()){
    std::istream &file = *file_ptr;
    //Node 0 finish reading data
    FP_FORMAT dataformat = FPformat<ScalarComplexType>::get();
    FPConv conv;
    
    conv.setHostFormat(dataformat);	
    FP_FORMAT fileformat = conv.setFileFormat(dformatbuf);    
    assert(conv.size(fileformat) == conv.size(dataformat));

    int dsize = conv.size(fileformat);

    std::string str;
    getline(file,str); assert(str == "BEGIN_DATA");

    static const int chunk = 32768; //32kb chunks
    assert(chunk % dsize == 0);
    int fdinchunk = chunk/dsize;
    char *rbuf = (char *)malloc(chunk * sizeof(char)); //leave room for auto null char
      
    char *dptr = (char *)mf;

    int off = 0;
    int nfd = 2*fsize;
    while(off < nfd){
      int grab = std::min(nfd-off, fdinchunk); //How many data elements to grab
      int grabchars = grab * dsize;

      file.read(rbuf,grabchars);
      int got = file.gcount();
      
      if(file.gcount() != grabchars)
	ERR.General("","","Only managed to read %d chars, needed %d\n",file.gcount(),grabchars);
      
      conv.file2host(dptr,rbuf,grab);

      off += grab;
      dptr += grabchars;
    }
    free(rbuf);

    file.ignore(1); //newline
    getline(file,str); assert(str == "END_DATA");
  }

  //Broadcast data
#ifdef USE_MPI
  ret = MPI_Bcast(mf, 2*fsize*sizeof(typename ScalarComplexType::value_type) , MPI_CHAR, head_mpi_rank, MPI_COMM_WORLD);
  if(ret != MPI_SUCCESS) ERR.General("A2AmesonField","read","Squirt data fail\n");
#endif
  
  //Every node do the checksum
  FPConv conv;
  FP_FORMAT dataformat = FPformat<ScalarComplexType>::get();
  unsigned int calc_cksum = conv.checksum((char*)mf, 2*fsize, dataformat);

  assert( calc_cksum == checksum );  
}



template<typename mf_Policies, template <typename> class A2AfieldL,  template <typename> class A2AfieldR>
void A2AmesonField<mf_Policies,A2AfieldL,A2AfieldR>::write(const std::string &filename, const std::vector<A2AmesonField<mf_Policies,A2AfieldL,A2AfieldR> > &mfs, FP_FORMAT fileformat){
  Float dtime = -dclock();
  if(!UniqueID()) printf("Writing meson-field vector of size %d to file %s\n",mfs.size(),filename.c_str());
  std::ofstream *file = !UniqueID() ? new std::ofstream(filename.c_str(),std::ofstream::out) : NULL;

  write(file,mfs,fileformat);

  if(!UniqueID()){
    file->close();
    delete file;
  }
    
  print_time("A2AmesonField","write Lt meson fields",dclock()+dtime);
}

template<typename mf_Policies, template <typename> class A2AfieldL,  template <typename> class A2AfieldR>
void A2AmesonField<mf_Policies,A2AfieldL,A2AfieldR>::write(std::ostream *file_ptr, const std::vector<A2AmesonField<mf_Policies,A2AfieldL,A2AfieldR> > &mfs, FP_FORMAT fileformat){
  if(!UniqueID()) assert(file_ptr != NULL);
  else assert(file_ptr == NULL);

  if(!UniqueID()){
    assert(!file_ptr->fail());
    file_ptr->exceptions ( std::ofstream::failbit | std::ofstream::badbit );

    (*file_ptr) << "BEGIN_MESONFIELD_VECTOR_HEADER\n";
    (*file_ptr) << "VECTOR_SIZE = " << int(mfs.size()) << '\n';
    (*file_ptr) << "END_MESONFIELD_VECTOR_HEADER\n";
    (*file_ptr) << "BEGIN_MESONFIELD_VECTOR_CONTENTS\n";
  }
    
  for(int i=0;i<mfs.size();i++)
    mfs[i].write(file_ptr,fileformat);

  if(!UniqueID())
    (*file_ptr) << "END_MESONFIELD_VECTOR_CONTENTS\n";
}

template<typename mf_Policies, template <typename> class A2AfieldL,  template <typename> class A2AfieldR>
void A2AmesonField<mf_Policies,A2AfieldL,A2AfieldR>::read(const std::string &filename, std::vector<A2AmesonField<mf_Policies,A2AfieldL,A2AfieldR> > &mfs){
  if(!UniqueID()) printf("Reading vector of meson fields from file %s\n",filename.c_str());
  std::ifstream *file = !UniqueID() ? new std::ifstream(filename.c_str()) : NULL;

  read(file,mfs);

  if(!UniqueID()){
    file->close();
    delete file;
  }
    
}

template<typename mf_Policies, template <typename> class A2AfieldL,  template <typename> class A2AfieldR>
void A2AmesonField<mf_Policies,A2AfieldL,A2AfieldR>::read(std::istream *file_ptr, std::vector<A2AmesonField<mf_Policies,A2AfieldL,A2AfieldR> > &mfs){
  if(!UniqueID()) assert(file_ptr != NULL);
  else assert(file_ptr == NULL);

  if(!UniqueID()){
    assert(!file_ptr->fail());
    file_ptr->exceptions ( std::ofstream::failbit | std::ofstream::badbit );
  }
#ifdef USE_MPI
  int my_mpi_rank = getMyMPIrank(); //Get this node's mpi rank  
  int head_mpi_rank = getHeadMPIrank(); //Broadcast to all nodes the mpi rank of the head node (UniqueID() == 0)
#else
  if(GJP.Xnodes()*GJP.Ynodes()*GJP.Znodes()*GJP.Tnodes()*GJP.Snodes() != 1) ERR.General("A2AmesonField","read(std::istream *, std::vector<A2AmesonField<mf_Policies,A2AfieldL,A2AfieldR> > &)","Parallel implementation requires MPI\n");
#endif
  
  int vsize;
  if(!UniqueID()){
    std::string str;
    
    getline(*file_ptr,str); assert(str == "BEGIN_MESONFIELD_VECTOR_HEADER");
    getline(*file_ptr,str); assert( sscanf(str.c_str(),"VECTOR_SIZE = %d",&vsize) == 1 );
    getline(*file_ptr,str); assert(str == "END_MESONFIELD_VECTOR_HEADER");
    getline(*file_ptr,str); assert(str == "BEGIN_MESONFIELD_VECTOR_CONTENTS");
  }
#ifdef USE_MPI
  int ret = MPI_Bcast(&vsize, 1, MPI_INT, head_mpi_rank, MPI_COMM_WORLD);
  if(ret != MPI_SUCCESS) ERR.General("A2AmesonField","read(vector)","Squirt 1 fail\n");
#endif
  mfs.resize(vsize);
  for(int i=0;i<vsize;i++)
    mfs[i].read(file_ptr);

  if(!UniqueID()){
    std::string str;
    getline(*file_ptr,str); assert(str == "END_MESONFIELD_VECTOR_CONTENTS");
  }
}

#endif
