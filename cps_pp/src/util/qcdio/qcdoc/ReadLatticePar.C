
// a slight modification on ReadLattice class (if not inheritance)
// to enable parallel reading/writing of "Gauge Connection Format" lattice data



// Read the format of Gauge Connection Format
// from QCDSP {load,unload}_lattice format

#include <config.h>
#include <math.h>
#include <util/gjp.h>
#include <sysfunc.h>
#include <util/ReadLatticePar.h>
#include <util/fpconv.h>
CPS_START_NAMESPACE

// GCFheaderPar class members
string elmSpacePar(string str)
{
  const int i0(str.find_first_not_of(" "));
  const int i1(str.find_last_not_of (" "));
  if(i1 - i0>0){ return(str.substr(i0,i1-i0+1)); }
  else         { return(str);  }
}

bool GCFheaderPar::add(string key_eq_value)
{
  const int eqp(key_eq_value.find("="));
  if( eqp  > 0  )
    {
      const string key( elmSpacePar( key_eq_value.substr(0,eqp) ) );
      const string val( elmSpacePar( key_eq_value.substr(eqp+1) ) );
      headerMap.insert(GCFHMapParT::value_type(key,val));
      return true;
    } 
  else 
    {
      return false;
    }
}


void GCFheaderPar::Show()
{
  for (GCFHMapParT::const_iterator iter = headerMap.begin(); 
       iter != headerMap.end(); ++iter) 
    {
      cout << iter->first << ":" << iter->second << endl;
    }
};


string GCFheaderPar::asString( const string key ) 
{
  GCFHMapParT::const_iterator n(headerMap.find(key));

  if (n == headerMap.end()) {
    cout << "header::asString key " << key << " not found. use Default." << endl;
    //exit(-1);
    prevFound = false;
    return string("");
  }
  else {
    prevFound = true;
    return ( n->second );
  }
}		  


int GCFheaderPar::asInt( const string key ) 
{
  GCFHMapParT::const_iterator n(headerMap.find(key));

  if (n == headerMap.end()) {
    cout << "header::asInt key "<<key<<" not found. use Default." << endl;
    //exit(-1);
    prevFound = false;
    return int(0);
  }

  else {
    prevFound = true;
    int tmp;
    sscanf((n->second).c_str() , "%d ", &tmp);
    return ( tmp );
  }
}		  

Float GCFheaderPar::asFloat( const string key ) 
{
  GCFHMapParT::const_iterator n(headerMap.find(key));

  if (n == headerMap.end()) {
    cout << "header::asFloat key " << key << " not found. use Default." << endl;
    //    exit(-1);
    prevFound = false;
    return Float(0.0);
  }
  else {
    prevFound = true;
    float tmp;
    sscanf((n->second).c_str() , "%f ", &tmp);
    return ( (Float) tmp );
  }
}


//////////////////////////////////////////
// ReadLatticeParallel class members

// checksum check
// modification added to enable parallel io
unsigned ReadLatticeParallel::calc_csum(char * fpoint)
{
  unsigned int csum(0);
  //  unsigned int *ip = (unsigned int*)fpoint;

  int size_ints ;
  const int size_matrices( do_arg.x_node_sites * do_arg.y_node_sites 
			   * do_arg.z_node_sites * do_arg.t_node_sites * 4); 
  if ( recon_row_3 ) {
    size_ints = 12*size_matrices;
  }
  else { 
    size_ints = 18*size_matrices;
  }

  csum = pc.globalSumUint(fpconv.checksum(fpoint,size_ints));

  if(pc.uniqueID() == 0) {
    unsigned long tmp;
    sscanf(hd.asString("CHECKSUM").c_str() , "%lx ", &tmp);
  
    if( tmp != csum ) {
      cout << "CheckSUM error !! Header:" 
	   << hd.asString("CHECKSUM") << " Host calc:"
	   <<hex << csum << dec << "\n";
      //exit(-13);
      error = 1;
    }
    else
      cout << "CheckSUM is ok\n";
  }

  if(pc.synchronize(error) != 0)  exit(-13);

  /*
  // Restore the byte ordering (will reverse again later before calculation)
  if(host_endian != file_endian)
    byterevn((type32*)fpoint,size_ints);
  */

  return csum ;

}


// modification added to enable parallel NFS
void ReadLatticeParallel::read( const char* file )
{
  // most codes coped from ReadLattice::read( )
  
  cout << endl << "Loading lattice..." << endl << endl;

  if ( allocated ) { delete[] lpoint; allocated=false; }

  // all open file and check error
  ifstream input(file);
  if ( !input.good() )
    {
      cout << "Could not open file:\n   "
	   << file
           << "\nfor input.\nNot even a little bit.\n";
      //exit(13);
      error = 1;
    }
  // first sync point
  // executed by all, sync and share error status information
  if(pc.synchronize(error) != 0)   exit(13);

  if (pc.uniqueID() == 0) { // commander, analyze file header
    string line;
    
    do {
      getline(input,line); // read one line
      hd.add(line);
    } while( line.find("END_HEADER") == string::npos);
    
    data_start = input.tellg();
  }
  pc.broadcastInt(&data_start);


  if(pc.uniqueID() == 0) {
    // Only implemented for 4D SU3 3x3
    //    int recon_row_3 = 0;
    recon_row_3 = 0;
    if(hd.asString("DATATYPE") != "4D_SU3_GAUGE_3x3"){
      cout << "ReadLatticeParallel: will reconstruct row 3 " 
	   << hd.asString("DATATYPE") << endl;
      recon_row_3 = 1;
    }
  }
  pc.broadcastInt(&recon_row_3);
  cout << "recon_row_3 = " << recon_row_3 << endl;


  int nx,ny,nz,nt;

  if(pc.uniqueID()==0) {
    //====================================
    // initialise the size of the lattice
    // int characters
    
    nx=hd.asInt("DIMENSION_1");
    ny=hd.asInt("DIMENSION_2");
    nz=hd.asInt("DIMENSION_3");
    nt=hd.asInt("DIMENSION_4");
  }
  pc.broadcastInt(&nx);
  pc.broadcastInt(&ny);
  pc.broadcastInt(&nz);
  pc.broadcastInt(&nt);

  if(GJP.Xnodes() * GJP.XnodeSites() != nx || GJP.Ynodes() * GJP.YnodeSites() != ny 
     || GJP.Znodes() * GJP.ZnodeSites() != nz || GJP.Tnodes() * GJP.TnodeSites() != nt) {
    cout << "GlobalJobParameter setting wrong!"<<endl;
    exit(-13);
  }

  // executed by all
  do_arg.x_nodes = pc.Xnodes(); 
  do_arg.y_nodes = pc.Ynodes();
  do_arg.z_nodes = pc.Znodes();
  do_arg.t_nodes = pc.Tnodes();
  do_arg.s_nodes = 1;

  do_arg.x_node_sites = nx/pc.Xnodes();
  do_arg.y_node_sites = ny/pc.Ynodes();
  do_arg.z_node_sites = nz/pc.Znodes();
  do_arg.t_node_sites = nt/pc.Tnodes();

//??    fermion's boundary condition or gauge ??
  if(pc.uniqueID() == 0) {
    // "" means not present in file
    do_arg.x_bc = (hd.asString("BOUNDARY_1") == "ANTIPERIODIC")
                 ? BND_CND_APRD : BND_CND_PRD;
    do_arg.y_bc = (hd.asString("BOUNDARY_2") == "ANTIPERIODIC")
                 ? BND_CND_APRD : BND_CND_PRD;
    do_arg.z_bc = (hd.asString("BOUNDARY_3") == "ANTIPERIODIC")
                 ? BND_CND_APRD : BND_CND_PRD;
    do_arg.t_bc = (hd.asString("BOUNDARY_4") == "ANTIPERIODIC")
                 ? BND_CND_APRD : BND_CND_PRD;
  }
  pc.broadcastInt((int*)&do_arg.x_bc);
  pc.broadcastInt((int*)&do_arg.y_bc);
  pc.broadcastInt((int*)&do_arg.z_bc);
  pc.broadcastInt((int*)&do_arg.t_bc);

  cout << "X bc:" << (do_arg.x_bc==BND_CND_PRD ? "PERIODIC":"ANTI-PERIODIC") << endl;
  cout << "Y bc:" << (do_arg.y_bc==BND_CND_PRD ? "PERIODIC":"ANTI-PERIODIC") << endl;
  cout << "Z bc:" << (do_arg.z_bc==BND_CND_PRD ? "PERIODIC":"ANTI-PERIODIC") << endl;
  cout << "T bc:" << (do_arg.t_bc==BND_CND_PRD ? "PERIODIC":"ANTI-PERIODIC") << endl;

  // check floating point format
  if(pc.uniqueID() == 0) {
    fpconv.setFileFormat(hd.asString("FLOATING_POINT").c_str());
  }
  pc.broadcastInt((int*)&fpconv.fileFormat);
  if(fpconv.fileFormat == FP_UNKNOWN) {
    cout << "Data file Floating Point format UNKNOWN" << endl;
    exit(13);
  }

  //  const int size_matrices(nx*ny*nz*nt*4);
  int size_matrices( do_arg.x_node_sites * do_arg.y_node_sites 
		     * do_arg.z_node_sites * do_arg.t_node_sites * 4); 

  int size_ints = recon_row_3 ? size_matrices*12 : size_matrices*18;
  int size_chars   ( size_ints * fpconv.fileFpSize()  );

  cout << "size_matrices = " << size_matrices << endl;


  lpoint = new Matrix[size_matrices]; allocated = true;

  cout << "ReadLatticeParallel::  gauge field allocated at " 
       << hex << lpoint << "  size " << dec << size_matrices << endl << endl;

  do_arg.start_conf_kind      = START_CONF_LOAD ;
  do_arg.start_conf_load_addr = lpoint;
  do_arg.start_seed_kind      = START_SEED_UNIFORM ;
  do_arg.start_seed_value     = 5296;
  //  do_arg.colors               = 3; 
  //  do_arg.verbose_level        = 3;

  
  if(pc.uniqueID() == 0)  hd.Show();

  // start reading data
  int size_per_mat = recon_row_3 ? 12*fpconv.fileFpSize() : 18*fpconv.fileFpSize();
  int size_per_site = size_per_mat * 4;

  char *fpoint = new char [size_chars];

  int xbegin = do_arg.x_node_sites * pc.Xcoor(), xend = do_arg.x_node_sites * (pc.Xcoor()+1);
  int ybegin = do_arg.y_node_sites * pc.Ycoor(), yend = do_arg.y_node_sites * (pc.Ycoor()+1);
  int zbegin = do_arg.z_node_sites * pc.Zcoor(), zend = do_arg.z_node_sites * (pc.Zcoor()+1);
  int tbegin = do_arg.t_node_sites * pc.Tcoor(), tend = do_arg.t_node_sites * (pc.Tcoor()+1);

  int tblk = nx*ny*nz*size_per_site;
  int zblk = nx*ny*size_per_site;
  int yblk = nx*size_per_site;

  cout << endl;
  cout << "Trying to read " << size_chars << " bytes"<<endl;

  // read in parallel manner, node 0 will assign & dispatch IO time slots
  pc.getIOTimeSlot();

  char *buf = fpoint;
  input.seekg(data_start,ios_base::beg);

  int jump = tbegin * tblk;
  for(int tr=tbegin;tr<tend;tr++) {
    jump += zbegin * zblk;
    for(int zr=zbegin;zr<zend;zr++) {
      jump += ybegin * yblk;
      for(int yr=ybegin;yr<yend;yr++) {
	jump += xbegin * size_per_site;
	input.seekg(jump,ios_base::cur);

	input.read(buf,(xend-xbegin) * size_per_site);
	buf += (xend-xbegin) * size_per_site;

	jump = (nx-xend) * size_per_site;  // "jump" restart from 0 and count
      }
      jump += (ny-yend) * yblk;
    }
    jump += (nz-zend) * zblk;
  }
  if ( !input.good() ) { cout << "blarg!\n"; error = 1; }

  pc.finishIOTimeSlot();
  //
  
  input.close();

  cout << "Actually read " << buf-fpoint << " bytes" << endl;

  if(pc.synchronize(error) != 0)  exit(-13);

  // STEP 1: checksum
  calc_csum(fpoint) ;


  /*
  // STEP 2: endian
  int fileBigEndian;
  if(pc.uniqueID() == 0) {
    if((hd.asString("FLOATING_POINT")=="IEEE32BIG")||
       (hd.asString("FLOATING_POINT")=="IEEE32")||
       (hd.asString("FLOATING_POINT")=="TIDSP32"))
      {
	cout<<"Data file Big Endian"<<endl ;
	fileBigEndian = 1;
      }
    else
      {
	cout << "Data file Little Endian" << endl;
	fileBigEndian = 0;
      }
  }
  pc.broadcastInt(&fileBigEndian);

  // by all
  int xendian=1;
  if( *(char *)&xendian == 1) 
    {
      cout << "Host Little Endian\n";
      if(fileBigEndian)
	{
	  cout<<"Lattice  needs byte reversal"<<endl ;
	  byterevn((type32*)fpoint,size_ints);
	}
    }
  else
    {
      cout << "Host Big Endian\n";
      if(!fileBigEndian)
	{
	  cout<<"Lattice needs byte reversal"<<endl ;
	  byterevn((type32*)fpoint,size_ints);
	}
    }
 

  // STEP 3: DSP => IEEE
  //====================================
  // check if dsp floating point format
  int fileIsTIDSP32 = 0;
  if(pc.uniqueID()==0) {
    if(hd.asString("FLOATING_POINT") == "TIDSP32") {
      fileIsTIDSP32 = 1;
    }
  }
  pc.broadcastInt(&fileIsTIDSP32);

  // by all
  if(fileIsTIDSP32) {
    cout << "Formating TIDSP32 => IEEE" << endl;
    ti2ieee((type32*)fpoint, size_ints);
  }


  ieee2cpsFloat(lpoint,fpoint,size_matrices,recon_row_3);
  */

  cout << "Transferring " << fpconv.name(fpconv.fileFormat) << " ==> " << fpconv.name(fpconv.hostFormat) << endl;
  
  if(recon_row_3) {
    cout << "Reconstructing row 3" << endl;
    for(int mat=0; mat<size_matrices; mat++) {
      char * src = fpoint + mat * size_per_mat;
      Float * rec = (Float*)&lpoint[mat];
      fpconv.file2host((char*)rec,src,12);
      // reconstruct the 3rd row
      rec[12] =  rec[2] * rec[10] - rec[3] * rec[11] - rec[4] * rec[8] + rec[5] * rec[9];
      rec[13] = -rec[2] * rec[11] - rec[3] * rec[10] + rec[4] * rec[9] + rec[5] * rec[8];
      rec[14] = -rec[0] * rec[10] + rec[1] * rec[11] + rec[4] * rec[6] - rec[5] * rec[7];
      rec[15] =  rec[0] * rec[11] + rec[1] * rec[10] - rec[4] * rec[7] - rec[5] * rec[6];
      rec[16] =  rec[0] * rec[ 8] - rec[1] * rec[ 9] - rec[2] * rec[6] + rec[3] * rec[7];
      rec[17] = -rec[0] * rec[ 9] - rec[1] * rec[ 8] + rec[2] * rec[7] + rec[3] * rec[6];
    }
  }
  else {
    fpconv.file2host((char*)lpoint,fpoint,size_ints);
  }

  delete[] fpoint;
  



  // plaq, linktrace in header
  if(pc.uniqueID()==0) {
    _plaq_inheader = hd.asFloat("PLAQUETTE");
    _linktrace_inheader = hd.asFloat("LINK_TRACE");
  }
  pc.broadcastFloat(&_plaq_inheader);
  pc.broadcastFloat(&_linktrace_inheader);    
};

// just add a extra level of sum (the global sum)
void ReadLatticeParallel::CheckPlaqLinktrace(Lattice &lattice, Float chkprec) 
{
  Float plaq = lattice.SumReTrPlaq() / 18.0 / GJP.VolSites() ;
  Float devplaq =   fabs(  (plaq - plaqInHeader()) / plaq ) ;
  printf("plaqerr::  calc: %f  header: %f dev.: %f\n",
         (float)plaq, (float)plaqInHeader(),(float)devplaq);

  if(devplaq > chkprec) {
    //       exit(-13);
    //    return;
  }


  Float linktrace(0);
  int is;
  Matrix *m =  lattice.GaugeField(); 
  for(is=0;is< GJP.VolNodeSites()*4; is++){
    linktrace += m->ReTr();
    m++;
  }
  linktrace = pc.globalSumFloat(linktrace) / (GJP.VolSites()*12.0);
  Float devlinktrace =   
    fabs(  (linktrace - linktraceInHeader()) / linktrace );

  printf("linktrace::  calc: %f  header: %f dev.: %f\n",
         (float) linktrace, (float)linktraceInHeader(),
         (float)devlinktrace);
  
  if(devlinktrace > chkprec) {
    //      exit(-13);
    return;
  }
}


CPS_END_NAMESPACE
