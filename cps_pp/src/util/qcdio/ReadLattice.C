

// Read the format of Gauge Connection Format
// from QCDSP {load,unload}_lattice format

#include <config.h>
#include <math.h>
#include <util/ReadLattice.h>
#include <util/Flconv.h>
CPS_START_NAMESPACE


// checksum check
unsigned ReadLattice::calc_csum(float *fpoint)
{
  unsigned int csum(0);
  unsigned int *ip = (unsigned int*)fpoint;

  int size_ints ;
  if ( recon_row_3 ) {
    size_ints = 4*12*
    hd.asInt("DIMENSION_1")*
    hd.asInt("DIMENSION_2")*
    hd.asInt("DIMENSION_3")*
    hd.asInt("DIMENSION_4");
  } else { 
    size_ints = 4*18*
    hd.asInt("DIMENSION_1")*
    hd.asInt("DIMENSION_2")*
    hd.asInt("DIMENSION_3")*
    hd.asInt("DIMENSION_4");
  }

  char end_check[4] = {0,0,0,1};
  unsigned long *lp = (unsigned long *)end_check;
  int host_endian;

  if ( *lp == 0x1 ) { 
    cout << "Host is big-endian\n";
    host_endian = 0;
  } else {
    cout << "Host is little-endian\n";
    host_endian = 1;
  }

  int BigEndian(0) ;

  if((hd.asString("FLOATING_POINT")=="IEEE32BIG")||
     (hd.asString("FLOATING_POINT")=="IEEE32")) 
    BigEndian = 1 ;

  if ( host_endian ) BigEndian = 1-BigEndian;

  // checksum check has to be done in BIG-endian
  // according to  NERC specification
  if(!BigEndian) byterevn((type32*)fpoint,size_ints);

  int i;
  for(i=0; i < size_ints;i++)  csum += ip[i];
  
  unsigned int tmp;
  sscanf(hd.asString("CHECKSUM").c_str() , "%x ", &tmp);
  
  if( tmp != csum ) {
    cout << "CheckSUM error !! " 
	 << hd.asString("CHECKSUM") << " "
	 <<hex << csum << "\n";
    exit(-13);
  }
  cout << "CheckSUM is ok\n";
  // Restore the byte ordering
  if ( host_endian ) BigEndian = 1-BigEndian;
  if(!BigEndian) byterevn((type32*)fpoint,size_ints);

  return csum ;
}

string elmSpace(string str)
{
  const int i0(str.find_first_not_of(" "));
  const int i1(str.find_last_not_of (" "));
  if(i1 - i0>0){ return(str.substr(i0,i1-i0+1)); }
  else         { return(str);  }
}

bool GCFheader::add(string key_eq_value)
{
  const int eqp(key_eq_value.find("="));
  if( eqp  > 0  )
    {
      const string key( elmSpace( key_eq_value.substr(0,eqp) ) );
      const string val( elmSpace( key_eq_value.substr(eqp+1) ) );
      headerMap.insert(GCFHMapT::value_type(key,val));
      return true;
    } 
  else 
    {
      return false;
    }
}


void GCFheader::Show()
{
  for (GCFHMapT::const_iterator iter = headerMap.begin(); 
       iter != headerMap.end(); ++iter) 
    {
      cout << iter->first << ":" << iter->second << endl;
    }
};


string GCFheader::asString( const string key ) 
{
  GCFHMapT::const_iterator n(headerMap.find(key));

  if (n == headerMap.end()) {
    cerr << "header::asString key not found " << key << endl;
    exit(-1);
  }
  else {
    return ( n->second );
  }
}		  


int GCFheader::asInt( const string key ) 
{
  GCFHMapT::const_iterator n(headerMap.find(key));

  if (n == headerMap.end()) {
    cerr << "header::asInt key not found " << key << endl;
    exit(-1);
  }

  else {
    int tmp;
    sscanf((n->second).c_str() , "%d ", &tmp);
    return ( tmp );
  }
}		  

Float GCFheader::asFloat( const string key ) 
{
  GCFHMapT::const_iterator n(headerMap.find(key));

  if (n == headerMap.end()) {
    cerr << "header::asFloat key not found " << key << endl;
    exit(-1);
  }
  else {
    float tmp;
    sscanf((n->second).c_str() , "%f ", &tmp);
    return ( (Float) tmp );
  }
}


ReadLattice::ReadLattice(const char* file )
{
    allocated =false;
//    read(file);
}

void ReadLattice::read( const char* file )
{
  if ( allocated ) { delete[] lpoint; allocated=false; }
  ifstream input(file);

  if ( input.bad() )
    {
      cerr << "Could not open file:\n   "
	   << file
           << "\nfor input.\nNot even a little bit.\n";
      exit(13);
    }

  string line;
  
  do {
      getline(input,line); // read one line
      hd.add(line);
  } while( line.find("END_HEADER") == string::npos);


  // Only implemented for 4D SU3 3x3
  int recon_row_3 = 0;
  if(hd.asString("DATATYPE") != "4D_SU3_GAUGE_3x3"){
    cout << "ReadLattice: reconstructing row 3 " 
	 << hd.asString("DATATYPE") << endl;
    recon_row_3 = 1;
  }


  //====================================
  // initialise the size of the lattice
  // int characters
  
  const int nx(hd.asInt("DIMENSION_1"));
  const int ny(hd.asInt("DIMENSION_2"));
  const int nz(hd.asInt("DIMENSION_3"));
  const int nt(hd.asInt("DIMENSION_4"));
  
  do_arg.x_nodes = 1;
  do_arg.y_nodes = 1;
  do_arg.z_nodes = 1;
  do_arg.t_nodes = 1;
  do_arg.s_nodes = 1;

  do_arg.x_node_sites = nx;
  do_arg.y_node_sites = ny;
  do_arg.z_node_sites = nz;
  do_arg.t_node_sites = nt;

//??    fermion's boundary condition or gauge ??
  /*
  do_arg.x_bc = (hd.asString("BOUNDARY_1") == "PERIODIC")
                 ? BND_CND_PRD : BND_CND_APRD;
  do_arg.y_bc = (hd.asString("BOUNDARY_2") == "PERIODIC")
                 ? BND_CND_PRD : BND_CND_APRD;
  do_arg.z_bc = (hd.asString("BOUNDARY_3") == "PERIODIC")
                 ? BND_CND_PRD : BND_CND_APRD;
  do_arg.t_bc = (hd.asString("BOUNDARY_4") == "PERIODIC")
                 ? BND_CND_PRD : BND_CND_APRD;
  */
  do_arg.x_bc = BND_CND_PRD;
  do_arg.y_bc = BND_CND_PRD;
  do_arg.z_bc = BND_CND_PRD;
  do_arg.t_bc = BND_CND_PRD;

  do_arg.start_conf_kind      = START_CONF_LOAD ;

  const int size_matrices(nx*ny*nz*nt*4);
  const int size_ints    ( size_matrices * ( 9 * 2 ) );
  const int size_chars   ( size_ints * sizeof(float)  );
  int size_float = size_matrices*6*2;

  lpoint = new Matrix[size_matrices]; allocated = true;

  cout << "ReadLattice::  gauge field allocated at " 
       << hex << lpoint << "  size " << dec << size_matrices << "\n";

  do_arg.start_conf_load_addr = lpoint;
  do_arg.start_seed_kind      = START_SEED_UNIFORM ;
  do_arg.start_seed_value     = 5296;
  do_arg.colors               = 3; 
  do_arg.verbose_level        = 3;

  
  hd.Show();

  float *fpoint = new float [size_matrices*2*9];
  if ( recon_row_3 ) { 
    bzero((char *)fpoint,size_chars);
    input.read((char*)fpoint,size_float*sizeof(float));
    if ( !input.good() ) { cout << "blarg!\n"; }
    calc_csum(fpoint) ;
    Copy(lpoint,fpoint,size_matrices,recon_row_3);
  } else { 
    input.read((char*)fpoint,size_chars);
    if ( !input.good() ) { cout << "blarg!\n"; }
    //check checksums
    calc_csum(fpoint) ;
    Copy(lpoint,fpoint,size_matrices,recon_row_3);
  }
  delete[] fpoint;
  

/*
  if ( recon_row_3 ) { 
    cout << "Reconstructing third row" << endl;
    float *fpoint= new float[size_float];
    bcopy(lpoint,fpoint,size_float*sizeof(float));
    ReconRow3(lpoint,fpoint,size_matrices);
    delete [] fpoint;
  }
*/

  // switch into correct byte ordering (ieee)
  int xend=1;
  if( *(char *)&xend == 1) 
    {
      //Little Endian
      cout << "Little Endian\n";
      
      if((hd.asString("FLOATING_POINT")=="IEEE32BIG")||
	 (hd.asString("FLOATING_POINT")=="TIDSP32"))
	{
	  cout<<"Lattice  needs byte reversal"<<endl ;
	  byterevn((type32*)lpoint,size_ints);
	}
    }
  else
    {
      //Big Endian
      cout << "Big Endian\n";
      if (hd.asString("FLOATING_POINT")=="IEEE32LITTLE")
	{
	  cout<<"Lattice needs byte reversal"<<endl ;
	  byterevn((type32*)lpoint,size_ints);
	}
    }
 

  //====================================
  // check if dsp floating point format
  if(hd.asString("FLOATING_POINT") == "TIDSP32") {
    type32 *lintp; 
    lintp = (type32*) lpoint;
    ti2ieee( lintp, size_ints );
  }

  // plaq, linktrace in header
  _plaq_inheader = hd.asFloat("PLAQUETTE");
  _linktrace_inheader = hd.asFloat("LINK_TRACE");

};

void ReadLattice::CheckPlaqLinktrace(Lattice &lattice, Float chkprec) 
{

  Float plaq = lattice.SumReTrPlaq() / 18.0 / GJP.VolSites() ;
  Float devplaq =   fabs(  (plaq - plaqInHeader()) / plaq ) ;
  printf("plaqerr::  calc: %f  header: %f dev.: %f\n",
         (float)plaq, (float)plaqInHeader(),(float)devplaq);

  if(devplaq > chkprec) {
       exit(-13);
  }

  Float linktrace(0);
  int is;
  Matrix *m =  lattice.GaugeField(); 
  for(is=0;is< GJP.VolSites()*4; is++){
    linktrace += m->ReTr();
    m++;
  }
  linktrace /= (GJP.VolSites()*12.0);
  Float devlinktrace =   
    fabs(  (linktrace - linktraceInHeader()) / linktrace );

  printf("linktrace::  calc: %f  header: %f dev.: %f\n",
         (float) linktrace, (float)linktraceInHeader(),
         (float)devlinktrace);
  
  if(devlinktrace > chkprec) {
      exit(-13);
  }
}


void ReadLattice::Copy(Matrix *lpoint,float *fpoint,int
size_matrices,int recon_row_3)
{
  int mat,j;
  float  *data;
  Float *rec;
  for ( mat = 0; mat <size_matrices; mat++ ) { 
    data = &fpoint[mat * 6*2 ];
    rec  = (Float *)&lpoint[mat];
    for (j = 0; j< 2*6; j++) rec[j] = (Float)data[j];
    if( recon_row_3){
      rec[12] =  rec[2] * rec[10] - rec[3] * rec[11] - rec[4] * rec[8] + rec[5] * rec[9];
      rec[13] = -rec[2] * rec[11] - rec[3] * rec[10] + rec[4] * rec[9] + rec[5] * rec[8];
      rec[14] = -rec[0] * rec[10] + rec[1] * rec[11] + rec[4] * rec[6] - rec[5] * rec[7];
      rec[15] =  rec[0] * rec[11] + rec[1] * rec[10] - rec[4] * rec[7] - rec[5] * rec[6];
      rec[16] =  rec[0] * rec[ 8] - rec[1] * rec[ 9] - rec[2] * rec[6] + rec[3] * rec[7];
      rec[17] = -rec[0] * rec[ 9] - rec[1] * rec[ 8] + rec[2] * rec[7] + rec[3] * rec[6];
    }
#if 0
    for(int i = 0;i<3;i++)
    for(int j = 0;j<3;j++){
      float sumf=0.0,sumi=0.0;
      for(int k = 0; k<3;k++){
        sumf += rec[(i*3+k)*2] * rec[(j*3+k)*2];
        sumf += rec[(i*3+k)*2+1] * rec[(j*3+k)*2+1];
        sumi += rec[(i*3+k)*2] * rec[(j*3+k)*2+1];
        sumi -= rec[(i*3+k)*2+1] * rec[(j*3+k)*2];
      }
      printf("%e %e ",sumf,sumi);
      if(j==2) printf("\n");
    }
#endif
  }
}
CPS_END_NAMESPACE
