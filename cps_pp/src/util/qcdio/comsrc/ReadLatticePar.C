
// a slight modification on ReadLattice class (if not inheritance)
// to enable parallel reading/writing of "Gauge Connection Format" lattice data



// Read the format of Gauge Connection Format
// from QCDSP {load,unload}_lattice format

#include <config.h>

#include <math.h>
#include <util/gjp.h>
#include <util/ReadLatticePar.h>
#include <util/fpconv.h>
#include <util/time.h>
CPS_START_NAMESPACE

#define PROFILE
void ReadLatticeParallel::read(Lattice & lat, const QioArg & rd_arg)
{
  // most codes coped from ReadLattice::read( ), modification added to enable parallel IO
  cout << endl << "Loading lattice..." << endl << endl;

#ifdef PROFILE
  struct timeval start,end;
  gettimeofday(&start,NULL);
#endif

  load_good = false;
  int error = 0;

  // all open file and check error
  ifstream input(rd_arg.FileName);
  if ( !input.good() )
    {
      cout << "Could not open file:\n   "
	   << rd_arg.FileName
           << "\nfor input.\nNot even a little bit.\n";
      error = 1;
    }
  // first sync point
  // executed by all, sync and share error status information
  if(synchronize(error) != 0)   return;


  int data_start;
  if (isRoot()) { // commander, analyze file header
    string line;
    do {
      getline(input,line); // read one line
      hd.add(line);
    } while( line.find("END_HEADER") == string::npos);
    
    data_start = input.tellg();
  }
  broadcastInt(&data_start);


  int recon_row_3 = 0;
  if(isRoot()) {
    if(hd.asString("DATATYPE") != "4D_SU3_GAUGE_3x3"){
      cout << "ReadLatticeParallel: will reconstruct row 3 " 
	   << hd.asString("DATATYPE") << endl;
      recon_row_3 = 1;
    }
  }
  broadcastInt(&recon_row_3);
  cout << "recon_row_3 = " << recon_row_3 << endl;


  // check dimensions, b.c, etc
  int nx = rd_arg.Xnodes() * rd_arg.XnodeSites();
  int ny = rd_arg.Ynodes() * rd_arg.YnodeSites();
  int nz = rd_arg.Znodes() * rd_arg.ZnodeSites();
  int nt = rd_arg.Tnodes() * rd_arg.TnodeSites();

  if(isRoot()) {
    int rdnx,rdny,rdnz,rdnt;

    //====================================
    // initialise the size of the lattice
    // int characters
    
    rdnx=hd.asInt("DIMENSION_1");
    rdny=hd.asInt("DIMENSION_2");
    rdnz=hd.asInt("DIMENSION_3");
    rdnt=hd.asInt("DIMENSION_4");

    cout << "Data dimensions: " << rdnx <<"x" << rdny <<"x" << rdnz <<"x" << rdnt << endl;
    if(rdnx != nx || rdny != ny || rdnz != nz || rdnt != nt) {
      cout << "Dimensions in file Contradicts GlobalJobParameter!"<<endl;
      error = 1;
    }

    int x_bc,y_bc,z_bc,t_bc;

    // "" means not present in file
    x_bc = (hd.asString("BOUNDARY_1") == "ANTIPERIODIC")
                 ? BND_CND_APRD : BND_CND_PRD;
    y_bc = (hd.asString("BOUNDARY_2") == "ANTIPERIODIC")
                 ? BND_CND_APRD : BND_CND_PRD;
    z_bc = (hd.asString("BOUNDARY_3") == "ANTIPERIODIC")
                 ? BND_CND_APRD : BND_CND_PRD;
    t_bc = (hd.asString("BOUNDARY_4") == "ANTIPERIODIC")
                 ? BND_CND_APRD : BND_CND_PRD;

    cout << "File specified:" << endl;
    cout << "X bc:" << (x_bc==BND_CND_PRD ? "PERIODIC":"ANTI-PERIODIC") << endl;
    cout << "Y bc:" << (y_bc==BND_CND_PRD ? "PERIODIC":"ANTI-PERIODIC") << endl;
    cout << "Z bc:" << (z_bc==BND_CND_PRD ? "PERIODIC":"ANTI-PERIODIC") << endl;
    cout << "T bc:" << (t_bc==BND_CND_PRD ? "PERIODIC":"ANTI-PERIODIC") << endl;
    
    if(x_bc != rd_arg.Xbc() || y_bc != rd_arg.Ybc() || z_bc != rd_arg.Zbc() || t_bc != rd_arg.Tbc()) {
      cout << "Boundary conditions in file DISAGREE with GlobalJobParameter!" << endl;
      error = 1;
    }
  }
  if(synchronize(error) != 0)  return;

  // check floating point format
  if(isRoot()) {
    fpconv.setFileFormat(hd.asString("FLOATING_POINT").c_str());
  }
  broadcastInt((int*)&fpconv.fileFormat);
  if(fpconv.fileFormat == FP_UNKNOWN) {
    cout << "Data file Floating Point format UNKNOWN" << endl;
    return;
  }


  int size_matrices( rd_arg.XnodeSites() * rd_arg.YnodeSites() 
		     * rd_arg.ZnodeSites() * rd_arg.TnodeSites() * 4); 

  cout << "size_matrices = " << size_matrices << endl;

  cout << "ReadLatticeParallel::  gauge field will be stored at " 
       << hex << rd_arg.StartConfLoadAddr << "  size " << dec << size_matrices << endl << endl;

  if(isRoot())  hd.Show();

  // start reading data
  //  int size_per_mat = recon_row_3 ? 12*fpconv.fileFpSize() : 18*fpconv.fileFpSize();
  //  int size_per_site = size_per_mat * 4;
  int Floats_per_site = recon_row_3? 4*12 : 4*18;
  int chars_per_site  = Floats_per_site * fpconv.fileFpSize();


  // TempBufAlloc is a Mem Allocator that prevents mem leak on "return"s
  TempBufAlloc fbuf(chars_per_site);  // buffer only stores one site


  int xbegin = rd_arg.XnodeSites() * rd_arg.Xcoor(), xend = rd_arg.XnodeSites() * (rd_arg.Xcoor()+1);
  int ybegin = rd_arg.YnodeSites() * rd_arg.Ycoor(), yend = rd_arg.YnodeSites() * (rd_arg.Ycoor()+1);
  int zbegin = rd_arg.ZnodeSites() * rd_arg.Zcoor(), zend = rd_arg.ZnodeSites() * (rd_arg.Zcoor()+1);
  int tbegin = rd_arg.TnodeSites() * rd_arg.Tcoor(), tend = rd_arg.TnodeSites() * (rd_arg.Tcoor()+1);

  int tblk = nx*ny*nz*chars_per_site;
  int zblk = nx*ny*chars_per_site;
  int yblk = nx*chars_per_site;

  // read in parallel manner, node 0 will assign & dispatch IO time slots
  int mat=0;
  unsigned int csum = 0;
  Matrix * lpoint = rd_arg.StartConfLoadAddr;

  setConcurIONumber(rd_arg.ConcurIONumber);
  getIOTimeSlot();

  input.seekg(data_start,ios_base::beg);

  int jump = tbegin * tblk;
  for(int tr=tbegin;tr<tend;tr++) {
    jump += zbegin * zblk;
    for(int zr=zbegin;zr<zend;zr++) {
      jump += ybegin * yblk;
      for(int yr=ybegin;yr<yend;yr++) {
	jump += xbegin * chars_per_site;
	input.seekg(jump,ios_base::cur);

	for(int xr=xbegin;xr<xend;xr++) {
	  input.read(fbuf,chars_per_site);
	  csum += fpconv.checksum(fbuf,Floats_per_site);
	  for(int i=0;i<4;i++) {
	    fpconv.file2host((char*)&lpoint[mat++],(char*)fbuf+chars_per_site/4*i,
			     Floats_per_site/4);
	  }
	}

	jump = (nx-xend) * chars_per_site;  // "jump" restart from 0 and count
      }
      jump += (ny-yend) * yblk;
    }
    jump += (nz-zend) * zblk;
  }
  if ( !input.good() ) { cout << "Input stream error!" << endl; error = 1; }

  finishIOTimeSlot();
  //
  
  input.close();

  if(synchronize(error) != 0)  return;

  // STEP 1: checksum
  if(rd_arg.Scoor() == 0)
    csum = globalSumUint(csum);
  else
    globalSumUint(0);

  if(isRoot()) {
    unsigned int cshead;
    sscanf(hd.asString("CHECKSUM").c_str() , "%lx ", &cshead);
  
    if( cshead != csum ) {
      cout << "CheckSUM error !! Header:" 
	   << hd.asString("CHECKSUM") << " Host calc:"
	   <<hex << csum << dec << "\n";
      error = 1;
    }
    else
      cout << "CheckSUM is ok\n";
  }

  if(synchronize(error) != 0) return;


  // STEP 2: reconstruct row 3
  if(recon_row_3) {
    cout << "Reconstructing row 3" << endl;
    for(int mat=0; mat<size_matrices; mat++) {
      Float * rec = (Float*)&lpoint[mat];
      // reconstruct the 3rd row
      rec[12] =  rec[2] * rec[10] - rec[3] * rec[11] - rec[4] * rec[8] + rec[5] * rec[9];
      rec[13] = -rec[2] * rec[11] - rec[3] * rec[10] + rec[4] * rec[9] + rec[5] * rec[8];
      rec[14] = -rec[0] * rec[10] + rec[1] * rec[11] + rec[4] * rec[6] - rec[5] * rec[7];
      rec[15] =  rec[0] * rec[11] + rec[1] * rec[10] - rec[4] * rec[7] - rec[5] * rec[6];
      rec[16] =  rec[0] * rec[ 8] - rec[1] * rec[ 9] - rec[2] * rec[6] + rec[3] * rec[7];
      rec[17] = -rec[0] * rec[ 9] - rec[1] * rec[ 8] + rec[2] * rec[7] + rec[3] * rec[6];
    }
  }

  // STEP 3: check plaq and linktrace
  Float plaq_inheader(0.0);
  Float linktrace_inheader(0.0);
  if(isRoot()) {
    plaq_inheader = hd.asFloat("PLAQUETTE");
    linktrace_inheader = hd.asFloat("LINK_TRACE");
  }

  if(lat.GaugeField() != lpoint) lat.GaugeField(lpoint);
  if(! CheckPlaqLinktrace(lat,rd_arg,plaq_inheader, linktrace_inheader))  return;

#ifdef PROFILE
  gettimeofday(&end,NULL);
  print_flops(cname,"read",0,&start,&end);
#endif

  load_good = true;

};

// just add a extra level of sum (the global sum)
bool ReadLatticeParallel::CheckPlaqLinktrace(Lattice &lat, const QioArg & rd_arg,
					     const Float plaq_inheader, const Float linktrace_inheader) 
{
  int error = 0;

  Float plaq = lat.SumReTrPlaq() / 18.0 / rd_arg.VolSites() ;
  Float devplaq(0.0);
  if(isRoot()) {
    devplaq =   fabs(  (plaq - plaq_inheader) / plaq ) ;
    printf("plaqerr::  calc: %f  header: %f dev.: %f\n",
	   (float)plaq, (float)plaq_inheader,(float)devplaq);
  }

  Float linktrace(0);
  int is;
  Matrix *m =  lat.GaugeField(); 
  for(is=0;is< rd_arg.VolNodeSites()*4; is++){
    linktrace += m->ReTr();
    m++;
  }
  if(rd_arg.Scoor() == 0) 
    linktrace = globalSumFloat(linktrace) / (rd_arg.VolSites()*12.0);
  else
    globalSumFloat(0.0);

  if(isRoot()) {
    Float devlinktrace =   
      fabs(  (linktrace - linktrace_inheader) / linktrace );

    printf("linktrace::  calc: %f  header: %f dev.: %f\n",
	   (float) linktrace, (float)linktrace_inheader,
	   (float)devlinktrace);
  
    Float chkprec = rd_arg.CheckPrecision;
    if(devplaq > chkprec || devlinktrace > chkprec) {
      cout << "Plaq trace and/or Link trace different from header" << endl;
      error = 1;
    }
  }

  if(synchronize(error) != 0) return false;

  return true;
}


CPS_END_NAMESPACE
