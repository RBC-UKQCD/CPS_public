#include <config.h>
#include <math.h>
#include <util/ReadLatticePar.h>

CPS_START_NAMESPACE
using namespace std;

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

  if (isRoot()) { // commander, analyze file header
    
    ifstream input(rd_arg.FileName);
    if ( !input.good() )
      {
	cout << "Could not open file:\n   "
	     << rd_arg.FileName
	     << "\nfor input.\nUSER: maybe you should kill the process!!\n";
	error = 1;
      }

    if(!error) {
      hd.read(input);
      input.close();
    }
  }
  if(synchronize(error) != 0) return;

  broadcastInt(&hd.data_start);
  broadcastInt(&hd.recon_row_3);
  cout << "recon_row_3 = " << hd.recon_row_3 << endl;


  // check all conditions between FILE and GJP
  int nx = rd_arg.Xnodes() * rd_arg.XnodeSites();
  int ny = rd_arg.Ynodes() * rd_arg.YnodeSites();
  int nz = rd_arg.Znodes() * rd_arg.ZnodeSites();
  int nt = rd_arg.Tnodes() * rd_arg.TnodeSites();

  if(isRoot()) {
    cout << "Data dimensions: " << hd.dimension[0] <<"x" << hd.dimension[1] <<"x" << hd.dimension[2] <<"x" << hd.dimension[3] << endl;
    if(hd.dimension[0] != nx || hd.dimension[1] != ny || hd.dimension[2] != nz || hd.dimension[3] != nt) {
      cout << "Dimensions in file DISAGREE with GlobalJobParameter!"<<endl;
      error = 1;
    }

    cout << "File specified:" << endl;
    cout << "X bc:" << (hd.boundary[0]==BND_CND_PRD ? "PERIODIC":"ANTI-PERIODIC") << endl;
    cout << "Y bc:" << (hd.boundary[1]==BND_CND_PRD ? "PERIODIC":"ANTI-PERIODIC") << endl;
    cout << "Z bc:" << (hd.boundary[2]==BND_CND_PRD ? "PERIODIC":"ANTI-PERIODIC") << endl;
    cout << "T bc:" << (hd.boundary[3]==BND_CND_PRD ? "PERIODIC":"ANTI-PERIODIC") << endl;
    
    if(hd.boundary[0] != rd_arg.Xbc() || hd.boundary[1] != rd_arg.Ybc() || hd.boundary[2] != rd_arg.Zbc() || hd.boundary[3] != rd_arg.Tbc()) {
      cout << "Boundary conditions in file DISAGREE with GlobalJobParameter!" << endl;
      error = 1;
    }
  }
  if(synchronize(error) != 0)  return;

  // see if file Floating Points is acceptable
  if(isRoot()) {
    fpconv.setFileFormat(hd.floating_point);
  }
  broadcastInt((int*)&fpconv.fileFormat);
  if(fpconv.fileFormat == FP_UNKNOWN) {
    cout << "Data file Floating Point format UNKNOWN" << endl;
    return;
  }

  cout << endl;
  if(isRoot())  hd.show();
  cout << endl;

  int data_per_site = hd.recon_row_3 ? 4*12 : 4*18;

  // read lattice data, using parallel style or serial (all on node 0) style
  unsigned int csum;


#if TARGET == QCDOC  // choice only applicable to QCDOC

  if(parIO()) {
    ParallelIO pario(rd_arg);
    if(! pario.load((char*)rd_arg.StartConfLoadAddr, data_per_site, sizeof(Matrix)*4,
		    hd, fpconv, 4, &csum))  return;  // failed to load
  }
  else {
    SerialIO serio(rd_arg);
    if(! serio.load((char*)rd_arg.StartConfLoadAddr, data_per_site, sizeof(Matrix)*4,
		    hd, fpconv, 4, &csum))  return;  // failed to load
  }

#else

  ParallelIO pario(rd_arg);
  if(! pario.load((char*)rd_arg.StartConfLoadAddr, data_per_site, sizeof(Matrix)*4,
		  hd, fpconv, 4, &csum))  return;  // failed to load
  
#endif

  //  cout << "loader finish, csum = " << hex << csum << dec << endl << endl;
  cout << "loader done" << endl << endl;


  // After reading...
  // STEP 1: checksum
  if(rd_arg.Scoor() == 0)
    csum = globalSumUint(csum);
  else
    globalSumUint(0);

  if(isRoot()) {
    if( hd.checksum != csum ) {
      cout << "CheckSUM error !! Header:" 
	   << hd.checksum << " Host calc:"
	   <<hex << csum << dec << "\n";
      error = 1;
    }
    else
      cout << "CheckSUM is ok\n";
  }

  if(synchronize(error) != 0) return;


  // STEP 2: reconstruct row 3
  int size_matrices( rd_arg.XnodeSites() * rd_arg.YnodeSites() 
		     * rd_arg.ZnodeSites() * rd_arg.TnodeSites() * 4); 

  Matrix * lpoint = rd_arg.StartConfLoadAddr;

  if(hd.recon_row_3) {
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
  if(lat.GaugeField() != lpoint) lat.GaugeField(lpoint);
  if(! CheckPlaqLinktrace(lat,rd_arg, hd.plaquette, hd.link_trace))  return;

#ifdef PROFILE
  gettimeofday(&end,NULL);
  print_flops(cname,"read",0,&start,&end);
#endif

  load_good = true;

};


bool ReadLatticeParallel::CheckPlaqLinktrace(Lattice &lat, const QioArg & rd_arg,
					     const Float plaq_inheader, const Float linktrace_inheader) 
{
  int error = 0;

  Float plaq = lat.SumReTrPlaq() / 18.0 / rd_arg.VolSites() ;
  Float devplaq(0.0);
  if(isRoot()) {
    devplaq =   fabs(  (plaq - plaq_inheader) / plaq ) ;
    cout << "plaquette::  calc: " << plaq << "  header: " << plaq_inheader
	 << "   rel.dev.: " << devplaq << endl;
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

    cout << "linktrace::  calc: " << linktrace << "  header: " << linktrace_inheader
	 << "   rel.dev.: " << devlinktrace << endl;
  
    Float chkprec = rd_arg.CheckPrecision;
    if(devplaq > chkprec || devlinktrace > chkprec) {
      cout << "Plaquette and/or Link trace different from header" << endl;
      error = 1;
    }
  }

  if(synchronize(error) != 0) return false;

  return true;
}

CPS_END_NAMESPACE
