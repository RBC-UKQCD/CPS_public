
#include <config.h>

#include <stdlib.h>
#include <math.h>
#include <util/latrngio.h>
#include <util/qioarg.h>
#include <util/intconv.h>
#include <util/random.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

CPS_START_NAMESPACE
using namespace std;

///////////////////////////////////////////////////////////////
// LatRngIO functions, common function for Read & Write
void LatRngIO::calcDim(const QioArg & io_arg) {
  for(int i=0;i<5;i++) local_n[i] = io_arg.NodeSites(i) / 2;
  if(local_n[4]==0)  local_n[4] = 1;
  for(int i=0;i<5;i++) {
    global_n[i] = io_arg.Nodes(i) * local_n[i];
    coor[i] = io_arg.Coor(i);
  }
}

int LatRngIO::uniqueSiteID5d(int local_id_5d) {
  int loc[5];
  for(int i=0; i<5;i++) {
    loc[i] = local_id_5d % local_n[i];
    local_id_5d /= local_n[i];
    loc[i] += local_n[i] * coor[i];
  }
  
  int global_id = loc[4];
  for(int i=3;i>=0;i--) {
    global_id = global_id * global_n[i+1] + loc[i];
  }
  return global_id;
}

int LatRngIO::uniqueSiteID4d(int local_id_4d) {
  int loc[4];
  for(int i=0; i<4;i++) {
    loc[i] = local_id_4d % local_n[i];
    local_id_4d /= local_n[i];
    loc[i] += local_n[i] * coor[i];
  }
  
  int global_id = loc[3];
  for(int i=2;i>=0;i--) {
    global_id = global_id * global_n[i+1] + loc[i];
  }
  return global_id;
}



///////////////////////////////////////////////////////////////
// Lat Random Generator Read  functions 

LatRngRead::LatRngRead() 
  : LatRngIO()
{ }

LatRngRead::~LatRngRead() 
{ }

void LatRngRead::read(UGrandomGenerator * ugran, UGrandomGenerator * ugran_4d,
		      const QioArg & rd_arg) { 
  cout << endl << "Loading RNGs..." << endl << endl;

  io_good = false;
  int error = 0;

  // prepare for pos-dep. checksum calc.
  calcDim(rd_arg);

  // all open file and check error
  ifstream input(rd_arg.FileName);
  if ( !input.good() )
    {
      cout << "Could not open file:\n   "
	   << rd_arg.FileName
           << "\nfor input.\n";
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


  if(isRoot()) {
    if(hd.asString("DATATYPE") != "LATTICE_RNG_5D_4D"){  // need both 5d & 4d
      cout << "Invalid RNG type: "
	   << hd.asString("DATATYPE") << endl;
      error = 1;
    }
  }
  if(synchronize(error) != 0)  return;


  // check dimensions, etc
  int nx = rd_arg.Xnodes() * rd_arg.XnodeSites();
  int ny = rd_arg.Ynodes() * rd_arg.YnodeSites();
  int nz = rd_arg.Znodes() * rd_arg.ZnodeSites();
  int nt = rd_arg.Tnodes() * rd_arg.TnodeSites();
  int ns = rd_arg.Snodes() * rd_arg.SnodeSites();

  if(isRoot()) {
    int rdnx,rdny,rdnz,rdnt,rdns;

    //====================================
    // initialise the size of the lattice
    // int characters
    
    rdnx=hd.asInt("DIMENSION_1");
    rdny=hd.asInt("DIMENSION_2");
    rdnz=hd.asInt("DIMENSION_3");
    rdnt=hd.asInt("DIMENSION_4");
    rdns=hd.asInt("DIMENSION_5");
    if(!hd.found())  rdns = 1;  // s-dim not found, use default value

    cout << "Lattice dimensions: " << rdnx <<"x" << rdny <<"x" << rdnz <<"x" << rdnt <<"x" << rdns<< endl;
    if(rdnx != nx || rdny != ny || rdnz != nz || rdnt != nt || rdns!=ns) {
      cout << "Dimensions in file DISAGREE with GlobalJobParameter!"<<endl;
      error = 1;
    }
  }

  if(synchronize(error) != 0)  return;

  // check floating point format
  if(isRoot()) {
    intconv.setFileFormat(hd.asString("INT_FORMAT").c_str());
  }
  broadcastInt((int*)&intconv.fileFormat);
  if(intconv.fileFormat == INT_UNKNOWN) {
    cout << "Data file Integer format UNKNOWN" << endl;
    return;
  }


  if(isRoot())  hd.Show();

  // the UGrandomGenerator has multiple type members, need a universal
  // format to make it cross-platform
  int size_rng_ints = ugran[0].RNGints();
  cout << "size_rng_ints = " << size_rng_ints << endl;
  int size_rng_chars = size_rng_ints * intconv.fileIntSize();

  // start reading data

  int xbegin = rd_arg.XnodeSites()/2 * rd_arg.Xcoor(), xend = rd_arg.XnodeSites()/2 * (rd_arg.Xcoor()+1);
  int ybegin = rd_arg.YnodeSites()/2 * rd_arg.Ycoor(), yend = rd_arg.YnodeSites()/2 * (rd_arg.Ycoor()+1);
  int zbegin = rd_arg.ZnodeSites()/2 * rd_arg.Zcoor(), zend = rd_arg.ZnodeSites()/2 * (rd_arg.Zcoor()+1);
  int tbegin = rd_arg.TnodeSites()/2 * rd_arg.Tcoor(), tend = rd_arg.TnodeSites()/2 * (rd_arg.Tcoor()+1);
  int sbegin = rd_arg.SnodeSites()/2 * rd_arg.Scoor(), send = rd_arg.SnodeSites()/2 * (rd_arg.Scoor()+1);
  if(rd_arg.SnodeSites()==1) send=sbegin+1;

  nx/=2; ny/=2; nz/=2; nt/=2; ns/=2;
  if(ns==0) ns = 1;

  int sblk = nx*ny*nz*nt*size_rng_chars;
  int tblk = nx*ny*nz*size_rng_chars;
  int zblk = nx*ny*size_rng_chars;
  int yblk = nx*size_rng_chars;

  
  // Step 1: read RNG to a buffer (Note: not all data in LatRanGen are stored)
  // and calc. checksum and pos_dep_csum at the same time

  // TempBufAlloc is a mem allocator that prevents mem leak on "return"s
  // see latrngio.h
  TempBufAlloc  filebuf(size_rng_chars);

  unsigned int csum = 0, pos_dep_csum = 0;
  int rngid=0;
  TempBufAlloc rng(size_rng_ints * sizeof(int));
  Float RandSum=0.0, Rand2Sum=0.0;

  // read in parallel manner, node 0 will assign & dispatch IO time slots
  setConcurIONumber(rd_arg.ConcurIONumber);
  getIOTimeSlot();

  input.seekg(data_start,ios_base::beg);

  // read 5d
  int jump = sbegin * sblk;
  for(int sr=sbegin;sr<send;sr++) {
    jump += tbegin * tblk;
    for(int tr=tbegin;tr<tend;tr++) {
      jump += zbegin * zblk;
      for(int zr=zbegin;zr<zend;zr++) {
	jump += ybegin * yblk;
	for(int yr=ybegin;yr<yend;yr++) {
	  jump += xbegin * size_rng_chars;
	  input.seekg(jump,ios_base::cur);

	  for(int xr=xbegin;xr<xend;xr++) {
	    input.read(filebuf,size_rng_chars);
	    csum += intconv.checksum(filebuf,size_rng_ints);
	    pos_dep_csum += intconv.posDepCsum(filebuf,size_rng_ints) * (uniqueSiteID5d(rngid)+1);
	    // load
	    intconv.file2host((char*)rng,filebuf,size_rng_ints);
	    ugran[rngid].load((int*)rng);
	    // next rand
	    Float rn = ugran[rngid].Grand();
	    RandSum += rn;
	    Rand2Sum += rn*rn;
	    // recover
	    ugran[rngid].load((int*)rng);
	    rngid++;
	  }

	  jump = (nx-xend) * size_rng_chars;  // "jump" restart from 0 and count
	}
	jump += (ny-yend) * yblk;
      }
      jump += (nz-zend) * zblk;
    }
    jump += (nt-tend) * tblk;
  }
  jump += (ns-send) * sblk;
  input.seekg(jump,ios_base::cur); // jump to the start of 4d data

  // read 4d
  // note: unlike write, all nodes in S dim. should read
  rngid = 0;
  jump = tbegin * tblk;
  for(int tr=tbegin;tr<tend;tr++) {
    jump += zbegin * zblk;
    for(int zr=zbegin;zr<zend;zr++) {
      jump += ybegin * yblk;
      for(int yr=ybegin;yr<yend;yr++) {
	jump += xbegin * size_rng_chars;
	input.seekg(jump,ios_base::cur);
	
	for(int xr=xbegin;xr<xend;xr++) {
	  input.read(filebuf,size_rng_chars);
	  if(rd_arg.Scoor() == 0) {
	    csum += intconv.checksum(filebuf,size_rng_ints);
	    pos_dep_csum += intconv.posDepCsum(filebuf,size_rng_ints) * (uniqueSiteID4d(rngid)+1);
	  }
	  // load
	  intconv.file2host((char*)rng,filebuf,size_rng_ints);
	  ugran_4d[rngid].load((int*)rng);
	  if(rd_arg.Scoor() == 0) {
	    // next rand
	    Float rn = ugran_4d[rngid].Grand();
	    RandSum += rn;
	    Rand2Sum += rn*rn;
	    // recover
	    ugran_4d[rngid].load((int*)rng);
	  }
	  rngid++;
	}
	
	jump = (nx-xend) * size_rng_chars;  // "jump" restart from 0 and count
      }
      jump += (ny-yend) * yblk;
    }
    jump += (nz-zend) * zblk;
  }

  finishIOTimeSlot();
  //

  if ( !input.good() ) { cout << "Input stream error!" << endl; error = 1; }
  input.close();  
  if(synchronize(error) != 0)  return;

  // Step 2.1: verify  checksum
  csum = globalSumUint(csum);

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

  // Step 2.2: verify position-dep. Checksum
  pos_dep_csum = globalSumUint(pos_dep_csum);
  if(isRoot()) {
    unsigned int cshead;
    sscanf(hd.asString("POS_DEP_CSUM").c_str() , "%lx ", &cshead);
  
    if( cshead != pos_dep_csum ) {
      cout << "Position Dependent CheckSUM error !! Header:" 
	   << hd.asString("POS_DEP_CSUM") << " Host_calc:"
	   <<hex << pos_dep_csum << dec << "\n";
      error = 1;
    }
    else
      cout << "Position Dependent CheckSUM is ok\n";
  }

  if(synchronize(error) != 0) return;


  // STEP 3: Verify Rand Average and Variance
  int total_rngs_4d = rd_arg.VolSites() / 16;
  int total_rngs_5d = total_rngs_4d * ns; 
  Float RandAvg = globalSumFloat(RandSum) / (total_rngs_5d + total_rngs_4d);
  Float RandVar = globalSumFloat(Rand2Sum) / (total_rngs_5d + total_rngs_4d)
                  - RandAvg * RandAvg;
  Float AvgInHeader, VarInHeader;
  if(isRoot()) {
    AvgInHeader = hd.asFloat("AVERAGE");
    VarInHeader = hd.asFloat("VARIANCE");
  }
  broadcastFloat(&AvgInHeader);
  broadcastFloat(&VarInHeader);

  cout << "Average::  calc: " << RandAvg << "  header: " << AvgInHeader 
       << "  dev.: " << fabs((RandAvg-AvgInHeader)/RandAvg) <<  endl;
  cout << "Variance:: calc: " << RandVar << "  header: " << VarInHeader 
       << "  dev.: " << fabs((RandVar-VarInHeader)/RandVar) << endl;

  io_good = true;

}


///////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// LatRngWrite functions

LatRngWrite::LatRngWrite() 
  : LatRngIO()
{ }

LatRngWrite::~LatRngWrite()
{ }



void LatRngWrite::write(UGrandomGenerator * ugran, UGrandomGenerator * ugran_4d,
			const QioArg & wt_arg) {
  cout << "LatRngWrite::write() start..." << endl;

  io_good = false;
  int error = 0;

  // prepare for pos-dep. checksum calc.
  calcDim(wt_arg);

  // some cross-platform issued to be discussed
  intconv.setFileFormat(wt_arg.FileIntFormat);
  cout << "Set File INT format : " << intconv.name(intconv.fileFormat) << endl;

  // size information
  int nx = wt_arg.Xnodes() * wt_arg.XnodeSites();
  int ny = wt_arg.Ynodes() * wt_arg.YnodeSites();
  int nz = wt_arg.Znodes() * wt_arg.ZnodeSites();
  int nt = wt_arg.Tnodes() * wt_arg.TnodeSites();
  int ns = wt_arg.Snodes() * wt_arg.SnodeSites();

  cout << "Lattice size = " << nx <<"x"<<ny<<"x"<<nz<<"x"<<nt<<"x"<<ns<<endl;

  // num of integers per RNG (NOTE: not all data in RNG are saved)
  int size_rng_ints = ugran[0].RNGints();
  cout << "size_rng_ints = " << size_rng_ints << endl;
  int size_rng_chars = size_rng_ints * intconv.fileIntSize();

  // TempBufAlloc is an mem allocator that prevents memory leak on "return"s
  TempBufAlloc  filebuf(size_rng_chars);
  TempBufAlloc  rng(size_rng_ints * sizeof(int));

  // start writing file
  ofstream output(wt_arg.FileName);
  if(!output.good())
    {
      cout << "Could not open file:\n   "
	   << wt_arg.FileName
           << "\nfor output.\n";
      error = 1;
    }

  if(synchronize(error) > 0)  return;

  // write header
  int csum_pos, pdcsum_pos;
  int avg_pos, var_pos;
  int data_start;

  if(isRoot()) {
    output << "BEGIN_HEADER" << endl;
    output << "HDR_VERSION = 1.0" << endl;
    output << "DATATYPE = LATTICE_RNG_5D_4D" << endl;
    output << "STORAGE_FORMAT = 1.0" << endl;
    for(int i=1;i<=5;i++){
      output << "DIMENSION_" << i << " = " << wt_arg.Nodes(i-1)*wt_arg.NodeSites(i-1) << "\n" ;
    }
    output << "CHECKSUM = ";
    csum_pos = output.tellp();
    output << hex << setw(8) << 0 << dec << endl;

    output << "POS_DEP_CSUM = ";
    pdcsum_pos = output.tellp();
    output << hex << setw(8) << 0 << dec << endl;

    output << "INT_FORMAT = " << intconv.name(intconv.fileFormat) << endl;
    
    output << "AVERAGE = ";
    avg_pos = output.tellp();
    output << setw(20) << left << ' ' << endl;  // 20 whitespaces, hopefully

    output << "VARIANCE = ";
    var_pos = output.tellp();
    output<< setw(20) << left << ' ' << endl;

    output << "END_HEADER" << endl;

    data_start = output.tellp();
  }
  
  broadcastInt(&data_start); // from 0 to all


  // start writing data
  int xbegin = wt_arg.XnodeSites()/2 * wt_arg.Xcoor(), xend = wt_arg.XnodeSites()/2 * (wt_arg.Xcoor()+1);
  int ybegin = wt_arg.YnodeSites()/2 * wt_arg.Ycoor(), yend = wt_arg.YnodeSites()/2 * (wt_arg.Ycoor()+1);
  int zbegin = wt_arg.ZnodeSites()/2 * wt_arg.Zcoor(), zend = wt_arg.ZnodeSites()/2 * (wt_arg.Zcoor()+1);
  int tbegin = wt_arg.TnodeSites()/2 * wt_arg.Tcoor(), tend = wt_arg.TnodeSites()/2 * (wt_arg.Tcoor()+1);
  int sbegin = wt_arg.SnodeSites()/2 * wt_arg.Scoor(), send = wt_arg.SnodeSites()/2 * (wt_arg.Scoor()+1);
  if(wt_arg.SnodeSites()==1) send=sbegin+1;

  nx/=2; ny/=2; nz/=2; nt/=2; ns/=2;
  if(ns==0) ns = 1;

  int sblk = nx*ny*nz*nt*size_rng_chars;
  int tblk = nx*ny*nz*size_rng_chars;
  int zblk = nx*ny*size_rng_chars;
  int yblk = nx*size_rng_chars;

  // write in parallel manner, node 0 will assign & dispatch IO time slots
  unsigned int csum=0, pos_dep_csum = 0;
  Float RandSum=0.0, Rand2Sum = 0.0;
  int rngid=0;

  setConcurIONumber(wt_arg.ConcurIONumber);
  getIOTimeSlot();

  output.seekp(data_start,ios_base::beg);

  // write 5d
  int jump = sbegin * sblk;
  for(int sr=sbegin;sr<send;sr++) {
    jump += tbegin * tblk;
    for(int tr=tbegin;tr<tend;tr++) {
      jump += zbegin * zblk;
      for(int zr=zbegin;zr<zend;zr++) {
	jump += ybegin * yblk;
	for(int yr=ybegin;yr<yend;yr++) {
	  jump += xbegin * size_rng_chars;
	  output.seekp(jump,ios_base::cur);

	  for(int xr=xbegin;xr<xend;xr++) {
//	    printf("%d %d %d %d %d rngid=%d ",xr,yr,zr,tr,sr,rngid);
	    // dump
	    ugran[rngid].store((int*)rng);
	    intconv.host2file(filebuf,(char*)rng,size_rng_ints);

	    // checksum used for verification
	    csum += intconv.checksum(filebuf,size_rng_ints);
	    pos_dep_csum += intconv.posDepCsum(filebuf,size_rng_ints) * (uniqueSiteID5d(rngid)+1);

	    // next rand
	    Float rn = ugran[rngid].Grand();
	    RandSum += rn;
	    Rand2Sum += rn*rn;
	    // restore
	    ugran[rngid].load((int*)rng);
#if 0
	    Float rn2 = ugran[rngid].Grand();
	    printf("rn=%0.15e rn2=%0.15e\n",rn,rn2);
	    ugran[rngid].load((int*)rng);
#endif

	    output.write(filebuf, size_rng_chars);
	    
	    rngid++;
	  }

	  jump = (nx-xend) * size_rng_chars;  // "jump" restart from 0 and count
	}
	jump += (ny-yend) * yblk;
      }
      jump += (nz-zend) * zblk;
    }
    jump += (nt-tend) * tblk;
  }
  jump += (ns-send) * sblk;
  output.seekp(jump,ios_base::cur); // jump to the start of 4d data


  // write 4d
  rngid=0;

  if(wt_arg.Scoor() == 0) { 
    // only nodes on S==0 face will write ugran_4d
    // and nodes w/ S>0 should have identical ugran_4d as node w/ S=0
    jump = tbegin * tblk;
    for(int tr=tbegin;tr<tend;tr++) {
      jump += zbegin * zblk;
      for(int zr=zbegin;zr<zend;zr++) {
	jump += ybegin * yblk;
	for(int yr=ybegin;yr<yend;yr++) {
	  jump += xbegin * size_rng_chars;
	  output.seekp(jump,ios_base::cur);
	  
	  for(int xr=xbegin;xr<xend;xr++) {
	    // dump
	    ugran_4d[rngid].store((int*)rng);
	    intconv.host2file(filebuf,(char*)rng,size_rng_ints);

	    // checksum used for verification
	    csum += intconv.checksum(filebuf,size_rng_ints);
	    pos_dep_csum += intconv.posDepCsum(filebuf,size_rng_ints) * (uniqueSiteID5d(rngid)+1);

#if 1
	    // next rand
	    Float rn = ugran_4d[rngid].Grand();
	    RandSum += rn;
	    Rand2Sum += rn*rn;
	    
	    // restore
	    ugran_4d[rngid].load((int*)rng);
#endif

	    output.write(filebuf, size_rng_chars);

	    rngid++;
	  }
	  
	  jump = (nx-xend) * size_rng_chars;  // "jump" restart from 0 and count
	}
	jump += (ny-yend) * yblk;
      }
      jump += (nz-zend) * zblk;
    }
  }  // end if (S==0)
  
  finishIOTimeSlot();
  //
  
  csum = globalSumUint(csum);
  pos_dep_csum = globalSumUint(pos_dep_csum);

  // fill in verification information
  if(isRoot()) {
    output.seekp(csum_pos,ios::beg);
    output << hex << setw(8) << csum << dec;
  
    output.seekp(pdcsum_pos, ios::beg);
    output <<hex << setw(8) << pos_dep_csum << dec;
  }

  int total_rngs_4d = wt_arg.VolSites() / 16;
  int total_rngs_5d = total_rngs_4d * ns; 
  Float RandAvg = globalSumFloat(RandSum) / (total_rngs_5d + total_rngs_4d);
  Float RandVar = globalSumFloat(Rand2Sum) / (total_rngs_5d + total_rngs_4d)
                  - RandAvg * RandAvg;

  char numstr[100];

  if(isRoot()) {
    output.seekp(avg_pos, ios::beg);
    sprintf(numstr,"%-20.10lf",RandAvg);
    output << numstr;

    output.seekp(var_pos, ios::beg);
    sprintf(numstr,"%-20.10lf",RandVar);
    output << numstr;
  }

  if ( !output.good() ) { cout << "Output stream error!" << endl; error = 1; }
  output.close();  
  if(synchronize(error) != 0)  return;

  io_good = true;
}


CPS_END_NAMESPACE
 
