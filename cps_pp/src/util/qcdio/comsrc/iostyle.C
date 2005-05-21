
// a slight modification on ReadLattice class (if not inheritance)
// to enable parallel reading/writing of "Gauge Connection Format" lattice data



// Read the format of Gauge Connection Format
// from QCDSP {load,unload}_lattice format

#include <config.h>

#include <math.h>
#include <util/gjp.h>
#include <util/iostyle.h>
#include <util/fpconv.h>
#include <util/time.h>

#include <iostream>
#include <fstream>
#include <string.h>

using namespace std;

CPS_START_NAMESPACE

#define PROFILE


/*********************************************************************/
/* ParallelIO functions ***********************************************/
/*********************************************************************/

// the last three pointers used to return information when loading Lattice Random Generators
int ParallelIO::load(char * data, const int data_per_site, const int site_mem,
		     const LatHeaderBase & hd, const DataConversion & dconv, 
		     const int dimension /* 4 or 5 */,
		     unsigned int * ptrcsum, unsigned int * ptrpdcsum,
		     Float * rand_sum, Float * rand_2_sum)  { 
  const char * fname = "load()";

  int error = 0;
  QioArg & rd_arg = qio_arg;

  // check dimensions, b.c, etc
  int nx = rd_arg.Xnodes() * rd_arg.XnodeSites();
  int ny = rd_arg.Ynodes() * rd_arg.YnodeSites();
  int nz = rd_arg.Znodes() * rd_arg.ZnodeSites();
  int nt = rd_arg.Tnodes() * rd_arg.TnodeSites();
  int ns = rd_arg.Snodes() * rd_arg.SnodeSites();

  int chars_per_site  = data_per_site * dconv.fileDataSize();

  int yblk = nx*chars_per_site;
  int zblk = ny * yblk;
  int tblk = nz * zblk;
  int sblk = nt * tblk;

  int xbegin = rd_arg.XnodeSites() * rd_arg.Xcoor(), xend = rd_arg.XnodeSites() * (rd_arg.Xcoor()+1);
  int ybegin = rd_arg.YnodeSites() * rd_arg.Ycoor(), yend = rd_arg.YnodeSites() * (rd_arg.Ycoor()+1);
  int zbegin = rd_arg.ZnodeSites() * rd_arg.Zcoor(), zend = rd_arg.ZnodeSites() * (rd_arg.Zcoor()+1);
  int tbegin = rd_arg.TnodeSites() * rd_arg.Tcoor(), tend = rd_arg.TnodeSites() * (rd_arg.Tcoor()+1);
  int sbegin = rd_arg.SnodeSites() * rd_arg.Scoor(), send = rd_arg.SnodeSites() * (rd_arg.Scoor()+1);


  // all open file and check error
  ifstream input(rd_arg.FileName);
  if ( !input.good() )
    {
      error = 1;
    }

  // executed by all, sync and share error status information
  if(synchronize(error) != 0)   
    ERR.FileR(cname, fname, rd_arg.FileName);

  // TempBufAlloc is a Mem Allocator that prevents mem leak on function exits
  TempBufAlloc fbuf(chars_per_site);  // buffer only stores one site
  
  // these two only needed when loading LatRng
  TempBufAlloc rng(data_per_site * dconv.hostDataSize());
  UGrandomGenerator * ugran = (UGrandomGenerator*)data;


  // read in parallel manner, node 0 will assign & dispatch IO time slots
  unsigned int csum = 0;
  unsigned int pdcsum = 0;
  Float RandSum = 0;
  Float Rand2Sum = 0;
  int siteid = 0;
  char * pd = data;

  VRB.Flow(cname, fname, "Parallel loading starting\n");
  setConcurIONumber(rd_arg.ConcurIONumber);
  //
  getIOTimeSlot();

  input.seekg(hd.dataStart(),ios_base::beg);

  int jump = 0;
  if(dimension == 5) jump = sbegin * sblk;

  for(int sr=sbegin; dimension==4 || sr<send; sr++) { // if 4-dim, has to enter once
    jump += tbegin * tblk;
    for(int tr=tbegin;tr<tend;tr++) {
      jump += zbegin * zblk;
      for(int zr=zbegin;zr<zend;zr++) {
	jump += ybegin * yblk;
	for(int yr=ybegin;yr<yend;yr++) {
	  jump += xbegin * chars_per_site;
	  input.seekg(jump,ios_base::cur);

	  for(int xr=xbegin;xr<xend;xr++) {
	    input.read(fbuf,chars_per_site);
	    csum += dconv.checksum(fbuf,data_per_site);
	    pdcsum += dconv.posDepCsum(fbuf, data_per_site, dimension,	rd_arg, siteid, 0);

	    if(hd.headerType() == LatHeaderBase::LATTICE_HEADER) {
	      for(int mat=0;mat<4;mat++) {
		dconv.file2host(pd, fbuf + chars_per_site/4*mat, data_per_site/4);
		pd += site_mem/4;
	      }
	    }
	    else { // LatHeaderBase::LATRNG_HEADER
	      // load
	      dconv.file2host(rng,fbuf,data_per_site);
	      ugran[siteid].load(rng.IntPtr());
	      // generate next rand for verification
	      Float rn = ugran[siteid].Grand();
	      RandSum += rn;
	      Rand2Sum += rn*rn;
	      // recover loading
	      ugran[siteid].load(rng.IntPtr());
	    }
	   
	    siteid++;
	  }
	  jump = (nx-xend) * chars_per_site;  // "jump" restart from 0 and count
	}
	jump += (ny-yend) * yblk;
      }
      jump += (nz-zend) * zblk;
      
      if(dimension == 4)
	cout << "Parallel loading: " << (int)((tr-tbegin+1) * 100.0 /(tend-tbegin)) << "% done." << endl;
    }
    
    if(dimension == 4) break;

    jump += (nt-tend) * tblk;

    cout << "Parallel loading: " << (int)((sr-sbegin+1) * 100.0 /(send-sbegin)) << "% done." << endl;
  }

  if ( !input.good() )    error = 1;

  VRB.Flow(cname,fname, "This Group Done!\n");

  finishIOTimeSlot();
  //
  
  input.close();


  if(synchronize(error) != 0)  
    ERR.General(cname, fname, "Loading failed\n");

  // This 3 lines differ from unloading part
  // these nodes participate in loading but not in summation
  if(dimension == 4 && rd_arg.Scoor()!=0) { 
    csum = pdcsum = 0;
    RandSum = Rand2Sum = 0;
  }

  if(ptrcsum) *ptrcsum = csum;
  if(ptrpdcsum) *ptrpdcsum = pdcsum;
  if(rand_sum) *rand_sum = RandSum;
  if(rand_2_sum) *rand_2_sum = Rand2Sum;

  return 1;
}

int ParallelIO::store(ostream & output,
		      char * data, const int data_per_site, const int site_mem,
		      LatHeaderBase & hd, const DataConversion & dconv,
		      const int dimension /* 4 or 5 */,
		      unsigned int * ptrcsum, unsigned int * ptrpdcsum,
		      Float * rand_sum, Float * rand_2_sum)  { 
  const char * fname = "store()";
  
  int error = 0;

  //  const int data_per_site = wt_arg.ReconRow3 ? 4*12 : 4*18;
  const int chars_per_site = data_per_site * dconv.fileDataSize();
  QioArg & wt_arg = qio_arg;

  int xbegin = wt_arg.XnodeSites() * wt_arg.Xcoor(), xend = wt_arg.XnodeSites() * (wt_arg.Xcoor()+1);
  int ybegin = wt_arg.YnodeSites() * wt_arg.Ycoor(), yend = wt_arg.YnodeSites() * (wt_arg.Ycoor()+1);
  int zbegin = wt_arg.ZnodeSites() * wt_arg.Zcoor(), zend = wt_arg.ZnodeSites() * (wt_arg.Zcoor()+1);
  int tbegin = wt_arg.TnodeSites() * wt_arg.Tcoor(), tend = wt_arg.TnodeSites() * (wt_arg.Tcoor()+1);
  int sbegin = wt_arg.SnodeSites() * wt_arg.Scoor(), send = wt_arg.SnodeSites() * (wt_arg.Scoor()+1);

  int nx = wt_arg.XnodeSites() * wt_arg.Xnodes();
  int ny = wt_arg.YnodeSites() * wt_arg.Ynodes();
  int nz = wt_arg.ZnodeSites() * wt_arg.Znodes();
  int nt = wt_arg.TnodeSites() * wt_arg.Tnodes();
  int ns = wt_arg.SnodeSites() * wt_arg.Snodes();


  int yblk = nx*chars_per_site;
  int zblk = yblk * ny;
  int tblk = zblk * nz;
  int sblk = tblk * nt;


  // TempBufAlloc is a mem allocator that prevents mem leak on function exits
  TempBufAlloc fbuf(chars_per_site);

  // these two only need when doing LatRng unloading
  TempBufAlloc rng(data_per_site * dconv.hostDataSize());
  UGrandomGenerator * ugran = (UGrandomGenerator*)data;

  // start parallel writing
  unsigned int csum = 0, pdcsum = 0;
  Float RandSum=0, Rand2Sum=0;
  const char * pd = data;
  int siteid=0;

  VRB.Flow(cname, fname, "Parallel unloading starting\n");
  setConcurIONumber(wt_arg.ConcurIONumber);
  getIOTimeSlot();

  if(dimension ==5 || wt_arg.Scoor() == 0) { // this line differs from read()
    output.seekp(hd.data_start,ios_base::beg);

    int jump=0;
    if(dimension==5) jump = sbegin * sblk;

    for(int sr=sbegin; dimension==4 || sr<send; sr++) {
      jump += tbegin * tblk;
      for(int tr=tbegin;tr<tend;tr++) {
	jump += zbegin * zblk;
	for(int zr=zbegin;zr<zend;zr++) {
	  jump += ybegin * yblk;
	  for(int yr=ybegin;yr<yend;yr++) {
	    jump += xbegin * chars_per_site;
	    output.seekp(jump, ios_base::cur);

	    for(int xr=xbegin;xr<xend;xr++) {
	      if(hd.headerType() == LatHeaderBase::LATTICE_HEADER) {
		for(int mat=0;mat<4;mat++) {
		  dconv.host2file(fbuf + chars_per_site/4*mat, pd, data_per_site/4);
		  pd += site_mem/4;
		}
	      }
	      else { //LatHeaderBase::LATRNG_HEADER
		// dump
		ugran[siteid].store(rng.IntPtr());
		dconv.host2file(fbuf,rng,data_per_site);
		// next rand
		Float rn = ugran[siteid].Grand();
		RandSum += rn;
		Rand2Sum += rn*rn;
		// recover
		ugran[siteid].load(rng.IntPtr());
#if 0
		Float rn2 = ugran[siteid].Grand();
		printf("rn=%0.15e rn2=%0.15e\n",rn,rn2);
		ugran[siteid].load(rng.IntPtr());
#endif
	      }

	      csum += dconv.checksum(fbuf,data_per_site);
	      pdcsum += dconv.posDepCsum(fbuf, data_per_site, dimension, wt_arg, siteid, 0);
	      output.write(fbuf,chars_per_site);

	      siteid++;
	    }	  
	  
	    jump = (nx-xend) * chars_per_site;  // jump restarted from 0
	  }
	  jump += (ny-yend) * yblk;
	}
	jump += (nz-zend) * zblk;

	if(dimension==4)
	  cout << "Parallel Unloading: " << (int)((tr-tbegin+1)*100.0/(tend-tbegin))
	       << "% done." << endl;
      }

      if(dimension == 4) break;

      jump += (nt-tend) * tblk;

      cout << "Parallel Unloading: " << (int)((sr-sbegin+1)*100.0/(send-sbegin))
	   << "% done." << endl;
    }
      
    if ( !output.good() )   error = 1;
  }

  VRB.Flow(cname, fname, "This Group done!\n");

  finishIOTimeSlot();
  //

  if(ptrcsum) *ptrcsum = csum;
  if(ptrpdcsum) *ptrpdcsum = pdcsum;
  if(rand_sum) *rand_sum = RandSum;
  if(rand_2_sum) *rand_2_sum = Rand2Sum;

  return 1;

}



#if TARGET == QCDOC

/*********************************************************************/
/* SerialIO functions ***********************************************/
/*********************************************************************/

int SerialIO::load(char * data, const int data_per_site, const int site_mem,
		   const LatHeaderBase & hd, const DataConversion & dconv, 
		   const int dimension /* 4 or 5 */,
		   unsigned int * ptrcsum, unsigned int * ptrpdcsum,
		   Float * rand_sum, Float * rand_2_sum) {
  const char * fname = "load()";

  // only node 0 is responsible for writing
  int error = 0;
  QioArg & rd_arg = qio_arg;

  // check dimensions, b.c, etc
  int nx = rd_arg.Xnodes() * rd_arg.XnodeSites();
  int ny = rd_arg.Ynodes() * rd_arg.YnodeSites();
  int nz = rd_arg.Znodes() * rd_arg.ZnodeSites();
  int nt = rd_arg.Tnodes() * rd_arg.TnodeSites();
  int ns = rd_arg.Snodes() * rd_arg.SnodeSites();

  //  int data_per_site = hd.recon_row_3? 4*12 : 4*18;
  int chars_per_site  = data_per_site * dconv.fileDataSize();

  ifstream input;
  if(isNode0()) {
    input.open(rd_arg.FileName);
    if ( !input.good() )
      {
	error = 1;
      }
  }

  // executed by all, sync and share error status information
  if(synchronize(error) != 0)   
    ERR.FileR(cname, fname, rd_arg.FileName);

  // TempBufAlloc is a Mem Allocator that prevents mem leak on function exits
  TempBufAlloc fbuf(chars_per_site);
  TempBufAlloc rng(data_per_site * dconv.hostDataSize());

  if(isNode0()) {
    input.seekg(hd.dataStart(),ios_base::beg);
    if(!input.good()) error = 1;
  }
  if(synchronize(error) != 0)   
    ERR.General(cname, fname, "Loading failed\n");
  
  int global_id = 0;
  unsigned int csum = 0;
  unsigned int pdcsum = 0;
  Float RandSum = 0;
  Float Rand2Sum = 0;
  UGrandomGenerator * ugran = (UGrandomGenerator*)data;

  VRB.Flow(cname, fname, "Serial loading <thru node 0> starting\n");
  for(int sc=0; dimension==4 || sc<ns; sc++) {
    for(int tc=0;tc<nt;tc++) {
      for(int zc=0;zc<nz;zc++) {
	for(int yc=0;yc<ny;yc++) {
	  for(int xnd=0; xnd<rd_arg.Xnodes(); xnd++) {
	    if(isNode0()) { // only node 0 reads
	      char * pd = data;

	      // only read file to lowest end of buffer,
	      // then shift (counter-shift) to final place
	      for(int xst=0; xst<rd_arg.XnodeSites(); xst++) {
		input.read(fbuf,chars_per_site);
		csum += dconv.checksum(fbuf,data_per_site);
		pdcsum += dconv.posDepCsum(fbuf, data_per_site, dimension, rd_arg,
					   -1, global_id);

		if(hd.headerType() == LatHeaderBase::LATTICE_HEADER) {
		  for(int mat=0;mat<4;mat++) {
		    dconv.file2host(pd, fbuf + chars_per_site/4*mat, data_per_site/4);
		    pd += site_mem/4;
		  }
		}
		else { //LatHeaderBase::LATRNG_HEADER
		  // load
		  dconv.file2host(rng,fbuf,data_per_site);
		  ugran[xst].load(rng.IntPtr());
		  // generate next rand for verification
		  Float rn = ugran[xst].Grand();
		  RandSum += rn;
		  Rand2Sum += rn*rn;
		  // recover loading
		  ugran[xst].load(rng.IntPtr());
		}
		
		global_id ++;

		if(!input.good()) {
		  error = 1;
		  break;
		}
	      }
	    } // endif(isNode0())
	  
	    // execute by all
	    if(synchronize(error) != 0)  
	      ERR.General(cname, fname, "Loading failed\n");

	    xShiftNode(data, site_mem * rd_arg.XnodeSites());
	  }

	  yShift(data, site_mem * rd_arg.XnodeSites());
	}

	zShift(data, site_mem * rd_arg.XnodeSites());
      }

      tShift(data, site_mem * rd_arg.XnodeSites());

      if(dimension==4)
	cout << "Serial loading: " << (int)((tc+1) * 100.0 / nt) << "% finished." << endl;
    }

    if(dimension==4) break;

    sShift(data, site_mem * rd_arg.XnodeSites());
    cout << "Serial loading: " << (int)((sc+1) * 100.0 / ns) << "% finished." << endl;
  }
  
  if(synchronize(error)!=0) 
    ERR.General(cname, fname, "Loading failed\n");

  if(isNode0()) 
    input.close();
  
  // spread (clone) lattice data along s-dim
  if(dimension==4) 
    sSpread(data, site_mem * rd_arg.VolNodeSites());
  
  VRB.Flow(cname, fname, "Loading done!\n");

  if(ptrcsum) *ptrcsum = csum;
  if(ptrpdcsum) *ptrpdcsum = pdcsum;
  if(rand_sum) *rand_sum = RandSum;
  if(rand_2_sum) *rand_2_sum = Rand2Sum;

  return 1;
}

int SerialIO::store(ostream & output,
		    char * data, const int data_per_site, const int site_mem,
		    LatHeaderBase & hd, const DataConversion & dconv,
		    const int dimension /* 4 or 5 */,
		    unsigned int * ptrcsum, unsigned int * ptrpdcsum,
		    Float * rand_sum, Float * rand_2_sum)  { 
  const char * fname = "store()";
  
  int error = 0;

  QioArg & wt_arg = qio_arg;

  //  const int data_per_site = wt_arg.ReconRow3 ? 4*12 : 4*18;
  const int chars_per_site = data_per_site * dconv.fileDataSize();

  int nx = wt_arg.XnodeSites() * wt_arg.Xnodes();
  int ny = wt_arg.YnodeSites() * wt_arg.Ynodes();
  int nz = wt_arg.ZnodeSites() * wt_arg.Znodes();
  int nt = wt_arg.TnodeSites() * wt_arg.Tnodes();
  int ns = wt_arg.SnodeSites() * wt_arg.Snodes();


  TempBufAlloc fbuf(chars_per_site);
  TempBufAlloc rng(data_per_site * dconv.hostDataSize());

  // start serial writing
  unsigned int csum = 0, pdcsum = 0;
  Float RandSum = 0, Rand2Sum = 0;
  UGrandomGenerator * ugran = (UGrandomGenerator*)data;
  int global_id = 0;

  output.seekp(hd.dataStart(), ios_base::beg);

  VRB.Flow(cname, fname, "Serial unloading <thru node 0> starting\n");
  for(int sc=0; dimension==4 || sc<ns; sc++) {
    for(int tc=0;tc<nt;tc++) {
      for(int zc=0;zc<nz;zc++) {
	for(int yc=0;yc<ny;yc++) {
	  for(int xnd=0; xnd<wt_arg.Xnodes(); xnd++) {
	    if(isNode0()) {
	      const char * pd = data;

	      for(int xst=0;xst<wt_arg.XnodeSites();xst++) {
		if(hd.headerType() == LatHeaderBase::LATTICE_HEADER) {
		  for(int mat=0;mat<4;mat++) {
		    dconv.host2file(fbuf + chars_per_site/4*mat, pd, data_per_site/4);
		    pd += site_mem/4;
		  }
		}
		else { //LatHeaderBase::LATRNG_HEADER
		  // dump
		  ugran[xst].store(rng.IntPtr());
		  dconv.host2file(fbuf,rng,data_per_site);
		  // next rand
		  Float rn = ugran[xst].Grand();
		  RandSum += rn;
		  Rand2Sum += rn*rn;
		  // recover
		  ugran[xst].load(rng.IntPtr());
#if 0
		  Float rn2 = ugran[xst].Grand();
		  printf("rn=%0.15e rn2=%0.15e\n",rn,rn2);
		  ugran[xst].load(rng.IntPtr());
#endif
		}

		csum += dconv.checksum(fbuf,data_per_site);
		pdcsum += dconv.posDepCsum(fbuf, data_per_site, dimension, wt_arg, 
					   -1, global_id);
		output.write(fbuf,chars_per_site);
		if(!output.good()) {
		  error = 1;
		  break;
		}

		global_id ++;
	      }
	    }

	    // run by all nodes
	    if(synchronize(error) != 0) 
	      ERR.General(cname, fname, "Unloading failed\n");
	    
	    xShiftNode(data, site_mem * wt_arg.XnodeSites());
	  }
	  
	  yShift(data, site_mem * wt_arg.XnodeSites());
	}
	
	zShift(data, site_mem * wt_arg.XnodeSites());
      }
      
      tShift(data, site_mem * wt_arg.XnodeSites());

      if(dimension==4)
	cout << "Serial unloading: " << (int)((tc+1)*100.0/nt) << "% done." <<endl;
    }
    
    if(dimension==4) break;

    sShift(data, site_mem * wt_arg.XnodeSites());
    cout << "Serial unloading: " << (int)((sc+1)*100.0/ns) << "% done." <<endl;
  }

  VRB.Flow(cname, fname, "Unloading done!\n");

  if(ptrcsum) *ptrcsum = csum;
  if(ptrpdcsum) *ptrpdcsum = pdcsum;
  if(rand_sum) *rand_sum = RandSum;
  if(rand_2_sum) *rand_2_sum = Rand2Sum;


  return 1;
}




// NOTE: !!!
// the x-shift is a little different from y-, z-, & t-shift.
// x-shift shift 1 NODE, while others shift a SITE
void SerialIO::xShiftNode(char * data, const int xblk, const int dir) const {

  if(qio_arg.Xnodes() <= 1) return;

  if(isRow0 ()) {  // x rotation only apply to nodes (0,0,0,x)
    //    const SCUDir pos_dir[] = { SCU_XP, SCU_YP, SCU_ZP, SCU_TP };
    //    const SCUDir neg_dir[] = { SCU_XM, SCU_YM, SCU_ZM, SCU_TM };
    
    TempBufAlloc  recvbuf (xblk);
    char * sendbuf = data;
    
    SCUDirArg send, recv;
    if(dir>0) {
      send.Init(sendbuf, SCU_XP, SCU_SEND, xblk);
      recv.Init(recvbuf, SCU_XM, SCU_REC,  xblk);
    }
    else {
      send.Init(sendbuf, SCU_XM, SCU_SEND, xblk);
      recv.Init(recvbuf, SCU_XP, SCU_REC,  xblk);
    }
    
    SCUTrans(&send);
    SCUTrans(&recv);
    SCUTransComplete();
    
    memcpy(sendbuf, recvbuf, xblk);
  }
}

void SerialIO::yShift(char * data, const int xblk, const int dir) const {
  int useSCU = 1;
  if(qio_arg.Ynodes() <= 1) useSCU = 0;

  if(isFace0()) {
    TempBufAlloc sendbuf(xblk);
    TempBufAlloc recvbuf(xblk);

    SCUDirArg send, recv;

    if(dir>0) {
      if(useSCU) {
	memcpy(sendbuf, data + xblk * (qio_arg.YnodeSites() - 1), xblk);
	send.Init(sendbuf, SCU_YP, SCU_SEND, xblk);
	recv.Init(recvbuf, SCU_YM, SCU_REC,  xblk);
      }	
      else {
	memcpy(recvbuf, data + xblk * (qio_arg.YnodeSites() - 1), xblk);
      }
    }
    else {
      if(useSCU) {
	memcpy(sendbuf, data, xblk);
	send.Init(sendbuf, SCU_YM, SCU_SEND, xblk);
	recv.Init(recvbuf, SCU_YP, SCU_REC,  xblk);
      }
      else {
	memcpy(recvbuf, data, xblk);
      }
    }

    if(useSCU) {
      SCUTrans(&send);
      SCUTrans(&recv);
    }

    // doing memory move at the same time
    if(dir>0) {
      for(int i=qio_arg.YnodeSites()-1;i>0;i--) 
	memcpy(data + i*xblk, data + (i-1)*xblk, xblk);
      //      memmove(data + xblk, data, xblk * (qio_arg.YnodeSites()-1));
    }
    else{
      for(int i=0;i<qio_arg.YnodeSites()-1;i++)
	memcpy(data+i*xblk, data+(i+1)*xblk, xblk);
      //      memmove(data, data+xblk, xblk * (qio_arg.YnodeSites()-1));
    }

    if(useSCU) {
      SCUTransComplete();
    }

    if(dir>0) {
      memcpy(data, recvbuf, xblk);
    }
    else{
      memcpy(data + xblk*(qio_arg.YnodeSites()-1), recvbuf, xblk);
    }
  }
}


void SerialIO::zShift(char * data, const int xblk, const int dir) const {
  int useSCU = 1;
  if(qio_arg.Znodes() <= 1) useSCU = 0;

  if(isCube0()) {
    int yblk = xblk * qio_arg.YnodeSites();

    TempBufAlloc sendbuf(xblk);
    TempBufAlloc recvbuf(xblk);

    SCUDirArg send, recv;

    char * loface = data;
    char * hiface = data + yblk * (qio_arg.ZnodeSites()-1);

    for(int yc=0;yc<qio_arg.YnodeSites();yc++) {
      if(dir>0) {
	if(useSCU) {
	  memcpy(sendbuf, hiface + yc*xblk, xblk);
	  send.Init(sendbuf, SCU_ZP, SCU_SEND, xblk);
	  recv.Init(recvbuf, SCU_ZM, SCU_REC,  xblk);
	}
	else {
	  memcpy(recvbuf, hiface + yc*xblk, xblk);
	}
      }
      else {
	if(useSCU) {
	  memcpy(sendbuf, loface + yc*xblk, xblk);
	  send.Init(sendbuf, SCU_ZM, SCU_SEND, xblk);
	  recv.Init(recvbuf, SCU_ZP, SCU_REC,  xblk);
	}
	else {
	  memcpy(recvbuf, loface + yc*xblk, xblk);
	}
      }

      if(useSCU) {
	SCUTrans(&send);
	SCUTrans(&recv);
      }

      if(dir > 0) {
	for(int zc=qio_arg.ZnodeSites()-1;zc>0;zc--) {
	  memcpy(data + zc*yblk + yc*xblk, data + (zc-1)*yblk + yc*xblk, xblk);
	}
      }
      else {
	for(int zc=0;zc<qio_arg.ZnodeSites()-1;zc++) {
	  memcpy(data + zc*yblk + yc*xblk, data + (zc+1)*yblk + yc*xblk, xblk);
	}
      }
	  
      if(useSCU) {
	SCUTransComplete();
      }

      if(dir>0) {
	memcpy(loface + yc*xblk, recvbuf, xblk);
      }
      else{
	memcpy(hiface + yc*xblk, recvbuf, xblk);
      }
    }
  }
}

void SerialIO::tShift(char * data, const int xblk, const int dir) const {
  int useSCU = 1;
  if(qio_arg.Tnodes() <= 1) useSCU = 0;

  if(isSdim0()) {
    int yblk = xblk * qio_arg.YnodeSites();
    int zblk = yblk * qio_arg.ZnodeSites();

    TempBufAlloc sendbuf(xblk);
    TempBufAlloc recvbuf(xblk);

    SCUDirArg send, recv;

    char * locube = data;
    char * hicube = data + zblk * (qio_arg.TnodeSites()-1);

    for(int zc=0;zc<qio_arg.ZnodeSites();zc++) {
      char * loface = locube + zc*yblk;
      char * hiface = hicube + zc*yblk;
      for(int yc=0;yc<qio_arg.YnodeSites();yc++) {
	if(dir>0) {
	  if(useSCU) {
	    memcpy(sendbuf, hiface + yc*xblk, xblk);
	    send.Init(sendbuf, SCU_TP, SCU_SEND, xblk);
	    recv.Init(recvbuf, SCU_TM, SCU_REC,  xblk);
	  }
	  else {
	    memcpy(recvbuf, hiface + yc*xblk, xblk);
	  }
	}
	else {
	  if(useSCU) {
	    memcpy(sendbuf, loface + yc*xblk, xblk);
	    send.Init(sendbuf, SCU_TM, SCU_SEND, xblk);
	    recv.Init(recvbuf, SCU_TP, SCU_REC,  xblk);
	  }
	  else {
	    memcpy(recvbuf, loface + yc*xblk, xblk);
	  }
	}

	if(useSCU) {
	  SCUTrans(&send);
	  SCUTrans(&recv);
	}

	if(dir > 0) {
	  for(int tc=qio_arg.TnodeSites()-1;tc>0;tc--) {
	    memcpy(data + tc*zblk + zc*yblk + yc*xblk, data + (tc-1)*zblk + zc*yblk + yc*xblk, xblk);
	  }
	}
	else {
	  for(int tc=0;tc<qio_arg.TnodeSites()-1;tc++) {
	    memcpy(data + tc*zblk + zc*yblk + yc*xblk, data + (tc+1)*zblk + zc*yblk + yc*xblk, xblk);
	  }
	}
	 
	if(useSCU) {
	  SCUTransComplete();
	}

	if(dir>0) {
	  memcpy(loface + yc*xblk, recvbuf, xblk);
	}
	else{
	  memcpy(hiface + yc*xblk, recvbuf, xblk);
	}
      }
    }
  }
}

void SerialIO::sShift(char * data, const int xblk, const int dir) const {
  if(qio_arg.Snodes() * qio_arg.SnodeSites() == 1) return;

  int useSCU = 1;
  if(qio_arg.Snodes() <= 1) useSCU = 0;

  int yblk = xblk * qio_arg.YnodeSites();
  int zblk = yblk * qio_arg.ZnodeSites();
  int tblk = zblk * qio_arg.TnodeSites();
  
  TempBufAlloc sendbuf(xblk);
  TempBufAlloc recvbuf(xblk);
  
  SCUDirArg send, recv;
  
  char * lohypcb = data;
  char * hihypcb = data + tblk * (qio_arg.SnodeSites()-1);
  
  for(int tc=0;tc<qio_arg.TnodeSites();tc++) {
    char * locube = lohypcb + tc*zblk;
    char * hicube = hihypcb + tc*zblk;
    
    for(int zc=0;zc<qio_arg.ZnodeSites();zc++) {
      char * loface = locube + zc*yblk;
      char * hiface = hicube + zc*yblk;
      for(int yc=0;yc<qio_arg.YnodeSites();yc++) {
	if(dir>0) {
	  if(useSCU) {
	    memcpy(sendbuf, hiface + yc*xblk, xblk);
	    send.Init(sendbuf, SCU_SP, SCU_SEND, xblk);
	    recv.Init(recvbuf, SCU_SM, SCU_REC,  xblk);
	  }
	  else {
	    memcpy(recvbuf, hiface + yc*xblk, xblk);
	  }
	}
	else {
	  if(useSCU) {
	    memcpy(sendbuf, loface + yc*xblk, xblk);
	    send.Init(sendbuf, SCU_SM, SCU_SEND, xblk);
	    recv.Init(recvbuf, SCU_SP, SCU_REC,  xblk);
	  }
	  else {
	    memcpy(recvbuf, loface + yc*xblk, xblk);
	  }	    
	}

	if(useSCU) {
	  SCUTrans(&send);
	  SCUTrans(&recv);
	}
	
	if(dir > 0) {
	  for(int sc=qio_arg.SnodeSites()-1;sc>0;sc--) {
	    memcpy(data + sc*tblk     + tc*zblk + zc*yblk + yc*xblk, 
		   data + (sc-1)*tblk + tc*zblk + zc*yblk + yc*xblk,
		   xblk);
	  }
	}
	else {
	  for(int sc=0;sc<qio_arg.SnodeSites()-1;sc++) {
	    memcpy(data + sc*tblk     + tc*zblk + zc*yblk + yc*xblk, 
		   data + (sc+1)*tblk + tc*zblk + zc*yblk + yc*xblk, 
		   xblk);
	  }
	}
	
	if(useSCU) {
	  SCUTransComplete();
	}
		
	if(dir>0) {
	  memcpy(loface + yc*xblk, recvbuf, xblk);
	}
	else{
	  memcpy(hiface + yc*xblk, recvbuf, xblk);
	}
      }
    }
  }
}


void SerialIO::sSpread(char * data, const int datablk) const {  

  const char *fname = "sSpread()";

  // clone the lattice from s==0 nodes to all s>0 nodes
  // only used in loading process
  if(qio_arg.Snodes() <= 1)  return;

  // eg. 8 nodes
  // step:  i   ii       iii       iv
  //        0   -->  1   -->   2   -->   3 
  //        |     
  //        \/    
  //        7   -->  6   -->   5   -->   4

  SCUDirArg  socket;

  VRB.Flow(cname, fname, "Spread on S dimension:: 0 ==> %d\n", qio_arg.Snodes()-1);

  if(qio_arg.Scoor() == 0) {
    socket.Init(data, SCU_SM, SCU_SEND, datablk);
    SCUTrans(&socket);
    SCUTransComplete();
  }
  else if(qio_arg.Scoor() == qio_arg.Snodes()-1) {
    socket.Init(data, SCU_SP, SCU_REC, datablk);
    SCUTrans(&socket);
    SCUTransComplete();
  }

  // spread simultaneously along +S and -S directions
  int sender[2] = { 0, qio_arg.Snodes()-1};
  int receiver[2] = { 1, sender[1]-1 };

  while(receiver[0] < receiver[1]) {  // send until two directions converge
    synchronize();

    VRB.Flow(cname, fname, "Spread on S dimension:: %d ==> %d\n", sender[0], receiver[0]);
    VRB.Flow(cname, fname, "Spread on S dimension:: %d ==> %d\n", sender[1], receiver[1]);

    if(qio_arg.Scoor() == sender[0]) {
      socket.Init(data, SCU_SP, SCU_SEND, datablk);
      SCUTrans(&socket);
      SCUTransComplete();
    }
    else if(qio_arg.Scoor() == receiver[0]) {
      socket.Init(data, SCU_SM, SCU_REC, datablk);
      SCUTrans(&socket);
      SCUTransComplete();
    }
    else if(qio_arg.Scoor() == sender[1]) {
      socket.Init(data, SCU_SM, SCU_SEND, datablk);
      SCUTrans(&socket);
      SCUTransComplete();
    }
    else if(qio_arg.Scoor() == receiver[1]) {
      socket.Init(data, SCU_SP, SCU_REC, datablk);
      SCUTrans(&socket);
      SCUTransComplete();
    }

    sender[0] = receiver[0];
    sender[1] = receiver[1];
    
    receiver[0]++;
    receiver[1]--;

  }

}



// Testing functions for *Shift() class

int SerialIO::backForthTest() {
  int error = 0;

  srand(1234);

  int datablk = 4*18*sizeof(Float) * qio_arg.VolNodeSites();
  int xblk = 4*18*sizeof(Float) * qio_arg.XnodeSites();
  char* data = new char[datablk];
  for(int i=0;i<datablk;i++)
    data[i] = rand() % 256 * uniqueID();

  TempBufAlloc buf(xblk);
  memcpy(buf, data, xblk);

  xShiftNode(data,xblk,-1);
  xShiftNode(data,xblk,1);
  if(memcmp(buf, data, xblk)) {
    error = 1;
    cout << "xShiftNode() error!!" << endl << endl;
  }

  xShiftNode(data,xblk,1);
  xShiftNode(data,xblk,-1);
  if(memcmp(buf, data, xblk)) {
    error = 1;
    cout << "xShiftNode() error!!" << endl << endl;
  }


  yShift(data,xblk,-1);
  yShift(data,xblk,1);
  if(memcmp(buf, data, xblk)) {
    error = 1;
    cout << "yShift() error!!" << endl << endl;
  }

  yShift(data,xblk,1);
  yShift(data,xblk,-1);
  if(memcmp(buf, data, xblk)) {
    error = 1;
    cout << "yShift() error!!" << endl << endl;
  }


  zShift(data,xblk,-1);
  zShift(data,xblk,1);
  if(memcmp(buf, data, xblk)) {
    error = 1;
    cout << "zShift() error!!" << endl << endl;
  }

  zShift(data,xblk,1);
  zShift(data,xblk,-1);
  if(memcmp(buf, data, xblk)) {
    error = 1;
    cout << "zShift() error!!" << endl << endl;
  }


  tShift(data,xblk,-1);
  tShift(data,xblk,1);
  if(memcmp(buf, data, xblk)) {
    error = 1;
    cout << "tShift() error!!" << endl << endl;
  }

  tShift(data,xblk,1);
  tShift(data,xblk,-1);
  if(memcmp(buf, data, xblk)) {
    error = 1;
    cout << "tShift() error!!" << endl << endl;
  }

  sShift(data,xblk,-1);
  sShift(data,xblk,1);
  if(memcmp(buf, data, xblk)) {
    error = 1;
    cout << "sShift() error!!" << endl << endl;
  }

  sShift(data,xblk,1);
  sShift(data,xblk,-1);
  if(memcmp(buf, data, xblk)) {
    error = 1;
    cout << "sShift() error!!" << endl << endl;
  }

  delete[] data;

  if(error) return 0;
  return 1;
}


int SerialIO::rotateTest() {
  int error = 0;

  srand(3456);

  int datablk = 4*18*sizeof(Float) * qio_arg.VolNodeSites();
  int xblk = 4*18*sizeof(Float) * qio_arg.XnodeSites();

  char * data = new char[datablk];
  for(int i=0;i<datablk;i++) 
    data[i] = rand() % 256 * uniqueID();

  TempBufAlloc buf(xblk);
  memcpy(buf, data, xblk);

  for(int i=0;i<qio_arg.Xnodes();i++) xShiftNode(data,xblk,-1);
  if(memcmp(buf, data, xblk)) { 
    error = 1;
    cout << "xShiftNode() error!!" << endl << endl;
  }
  for(int i=0;i<qio_arg.Xnodes();i++) xShiftNode(data,xblk,1);
  if(memcmp(buf, data, xblk)) { 
    error = 1;
    cout << "xShiftNode() error!!" << endl << endl;
  }

  for(int i=0;i<qio_arg.Ynodes()*qio_arg.YnodeSites();i++) yShift(data,xblk,-1);
  if(memcmp(buf, data, xblk)) { 
    error = 1;
    cout << "yShift() error!!" << endl << endl;
  }
  for(int i=0;i<qio_arg.Ynodes()*qio_arg.YnodeSites();i++) yShift(data,xblk,1);
  if(memcmp(buf, data, xblk)) { 
    error = 1;
    cout << "yShift() error!!" << endl << endl;
  }

  for(int i=0;i<qio_arg.Znodes()*qio_arg.ZnodeSites();i++) zShift(data,xblk,-1);
  if(memcmp(buf, data, xblk)) { 
    error = 1;
    cout << "zShift() error!!" << endl << endl;
  }
  for(int i=0;i<qio_arg.Znodes()*qio_arg.ZnodeSites();i++) zShift(data,xblk,1);
  if(memcmp(buf, data, xblk)) { 
    error = 1;
    cout << "zShift() error!!" << endl << endl;
  }

  for(int i=0;i<qio_arg.Tnodes()*qio_arg.TnodeSites();i++) tShift(data,xblk,-1);
  if(memcmp(buf, data, xblk)) { 
    error = 1;
    cout << "tShift() error!!" << endl << endl;
  }
  for(int i=0;i<qio_arg.Tnodes()*qio_arg.TnodeSites();i++) tShift(data,xblk,1);
  if(memcmp(buf, data, xblk)) { 
    error = 1;
    cout << "tShift() error!!" << endl << endl;
  }

  for(int i=0;i<qio_arg.Snodes()*qio_arg.SnodeSites();i++) sShift(data,xblk,-1);
  if(memcmp(buf, data, xblk)) { 
    error = 1;
    cout << "sShift() error!!" << endl << endl;
  }
  for(int i=0;i<qio_arg.Snodes()*qio_arg.SnodeSites();i++) sShift(data,xblk,1);
  if(memcmp(buf, data, xblk)) { 
    error = 1;
    cout << "sShift() error!!" << endl << endl;
  }

  delete[] data;

  if(error) return 0;
  return 1;

}


#endif // TARGET == QCDOC



CPS_END_NAMESPACE
