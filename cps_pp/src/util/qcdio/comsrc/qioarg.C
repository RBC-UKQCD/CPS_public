#include <iostream>
#include <unistd.h>
#include <util/gjp.h>
#include <sys/time.h>
#include <iomanip>
#include <cstring>
using namespace std;

#include <util/qioarg.h>

#if TARGET == QCDOC
#include <util/gsum64ext.h>
#endif

CPS_START_NAMESPACE


/////////////////////////////////////////////////////////////////
// QioArg members////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void QioArg::init(const char * file, const int concur_io_number, const Float chk_prec,
		  const FP_FORMAT file_format, const INT_FORMAT file_int_format,
		  const int recon_row_3) {
  for(int dir=0;dir<5;dir++) {
    nodes[dir] = GJP.Nodes(dir);
    node_sites[dir] = GJP.NodeSites(dir);
    coor[dir] = GJP.NodeCoor(dir);
  }

  for(int dir=0;dir<4;dir++) 
    bc[dir] = GJP.Bc(dir);

  StartConfLoadAddr = GJP.StartConfLoadAddr();

  // user set params
  ConcurIONumber = concur_io_number;
  strcpy(FileName, file);
  CheckPrecision = chk_prec;
  FileFpFormat = file_format;
  FileIntFormat = file_int_format;
  ReconRow3 = recon_row_3;
}


/////////////////////////////////////////////////////////////////////
// QioControl members////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
QioControl::QioControl() 
  : num_concur_io(8)
{
  cout << "I am on a " << GJP.Xnodes() << "x"<< GJP.Ynodes() << "x"<< GJP.Znodes() 
       << "x"<< GJP.Tnodes() <<"x"<<GJP.Snodes() << " machine" << endl;
  cout << "My pos is (" << GJP.XnodeCoor() << ","<< GJP.YnodeCoor() << ","<< GJP.ZnodeCoor() << ","
       << GJP.TnodeCoor() << "," << GJP.SnodeCoor() << ")" << endl;

  unique_id = GJP.XnodeCoor() + GJP.Xnodes() * (GJP.YnodeCoor() + GJP.Ynodes() * (GJP.ZnodeCoor() + GJP.Znodes() * (GJP.TnodeCoor() + GJP.Tnodes() * GJP.SnodeCoor() ) ) );
  cout << "My UniqueID() = " << unique_id << endl;

  number_nodes = GJP.Xnodes() * GJP.Ynodes() * GJP.Znodes() * GJP.Tnodes() * GJP.Snodes();
  cout << "Total number of nodes = " << number_nodes << endl;

}

QioControl::~QioControl() {

}


// The following are NEW functions added to QioControl class 
// to enable message passing between parallel processors, based on QMP calls
// (some are pretty useful...)

int QioControl::synchronize(const int errorStatus)  const {
  int error = errorStatus;
 
  if(NumNodes()>1) {
    error = globalSumInt(error);
    if(error > 0) {
      cout << "Totally " << error << " nodes reported error!" << endl;
    }
  }
  return error;
}

void QioControl::broadcastInt(int * data, int fromID)  const {
  if(NumNodes() > 1) {
    if(unique_id != fromID) {
      *data = 0;
    }

    *data = globalSumInt(*data);
  }
}

void QioControl::broadcastFloat(Float * data, int fromID) const {
  if(NumNodes() > 1) {
    if(unique_id != fromID) {
      * data = 0;
    }

    *data = globalSumFloat(*data);
  }
}

int QioControl::globalSumInt(const int data) const{
#if TARGET == QCDOC
  Gsum64Ext  gsum;
  return gsum.Sum(data);
#else
  return data;
#endif
}

unsigned int QioControl::globalSumUint(const unsigned int data) const{
#if TARGET == QCDOC
  Gsum64Ext  gsum;
  return gsum.Sum(data);
#else
  return data;
#endif
}

Float QioControl::globalSumFloat(const Float data) const {
#if TARGET == QCDOC
  Gsum64Ext  gsum;
  return gsum.Sum(data);
#else
  return data;
#endif
}

// IO control pattern:  two broadcast to set id range who got control
//                      read/write
//                      one sync to indicate finish, ret<0 means all finished
int QioControl::getIOTimeSlot(int block) const {
  if(!block) {
    cout << "non-blocking mode not implemented! sorry!" << endl;
    cout << "using blocking mode..." << endl;
  }

  if(NumNodes() > 1) {
    // using intelligent commander(node-0), dumb server(others) mode
    if(unique_id == 0) {
      return  IOCommander(0);
    }
    else {
      int firstID, lastID;
      while(1) {
	broadcastInt(&firstID);
	broadcastInt(&lastID);
	if(unique_id >= firstID && unique_id <= lastID) // got time slot
	  return 1;
	
	synchronize();
      }
    }
  }

  return 1;
}

int QioControl::finishIOTimeSlot(int block) const {
  if(!block) {
    cout << "non-blocking mode not implemented! sorry!" << endl;
    cout << "using blocking mode..." << endl;
  }
  
  if(NumNodes() > 1) {
    if(unique_id == 0) {
      return IOCommander(1);
    }
    else {
      if(synchronize()<0) return 0; // io finished
      while(1) {
	int dummy;
	broadcastInt(&dummy);
	broadcastInt(&dummy);  
	if(synchronize()<0) break;
      }
    }
  }

  return 0;
}


int QioControl::IOCommander(int caller, int block) const {
  if(!block) {
    cout << "non-blocking mode not implemented! sorry!" << endl;
    cout << "using blocking mode..." << endl;
  }

  int totalnodes = NumNodes();
  int batches = totalnodes / num_concur_io;
  if(num_concur_io * batches < totalnodes)  batches ++;

  int firstID, lastID;

  if(caller == 0)  { // let node 0 finish its task first (w/ the first batch)
    firstID = 0;
    lastID = num_concur_io-1;
    if(lastID > totalnodes-1)  lastID = totalnodes-1;

    broadcastInt(&firstID);
    broadcastInt(&lastID);
    cout << "Parallel IO: Group 1, Node " << firstID << " thru Node " << lastID << endl;
    return 1;
  }
  else { // now node 0 finished his own io, can control others
    if(batches==1) {
      synchronize(-1);  // io finished
      return 0;
    }

    for(int i=1;i<batches;i++) {
      synchronize(0);  // one batch done, but still more
      firstID = i * num_concur_io;
      lastID = (i+1) * num_concur_io - 1;
      if(lastID > totalnodes-1)  lastID = totalnodes-1;

      broadcastInt(&firstID);
      broadcastInt(&lastID);
      cout << "Parallel IO: Group " << i+1 << ", Node "<<firstID<<" thru Node "<<lastID<<endl;
    }

    synchronize(-1);  // io finished
    return 0;
  }
}




CPS_END_NAMESPACE
