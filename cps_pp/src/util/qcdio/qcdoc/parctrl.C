#include <sysfunc.h>
#include <iostream>
using namespace std;

#include <util/parctrl.h>

CPS_START_NAMESPACE

ParallelControl::ParallelControl() 
  : num_concur_io(8)
{
  // TODO: set nodes[4], set coor[4], set unique_id
  // waiting for Bob's reply...

  cout << "I am on a " << SizeX() << "x"<< SizeY() << "x"<< SizeZ() << "x"<< SizeT() <<"x"
       <<SizeS()<<"x"<<SizeW()<< " lattice" << endl;
  cout << "My pos is (" << CoorX() << ","<< CoorY() << ","<< CoorZ() << ","<< CoorT() << ","
       << CoorS() << "," << CoorW() << ")" << endl;
  cout << "My UniqueID() = " << UniqueID() << endl;

  // naive version, just check dim_s = 1, dim_w = 1, and other 4 directions are used
  if(SizeS()!= 1 || SizeW()!=1) {
    cout << "Dimension > 4 not implemented!" << endl;
    exit(13);
  }
  nodes[0] = SizeX();
  nodes[1] = SizeY();
  nodes[2] = SizeZ();
  nodes[3] = SizeT();

  coor[0] = CoorX();
  coor[1] = CoorY();
  coor[2] = CoorZ();
  coor[3] = CoorT();

  unique_id = CoorX() + SizeX() * (CoorY() + SizeY() * (CoorZ() + SizeZ() * CoorT()));

}

ParallelControl::~ParallelControl() {

}


// The following are NEW functions added to ParallelControl class 
// to enable message passing between parallel processors, based on QMP calls
// (some are pretty useful...)

int ParallelControl::synchronize(const int errorStatus)  const {
  int error = errorStatus;
 
  if(NumNodes()>1) {
    QMP_sum_int(&error);
    if(error > 0) {
      cout << "Totally " << error << " nodes reported error!" << endl;
    }
  }
  return error;
}

void ParallelControl::broadcastInt(int * data, int fromID)  const {
  if(NumNodes() > 1) {
    if(unique_id != fromID) {
      *data = 0;
    }
    QMP_sum_int(data);
  }
}

void ParallelControl::broadcastFloat(Float * data, int fromID) const {
  if(NumNodes() > 1) {
    if(unique_id != fromID) {
      * data = 0;
    }
    QMP_float_t  buf = *data;
    QMP_sum_float(&buf);
    *data = buf;
  }
}

unsigned int ParallelControl::globalSumUint(const unsigned int data) const{
  if(NumNodes() > 1) {
    // i cannot find global sum on uint, so sum on 2 halves respectively
    // on 32-bit int machine, this proc may break down when more then 32768 nodes
    int halfword[2];
    int shift = sizeof(unsigned int) / 2 * 8;
    int mask = (0x1<<shift)-0x1;
    halfword[0] = data & mask;             // right bits
    halfword[1] = data >> shift;           // left bits

    QMP_sum_int(&halfword[0]);
    QMP_sum_int(&halfword[1]);
    
    unsigned int result;
    result = halfword[1];
    result = (result << shift) + halfword[0];
    
    return result;
  }
  else
    return data;
}

Float ParallelControl::globalSumFloat(const Float data) const {
  QMP_float_t buf = data;
  
  if(NumNodes() > 1) {
    QMP_sum_float(&buf);
  }

  return buf;
}

// IO control pattern:  two broadcast to set id range who got control
//                      read/write
//                      one sync to indicate finish, ret<0 means all finished
int ParallelControl::getIOTimeSlot(int block) const {
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

int ParallelControl::finishIOTimeSlot(int block) const {
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


int ParallelControl::IOCommander(int caller, int block) const {
  if(!block) {
    cout << "non-blocking mode not implemented! sorry!" << endl;
    cout << "using blocking mode..." << endl;
  }

  int totalnodes = Xnodes() * Ynodes() * Znodes() * Tnodes();
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
