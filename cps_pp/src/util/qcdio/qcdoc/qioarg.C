#include <iostream>
#include <util/gjp.h>
using namespace std;

#include <util/qioarg.h>

CPS_START_NAMESPACE


/////////////////////////////////////////////////////////////////
// QioArg members////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void QioArg::init(const char * file, const int concur_io_number, const Float chk_prec,
		  const FP_FORMAT file_format, const int recon_row_3) {
  Xnodes = GJP.Xnodes();  Ynodes = GJP.Ynodes();  Znodes = GJP.Znodes();  Tnodes = GJP.Tnodes(); 
  Snodes = 1;

  XnodeSites = GJP.XnodeSites();  YnodeSites = GJP.YnodeSites(); 
  ZnodeSites = GJP.ZnodeSites();  TnodeSites = GJP.TnodeSites(); 
  SnodeSites = 1;

  Xcoor = GJP.XnodeCoor(); Ycoor = GJP.YnodeCoor(); Zcoor = GJP.ZnodeCoor(); Tcoor = GJP.TnodeCoor();
  Scoor = 0;

//  Xbc = GJP.Xbc(); Ybc = GJP.Ybc(); Zbc = GJP.Zbc(); Tbc = GJP.Tbc();
  Xbc = Ybc = Zbc = Tbc = BND_CND_PRD;

  StartConfLoadAddr = GJP.StartConfLoadAddr();

  // user set params
  ConcurIONumber = concur_io_number;
  strcpy(FileName, file);
  CheckPrecision = chk_prec;
  FileFpFormat = file_format;
  ReconRow3 = recon_row_3;
}


/////////////////////////////////////////////////////////////////////
// QioControl members////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
QioControl::QioControl() 
  : num_concur_io(8)
{
  // TODO: set nodes[4], set coor[4], set unique_id
  // waiting for Bob's reply...

  cout << "I am on a " << GJP.Xnodes() << "x"<< GJP.Ynodes() << "x"<< GJP.Znodes() 
       << "x"<< GJP.Tnodes() <<"x"<<GJP.Snodes() << " lattice" << endl;
  cout << "My pos is (" << GJP.XnodeCoor() << ","<< GJP.YnodeCoor() << ","<< GJP.ZnodeCoor() << ","
       << GJP.TnodeCoor() << "," << GJP.SnodeCoor() << ")" << endl;

  // naive version, just check dim_s = 1, dim_w = 1, and other 4 directions are used
  if(GJP.Snodes()!= 1) {
    cout << "Dimension > 4 not implemented!" << endl;
    exit(13);
  }
  nodes[0] = GJP.Xnodes();
  nodes[1] = GJP.Ynodes();
  nodes[2] = GJP.Znodes();
  nodes[3] = GJP.Tnodes();

  coor[0] = GJP.XnodeCoor();
  coor[1] = GJP.YnodeCoor();
  coor[2] = GJP.ZnodeCoor();
  coor[3] = GJP.TnodeCoor();

  unique_id = GJP.XnodeCoor() + GJP.Xnodes() * (GJP.YnodeCoor() + GJP.Ynodes() * 
						(GJP.ZnodeCoor() + GJP.Znodes() * GJP.TnodeCoor()));
  cout << "My UniqueID() = " << unique_id << endl;

}

QioControl::~QioControl() {

}


// The following are NEW functions added to QioControl class 
// to enable message passing between parallel processors, based on QMP calls
// (some are pretty useful...)

int QioControl::synchronize(const int errorStatus)  const {
  int error = errorStatus;
 
  if(NumNodes()>1) {
    //    QMP_sum_int(&error);
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
    //    QMP_sum_int(data);
    *data = globalSumInt(*data);
  }
}

void QioControl::broadcastFloat(Float * data, int fromID) const {
  if(NumNodes() > 1) {
    if(unique_id != fromID) {
      * data = 0;
    }
    //    QMP_float_t  buf = *data;
    //    QMP_sum_float(&buf);
    *data = globalSumFloat(*data);
    //    *data = buf;
  }
}

int QioControl::globalSumInt(const int data) const{
  Gsum64Ext  gsum;
  return gsum.Sum(data);
}

unsigned int QioControl::globalSumUint(const unsigned int data) const{
  Gsum64Ext  gsum;
  return gsum.Sum(data);

  /*
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
  */
}

Float QioControl::globalSumFloat(const Float data) const {
  Gsum64Ext  gsum;
  return gsum.Sum(data);

  /*
  QMP_float_t buf = data;
  
  if(NumNodes() > 1) {
    QMP_sum_float(&buf);
  }

  return buf;
  */
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
