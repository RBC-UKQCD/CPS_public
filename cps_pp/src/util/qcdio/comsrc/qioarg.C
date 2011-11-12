#include <config.h>
#include <iostream>
#include <unistd.h>
#include <util/gjp.h>
#include <sys/time.h>
#include <iomanip>
#include <cstring>
#include <comms/glb.h>
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

// Make it all periodic as NERSC header specifies gauge boundary condition, 
// 04/03/05 CJ
  for(int dir=0;dir<4;dir++) 
    bc[dir] = BND_CND_PRD;
//    bc[dir] = GJP.Bc(dir);

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
  : num_concur_io(0), do_log(0), cname("QioControl"), io_good(false)
{
  //  cout << "I am on a " << GJP.Xnodes() << "x"<< GJP.Ynodes() << "x"<< GJP.Znodes() 
  //       << "x"<< GJP.Tnodes() <<"x"<<GJP.Snodes() << " machine" << endl;
  //  cout << "My pos is (" << GJP.XnodeCoor() << ","<< GJP.YnodeCoor() << ","<< GJP.ZnodeCoor() << ","
  //       << GJP.TnodeCoor() << "," << GJP.SnodeCoor() << ")" << endl;

  unique_id = GJP.XnodeCoor() + GJP.Xnodes() * (GJP.YnodeCoor() + GJP.Ynodes() * (GJP.ZnodeCoor() + GJP.Znodes() * (GJP.TnodeCoor() + GJP.Tnodes() * GJP.SnodeCoor() ) ) );
  //  cout << "My UniqueID() = " << unique_id << endl;

  number_nodes = GJP.Xnodes() * GJP.Ynodes() * GJP.Znodes() * GJP.Tnodes() * GJP.Snodes();
  //  cout << "Total number of nodes = " << number_nodes << endl;

}

QioControl::~QioControl() {

}

//#ifndef USE_QMP
#if (!defined USE_QMP ) &&  (TARGET != NOARCH )
// The following are NEW functions added to QioControl class 
// to enable message passing between parallel processors, based on QMP calls
// (some are pretty useful...)

int QioControl::synchronize(const int errorStatus)  const {
  const char * fname = "synchronize()";
  int error = errorStatus;
 
  if(NumNodes()>1) {
    error = globalSumInt(error);
    if(error > 0) {
      VRB.Flow(cname,fname,"Totally %d nodes reported error!\n",error);
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

int QioControl::round(const Float fdata) const{
  int ndata = (int)fdata;
  if(fdata - ndata >= 0.5) ndata++;
  if(fdata - ndata < -0.5) ndata--;
  return ndata;
}

int QioControl::globalSumInt(const int data) const{
#ifdef PARALLEL
  //  Gsum64Ext  gsum;
  //  return gsum.Sum(data);
  int hfbits = sizeof(unsigned int) * 8 / 2;
  unsigned int mask = (1 << hfbits) - 1;

  int sumd = data;
  int hi = sumd >> hfbits;
  int lo = sumd & mask;
  hi = round(globalSumFloat(hi));
  lo = round(globalSumFloat(lo));
  sumd = (hi<<hfbits)+lo;
  return sumd;
#else
  return data;
#endif
}

unsigned int QioControl::globalSumUint(const unsigned int data) const{
#ifdef  PARALLEL 
  //  Gsum64Ext  gsum;
  //  return gsum.Sum(data);
  int hfbits = sizeof(unsigned int) * 8 / 2;
  unsigned int mask = (1 << hfbits) - 1;

  unsigned int sumd = data;
  unsigned int hi = sumd >> hfbits;
  unsigned int lo = sumd & mask;
  hi = round(globalSumFloat(hi));
  lo = round(globalSumFloat(lo));
  sumd = (hi<<hfbits)+lo;
  return sumd;
#else
  return data;
#endif
}

Float QioControl::globalSumFloat(const Float data) const {
#ifdef  PARALLEL
  //  Gsum64Ext  gsum;
  //  return gsum.Sum(data);
  Float sumdata = data;
  glb_sum_five(&sumdata);
  return sumdata;
#else
  return data;
#endif
}

int QioControl::globalMinInt(const int data) const{
#ifdef  PARALLEL
  Float fdata = data;
  glb_min(&fdata);
  int res = round(fdata);
  return res;
#else
  return data;
#endif
}
#endif

// IO control pattern:  two broadcast to set id range who got control
//                      read/write
//                      one sync to indicate finish, ret<0 means all finished
int QioControl::getIOTimeSlot() const {
  const char * fname = "getIOTimeSlot()";

//  printf("Node %d: GetIOTimeSlot()\n",UniqueID());
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
	if(unique_id >= firstID && unique_id <= lastID){ // got time slot
//          printf("Node %d: Got time slot!\n",UniqueID());
	  return 1;
        }
	
	synchronize();
      }
    }
  }

  return 1;
}

int QioControl::finishIOTimeSlot() const {
//  printf("Node %d: finishIOTimeSlot()\n",UniqueID());
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


int QioControl::IOCommander(int caller) const {
  const char * fname = "IOCommander()";
  int totalnodes = NumNodes();
  int do_concur_io = num_concur_io;
  if(do_concur_io <= 0) do_concur_io = totalnodes;

  int batches = totalnodes / do_concur_io;
  if(do_concur_io * batches < totalnodes)  batches ++;

  int firstID, lastID;

  if(caller == 0)  { // let node 0 finish its task first (w/ the first batch)
    firstID = 0;
    lastID = do_concur_io-1;
    if(lastID > totalnodes-1)  lastID = totalnodes-1;
    printf("Node %d: IOCommander(%d) batches=%d firstID=%d lastID=%d\n",UniqueID(),caller, batches, firstID,lastID);

    broadcastInt(&firstID);
    broadcastInt(&lastID);
    VRB.Flow(cname, fname, "Parallel IO: Group 1, Node %d thru Node %d\n",firstID,lastID);
    return 1;
  }
  else { // now node 0 finished his own io, can control others
    if(batches==1) {
      synchronize(-1);  // io finished
      return 0;
    }

    for(int i=1;i<batches;i++) {
      synchronize(0);  // one batch done, but still more
      firstID = i * do_concur_io;
      lastID = (i+1) * do_concur_io - 1;
//    printf("Node %d: IOCommander(%d) firstID=%d lastID=%d\n",UniqueID(),caller, firstID,lastID);
      if(lastID > totalnodes-1)  lastID = totalnodes-1;

      broadcastInt(&firstID);
      broadcastInt(&lastID);
      VRB.Flow(cname,fname,"Parallel IO: Group %d, Node %d thru Node %d\n",i+1,firstID,lastID);
    }

    synchronize(-1);  // io finished
    return 0;
  }
}


void QioControl::buildNodesList(int * active_num, int * active_node_list, int this_active) const {
  *active_num = globalSumInt(this_active?1:0);
  for(int i=0;i < *active_num; i++) {
    int sendid;
    if(this_active)  sendid = uniqueID();
    else             sendid = NumNodes(); // > all possible uniqueID();
    active_node_list[i] = globalMinInt(sendid);
    if(active_node_list[i] == uniqueID()) this_active = 0;  // exclude the nodes already in list
  }
}

int QioControl::syncError(int this_error) const {
  const char * fname = "testError()";
  TempBufAlloc nodes_list_buf(NumNodes()*sizeof(int));
  int * nodes_list = nodes_list_buf.IntPtr();

  int error_nodes;
  buildNodesList(&error_nodes, nodes_list, this_error);
  if(error_nodes>0) {
    VRB.Flow(cname, fname,
	     "%d nodes report error! They are (if more than 10 nodes, only list first 10 ids):\n",
	     error_nodes);
    for(int i=0;i<10 && i<error_nodes; i++)      VRB.Flow(cname,fname,"error id = %d\n",nodes_list[i]);
    VRB.Flow(cname,fname,"####\n");
  }
  return error_nodes;
}

void QioControl::setLogDir(const char * LogDir) {
  do_log = 1;
  logging = 0;
  strcpy(log_dir,LogDir);
}


void QioControl::startLogging(const char * action) {
  const char * fname = "startLogging()";

  int error = 0;

  if(!do_log) return;
  logs.clear();

  char logname[256];
  sprintf(logname,    "%s/qcdio.log.%d",log_dir,uniqueID());
  //  sprintf(oldlogname, "%s/%d.log.old",log_dir,uniqueID());

  /*  
  // copy logs to oldlogs
  logs.open(oldlogname);
  if(!logs.is_open()) {
    cout << "LOG file [" << oldlogname << "] open failed!" << endl;
    logging = 0;
    return;
  }

  ifstream prevlogs(logname);
  if(prevlogs.is_open()) {
    logs << prevlogs.rdbuf();
    prevlogs.close();
  }
  logs.close();

  logs.clear();
  prevlogs.clear();

  */

  // clear new logs, copy oldlogs to new logs thus to
  // set file pointer to the end so that we can append new logs
  //  VRB.Flow(cname,fname,"Try open file %s",logname);
  logs.open(logname, ios_base::in | ios_base::out | ios_base::ate);
  if(!logs.good()) { // file doesn't exist?
    logs.clear();
    logs.open(logname, ios_base::out | ios_base::trunc); 
    if(!logs.good()) {
      logs.clear();
      VRB.Flow(cname,fname,"LOG file [%s] open failed!\n",logname);
      logging = 0;
      error = 1;
    }
  }
  if(syncError(error)>0) {
    ERR.FileA(cname,fname,logname); 
  }

  /*
  prevlogs.open(oldlogname);
  if(prevlogs.is_open()) {
    logs << prevlogs.rdbuf();
    prevlogs.close();
  }

  logs.clear(); // if prevlogs is empty, the logs may have a error bit set
  */
  
  /*
  char logfile[200];
  strcpy(logfile,log_dir);
  strcat(logfile,"/qcdio.log");
  logs = Fopen(ADD_ID, logfile, "a");
  if(!logs) error = 1;
  if(testError(error) > 0) {
    ERR.FileA(cname,fname,logfile);
  }
  */

  //  cout << "start logging..." << endl;

  // start logging  
  struct timeval tp;
  gettimeofday(&tp,NULL);
  log_start = tp.tv_sec;
  char logtime[100];
  strcpy(logtime,ctime(&log_start));
  logtime[strlen(logtime)-1] = '\0';  // cut the last '\n'

  logs << "LOG<" << uniqueID() << ">["<< logtime << "] ";
  if(action) logs << action;
  logs<<" : \t";
  log_point = logs.tellp();
  logs << "Processing" << endl << flush;

  logging = 1;
}

void QioControl::log(const char * short_note) {
  const char * fname = "log()";
  int error = 0;

  if(!do_log || !logging) return;
  //  if(!logs.is_open() || !logs.good())  error = 1;
  if(syncError(error)>0) {
    ERR.Hardware(cname,fname,"Wrinting to file qcdio.log.* failed");
  }

  //  cout << "continue logging..." << endl;
  struct timeval tp;
  gettimeofday(&tp,NULL);
  time_t tm_elapse = tp.tv_sec - log_start;

  logs.seekp(log_point);
  logs << tm_elapse;
  if(short_note) logs << "(" << short_note << ")";
  logs<<"\t";
  log_point = logs.tellp();
  logs<<"Processing" << endl << flush;
}

void QioControl::finishLogging(const char * ending_word) {
  const char * fname = "finishLogging()";
  int error = 0;

  if(!do_log || !logging) return;
  //  if(!logs.is_open() || !logs.good()) error=1;
  if(syncError(error)>0) {
    ERR.Hardware(cname,fname,"Closing file qcdio.log.* failed");
  }

  //  cout << "finish logging..." << endl;

  struct timeval tp;
  gettimeofday(&tp,NULL);
  time_t tm_elapse = tp.tv_sec - log_start;

  logs.seekp(log_point);
  if(ending_word) logs << ending_word;
  logs<< "["<<tm_elapse << " sec]";
  logs<< "          " << endl << flush; // erase any chars not overwritten
  logs.close();

  logging = 0;
}






CPS_END_NAMESPACE
