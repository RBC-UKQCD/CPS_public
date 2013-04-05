#ifndef __QIO_ARG__
#define __QIO_ARG__

#include <iostream>
#include <fstream>
#include <string>
#include <map>

#include <util/lattice.h>
#include <util/gjp.h>
#include <util/fpconv.h>
#include <util/intconv.h>


CPS_START_NAMESPACE
using namespace std;


class QioArg {
 public:
  QioArg() { /* should init before use */ }   // required for dec. of a indep. object

 private:
  // added s axis for future expansion
  //  int Xnodes,     Ynodes,     Znodes,     Tnodes,     Snodes;
  int nodes[5];
  //  int XnodeSites, YnodeSites, ZnodeSites, TnodeSites, SnodeSites;
  int node_sites[5];
  //  int Xcoor,      Ycoor,      Zcoor,      Tcoor,      Scoor; 
  int coor[5];
  //  BndCndType  Xbc, Ybc, Zbc, Tbc;
  BndCndType bc[4];

 public:
  inline int Xnodes() const  { return nodes[0]; }
  inline int Ynodes() const  { return nodes[1]; }
  inline int Znodes() const  { return nodes[2]; }
  inline int Tnodes() const  { return nodes[3]; }
  inline int Snodes() const  { return nodes[4]; }

  inline int XnodeSites() const  { return node_sites[0]; }
  inline int YnodeSites() const  { return node_sites[1]; }
  inline int ZnodeSites() const  { return node_sites[2]; }
  inline int TnodeSites() const  { return node_sites[3]; }
  inline int SnodeSites() const  { return node_sites[4]; }

  inline int Xcoor() const  { return coor[0]; }
  inline int Ycoor() const  { return coor[1]; }
  inline int Zcoor() const  { return coor[2]; }
  inline int Tcoor() const  { return coor[3]; }
  inline int Scoor() const  { return coor[4]; }

  inline int Xbc() const  { return bc[0]; }
  inline int Ybc() const  { return bc[1]; }
  inline int Zbc() const  { return bc[2]; }
  inline int Tbc() const  { return bc[3]; }

  inline BndCndType Bc(int dir) const { return bc[dir]; }
  inline int Nodes(int dir) const { return nodes[dir]; }
  inline int NodeSites(int dir) const { return node_sites[dir]; }
  inline int Coor(int dir) const { return coor[dir]; }

  inline void cutHalf() {  // used in loading/unloading LatRng data
    for(int d=0;d<5;d++)  node_sites[d] /= 2;
    if(node_sites[4]==0)  node_sites[4] = 1;
  }

 public:
  Matrix * StartConfLoadAddr;
  Float * StartU1ConfLoadAddr;

  char FileName[256];
//  string FileName;
  int ConcurIONumber;

  Float CheckPrecision;  // used in read

  enum FP_FORMAT   FileFpFormat;  // used in NERSC unloading
  enum INT_FORMAT  FileIntFormat; // used in RNG unloading
  int ReconRow3;

  void init(const char * file, const int concur_io_number, 
	    const Float chk_prec,
	    const FP_FORMAT file_format, const INT_FORMAT file_int_format,
	    const int recon_row_3);


 public:
  inline int VolNodeSites() const 
    { return node_sites[0] * node_sites[1] * node_sites[2] * node_sites[3]; }
  inline int VolSites()     const 
    { return nodes[0] * nodes[1] * nodes[2] * nodes[3] * VolNodeSites(); }

 public:
  QioArg(const char * file) {  
    init(file, 0, 0.01, FP_AUTOMATIC, INT_AUTOMATIC, 1); 
  }
  QioArg(const char * file, int concur_io_number) {  
    init(file, concur_io_number, 0.01, FP_AUTOMATIC, INT_AUTOMATIC, 1); 
  }
  QioArg(const char * file, const Float chkprec) {  
    // used in ReadLatticePar
    init(file, 0, chkprec, FP_AUTOMATIC, INT_AUTOMATIC, 1);
  }

  QioArg(const char * file, const FP_FORMAT dataformat, const int recon_row_3){
    // used in WriteLatticePar
    init(file, 0, 0.01, dataformat, INT_AUTOMATIC, recon_row_3);
  }
  
  QioArg(const char * file, const INT_FORMAT dataintformat) {
    // used in LatRngIO
    init(file, 0, 0.01, FP_AUTOMATIC, dataintformat, 1);
  }

};



class QioControl {
 private:
  char * cname;
  int unique_id;  // = my_data_pos_in_file = (((s*dim_t+t)*dim_z+z)*dim_y+y)*dim_x+x
  int number_nodes;
  int num_concur_io; // number of nodes to excecute concurrent io
  int IOCommander(int caller) const; // caller 0 <get IO>;   
                                                 // caller 1 <finish IO>
  void buildNodesList(int * active_num, int * active_node_list, int this_active) const;

 public:
  QioControl();
  virtual ~QioControl();

  int synchronize(const int errorStatus=0) const; // sync processors and share error status
  void broadcastInt(int * data, int fromID = 0) const;
  void broadcastFloat(Float * data, int fromID = 0) const;
  int globalSumInt(const int data)  const;
  unsigned int globalSumUint(const unsigned int data) const;
  Float globalSumFloat(const Float data) const;
  int round(const Float fdata) const;
  int globalMinInt(const int data) const;

  int getIOTimeSlot() const;
  int finishIOTimeSlot() const;

  inline void setConcurIONumber(int set_concur) { num_concur_io = set_concur; }

  inline bool isRoot() { return (uniqueID() == 0); }
  inline int uniqueID() const { return unique_id; }
  inline int NumNodes()  const { return number_nodes; }

  int syncError(int this_error) const;

 private:
  int do_log;
  int logging;
  fstream logs;
  long log_point;
  time_t log_start;
  char log_dir[200];

 public:
  void setLogDir(const char * LogDir);
  void startLogging(const char * action=0);
  void finishLogging(const char * ending_word=0);
  void log(const char * short_note=0);

 protected:
  FPConv fpconv;
  bool io_good;
 public:
  inline bool good() const { return io_good; }

 public:
  void SimQCDSP(int sim) { fpconv.SimQCDSP(sim); }
};



// a memory allocator that prevents memory leak
class TempBufAlloc { 
 private:
  char * buf;
 public:
  TempBufAlloc(const size_t size_chars) {  buf = (char*)fmalloc(size_chars); } // fast memory
  virtual ~TempBufAlloc() {  ffree(buf); }

  inline operator char*() { return buf; }
  inline char * CharPtr() { return buf; }
  //  operator int*() {  return (int*)buf; }
  inline int *IntPtr() { return (int*)buf; }
  inline Float *FPtr() { return (Float*)buf; }
};


CPS_END_NAMESPACE

#endif
