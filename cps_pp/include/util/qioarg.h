#ifndef __QIO_ARG__
#define __QIO_ARG__

#include <iostream>
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
  char FileName[256];
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
    init(file, 8, 0.01, FP_AUTOMATIC, INT_AUTOMATIC, 1); 
  }
  QioArg(const char * file, const Float chkprec) {  
    // used in ReadLatticePar
    init(file, 8, chkprec, FP_AUTOMATIC, INT_AUTOMATIC, 1);
  }

  QioArg(const char * file, const FP_FORMAT dataformat, const int recon_row_3){
    // used in WriteLatticePar
    init(file, 8, 0.01, dataformat, INT_AUTOMATIC, recon_row_3);
  }
  
  QioArg(const char * file, const INT_FORMAT dataintformat) {
    // used in LatRngIO
    init(file, 8, 0.01, FP_AUTOMATIC, dataintformat, 1);
  }

};



class QioControl {
 private:
  int unique_id;  // = my_data_pos_in_file = (((s*dim_t+t)*dim_z+z)*dim_y+y)*dim_x+x
  int number_nodes;
  int num_concur_io; // number of nodes to excecute concurrent io
  int IOCommander(int caller,int block=1) const; // caller 0 <get IO>;   
                                                 // caller 1 <finish IO>
 public:
  QioControl();
  virtual ~QioControl();

  int synchronize(const int errorStatus=0) const; // sync processors and share error status
  void broadcastInt(int * data, int fromID = 0) const;
  void broadcastFloat(Float * data, int fromID = 0) const;
  int globalSumInt(const int data)  const;
  unsigned int globalSumUint(const unsigned int data) const;
  Float globalSumFloat(const Float data) const;

  int getIOTimeSlot(int block=1) const;
  int finishIOTimeSlot(int block=1) const;

  inline void setConcurIONumber(int set_concur)   { 
    num_concur_io = (set_concur <= 0? NumNodes() : set_concur);
  }

  inline bool isRoot() { return (uniqueID() == 0); }
  inline int uniqueID() const { return unique_id; }
  inline int NumNodes()  const { return number_nodes; }
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
};


CPS_END_NAMESPACE

#endif
