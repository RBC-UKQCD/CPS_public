#ifndef __QIO_ARG__
#define __QIO_ARG__

#include <iostream>

//#include <qmp.h>
#include <util/lattice.h>

#include <util/fpconv.h>
#if TARGET == QCDOC
#include <util/gsum64ext.h>
#endif


CPS_START_NAMESPACE
using namespace std;


class QioArg {
 public:
  // added s axis for future expansion
  int Xnodes,     Ynodes,     Znodes,     Tnodes,     Snodes;
  int XnodeSites, YnodeSites, ZnodeSites, TnodeSites, SnodeSites;
  int Xcoor,      Ycoor,      Zcoor,      Tcoor,      Scoor; 
  BndCndType  Xbc, Ybc, Zbc, Tbc;
  Matrix * StartConfLoadAddr;
  char FileName[256];
  int ConcurIONumber;

  Float CheckPrecision;  // used in read

  enum FP_FORMAT   FileFpFormat;  // used in write
  int ReconRow3;

  void init(const char * file, const int concur_io_number, const Float chk_prec,
	    const FP_FORMAT file_format, const int recon_row_3);

 public:
  inline int VolNodeSites() const { return XnodeSites * YnodeSites * ZnodeSites * TnodeSites; }
  inline int VolSites()     const { return Xnodes * Ynodes * Znodes * Tnodes * VolNodeSites(); }
  inline BndCndType Bc(int dir) const { 
    BndCndType bc[4] = {Xbc, Ybc, Zbc, Tbc};
    return bc[dir];
  }
  inline int Nodes(int dir) const {
    int nodes[5] = {Xnodes, Ynodes, Znodes, Tnodes, Snodes};
    return nodes[dir];
  }
  inline int NodeSites(int dir) const {
    int nodesites[5] = {XnodeSites, YnodeSites, ZnodeSites, TnodeSites, SnodeSites};
    return nodesites[dir];
  }

 public:
  QioArg(const char * file) {  
    init(file, 8, 0.01, FP_AUTOMATIC, 1); 
  }
  QioArg(const char * file, const Float chkprec) {   // used in ReadLattice
    init(file, 8, chkprec, FP_AUTOMATIC, 1);
  }

  QioArg(const char * file, const FP_FORMAT dataformat, const int recon_row_3){ // used in WriteLattice
    init(file, 8, 0.01, dataformat, recon_row_3);
  }

};



class QioControl {
 private:
  int unique_id;  // = my_data_pos_in_file = ((t*dim_z+z)*dim_y+y)*dim_x+x
  int coor[4];
  int nodes[4];
  int num_concur_io; // number of nodes to excecute concurrent io
  int IOCommander(int caller,int block=1) const; // caller 0 <get IO>;   
                                                 // caller 1 <finish IO>
  /*
 private:  // information from GJP
  inline int Xnodes() const { return nodes[0]; }
  inline int Ynodes() const { return nodes[1]; }
  inline int Znodes() const { return nodes[2]; }
  inline int Tnodes() const { return nodes[3]; }
  inline int Xcoor()  const { return coor[0]; }
  inline int Ycoor()  const { return coor[1]; }
  inline int Zcoor()  const { return coor[2]; }
  inline int Tcoor()  const { return coor[3]; }

  */

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
  inline int NumNodes()  const { return nodes[0]*nodes[1]*nodes[2]*nodes[3]; }
};


CPS_END_NAMESPACE

#endif
