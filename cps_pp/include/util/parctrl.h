#ifndef __PAR_CTRL__
#define __PAR_CTRL__

#include <iostream>

#include <qmp.h>
#include <util/lattice.h>


CPS_START_NAMESPACE
using namespace std;

class ParallelControl {
 protected:
  int unique_id;  // = my_data_pos_in_file = ((t*dim_z+z)*dim_y+y)*dim_x+x
  int coor[4];
  int nodes[4];
  int num_concur_io; // number of nodes to excecute concurrent io
  int IOCommander(int caller,int block=1) const; // caller 0 <get IO>;   
                                                 // caller 1 <finish IO>

 public:
  ParallelControl();
  virtual ~ParallelControl();

  int synchronize(const int errorStatus=0) const; // sync processors and share error status
  void broadcastInt(int * data, int fromID = 0) const;
  void broadcastFloat(Float * data, int fromID = 0) const;
  unsigned int globalSumUint(const unsigned int data) const;
  Float globalSumFloat(const Float data) const;
  int getIOTimeSlot(int block=1) const;
  int finishIOTimeSlot(int block=1) const;

 public:  // trivial functions
  inline int uniqueID() const { return unique_id; }
  inline void setConcurIONumber(int set_concur) 
    { 
      num_concur_io = (set_concur <= 0? NumNodes() : set_concur);
    }
  inline int Xnodes() const { return nodes[0]; }
  inline int Ynodes() const { return nodes[1]; }
  inline int Znodes() const { return nodes[2]; }
  inline int Tnodes() const { return nodes[3]; }
  inline int Xcoor()  const { return coor[0]; }
  inline int Ycoor()  const { return coor[1]; }
  inline int Zcoor()  const { return coor[2]; }
  inline int Tcoor()  const { return coor[3]; }
  inline int NumNodes()  const { return Xnodes()*Ynodes()*Znodes()*Tnodes(); }
};

CPS_END_NAMESPACE

#endif
