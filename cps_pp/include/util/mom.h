#include<config.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------
//
// mom.h
//
// Header file for the Mom class.
//
// Mom calculates the phase factor for each
// lattice site given a number of momenta and the source parameters
//
//------------------------------------------------------------------


#ifndef INCLUDED_MOM_H
#define INCLUDED_MOM_H

CPS_END_NAMESPACE
#include <math.h>    // for cos and sin
#include<util/data_types.h>
CPS_START_NAMESPACE

struct MomArg {
  int no_of_momenta; // number of different momenta    
  int max_p1;        // max. number of momentum units into 1-direction
  int max_p2;        // max. number of momentum units into 2-direction
  int max_p3;        // max. number of momentum units into 3-direction
  int deg;           // control flag: average over degenerate momenta on/off
  int dir;           // propagation direction
  int src_begin[4];  // source 
  int src_end[4];

};



class Mom
{
 private:
    char *cname;
    MomArg *mom_arg; // argument structure for momentum states

    Float PI;
    int dir,i,j,k;      // propagation direction and 3 orthogonal
    int no_of_mom;      // number of momenta to be calculated
    int deg;            // calculate degenerate momenta separately/together
    int max_p1,max_p2,max_p3; // maximal momentum in orthogonal directions

    int nx[4];                // local lattice extent
    int glb_L[4];             // global lattice extent
    int glb_sour_center[4];   // global source location
    Complex *mom_fact;

 public:
    Mom();
    virtual ~Mom();
    void Init(MomArg *arg);
    void run(void);
        
    // returns the complex phase factor for momentum "imom" at site "s"
    Complex fact(int imom, int *s);
};

//------------------------------------------------------------------
// External declarations.
//------------------------------------------------------------------
extern Mom MOM;


#endif




CPS_END_NAMESPACE
