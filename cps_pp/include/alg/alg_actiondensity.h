#ifndef __ALG__ACTIONDENSITY__
#define __ALG__ACTIONDENSITY__
#include <config.h>

#include <util/lattice.h>
#include <util/smalloc.h>
#include <util/pmalloc.h> 
#include "alg_base.h" 
#include "common_arg.h"
#include "no_arg.h"

CPS_START_NAMESPACE

class AlgActionDensity : public Alg
{
private:

  const char *cname;

  void ZeroReal(Matrix& m);
  Complex ActionDensity(Matrix clovers[]);
  void CloverLeaf(Lattice& lattice, Matrix& pl,  int* pos, int mu, int nu);
  inline Float * GsiteOffset(Float * p, const int *x, const int *g_dir_offset);
  void PathOrdProd(Matrix & mat, int* x, int* dirs, int n, Float *gfield, int *g_dir_offset);
  void CloverLeaf(Matrix& pl, int* pos, int mu, int nu, Float *gfield, int *g_dir_offset);



public:


  AlgActionDensity(Lattice&    latt, 
             CommonArg *c_arg ):
    Alg(latt,c_arg),
    cname("AlgActionDensity")
  {;}
  
  virtual ~AlgActionDensity() {;}

  //If result!=NULL the result will be copied to that memory address as well as saving to disk
  void run(Float *result = NULL);
  void smartrun(Float *result = NULL);

};

CPS_END_NAMESPACE
#endif 
