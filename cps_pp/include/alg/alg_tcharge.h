#ifndef __ALG__TCHARGE__
#define __ALG__TCHARGE__
#include <config.h>

#include <util/lattice.h>
#include <util/smalloc.h>
#include <util/pmalloc.h> 
#include "alg_base.h" 
#include "common_arg.h"
#include "no_arg.h"

CPS_START_NAMESPACE

/*!
  Constructs the topological charge from 
  
  1) the 1x1 plaquette loops about a point
  2) the 1x2 rectangle loops about a point

  expanding about the classical continuum limit. Two O(a^2)
  improved combinations are also quoted: one from combining the topological
  charge from 1) and 2) so that the O(a^2) term vanishes, and one from
  first constructing an O(a^2) improved definition of the field-strength
  tensor and from this the topological charge.
*/
class AlgTcharge : public Alg
{
private:

  const char *cname;
 
public:

  Float charge_2  ;
  Float charge     ;
  Float charge_clov;
  Float charge_rect;

public:


  AlgTcharge(Lattice&    latt, 
             CommonArg *c_arg ):
    Alg(latt,c_arg),
    cname("AlgTcharge")
  {;}
  
  virtual ~AlgTcharge() {;}
  
  void run();
  void smartrun();

};

CPS_END_NAMESPACE
#endif 
