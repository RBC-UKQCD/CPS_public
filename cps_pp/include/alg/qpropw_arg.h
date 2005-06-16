/*  qpropw_arg.h 

    Kostas Orginos (February 2002)

*/

/*  The structure type QPropWArg holds the parameters specific to
    QPropW class   */

#ifndef INCLUDED_QPropW_ARG_H
#define INCLUDED_QPropW_ARG_H

#include <alg/cg_arg.h>
#include <util/vector.h>

CPS_START_NAMESPACE

enum SourceType {
  POINT      = 0 ,  
  VOLUME     = 1 , 
  WALL       = 2 ,
  BOX        = 3 ,
  RANDVOLUME = 4 , 
  RANDWALL   = 5 , 
  RANDSLAB   = 6 , 
  MESSEQ     = 7 , 
  PROT_U_SEQ = 8 , 
  PROT_D_SEQ = 9 ,
  UNDEF      = 10,
  DERIV      = 11 };

enum RandomType {
  GAUSS  = 0, 
  UONE   = 1, 
  ZTWO   = 2, 
  NORAND = 3 } ;

class QPropWArg {
public:

  //! CG arguments for quark propagator inversion
  CgArg cg;		

  //! filename from which propagator may be loaded
  char file[100];

  //! source location
  int x,y,z,t;

  //! width of slab for random source
  int slab_width;

  // box source size
  int box_start;
  int box_end;

  //! Gauge Fixing flags
  int gauge_fix_src;
  int gauge_fix_snk;

  //! random number generator type
  RandomType rng;

  //! random number seed
  int seed;

  //! should midpoint propagator be stored?
  int store_midprop;

  //! should propagator be saved to disk?
  int save_prop;

  //!  should (1+gamma_t)/2 projected sources be used? (good far baryons)
  int do_half_fermion;

  QPropWArg() :
    x(0),y(0),z(0),t(0),
    slab_width(1),
    box_start(0),box_end(0),
    gauge_fix_src(0), 
    gauge_fix_snk(0), 
    rng(NORAND),
	seed(1111),
    store_midprop(0),
    save_prop(0),
    do_half_fermion(0)
    { 
      file[0]=(char)0; // empty string 
    };
  
};
CPS_END_NAMESPACE

#endif /* !INCLUDED_QPropW_ARG_H */
