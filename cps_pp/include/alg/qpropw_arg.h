/*  qpropw_arg.h 

    Kostas Orginos (February 2002)

*/

/*  The structure type QPropWArg holds the parameters specific to
    QPropW class   */

#ifndef INCLUDED_QPropW_ARG_H
#define INCLUDED_QPropW_ARG_H

#include <alg/cg_arg.h>
#include <util/vector.h>


enum SourceType {POINT      = 0 ,  
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

enum RandomType {GAUSS  = 0, 
		 UONE   = 1, 
		 ZTWO   = 2, 
		 NORAND = 3 } ;

class QPropWArg 
{
public:

  /*! 
    The conjugate gradient argument for 
    the quark propagator inversion 
  */
  CgArg cg;		

  //! Propagator filename
  /*!
    100 characters i.e 100 words should be enough I guess
  */
  char file[100]    ; 

  // source location:
  int x;

  int y;
 
  int z;

  int t;

  //! end of time slice range for random src
  int tEnd;

  // box source size
  int bstart;
  int bend;

  //Gauge Fixing flags

  int GaugeFixSrc ;
  int GaugeFixSnk ;

  //! Random number generator type
  RandomType rnd ;

  /*!
    store mid point prop ?
  */
  int StoreMidpointProp;

  /*!
    flag for saving prop to disk
  */
  int SaveProp;

  /*!
    flag for doing the (1+gamma_t)/2 projected sources
    this saves a factor of 4 for the baryon matrix elemets
  */
  int DoHalfFermion ;

  QPropWArg() :
    x(0),y(0),z(0),t(0),
    tEnd(0),
    bstart(0),bend(0),
    GaugeFixSrc(0), 
    GaugeFixSnk(0), 
    rnd(NORAND),
    StoreMidpointProp(0),
    SaveProp(0),
    DoHalfFermion(0)
    { 
      file[0]=(char)0 ; // empty string 
    } ;
  
};

#endif /* !INCLUDED_QPropW_ARG_H */
