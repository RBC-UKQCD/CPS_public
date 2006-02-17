/*!
  qpropwGFL.h
   
  chris dawson
*/ 

#ifndef _QPROPWGFL_CD_
#define _QPROPWGFL_CD_

#include <config.h>
#include <alg/qpropw.h>
#include <util/error.h>

CPS_START_NAMESPACE

/*!
  Class to calculate a quark propagator fixed
  into landau gauge on source *and* sink
*/
class QPropWGFPointSrc: public QPropW
{
public:

  /*!
    Just set the CommonArg
   */
  QPropWGFPointSrc( Lattice  &  lat ,
                    CommonArg* c_arg ):
    QPropW( lat, c_arg )
  {;}


  /*!
    Sets the CommonArg and the QPropWArg
    does not calculate the propagator
  */
  QPropWGFPointSrc( Lattice  &  lat ,
                    QPropWArg*  arg ,
                    CommonArg* c_arg ):
    QPropW( lat, arg, c_arg )
  { ; }

  void SetSource( FermionVectorTp& src,
                  int spin,
                  int colour )
  {
    src.SetGFPointSource(AlgLattice(),
                         colour,
                         spin,
                         qp_arg.x,
                         qp_arg.y,
                         qp_arg.z,
                         qp_arg.t);
  }
};


/*!
  A disconnected propagator with a fixed momentum on one leg .
*/

class QPropWLandauGaugeVolumeSrc: public QPropWPointSrc
{
private:

  // the fixed momenta
  int p[4];

public:
  
  QPropWLandauGaugeVolumeSrc( Lattice  &  lat ,
                              CommonArg* c_arg ):
    QPropWPointSrc(lat,c_arg)
  {
    ERR.NotImplemented("QPropWLandauGaugeVolumeSrc",
                       "QPropWLandauGaugeVolumeSrc");
  }

  // This ctor does not automatically run the inverter
  // (unlike QPropWPointSrc)
  QPropWLandauGaugeVolumeSrc( Lattice  &  lat ,
                              QPropWArg*  arg ,
                              CommonArg* c_arg ):
    QPropWPointSrc(lat, c_arg )
  { qp_arg=*arg; }
  
  void SetMomenta( const int* pval )
  {
    int i;
    for (i=0;i<4;i++) { p[i] = pval[i]; }
  }

  void SetSource(  FermionVectorTp& src,
                   int spin,
                   int colour )
  {
    src.SetLandauGaugeMomentaSource(AlgLattice(), 
                                    colour      , 
                                    spin        ,
                                    p            ); 
  }

};


CPS_END_NAMESPACE

#endif /* _QPROPWGFL_CD_ */






