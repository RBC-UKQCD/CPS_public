/*!
  alg_smear.h
  
  chris dawson
*/
#ifndef _ALG_SMEAR_CD_
#define _ALG_SMEAR_CD_

#include <alg/alg_base.h>


/*!
  Base class for smearing a lattice, with or without su3
  projection of the smeared link (flag to constructor).

  Also allows for smearing on hyperplanes orthoganal to a
  specified direction. On the base class level this will
  only excludes the links in that direction from being
  smeared, it is the responsibility of the derived classes
  to exclude links from being smeared with paths that move in
  that direction.
*/
CPS_START_NAMESPACE
class AlgSmear:public Alg
{
private:
  
  int bool_su3_proj;
  int orthog       ;
  Matrix* lat_back ;
  
public:
  
  AlgSmear( Lattice&   lat,
            CommonArg* ca ,
            int _bool_su3_proj );
  
  virtual ~AlgSmear();
  
  void run();
  
  // smear only on hyperplanes orthogonal to this direction
  void set_orthog( int i ) { orthog=i; }
  int  get_orthog() const  { return orthog; }

  virtual void smear_link( Matrix& link,
                           int*    pos,
                           int      mu )=0;
};

/*!
  APE smearing 
*/

class AlgApeSmear:public AlgSmear
{
private:
  
  //! the smearing coefficient
  Float c;

public:

  AlgApeSmear(Lattice&   lat,
              CommonArg* ca  ):
    AlgSmear(lat,ca,1)
  {;}
  
  ~AlgApeSmear()
  {;}

  void set_coef( Float x ) { c = x ; }
  
  void smear_link( Matrix& mat,
                   int*    pos,
                   int      mu );
};


/*!
  Performs the smearing needed to define an
  improved kinetic term vertex. Because of this the 
  resulting smoothed configuration need not be 
  unitary.
*/
class AlgKineticSmear:public AlgSmear
{
private:
  
  Float _coef[5];
  
public:
  
  AlgKineticSmear(Lattice&   lat,
                  CommonArg* ca  ):
    AlgSmear(lat,ca,0)
  {
    int i;
    for (i=0;i<5;++i) { _coef[i]=0.0; }
  }
  
  ~AlgKineticSmear()
  {;}
  
  
  void run();

  void smear_link( Matrix& link,
                   int*    pos,
                   int      mu );

  /*
    Setting the coefficients in the smeared action
  */
  void single_link( Float coef ) { _coef[0] = coef; }
  void three_link ( Float coef ) { _coef[1] = coef; }
  void five_link  ( Float coef ) { _coef[2] = coef; }
  void seven_link ( Float coef ) { _coef[3] = coef; }
  void lepage     ( Float coef ) { _coef[4] = coef; }
};



/*!
  Hyp smearing: taking the definition from hep-lat/0103029
*/
class  AlgHypSmear:public AlgSmear
{
private:
  
  Float c1, c2, c3;
  
private:

  void get_vtilde( Matrix& , int*, int, int );
  void get_vbar  ( Matrix& , int*, int, int, int );

public:

  AlgHypSmear(Lattice&   lat,
              CommonArg* ca  ):
    AlgSmear(lat,ca,1)
  {;}

  ~AlgHypSmear()
  {;}

  void run();
  
  void smear_link( Matrix& link,
                   int*    pos,
                   int      mu );
  
  
  void set_c1( Float x ) { c1=x; }
  void set_c2( Float x ) { c2=x; }
  void set_c3( Float x ) { c3=x; }

};


CPS_END_NAMESPACE
#endif /* _ALG_SMEAR_CD_ */

