//------------------------------------------------------------------
/*!\file
  \brief Definition of smearing classes.

  AlgSmear, AlgApeSmear, AlgKineticSmear and AlgHypSmear classes.
  
  $Id: alg_smear.h,v 1.4 2004-09-02 16:53:01 zs Exp $
*/
//------------------------------------------------------------------
/*
  alg_smear.h
  
  chris dawson
*/
#ifndef _ALG_SMEAR_CD_
#define _ALG_SMEAR_CD_

#include <config.h>
#include <alg/alg_base.h>
CPS_START_NAMESPACE

//! Base class for smearing a lattice.
/*!
  Base class for smearing a lattice, with or without SU(3)
  projection of the smeared links (controlled with a flag to constructor).

  Also allows for restricting smearing to hyperplanes orthogonal to a
  specified direction. On the base class level this will
  only excludes the links in that direction from being
  smeared; it is the responsibility of the derived classes
  to exclude links from being smeared with paths that move in
  that direction.

  \ingroup alg
*/

class AlgSmear:public Alg
{
  private:
  
    int bool_su3_proj;
    int orthog       ;
    Matrix* lat_back ;
    char *cname;

    void sub( Matrix& x, Matrix& y , int ind );

  public:
  
    AlgSmear( Lattice&   lat,
	      CommonArg* ca ,
	      int su3_proj );
  
    virtual ~AlgSmear();

    //! Perform the smearing
    void run();
  
    //! Restrict the smearing to some hyperplanes.
    /*!
      Smear only on hyperplanes orthogonal to this direction.
      \param i The direction.
    */
    void set_orthog( int i ) { orthog=i; }

    //! Get the direction orthogonal to which smearing is done.
    /*!
      \return The direction normal to the smeared hyperplanes.
    */
    int  get_orthog() const  { return orthog; }

  protected:

    //! Smear a  link.
    /*!
      \param link The link to smear.
      \param pos The coordinates of the link.
      \param mu The direction of the link.
    */
    virtual void smear_link( Matrix& link, 
			     int*    pos,
			     int      mu )=0;

    //! Calculates the three-link staple around a link.
    void three_staple(Lattice& latt, Matrix& link, int *pos, int u, int orth);

    //! Projects a matrix on to the SU(3)manifold.
    int su3_proj( Matrix& x );
    
};


//! Performs APE smearing 
/*!
  Performs APE smearing with SU(3) projection of the smeared links

  Also allows for restricting smearing to hyperplanes orthogonal to a
  specified direction.

  \ingroup alg  
*/

class AlgApeSmear:public AlgSmear
{
private:
  
  //! the smearing coefficient
  Float c;
  char *cname;
  
public:

  /*!
    \param lat The Lattice object containing the gauge field with which
    smearing is done.
    \param ca Container for generic parameters. .
  */
  
  AlgApeSmear(Lattice&   lat,
              CommonArg* ca  ):
    AlgSmear(lat,ca,1)
      {
	  cname = "AlgApeSmear";
      }
  
  ~AlgApeSmear()
  {;}

  //! Set the smearing coefficient
  /*! \param c The smearing coefficient. */
  void set_coef( Float x ) { c = x ; }

  protected:
  
  void smear_link( Matrix& mat,
                   int*    pos,
                   int      mu );
};


//! Kinetic smearing
/*!
  Performs the smearing needed to define an improved kinetic term vertex.
  Because of this the resulting smoothed configuration need not be unitary.

  Allows for restricting smearing to hyperplanes orthogonal to a
  specified direction.
  
  The default smearing coefficients are all zero.

  \ingroup alg
*/
class AlgKineticSmear:public AlgSmear
{

  private:
  
    Float _coef[5];
    char *cname;

    void five_staple(Lattice& latt, Matrix& link, int *pos, int u, int orth);
    void seven_staple(Lattice& latt, Matrix& link, int *pos, int u, int orth);
    void lepage_staple(Lattice& latt,  Matrix& link, int *pos, int u, int orth);

  public:
  
  /*!
    \param lat The Lattice object containing the gauge field with which
    smearing is done.
    \param ca Container for generic parameters. .
    \post The  smearing coefficients are all set to zero.
  */
  AlgKineticSmear(Lattice&   lat,
                  CommonArg* ca  ):
    AlgSmear(lat,ca,0)
  {
    for (int i=0;i<5;++i)  _coef[i]=0.0;
    cname = "AlgKineticSmear";
  }
  
  ~AlgKineticSmear()
  {;}
  
  //! Do the smearing  
  void run();

  //!  Sets the single link coefficient/
  void single_link( Float coef ) { _coef[0] = coef; }
  //!  Sets the three-link staple coefficient.
  void three_link ( Float coef ) { _coef[1] = coef; }
  //!  Sets the five-link staple coefficient.
  void five_link  ( Float coef ) { _coef[2] = coef; }
  //!  Sets the seven-link staple coefficient.
  void seven_link ( Float coef ) { _coef[3] = coef; }
  //!  Sets the Lepage-term staple coefficient.
  void lepage     ( Float coef ) { _coef[4] = coef; }

  protected:
  
  void smear_link( Matrix& link,
                   int*    pos,
                   int      mu );

};


//! Performs HYP smearing
/*!
  Hyp smearing: taking the definition from hep-lat/0103029
  (http://xxx.lanl.gov/abs/hep-lat/0103029)

  Allows for restricting smearing to hyperplanes orthogonal to a
  specified direction.

  \ingroup alg  
*/
class AlgHypSmear:public AlgSmear
{
private:
  
  Float c1, c2, c3;
  char *cname;

  const Matrix GetLink( Lattice& lat, const int* x, int mu );
  void get_vtilde( Matrix& , int*, int, int );
  void get_vbar  ( Matrix& , int*, int, int, int );

public:
  /*!
    \param lat The Lattice object containing the gauge field with which
    the smearing is done. 
    \param ca Container for generic parameters. .
  */
  AlgHypSmear(Lattice&   lat, CommonArg* ca  ): AlgSmear(lat,ca,1)
      {
	  cname = "AlgHypSmear";
      }

  ~AlgHypSmear()
  {;}

  //! Do the smearing
  void run();

  //!  Sets a smearing coefficient
  void set_c1( Float x ) { c1=x; }
  //!  Sets a smearing coefficient
  void set_c2( Float x ) { c2=x; }
  //!  Sets a smearing coefficient  
  void set_c3( Float x ) { c3=x; }
  
  protected:
  
  void smear_link( Matrix& link,
                   int*    pos,
                   int      mu );

};


CPS_END_NAMESPACE
#endif /* _ALG_SMEAR_CD_ */

