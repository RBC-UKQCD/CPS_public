//------------------------------------------------------------------
/*!\file
  \brief Definition of smearing classes.

  AlgSmear, AlgApeSmear, AlgKineticSmear and AlgHypSmear classes.
  
  AlgOlegSmear added by Yasumichi 6/8/07.

  $Id: alg_smear.h,v 1.9 2012-08-10 21:34:27 chulwoo Exp $
*/
//------------------------------------------------------------------
#ifndef INCLUDED_ALG_SMEAR_H
#define INCLUDED_ALG_SMEAR_H
#include <config.h>
#include <alg/alg_base.h>
#include <alg/ape_smear_arg.h>
#include <alg/kinetic_smear_arg.h>
#include <alg/hyp_smear_arg.h>


CPS_START_NAMESPACE

//! Base class for smearing a lattice.
/*!
  Base class for smearing a lattice, with or without SU(3)
  projection of the smeared links (controlled with a flag to constructor).

  Also allows for restricting smearing to hyperplanes orthogonal to a
  specified direction. On the base class level this will
  only exclude the links in that direction from being
  smeared; it is the responsibility of the derived classes
  to exclude links from being smeared with paths that move in
  that direction.

  \ingroup alg
*/

class AlgSmear:public Alg
{
private:
  char *cname;

  int bool_su3_proj;
  Float tolerance  ;
  int orthog       ;
  Matrix* lat_back ;
    
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

  //! set tolerance for the SU(3) projection 
  void set_tol( Float x ) { tolerance = x; }
  
  //! get tolerance for the SU(3) projection
  Float get_tol() const { return tolerance; }

  // YA
  int ifSu3Proj() const {return bool_su3_proj;}

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
  
};

//! Calculates the three-link staple around a link.
void three_staple(Lattice& latt, Matrix& link, int *pos, int u, int orth);

//! Calculates the five-link staple around a link.
void five_staple  (Lattice& latt, Matrix& link, int *pos, int u, int orth);

//! Calculates the seven-link staple around a link.
void seven_staple (Lattice& latt, Matrix& link, int *pos, int u, int orth);

//! Calculates the lepage staple around a link.
void lepage_staple(Lattice& latt, Matrix& link, int *pos, int u, int orth);
  
//! Projects a matrix on to the SU(3) manifold.
int su3_proj( Matrix& x , Float tolerance );


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
  
  AlgApeSmear(Lattice&     lat,
              CommonArg*   ca ,
              ApeSmearArg* asa,
  	      int 	 in_bool_su3_proj=1);

  ~AlgApeSmear()
  {;}
  
  /*!
    If an output file is specified in the CommonArg argument, then
    the smearing coefficients are written to the file.
  */
  void run();
  void smartrun();

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
  
  char *cname;
  Float _coef[5];
    
public:
  
  /*!
    \param lat The Lattice object containing the gauge field with which
    smearing is done.
    \param ca Container for generic parameters. .
    \post The  smearing coefficients are all set to zero.
  */
  AlgKineticSmear(Lattice&   lat,
                  CommonArg* ca ,
                  KineticSmearArg* ksa );
  
  ~AlgKineticSmear()
  {;}
  
  //! Do the smearing  
  void run();

  void single_link( Float coef ) { _coef[0] = coef; }
  void three_link ( Float coef ) { _coef[1] = coef; }
  void five_link  ( Float coef ) { _coef[2] = coef; }
  void seven_link ( Float coef ) { _coef[3] = coef; }
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

  Note: the current implementation is aggressively inefficient

  \ingroup alg  
*/

class AlgHypSmear:public AlgSmear
{
private:

  char *cname;  
  Float c1, c2, c3;
  
  const Matrix GetLink( Lattice& lat, const int* x, int mu );
  void get_vtilde( Matrix& , int*, int, int );
  void get_vbar  ( Matrix& , int*, int, int, int );
  
public:
  /*!
    \param lat The Lattice object containing the gauge field with which
    the smearing is done. 
    \param ca Container for generic parameters. .
  */
  AlgHypSmear(Lattice&     lat, 
              CommonArg*    ca,
              HypSmearArg* hsa );
  
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

class AlgOlegSmear:public AlgApeSmear {

private:
  
  char *cname;
  Matrix* lat_back2 ;

public:
  AlgOlegSmear(Lattice&     lat,
	       CommonArg*   ca,
	       ApeSmearArg* asa);
  ~AlgOlegSmear();

  Float coef;
  void run();
};

CPS_END_NAMESPACE
#endif /* _ALG_SMEAR_CD_ */

