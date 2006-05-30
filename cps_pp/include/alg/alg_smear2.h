//------------------------------------------------------------------
/*!\file
  \brief Definition of smearing classes.

  AlgSmear, AlgApeSmear, AlgKineticSmear and AlgHypSmear classes.
  
  $Id: alg_smear2.h,v 1.2 2006-05-30 20:32:27 chulwoo Exp $
*/
//------------------------------------------------------------------
#ifndef _ALG_SMEAR_CD_
#define _ALG_SMEAR_CD_
#include <config.h>
#include <alg/alg_base.h>
#include <alg/ape_smear_arg.h>
//#include <alg/kinetic_smear_arg.h>
//#include <alg/hyp_smear_arg.h>


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

class AlgSmear2:public Alg
{
private:
  char *cname;

  int bool_su3_proj;
  Float tolerance  ;
  int orthog       ;
  Matrix* lat_back ;
    
public:
  
  AlgSmear2( Lattice&   lat,
            CommonArg* ca ,
            int su3_proj );
  
  virtual ~AlgSmear2();
  
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

protected:
  
  //! Smear a  link.
  /*!
    \param link The link to smear.
    \param pos The coordinates of the link.
    \param mu The direction of the link.
  */
  virtual void smear_link2( Matrix& link, 
                           int*    pos,
                           int      mu )=0;
  
  
};

//! Calculates the three-link staple around a link.
void three_staple2(Lattice& latt, Matrix& link, int *pos, int u, int orth);

//! Calculates the five-link staple around a link.
//void five_staple  (Lattice& latt, Matrix& link, int *pos, int u, int orth);

//! Calculates the seven-link staple around a link.
//void seven_staple (Lattice& latt, Matrix& link, int *pos, int u, int orth);

//! Calculates the lepage staple around a link.
//void lepage_staple(Lattice& latt, Matrix& link, int *pos, int u, int orth);
  
//! Projects a matrix on to the SU(3) manifold.
int su3_proj( Matrix& x , Float tolerance );


//! Performs APE smearing 
/*!
  Performs APE smearing with SU(3) projection of the smeared links

  Also allows for restricting smearing to hyperplanes orthogonal to a
  specified direction.

  \ingroup alg  
*/
class AlgApeSmear2:public AlgSmear2
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
  
  AlgApeSmear2(Lattice&     lat,
              CommonArg*   ca ,
              ApeSmearArg* asa );
  
  ~AlgApeSmear2()
  {;}
  
  /*!
    If an output file is specified in the CommonArg argument, then
    the smearing coefficients are written to the file.
  */
  void run();

  void set_coef( Float x ) { c = x ; }

protected:
  
  void smear_link2( Matrix& mat,
                   int*    pos,
                   int      mu );
};


CPS_END_NAMESPACE
#endif /* _ALG_SMEAR_CD_ */

