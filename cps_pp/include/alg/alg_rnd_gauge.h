/*!\file
  \brief Definition of AlgRandomGauge and AlgRotateGauge classes.

  $Id: alg_rnd_gauge.h,v 1.2 2004-09-02 16:56:31 zs Exp $
*/
//----------------------------------------------------------------------

#ifndef __ROTATEGAUGE__CD
#define __ROTATEGAUGE__CD

#include <util/lattice.h>
#include <alg/alg_base.h>
#include <alg/common_arg.h>
#include <util/vector.h>

CPS_START_NAMESPACE



//! Produces a ramdom gauge tranfromation field
/*!
  Calculates the gauge transformation matrices for a random local gauge
  transformation and places them in the matrices pointed to by
  Lattice::fix_gauge_ptr

  First a random choice is made of of one of three su(2) subgroups, which is
  them filled with random elements, controlled by an angle theta ( As might
  be expected 0 and 2pi correspond to no rotation )

    \ingroup alg
*/
class AlgRandomGauge: public Alg
{
  private:

    char *cname;

    //! rotation angle, defaults to 1 radian
    Float theta;

    // parint out matrix (for debugging)
    void printMatrix( const Matrix& x );
  
    void RandMatrix( Float theta , Matrix& x , int sub );

    void UpperRandMatrix( Float theta , Matrix& x ) { 
	RandMatrix( theta, x, 0 ); 
    }

    void LowerRandMatrix( Float theta , Matrix& x ) { 
	RandMatrix( theta, x, 1 ); 
    }

    void MixedRandMatrix( Float theta , Matrix& x ) { 
	RandMatrix( theta, x, 2 ); 
    }

    int NX, NY, NZ, NT;
    int index(int, int, int, int);
	
  public:

    AlgRandomGauge ( Lattice& latt , CommonArg *c_arg );


    /*!
      \note The destructor does not free the memory allocated.
      Use AlgRandomGauge::free for this or Lattice::FixGaugeFree
    */
    ~AlgRandomGauge() {;}

    //! Set the rotation angle.
    /*!
      \param t The rotation angle.
    */
    void set_theta( Float t ) { theta=t; } 
  
    //! Run this algorithm.
    void run();

    //! Free memory allocated by this class.
    /*!
      You might want to keep the gauge fixing matrices
      longer than the scope of this class so the memory
      for them isn't automaticaly freed in the destructor.
      (ugly, but follows convention of gauge-fixing code)
    */
    void free() { AlgLattice().FixGaugeFree(); }
};


//! Application of a gauge transformation.
/*!
  Rotates the gauge field by the gauge transformation stored at
  Lattice::fix_gauge_ptr

  The gauge fixing matrices need to be stored LANDAU style; one matrix
  for each site.

  \ingroup alg
*/
class AlgRotateGauge : public Alg
{
private:

    char  *cname;
    int NX, NY, NZ, NT;
    int index(int, int, int, int, int);    
  
public:

  /*!
    \param latt The lattice object containing the gauge field and the gauge
    transformation field.
    \param c_arg Generic algorithm arguments.
  */
  AlgRotateGauge( Lattice& latt, CommonArg *c_arg ):
    Alg(latt,c_arg),
    cname("RotateGauge")
  {;}
  
  ~AlgRotateGauge() {;}

  // Apply the gauge transformation.
  void run();
    

private:

  /*!
    gets gauge transf. matrix (from neighbouring node
    if needed).
  */
  Matrix GetMat(const Matrix *m_offset, const int *x, int nu,
		const int *node_sites, const int *m_dir_offset);

};

CPS_END_NAMESPACE
#endif




