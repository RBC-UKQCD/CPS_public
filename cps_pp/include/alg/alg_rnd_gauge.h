#ifndef __ROTATEGAUGE__CD
#define __ROTATEGAUGE__CD

#include <util/lattice.h>
#include <util/smalloc.h>
#include <util/pmalloc.h>
#include <alg/alg_base.h>
#include <alg/common_arg.h>
#include <alg/no_arg.h>
#include <util/random.h>
#include <util/vector.h>

CPS_START_NAMESPACE


//===================================
// print out matrix (for debugging)
//===================================

void printMatrix( const Matrix& x );

/*!
  Fill one of three su(2) subgroups with

  cos(theta) -i nz sin(theta)  , -(inx+ny)sin(theta)
  (-inx+ny)sin(theta)          , cos(theta) + i nz sin (theta) ..,

  where n* are the components of a 3d unit vector, and theta
  is a specified argument.
  
  other co-ord 1 on diagonal, zero otherwise.
*/

void RandMatrix( Float theta , Matrix& x , int sub );

inline void UpperRandMatrix( Float theta , Matrix& x ) { 
  RandMatrix( theta, x, 0 ); 
}

inline void LowerRandMatrix( Float theta , Matrix& x ) { 
  RandMatrix( theta, x, 1 ); 
}

inline void MixedRandMatrix( Float theta , Matrix& x ) { 
  RandMatrix( theta, x, 2 ); 
}
     
/*!  
  Calculates the gauge transformation matrices for a random local gauge
  transformation and places them in the matrices pointed to by
 
  Lattice::fix_gauge_ptr

  First a random choice is made of of one of three su(2) subgroups, which is
  them filled with random elements, controlled by an angle \theta ( As might
  be expected 0 and 2\pi correspond to no rotation )
*/
class AlgRandomGauge:public Alg
{
private:

  char *cname;

  //! rotation angle, defaults to 1 radian
  Float theta;

public:

  AlgRandomGauge ( Lattice& latt , CommonArg *c_arg );

  ~AlgRandomGauge() {;}


  void set_theta( Float t ) { theta=t; } 
  

  void run();

  /*!
    You might want to keep the gauge fixing matrices
    longer than the scope of this class so the memory
    for them isn't automaticaly freed in the destructor.
    (ugly, but follows convention of gauge-fixing code)
  */
  void free() { AlgLattice().FixGaugeFree(); }
};


/*!
  rotates the gauge field by the
  gauge tranformation stored in 

  Lattice::fix_gauge_ptr

  The gauge fixing matrices need to be
  stored LANDAU style: one matrix
  for each site.

*/
class AlgRotateGauge : public Alg
{
private:

  char  *cname;
  
public:
  
  AlgRotateGauge( Lattice& latt,
               CommonArg *c_arg ):
    Alg(latt,c_arg),
    cname("RotateGauge")
  {;}
  
  ~AlgRotateGauge() {;}
  
  void run();
    

private:

  /*!
    gets gauge transf. matrix (from neighbouring node
    if needed).
  */
  Matrix GetMat(Matrix *m_offset, const int *x, int nu,
		int *node_sites, int *m_dir_offset);

};

CPS_END_NAMESPACE
#endif




