#include<config.h>
CPS_START_NAMESPACE
#ifndef INCLUDED_W_FERM_VEC
#define INCLUDED_W_FERM_VEC

//---------------------------------------------------------------------------
// class FermionVector 
//---------------------------------------------------------------------------
class FermionVector : public WspectGinfo {

private:
  static char *d_class_name;  

  // data -- local fermionic field and its size in units of Float
  //Storage order: d_fermion_p[T][Z][Y][X][DIRACS][COLORS][COMPLEXES]
  IFloat *d_fermion_p;
  int    d_size; 
  // flag to determine if memory is allocated
  static int allocated; 

  static Vector m_tmp1; 

  //location of the box source/sink
  //  int box_b[LORENTZs];//start coordniates
  //  int box_e[LORENTZs];//end coordinates

  // helper member function
  void zeroOut()  const;
      
public:
  // CTOR
  FermionVector();
  // secondary CTOR; use an existing array
  // This is dangerous, since user has to check if the array has the
  // right size for the fermion vector
  FermionVector(Float *);
  // DTOR
  ~FermionVector();

  // Copy Constructor
  FermionVector(const FermionVector& fv)
    {*this = fv;}
  
  // Assignment Operator
  FermionVector& operator=(const FermionVector& rhs);
  FermionVector& operator+=(const FermionVector &f1);
  
  // ACCESSORS
  Float * data(void)          const {return (Float *)d_fermion_p;}
  // Float * proj_data(void)     const {return (Float *)d_proj_fermion_p;}
  void    print(char *file)   const;
  void    printWaveFunc(char *file) const;
  // MANIPULATORS

  // Caller's responsibility:  check argument color and spin.

  void setPointSrc(int color, int spin, 
		   const int *point_global_coord) const;
  // point_global_coord = 0 means the source point is off-node.

  // removed overloaded in w_ferm_vec.C, not implemented any longer
  void setPointSrc(const int *point_global_coord) const;

  void setSource(int color, int spin, IFloat rs_fac,
		     const IFloat *src_matrix,
		     int          src_plane_dir,
		     const int *  src_local_origin ,
		     const int *  src_local_end  ) const;
  // src_matrix is a 3x3 matrix for each (local) site
  // if src_matrix = unit matrix at origin --> same as point source
  // if src_matrix = gauge_fixed_matrix --> same as setSmearedSource(Box,Wall)
  // more general cases are possible where
  // a) a derivative has been applied to (Point,Box or Wall) --> source_slice
  // b) a gauge invariant smearing has been defined on point_source --> source_slice
  // added by Thomas and Xiaodong, April 3rd.
  



  void setSmearedSrc(int color, int spin, 
		     const IFloat *src_matrix,
		     int          src_plane_dir,
		     const int *  src_local_origin = 0,
		     const int *  src_local_end = 0) const;
  // src_matrix = 0 means the source region is entirely off-node:
  //         the routine simply zero out the FermionVector.


  // removed overloaded function from w_ferm_vec.C
  // not implemented anylonger
  void setSmearedSrc(const IFloat *src_matrix,
  		     int          src_plane_dir,
  		     const int *  src_local_origin = 0,
  	             const int *  src_local_end = 0) const;
  



  // added by Thomas and Xiaodong: setJacobiSrc, doSourceOperator

  void setJacobiSrc(const int *point_global_coord, const Lattice &lat, Float epsi, int n_iter) const;

  // point_global_coord = 0 means the source point is off-node.
  // src_matrix should contain the source vector (ncolour x n_local_sites)
  // f_c(x). This vector should be set up in the cstor of WspectQuark

  const Vector *GetFermion(const int *site, int dirac) const;

  // removed, not implemented after March 30, 2000 (T&X)
  void projectSource(const int Cpick, const int Dpick) const;

  //----------------------------------------------------------------------
  // Coulomb gauge fix fermion solution vectors. added by mflin 01/31/06
  //----------------------------------------------------------------------
  void gaugeFixSink(Lattice &lat, int dir) const;

  //-------------------------------------------------------------------
  // Sum over hyperplanes. mflin 02/01/06
  //--------------------------------------------------------------------
  void sumOverHyperPlane(int dir, int box_b[], int box_e[]); //all hyperplanes in direction dir
  void sumOverHyperPlane(int dir, int hp, int box_b[],int box_e[]); //one hyperplanevoid

  void sumOverHyperPlaneStride(int dir, int box_b[], int box_e[]); //all hyperplanes in direction dir
  void sumOverHyperPlaneStride(int dir, int hp, int box_b[],int box_e[]); //one hyperplanevoid
  
  void sumOverHyperPlaneZeroMom(int dir, int box_b[], int box_e[]);


};

#endif // ! _INCLUDED_W_FERM_VEC

CPS_END_NAMESPACE
