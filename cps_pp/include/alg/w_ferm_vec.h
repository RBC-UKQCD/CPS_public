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

  static Vector m_tmp1; 
  //DRAM buffer for scu transfer in GetFermion, has to be static, since it
  //will be returned by GetFermion by "retrun &recv", 
  //will be deallocated if not static 

  // not implemented
  FermionVector(const FermionVector&);
  FermionVector& operator=(const FermionVector&);

  // helper member function
  void zeroOut()  const;
      
public:
  // CTOR
  FermionVector();
  // DTOR
  ~FermionVector();

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
};

#endif // ! _INCLUDED_W_FERM_VEC

CPS_END_NAMESPACE
