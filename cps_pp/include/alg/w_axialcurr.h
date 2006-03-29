#include<config.h>
CPS_START_NAMESPACE
#ifndef INCLUDED_W_AXIALCURR
#define INCLUDED_W_AXIALCURR

CPS_END_NAMESPACE
// #include <util/data_types.h>        // Float, Complex
// #include <alg/w_spect_arg.h>        // SourceKind
// #include <alg/w_ginfo.h>
// #include <alg/w_hyper_rect.h>
// #include <alg/w_gamma_mat.h>
CPS_START_NAMESPACE

//---------------------------------------------------------------------------
// Forward Declarations  -- classes defined in other translation units
//---------------------------------------------------------------------------
// class Lattice;                      // defined in util/include/lattice.h
// class CgArg;                        // defined in  alg/include/cg_arg.h
// class CommonArg;                    // defined in  alg/include/common_arg.h
// class AlgWspect;                    // defined in  alg/include/alg_w_spect.h
// class WspectGinfo;
// class WspectHyperRectangle;
// class WspectQuark;
// class WspectMomenta;

//---------------------------------------------------------------------------
// class WspectAxialCurrent
// This class is used to calculate <A_0 P> correlators for DWF lattice 
// where A_0 is the conserved axial current or the local axial current
// and P is \bar\psi \gamma^5 \psi
//---------------------------------------------------------------------------
class WspectAxialCurrent : public WspectGinfo {

public:
  // CTORs
  //-------------------------------------------------------------------------
  WspectAxialCurrent(Lattice & latt,
		     WspectArg & arg,
	             const WspectHyperRectangle & whr ,
		     char * ap_corr_outfile);
      
  // DTOR
  //-------------------------------------------------------------------------
  ~WspectAxialCurrent();

  // ACCESSOR
  //-------------------------------------------------------------------------
  void print() const;  


  // FOR THE PURPOSE OF DEBUGGING 
  //-------------------------------------------------------------------------
  void dumpData(char *filename) const;  

  // This member function does the <A_0 P> calculation for DWF lattice. 
  // It calls function measureConserved for conserved axial current A0,
  // and function measureLocal for local axial current A0.
  void measureAll(Vector *sol_5d); 


  // This member function does the global sum of all the data 
  // called after all the measurements are done
  void doSum();


  // PRIVATEs
  //-------------------------------------------------------------------------
private:
  static char *  d_class_name;
  char *          ap_filename;

  Float *           d_local_p;        // Complex[global_slice]
  Float *      d_local_p_wall;        // Complex[global_slice]
  Float *       d_conserved_p;        // Complex[global_slice]
  WspectArg &            warg;

  // references
  Lattice & d_lat; 
  const WspectHyperRectangle & d_whr;  

  // The following data members are here for computation efficiency.
  int                          ls_glb;
  int                          prop_dir;  
  int                          lclMin[LORENTZs]; 
  int                          lclMax[LORENTZs];
  int                          glb_walls;
  

  // The following two 4-d fields are used to store 
  // results obtained by applying proper sink to 5d propagator
  IFloat *       d_data_p1;      // 4d field 
  IFloat *       d_data_p2;      // 4d field

  // temporary spinor fields
  Float * v1_g5, * v1_next_g5, * tmp_p1, * tmp_p2;

  // Array with entries that point to 8 non-member functions.
  // These functions are called as follows:
  // sproj_tr[SprojType mu]();
  // For the various SprojTypes see enum.h
  // These functions return a color matrix in f constructed from
  // the spinors v, w using: 
  // f_(i,j) = Tr_spin[ (1 +/- gamma_mu) v_i w^dag_j ]
  //
  // num_blk is the number of spinors v, w. The routines 
  // accumulate the sum over spinors in f.
  //
  // v_stride and w_stride are the number of Floats between spinors
  // (a stride = 0 means that the spinors are consecutive in memory)
    void (*sproj_tr[8])(IFloat *f, 
		      IFloat *v, 
		      IFloat *w, 
		      int num_blk, 
		      int v_stride,
		      int w_stride) ;

  // This member function does the <A_0 P> calculation for DWF lattice
  // where A_0 is the conserved axial current and P is \bar\psi \gamma^5 \psi
  void measureConserved(Vector * );

  // A generic function to do the mesonic contraction for <A_0 P>
  // in is the pointer to the 4d fermion vector.
  // out is the resulting meson correlator. 
  // ---mflin,  March 2006
  void measureLocal(const Float *in, Float *out);

  // This member function does the <A_0 P> calculation for DWF lattice with
  // a point sink
  // where A_0 is the local axial current and P is \bar\psi \gamma^5 \psi
  void measureLocalPoint(Vector *ferm_vec_5d);

  // This member function does the <A_0 P> calculation for DWF lattice with
  // a wall sink
  // where A_0 is the local axial current and P is \bar\psi \gamma^5 \psi
  void measureLocalWall(Vector *ferm_vec_5d);

};

//-------------------------end of declaration of class WspectAxialCurrent
#endif // ! _INCLUDED_W_AXIALCURR

CPS_END_NAMESPACE
