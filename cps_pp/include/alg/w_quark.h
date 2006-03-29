#include<config.h>
CPS_START_NAMESPACE
#ifndef INCLUDED_W_QUARK
#define INCLUDED_W_QUARK
//---------------------------------------------------------------------------
// class WspectQuark
//---------------------------------------------------------------------------

CPS_END_NAMESPACE
#include <alg/w_fuzzing.h>
#include <alg/w_field.h>
#include <alg/w_ferm_vec.h>
CPS_START_NAMESPACE

class WspectQuark : public WspectGinfo {

public:
  // CTOR
  // pbp_outfile, mid_point_outfile and ap_corr_outfile are from phys_v4.0.0
  // DEVOperatorKind src_op_kind=UNIT, WspectFuzzing *srcfuzz_p=0 and
  // WspectFuzzing *srcfuzz_p=0 are included from phys_v3.11.4.xiaodong
  //-------------------------------------------------------------------------
  WspectQuark(Lattice &, 
	      char *outfile,
	      char *pbp_outfile,
	      char *mid_point_outfile,
	      char *ap_corr_outfile,       // output file for <A_0 P> correlators   
	      WspectArg & warg,      
	      CgArg &cg,// const WspectArg actually
	      const WspectHyperRectangle & whr,
	      DEVOperatorKind src_op_kind=UNIT,
	      WspectFuzzing *srcfuzz_p=0,
	      WspectField *fld_p=0);

  
  // DTOR
  //-------------------------------------------------------------------------
  ~WspectQuark();
  
  // ACCESSORS
  // const IFloat* Data_SP() is from phys_v4.0.0
  // int dataSize(), const IFloat* SourceSlice(), DEVOperatorKind srcOpKind()
  // are from phys_v3.11.9.aniso
  // all others are common
  //-------------------------------------------------------------------------
  const IFloat* Data()             const     { return d_data_p; }
  const IFloat* Data_SP1()      const     { return d_data_mid_point_1; }
  const IFloat* Data_SP2()      const     { return d_data_mid_point_2; }
  int dataSize()                  const          {return d_size;}
  const IFloat* SourceSlice()      const     { return source_matrix_p; }
  DEVOperatorKind srcOpKind() const { return src_op_kind;}

  operator const IFloat*()  const     { return d_data_p; }  // conversion

  // ACCESSOR  -- twice of the center of the source, for non-zero mom proj
  // left as in phys_v3.11.4.xiaodong (instead of SourceCenter as in phys_v4.0.0)
  //-------------------------------------------------------------------------
  const int* SourceCenter2() const    { return d_source_center2;}

  // CLONE   
  //-------------------------------------------------------------------------
  void Clone(void *dst)    const     {
    moveMem(dst, d_data_p, d_size * sizeof(*d_data_p));
  }

  // Static member functions
  //-------------------------------------------------------------------------
  static int weightSrcDirac()        {return d_weight_Dy;}
  static int weightSrcColor()        {return d_weight_Cy;}

  
  // FOR THE PURPOSE OF DEBUGGING
  //-------------------------------------------------------------------------
  void dumpData(char *filename) const;  
  void dumpSource(char *filename, FermionVector &source) const;  
  void CheckSU3(const Float *gauge) const;
  //void DisplayMatrix(const Float *matrix, const char *text) const;

private:
  // Static data members
  //-------------------------------------------------------------------------
  static char *d_class_name;
  static int   d_weight_Dy;          // to ease the access to
  static int   d_weight_Cy;          // the quark propagators
  static Matrix m_tmp1;              //static buffer for internode communication
  static Vector v_tmp1;              //static buffer for commu.
  static Vector v_tmp2[DIRACs];      //static buffer for commu(sink_dev)
  IFloat *      d_data_p;             
  IFloat *      d_data_mid_point_1;    // propagator with mid_point sink
  IFloat *      d_data_mid_point_2;    // second propagator with mid_point sink when the sink is not Ls/2
  int          d_size;  
  int          midplane;              //location of the midpoint plane          
  
  //constant references
  const Lattice &d_lat;  
  const int prop_direction;
  const int src_plane_position; //local coord
  const DEVOperatorKind src_op_kind;

  // source matrix - colour*colour matrix at each point on the
  // source_slice (typically [x,y.z], but could also be [x,y,t]
  // if propagation is along the z-direction.
  IFloat *source_matrix_p;
  int    source_matrix_size;
  
  // twice the center point of the source, [0,0,0,0] for point source.
  // left as in phys_v3.11.4.xiaodong rather than d_source_center
  int          d_source_center2[LORENTZs];  


  
  // below many additions follow which are inherited from
  // phys_3.11.4.xiaodong
  // they are here to accommodate more general sources and extended operators

  //fuzzing object pointers
  WspectFuzzing *src_fuzz_p;
  WspectFuzzing *sink_fuzz_p;
  //field object pointers
  WspectField *fld_p;

  // not implemented
  WspectQuark(const WspectQuark&);
  WspectQuark& operator=(const WspectQuark&);

  // added by T&X
  // these are all functions to modify the source
  // using the new pointer Float *source_slice
  // these functions are implemented in w_quark.C


  void copyGaugefixToSource(const Float *src_matrix, 
			    const WspectHyperRectangle &whr ) const;

  //convert derivative operator DEV_I into actual polarisation direction
  //because dev_dir!=prop_dir
  int devDirToPolDir(int dev_dir, int prop_dir) const;
 
  void doSourceOperator(Lattice &lat, const WspectHyperRectangle &whr,
			DEVOperatorKind src_op_kind) const;

  void doBEOperator(DEVOperatorKind src_op, const WspectHyperRectangle &whr) const;
  void symmetricDerivative(Lattice &lat, const WspectHyperRectangle &whr, int polarisation) const;

  void JacobiSource(Lattice &lat, const WspectHyperRectangle &whr, Float r) const;
  const Matrix *GetSource(const int *site, const WspectHyperRectangle &whr ) const;

  // used to calculated the Dagger in a very simplistic way
  //void Dagger(Float * dagger, const Float * matrix) const;
  void AddToSum(Float *sum, Float weight) const;
  void Equal(Float *orig) const;
  //void DisplayAllSource(const WspectHyperRectangle &whr) const;

  //sinkOperators 
 public:
 
  //do sink operator on a single sink time slice
  void doSinkOperator(int lclWall, DEVOperatorKind sink_op_kind, Float *prop_out, Float *prop_tmp, WspectFuzzing *sink_fuzz_p);
 private:
  //isFullProp=0/1 indicates whether the prop_in is a full propogator or a partial propogator
  //                         on a single time slice(eg. after one derivative operator)
  void sink_deriv(int lclWall, int isFullProp, const Float *prop_in, Float *prop_out, int dir) const ;

  const Vector *getPropData(int isFullProp, const IFloat *prop_data_p, const int *site, int Dy, int Cy, int Dx) const;
};


#endif // ! _INCLUDED_W_QUARK

/*SourceOperatorKind
 *SinkOperatorKind
 *DEV1, DEV2, DEV3 
 *Depending on proppagation direction, the actual direction of derivatives are:
 *   Note: prop_dir are 0(X),1(Y),2(Z),3(T)
 *      propdir-> 0  1  2  3
 *       DEV1     Y  X  X  X
 *       DEV2     Z  Z  Y  Y
 *       DEV3     T  T  T  Z
 */

CPS_END_NAMESPACE
