//------------------------------------------------------------------
//
// qpropw.h
//
// Kostas Orginos  (February 2002)
//
// Header file for the QPropW and WilsonMatrix classes for 
// Wilson-like quarks.
//
// For now this is specific to three colors. The constructor
// will exit if the number of colors is not equal to three.
//
// This file contains the declarations of the QPropW class.
//
//------------------------------------------------------------------

#ifndef INCLUDED_QPROPW_H
#define INCLUDED_QPROPW_H

#include <math.h>
#include <string.h>
#include <util/gjp.h>
#include <util/lattice.h>
#include <util/random.h>
#include <util/rcomplex.h>
#include <util/vector.h>
#include <util/wilson.h>
#include <alg/alg_base.h>
#include <alg/common_arg.h>
#include <alg/fermion_vector.h>
#include <alg/wilson_matrix.h>
#include <alg/qpropw_arg.h>
#include <alg/diquark.h>

#define MIDPROP 1
#define PROP  0

CPS_START_NAMESPACE

//! The Quark Propagator (Wilson type) class.

class QPropW : public Alg  {

  //! pointer to 4d prop at s-mid-point
  WilsonMatrix* midprop;

protected:

  // Hueywen: move prop from private to protected
  //! pointer to 4d prop
  WilsonMatrix* prop;
  
  QPropWArg qp_arg; 
  
  //! The class name
  char *cname; 
  
public:

  QPropW(Lattice& lat, CommonArg* c_arg);

  QPropW(Lattice& lat, QPropWArg* arg, CommonArg* c_arg);

  //! copy constructor
  QPropW(const QPropW& rhs); 

  //! averaging constructor
  QPropW(QPropW& prop1, QPropW& prop2 ); 

  Float Mass()      const { return qp_arg.cg.mass; }    
  int SourceTime()  const { return qp_arg.t; }       
#if 0
  int PointSrcX()   const { return qp_arg.x; }       
  int PointSrcY()   const { return qp_arg.y; }        
  int PointSrcZ()   const { return qp_arg.z; }
#endif

  //! Is sink gauge fixed?
  int GFixedSnk() const { return qp_arg.gauge_fix_snk; } 
  //! Is source gauge fixed?


  int StoreMidprop() const { return  qp_arg.store_midprop; }

  //! Returns the half fermion flag
  int DoHalfFermion() const { return qp_arg.do_half_fermion; }

  void Run();
  
  void CG(FermionVectorTp&, FermionVectorTp&, FermionVectorTp&, int&, Float&);
   // HueyWen: addtional CG definition
  void CG(Lattice &lat, CgArg *arg, FermionVectorTp& source,
        FermionVectorTp& sol , int& iter, Float& true_res);
  void FixSol(FermionVectorTp& sol);
  void LoadRow(int spin, int color, FermionVectorTp&, FermionVectorTp&);
  void SetFileName(char *nm);
  void ShiftPropForward(int n);
  void ShiftPropBackward(int n);

  //! Allocates memory for prop or prop_mid
  void Allocate(int);

  //! Deallocates memory of prop or prop_mid 
  void Delete(int);   

  //! Sets the QPropW arguments
  void SetArgs(QPropWArg& arg){ qp_arg = arg; }

  //! computes .5(prop+Q)
  void Average(QPropW& Q); 

  //! Comunicate Wilson Matrices...
  WilsonMatrix& GetMatrix(const int *, WilsonMatrix&) const;

  virtual void RestoreQProp(char*, int mid);
  virtual void SaveQProp(char*, int mid);

  virtual void SetSource(FermionVectorTp& src, int spin, int color);

  virtual SourceType SrcType() { return UNDEF; }

  virtual Complex& rand_src(int i) const;

  /*! This is a better name for the WallWallProp */
  WilsonMatrix WallSinkProp(int t_sink);

  // This is for compatibility with the old alg_threept code 
#define   WallWallProp WallSinkProp

  // operator functions
  
  QPropW& operator=(const QPropW& rhs);

  /*! Returns the prop */
  WilsonMatrix& operator[](int i){ return prop[i]; }

  /*! Returns the midpoint prop */
  WilsonMatrix& operator()(int i){ return midprop[i]; }
  
  // DESTRUCTORS
  virtual ~QPropW();
  int GFixedSrc() const { return qp_arg.gauge_fix_src; } 
  int PointSrcX()   const { return qp_arg.x; }       
  int PointSrcY()   const { return qp_arg.y; }        
  int PointSrcZ()   const { return qp_arg.z; }

};
  


class QPropWWallSrc : public QPropW {

public:
  
  // CONSTRUCTORS
  QPropWWallSrc(Lattice& lat, CommonArg* c_arg);
  // Most likely it is not needed since the base class copy constructor 
  // will be used. CHECK it!!!
  //    QPropWWallSrc(const QPropWWallSrc& rhs);
  QPropWWallSrc(Lattice& lat, QPropWArg* arg, CommonArg* c_arg);
  QPropWWallSrc(QPropWWallSrc& prop1, QPropWWallSrc& prop2 );
  
  void SetSource(FermionVectorTp& src, int spin, int color);
  SourceType SrcType() { return WALL; }
};

class QPropWMomSrc : public QPropWWallSrc {

  ThreeMom mom;
  
public:
  
  // CONSTRUCTORS
  QPropWMomSrc(Lattice& lat, CommonArg* c_arg);
  QPropWMomSrc(const QPropWMomSrc& rhs);
  QPropWMomSrc(Lattice& lat, QPropWArg* arg, int *p, CommonArg* c_arg);
  
  void SetSource(FermionVectorTp& src, int spin, int color);
  
  ThreeMom Mom() { return mom; } 
};

class QPropWVolSrc : public QPropW {

public:
  
  // CONSTRUCTORS
  QPropWVolSrc(Lattice& lat, CommonArg* c_arg);
  // Most likely it is not needed since the base class copy constructor 
  // will be used. CHECK it!!!
  //    QPropWVolSrc(const QPropWVolSrc& rhs);
  QPropWVolSrc(Lattice& lat, QPropWArg* arg, CommonArg* c_arg);
  
  void SetSource(FermionVectorTp& src, int spin, int color);
  SourceType SrcType() { return VOLUME; }
};

//HueyWen and Peter for general smearing function
class QPropWGFLfuncSrc : public QPropW
{

  Float (*func)(int gx,int gy,int gz,int gt);

  public:

    // CONSTRUCTORS
  QPropWGFLfuncSrc(Lattice &lat,
                   CgArg *arg,
                   CommonArg *common_arg,
                   Float (*fn) (int,int,int,int)
                  );

    // DESTRUCTOR
  ~QPropWGFLfuncSrc();
  SourceType SrcType() { return FL_FUNC; }

};


class QPropWPointSrc : public QPropW {

public:
  
  // CONSTRUCTORS
  QPropWPointSrc(Lattice& lat, CommonArg* c_arg);
  // Most likely it is not needed since the base class copy constructor 
  // will be used. CHECK it!!!
  //    QPropWPointSrc(const QPropWPointSrc& rhs);
  QPropWPointSrc(Lattice& lat, QPropWArg* arg, CommonArg* c_arg);
  
  void SetSource(FermionVectorTp& src, int spin, int color);
  SourceType SrcType() { return POINT; }
};

class QPropWBoxSrc : public QPropW {

protected:
  QPropWBoxArg box_arg;
  
public:
  
  // CONSTRUCTORS
  QPropWBoxSrc(Lattice& lat, CommonArg* c_arg);
  // Most likely it is not needed since the base class copy constructor 
  // will be used. CHECK it!!!
  //    QPropWBoxSrc(const QPropWBoxSrc& rhs);
  QPropWBoxSrc(Lattice& lat, QPropWArg* arg, CommonArg* c_arg);
  
  void SetSource(FermionVectorTp& src, int spin, int color);
  SourceType SrcType(){ return BOX; }
//  int BoxSrcStart() const { return qp_arg.box_start; } 
//  int BoxSrcEnd()   const { return qp_arg.box_end; }     
};

class QPropWRand : public QPropW {
  
protected:

  QPropWRandArg rand_arg;
  Float* rsrc;
  
public:
  
  // CONSTRUCTORS
  QPropWRand(Lattice& lat, CommonArg* c_arg);
  
  QPropWRand(Lattice& lat, QPropWArg* , QPropWRandArg *, CommonArg*);

  QPropWRand(const QPropWRand& rhs);  //copy constructor

  //MANIPULATORS
  Complex& rand_src(int i) const;

  // operator functions
  QPropWRand& operator=(const QPropWRand& rhs);

  void ShiftPropForward(int n);
  void ShiftPropBackward(int n);

  //Free prop
  void AllocateRsrc(); // Allocates memory for rsrc
  void DeleteRsrc();    // Deallocates memory for rsrc
  void RestoreQProp(char*, int mid);
  void SaveQProp(char*, int mid);

  RandomType Rng() const { return rand_arg.rng; }
  
  // DESTRUCTOR
  virtual ~QPropWRand();

};


class QPropWRandWallSrc : public QPropWRand {
  
public:
  QPropWRandWallSrc(Lattice& lat, CommonArg* c_arg);
  QPropWRandWallSrc(Lattice& lat, QPropWArg* , QPropWRandArg *, CommonArg*);
  
  void SetSource(FermionVectorTp& src, int spin, int color);
  SourceType SrcType() { return RANDWALL; }

};

class QPropWRandVolSrc : public QPropWRand {

public:
  QPropWRandVolSrc(Lattice& lat, CommonArg* c_arg);
  QPropWRandVolSrc(Lattice& lat, QPropWArg* , QPropWRandArg *, CommonArg*);
   
  void SetSource(FermionVectorTp& src, int spin, int color);
  SourceType SrcType() { return RANDVOLUME; }
};

class QPropWRandSlabSrc : public QPropWRand {

protected:
  int slab_width;

public:
  QPropWRandSlabSrc(Lattice& lat, CommonArg* c_arg);
  QPropWRandSlabSrc(Lattice& lat, QPropWArg* , QPropWSlabArg *, CommonArg*);
  QPropWRandSlabSrc(Lattice& lat, QPropWArg* arg, Float* src, CommonArg* );
  
  void SetSource(FermionVectorTp& src, int spin, int color);
  SourceType SrcType() { return RANDSLAB; }

};

class QPropWSeq : public QPropW {

protected:
  QPropW& quark;
  ThreeMom mom;
  Float quark_mass;

public:
  
  // CONSTRUCTORS
  
  QPropWSeq(Lattice& lat, QPropW& q,  int *p, QPropWArg*, CommonArg*);
  
  //QPropWSeq(const QPropWSeq& rhs);  //copy constructor not needed
  
  ThreeMom Mom() { return mom; }

  //! Returns the quark mass of the source propagator  
  Float SourceMass(){ return quark_mass; }
};

class QPropWSeqMesSrc : public QPropWSeq {
  /*!
    Only does a single gamma matrix insersion at the sink
    This is only used for pseudoscalars and vectors only anyway.
    should be enough for the Nucleon Decay project
  */
  int gamma;
  
public:
  
  // CONSTRUCTORS
  
  QPropWSeqMesSrc(Lattice& lat, QPropW& quark,  int *p, 
				  int g, QPropWArg*, CommonArg*);
  
  void SetSource(FermionVectorTp& src, int spin, int color);
  SourceType SrcType() { return MESSEQ; }
};


class QPropWSeqBar : public QPropWSeq {

protected:
  
  ProjectType proj;
  
public:
  
  // CONSTRUCTORS
  
  QPropWSeqBar(Lattice& lat, QPropW& quark,  int *p, 
			   ProjectType pp, QPropWArg*, CommonArg*);

  ProjectType Projection() { return proj; }

};

class QPropWSeqProtUSrc : public QPropWSeqBar {
  
public:
  
  // CONSTRUCTORS

  QPropWSeqProtUSrc(Lattice& lat, QPropW& quark,  int *p, 
		ProjectType pp, QPropWArg*, CommonArg*);

  //QPropWSeq(const QPropWSeq& rhs);  //copy constructor not needed

  void SetSource(FermionVectorTp& src, int spin, int color);
  SourceType SrcType() { return PROT_U_SEQ; }
};

class QPropWSeqProtDSrc : public QPropWSeqBar {

public:
  
  // CONSTRUCTORS

  QPropWSeqProtDSrc(Lattice& lat, QPropW& quark,  int *p, 
					ProjectType pp, QPropWArg*, CommonArg*);

  //QPropWSeq(const QPropWSeq& rhs);  //copy constructor not needed
  
  void SetSource(FermionVectorTp& src, int spin, int color);
  SourceType SrcType() { return PROT_D_SEQ; }
};

// === for exponential smeared source ==================================

// exponential smeared source
class QPropWExpSrc : public QPropW
{

  protected:
  // exponential smeared source
  // phi(|x-y|) = exp_A*exp(-exp_B*|x-y|) for |x-y| <= exp_C
  //            = 0                       for |x-y| >  exp_C
//  Float exp_A,exp_B,exp_C ;

  QPropWExpArg exp_arg;

 public:
  
  // CONSTRUCTORS
  QPropWExpSrc(Lattice& lat, CommonArg* c_arg) ;

  QPropWExpSrc(Lattice& lat, QPropWArg* arg, QPropWExpArg *exp_arg, 
    CommonArg* c_arg);
  QPropWExpSrc(QPropWExpSrc& prop1, QPropWExpSrc& prop2 );
  QPropWExpSrc(QPropWExpSrc* prop1, QPropWExpSrc* prop2 );
  QPropWExpSrc(QPropW& prop1) ;

  void SetSource(FermionVectorTp& src, int spin, int color);
  SourceType SrcType(){return EXP ; }
};

// =====================================================================

class QPropWFactory {
 public:
  static QPropW * Create (Lattice &lat,SourceType type, QPropWArg *arg, 
       CommonArg *c_arg, void *e_arg);
  static void Destroy(QPropW *qp);
};

CPS_END_NAMESPACE

#endif
