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

#define PROP_MID 1
#define PROP  0

CPS_START_NAMESPACE

//! The QPropW class.

class QPropW : public Alg 
{
  //! pointer to 4d prop 
  WilsonMatrix* prop;

  //! pointer to 4d prop at s-mid-point
  WilsonMatrix* prop_mid;

protected:
  
  /*!
    Copy of the arguments used for construction
    could be valuable in checking the restore from disk
  */
  QPropWArg Arg ; 
  
  //! The class name
  char *cname ; 
  
public:

  QPropW(Lattice& lat, CommonArg* c_arg);

  QPropW(Lattice& lat, QPropWArg* arg, CommonArg* c_arg);

  //! copy constructor
  QPropW(const QPropW& rhs); 

  //! averaging constructor
  QPropW(QPropW& prop1, QPropW& prop2 ); 

  Float Mass() const { return Arg.cg.mass ;}    
  int SourceTime() const { return Arg.t ;}       
  int BoxSrcStart() const { return Arg.bstart ;} 
  int BoxSrcEnd() const { return Arg.bend ;}     
  int PointSrcX() const { return Arg.x ;}       
  int PointSrcY() const { return Arg.y ;}        
  int PointSrcZ() const { return Arg.z ;}

  //! Is sink gauge fixed?
  int GFixedSink() const { return Arg.GaugeFixSnk;} 
  //! Is source gauge fixed?
  int GFixedSource() const { return Arg.GaugeFixSrc;} 

  RandomType Rnd() const { return Arg.rnd ;}

  int MidpointProp() const { return  Arg.StoreMidpointProp ;}

  //! Returns the half fermion flag
  int  HalfFermion() const { return Arg.DoHalfFermion ;}

  void run();
  
  void CG(FermionVectorTp&, FermionVectorTp&, FermionVectorTp&, int&, Float&);
  void FixSol(FermionVectorTp& sol);
  void LoadRow(int spin, int color, FermionVectorTp&, FermionVectorTp&);
  void setFileName(char *nm) ;
  void ShiftPropForward(int n);
  void ShiftPropBackward(int n);

  //! Allocates memory for prop or prop_mid
  void Allocate(int) ;

  //! Deallocates memory of prop or prop_mid 
  void Delete(int);   

  //! Sets the QPropW arguments
  void setArgs(QPropWArg& qargs){ Arg=qargs; }

  /*! computes .5(prop+Q) */
  void Average(QPropW& Q) ; 

  //! Comunicate Wilson Matrices...
  WilsonMatrix& GetMatrix(const int *, WilsonMatrix&) const;

  virtual void RestoreQProp(char*, int mid);
  virtual void SaveQProp(char*, int mid);

  virtual void SetSource(FermionVectorTp& src, int spin, int color);

  virtual SourceType SrcType() { return UNDEF; }

  virtual Complex& rand_src(int i) const ;

  /*! This is a better name for the WallWallProp */
  WilsonMatrix WallSinkProp(int t_sink);

  // This is for compatibility with the old alg_threept code 
#define   WallWallProp WallSinkProp

  // operator functions
  
  QPropW& operator=(const QPropW& rhs);

  /*! Returns the prop */
  WilsonMatrix& operator[](int i){return prop[i];}

  /*! Returns the midpoint prop */
  WilsonMatrix& operator()(int i){return prop_mid[i];}
  
  // DESTRUCTORS
  virtual ~QPropW();
  
  
};
  


class QPropWWallSrc : public QPropW
{
  
public:
  
  // CONSTRUCTORS
  QPropWWallSrc(Lattice& lat, CommonArg* c_arg) ;
  // Most likely it is not needed since the base class copy constructor 
  // will be used. CHECK it!!!
  //    QPropWWallSrc(const QPropWWallSrc& rhs);
  QPropWWallSrc(Lattice& lat, QPropWArg* arg, CommonArg* c_arg);
  QPropWWallSrc(QPropWWallSrc& prop1, QPropWWallSrc& prop2 );
  QPropWWallSrc(QPropWWallSrc* prop1, QPropWWallSrc* prop2 );
  
  //Not used
  //QPropWWallSrc& Sum(QPropWWallSrc* prop1, QPropWWallSrc* prop2 );
  
  void SetSource(FermionVectorTp& src, int spin, int color);
  SourceType SrcType(){return WALL ; }
};

class QPropWMomSrc : public QPropWWallSrc
{

  ThreeMom mom ;
  
public:
  
  // CONSTRUCTORS
  QPropWMomSrc(Lattice& lat, CommonArg* c_arg) ;
  QPropWMomSrc(const QPropWMomSrc& rhs);
  QPropWMomSrc(Lattice& lat, QPropWArg* arg, int *p, CommonArg* c_arg);
  
  void SetSource(FermionVectorTp& src, int spin, int color);
  
  ThreeMom Mom(){return mom;} 
};

class QPropWVolSrc : public QPropW
{
  
public:
  
  // CONSTRUCTORS
  QPropWVolSrc(Lattice& lat, CommonArg* c_arg) ;
  // Most likely it is not needed since the base class copy constructor 
  // will be used. CHECK it!!!
  //    QPropWVolSrc(const QPropWVolSrc& rhs);
  QPropWVolSrc(Lattice& lat, QPropWArg* arg, CommonArg* c_arg);
  
  void SetSource(FermionVectorTp& src, int spin, int color);
  SourceType SrcType(){return VOLUME ; }
};

class QPropWPointSrc : public QPropW
{
  
public:
  
  // CONSTRUCTORS
  QPropWPointSrc(Lattice& lat, CommonArg* c_arg) ;
  // Most likely it is not needed since the base class copy constructor 
  // will be used. CHECK it!!!
  //    QPropWPointSrc(const QPropWPointSrc& rhs);
  QPropWPointSrc(Lattice& lat, QPropWArg* arg, CommonArg* c_arg);
  
  void SetSource(FermionVectorTp& src, int spin, int color);
  SourceType SrcType(){return POINT ;}
};

class QPropWBoxSrc : public QPropW
{

  
  
public:
  
  // CONSTRUCTORS
  QPropWBoxSrc(Lattice& lat, CommonArg* c_arg) ;
  // Most likely it is not needed since the base class copy constructor 
  // will be used. CHECK it!!!
  //    QPropWBoxSrc(const QPropWBoxSrc& rhs);
  QPropWBoxSrc(Lattice& lat, QPropWArg* arg, CommonArg* c_arg);
  
  void SetSource(FermionVectorTp& src, int spin, int color);
  SourceType SrcType(){return BOX ;}
};

class QPropWRand : public QPropW
{
  
protected:

  Float* rsrc;
  
public:
  
  // CONSTRUCTORS
  QPropWRand(Lattice& lat, CommonArg* c_arg);
  
  QPropWRand(Lattice& lat, QPropWArg* arg,  CommonArg*);

  QPropWRand(const QPropWRand& rhs);  //copy constructor

    
    
  //MANIPULATORS
  Complex& rand_src(int i) const ;

  // operator functions
  QPropWRand& operator=(const QPropWRand& rhs);

    

  void ShiftPropForward(int n);
  void ShiftPropBackward(int n);

  //Free prop
  void AllocateRsrc() ; // Allocates memory for rsrc
  void DeleteRsrc();   // Deallocates memory for rsrc
  void RestoreQProp(char*, int mid);
  void SaveQProp(char*, int mid);
  
  // DESTRUCTOR
  virtual ~QPropWRand();

};


class QPropWRandWallSrc : public QPropWRand
{
  
public:
  QPropWRandWallSrc(Lattice& lat, CommonArg* c_arg) ;
  QPropWRandWallSrc(Lattice& lat, QPropWArg* arg, CommonArg* c_arg);
  
  void SetSource(FermionVectorTp& src, int spin, int color);
  SourceType SrcType(){return RANDWALL; }

} ;

class QPropWRandVolSrc : public QPropWRand
{

public:
  QPropWRandVolSrc(Lattice& lat, CommonArg* c_arg) ;
  QPropWRandVolSrc(Lattice& lat, QPropWArg* arg, CommonArg* c_arg);
   
  void SetSource(FermionVectorTp& src, int spin, int color);
  SourceType SrcType(){return RANDVOLUME ; }
} ;

class QPropWRandSlabSrc : public QPropWRand
{

public:
  QPropWRandSlabSrc(Lattice& lat, CommonArg* c_arg) ;
  QPropWRandSlabSrc(Lattice& lat, QPropWArg* arg, CommonArg* c_arg);
  QPropWRandSlabSrc(Lattice& lat, QPropWArg* arg, Float* src, CommonArg* );
  
  void SetSource(FermionVectorTp& src, int spin, int color);
  SourceType SrcType(){return RANDSLAB; }

} ;

class QPropWSeq : public QPropW
{

protected:
  QPropW&  quark;
  ThreeMom mom ;
  Float quark_mass ;

public:
  
  // CONSTRUCTORS
  
  QPropWSeq(Lattice& lat, QPropW& q,  int *p, QPropWArg*, CommonArg*);
  
  //QPropWSeq(const QPropWSeq& rhs);  //copy constructor not needed
  
  ThreeMom Mom(){return mom;}

  /*!
    Returns the quark mass of the source propagator  
   */
  Float SourceMass(){ return quark_mass ;}
} ;

class QPropWSeqMesSrc : public QPropWSeq
{
  /*!
    Only does a single gamma matrix insersion at the sink
    This is only used for pseudoscalars and vectors only anyway.
    should be enough for the Nucleon Decay project
  */
  int gamma ;
  
public:
  
  // CONSTRUCTORS
  
  QPropWSeqMesSrc(Lattice& lat, QPropW& quark,  int *p, 
		  int g, QPropWArg*, CommonArg*);
  
  void SetSource(FermionVectorTp& src, int spin, int color);
  SourceType SrcType(){return MESSEQ; }
} ;


class QPropWSeqBar : public QPropWSeq
{

protected:
  
  ProjectType P ;
  
public:
  
  // CONSTRUCTORS
  
  QPropWSeqBar(Lattice& lat, QPropW& quark,  int *p, 
	       ProjectType pp, QPropWArg*, CommonArg*);

  ProjectType Projection(){return P;}

} ;

class QPropWSeqProtUSrc : public QPropWSeqBar
{
  
public:
  
  // CONSTRUCTORS

  QPropWSeqProtUSrc(Lattice& lat, QPropW& quark,  int *p, 
		ProjectType pp, QPropWArg*, CommonArg*);

  //QPropWSeq(const QPropWSeq& rhs);  //copy constructor not needed

  void SetSource(FermionVectorTp& src, int spin, int color);
  SourceType SrcType(){return PROT_U_SEQ;}
} ;

class QPropWSeqProtDSrc : public QPropWSeqBar
{

public:
  
  // CONSTRUCTORS

  QPropWSeqProtDSrc(Lattice& lat, QPropW& quark,  int *p, 
		ProjectType pp, QPropWArg*, CommonArg*);

  //QPropWSeq(const QPropWSeq& rhs);  //copy constructor not needed
  
  void SetSource(FermionVectorTp& src, int spin, int color);
  SourceType SrcType(){return PROT_D_SEQ;}
} ;
CPS_END_NAMESPACE

#endif
