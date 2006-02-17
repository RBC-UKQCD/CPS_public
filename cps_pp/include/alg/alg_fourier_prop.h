//------------------------------------------------------------------
//
// alg_fourier_prop.h
//
// Header file for all alg classes relevant to Wilson-type fermion
// three point functions. The type of glue or fermion is given as
// an argument of type Lattice& to the constructor. If the 
// fermion type is not F_CLASS_WILSON or F_CLASS_CLOVER or 
// F_CLASS_DWF the constructors exit with a general error.
//
//------------------------------------------------------------------

#ifndef INCLUDED_ALG_FPROP_H
#define INCLUDED_ALG_FPROP_H

#include <config.h>
#include <util/lattice.h>
#include <util/smalloc.h>
#include <util/pmalloc.h>
#include <alg/alg_base.h>
#include <alg/common_arg.h>
#include <alg/fourierprop_arg.h>
#include <alg/wilson_matrix.h>
#include <alg/qpropw.h>
#include <alg/qpropwGFL.h>
#include <stdio.h>

CPS_START_NAMESPACE

class AlgFourierProp : public Alg
{
private:

  char* cname;

  //! Node checkerboard size of the fermion field
  int f_size;

protected:
  
  //! store the boundary conditions
  int bc[4];

  /*!
    The argument structure for the
    fourier transform of the quark propagator
  */
  
  FourierPropArg* alg_FourierProp_arg;
  
  /*!
    propagator for one momenta
  */
  WilsonMatrix momprop;

  /*!
    pointer to current configuration space
    propagagtor
  */
  QPropW *prop;

public:
  
  AlgFourierProp( Lattice&         latt, 
                  CommonArg*      c_arg, 
                  FourierPropArg*   arg  );
  
  virtual ~AlgFourierProp();
  
  virtual void run();

  /*!
    fourier transform for all chosen
    momenta and output to file
  */
  virtual void ft(FILE*,MomentaList&);

  /*!
    Calculate a single momentum ft
  */
  void calcmom(const FourMom&);

  //! print out file header
  void print_header( FILE* , const char* title );
};


/*!
  Class for calculating the needed ft's for 
  all the dS=1 contractions

*/

class AlgFourierPropDis:public AlgFourierProp
{
private:

  char* cname;

  /*! do the disconnected */
  bool discon;

public:
  
  AlgFourierPropDis(Lattice & latt, CommonArg* c_arg, FourierPropArg* arg):
    AlgFourierProp(latt, c_arg, arg),
      discon(true)
  {;}
  
  virtual ~AlgFourierPropDis()
  {;}
  

  void set_discon( bool val ) { discon=val; } 

  void run();
  
};

CPS_END_NAMESPACE

#endif
