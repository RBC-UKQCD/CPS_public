#ifndef INCLUDED_FPROP_ARG_H
#define INCLUDED_FPROP_ARG_H

#include <config.h>
#include <alg/cg_arg.h>
#include <alg/FourMom.h>

CPS_START_NAMESPACE

/*!  
  The structure type FourierPropArg holds the parameters specific to
  fourier transforming a point source quark propagator 
*/
class  FourierPropArg 
{
public:

  /*! 
    The conjugate gradient argument for 
    the quark propagator inversion 
  */
  CgArg cg;		      

  /*! 
    position of the src in global co-ordinates
  */
  int x_src, y_src, z_src, t_src;
  
  //! momenta list
  
  MomentaList plist;
  
 
  /*! 
    'special' momenta directions
    These are the fixed momenta
    used in the delta S=1 renormalisation
  */
  MomentaList smom;

  /*!
    ft momenta smom.size() arrray
    of MomentaLists
  */
  MomentaList* smom_comp; 
  
  /*!
    points to char* for containing the name of the 
    output file
  */
  const void *results;

  /*!
    used for the dissconnected calculation
  */
  const void *results2;
  const void *traces;
 

public:

  FourierPropArg():
    x_src(0),
    y_src(0),
    z_src(0),
    t_src(0)
  {;}
  
  
};

CPS_END_NAMESPACE

#endif /* !INCLUDED_FPROP_ARG_H */
