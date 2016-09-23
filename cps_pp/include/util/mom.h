#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Declaration of Mom class.

  $Id: mom.h,v 1.3 2004/08/18 11:57:37 zs Exp $
*/
//------------------------------------------------------------------
//
// mom.h
//
// Header file for the Mom class.
//
// Mom calculates the phase factor for each
// lattice site given a number of momenta and the source parameters
//
//------------------------------------------------------------------


#ifndef INCLUDED_MOM_H
#define INCLUDED_MOM_H

CPS_END_NAMESPACE
#include <math.h>    // for cos and sin
#include <util/data_types.h>
CPS_START_NAMESPACE

//! Structure for arguments to the ::Mom class constructor.

struct MomArg {
    int no_of_momenta; /*!< The number of different momenta.
			 Why do we need to specify this? */
    int max_p1;        /*!< The maximum momentum in lattice units in the 1st direction */   
  int max_p2;          /*!< The maximum momentum in lattice units in the 1st direction */   
  int max_p3;          /*!< The maximum momentum in lattice units in the 1st direction */   
    int deg;           /*!< Whether to average over degenerate momenta. */
    int dir;           /*!< The propagation direction.*/
    int src_begin[4];  /*!< Source coordinates */ 
  int src_end[4];      /*!< Source coordinates */

};


//! A class to compute momentum phases.
/*!
  The phases \e exp(ip.x) are computed where \e x is the lattice site
  coordinates relative to some source point and \e p is a momentum.
  The momentum is zero in a chosen direction of propagation.
  The values of the other momentum components \f$p_\mu\f$
  range from 0 to a chosen maximum in lattice units \e i.e. in units of
  \f$ 2\pi/L_\mu \f$ where \e L is the size of the lattice in the \f$\mu\f$
  direction. If desired, momenta of the same magnitude can be averaged
  (this should only be done if the lattice sizes are equal).

  Each different momentum used is assigned a number. These numbers and the
  corresponding momenutm components are written to a file called
  \c mom_table.log.
    
 */
class Mom
{
 private:
    char *cname;
    MomArg *mom_arg; // argument structure for momentum states

    Float PI;
    int dir,i,j,k;      // propagation direction and 3 orthogonal
    int no_of_mom;      // number of momenta to be calculated
    int deg;            // calculate degenerate momenta separately/together
    int max_p1,max_p2,max_p3; // maximal momentum in orthogonal directions

    int nx[4];                // local lattice extent
    int glb_L[4];             // global lattice extent
    int glb_sour_center[4];   // global source location
    Complex *mom_fact;

 public:
    Mom();
    /*!
      \todo Why is the destructor virtual? Is this class used polymorphically?
    */
    virtual ~Mom();
    //! Initialisation of parameters.
    void Init(MomArg *arg);
    //! Calculation of the phases
    void run(void);
        
    // returns the complex phase factor for momentum "imom" at site "s"
    //! Get a particular phase
    Complex fact(int imom, int *s);
};

//------------------------------------------------------------------
// External declarations.
//------------------------------------------------------------------
extern Mom MOM;


#endif





CPS_END_NAMESPACE
