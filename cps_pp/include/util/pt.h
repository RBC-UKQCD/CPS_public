#include<config.h>
#ifdef PARALLEL
#include<comms/sysfunc.h>
#endif
#include<util/pt_int.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Definition of the parallel transport classes.

  $Id: pt.h,v 1.11 2004-08-08 05:05:29 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-08-08 05:05:29 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/pt.h,v 1.11 2004-08-08 05:05:29 chulwoo Exp $
//  $Id: pt.h,v 1.11 2004-08-08 05:05:29 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: pt.h,v $
//  $Revision: 1.11 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/pt.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------


#ifndef INCLUDED_PT_H
#define INCLUDED_PT_H

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/vector.h>
#include <comms/scu.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
//! A class implementing parallel transports.
/*!
  These are operations of the form
  \f$ u(x) = U_\mu(x) v(x+\mu) \f$
  where \e u and \e v are fermionic vectors or 3x3 matrices which are
  defined as a  lattice field on all  lattice sites.

  The class is designed to perform a number of such parallel transports,
  each with an independent input field, output field and direction.
  
  This should used like an abstract base class.
 */
//------------------------------------------------------------------
class ParTrans
{
 private:
  char *cname;                     // Class name.


 protected:
  Lattice& lat;                    //!< The lattice..
  Matrix *gauge_field;             //!< Pointer to the gauge field.
  static int node_sites[5];
  static int bc[4];

 public:

  static int scope_lock;           // lock that forbids more than
                                   // one ParTrans object to be on
                                   // scope at any time.
  static int PTflops;              //! Count the flops 

  ParTrans(Lattice& latt);         

  virtual ~ParTrans();

  static void BondCond(Lattice& lat, Matrix *u_base);

};


//------------------------------------------------------------------

//! A class describing parallel transports for all sorts of staggered fermions.
/*!
  These are operations of the form
  \f$ u(x) = U_\mu(x) v(x+\mu) \f$
  where \e u and \e v are fermionic vectors or 3x3 matrices.

  This class just reimplements ParTrans; the derived class should be used.
*/
//------------------------------------------------------------------
class ParTransStagTypes : public ParTrans
{
 private:
  char *cname;    // Class name.

 public:
  ParTransStagTypes(Lattice& latt);            // Lattice object.

  virtual ~ParTransStagTypes();

 protected:
  enum
    {
      VECT_LEN=6,          //!< Number of floats in  a Vector
      MATRIX_SIZE=18,  //!< Number of floats in  a Matrix
      SITE_LEN=72      //!< Number of floats in four Matrix's
    };

};

struct gauge_agg{
  int src;
  int dest;
  IFloat mat[18];
};

struct hop_pointer {
  int src;
  int dest;
};

//------------------------------------------------------------------
//! A class describing the Parallel Transport operator for staggered fermions.
/*!
  These are operations of the form
  \f$ u(x) = U_\mu(x) v(x+\mu) \f$
  where \e u and \e v are fermionic vectors or 3x3 matrices.

  The gauge field must be in staggered ::STAG order.
*/
//------------------------------------------------------------------
class ParTransAsqtad : public ParTransStagTypes
{
  private:

    char *cname;         // Class name.

    int f_size_cb;       //The node checkerbrd. size of the ferm. field

    Vector *frm_tmp;     // Temporary fermion field

    CnvFrmType converted;

    int Offset(int dir, int hop);
    void setHopPointer();

  public:

    /*!
      \param latt The lattice containing the gauge field on which this
      operation is defined.
    */
    ParTransAsqtad(Lattice& latt);            // Lattice object.

    //! Parallel transports fields of 3x3 matrices.
    /*!
      \param n The number of fields.
      \param vout Array of the transported fields.
      \param vin Array of the fields to be transported.
      \param dir Array of directions in which the transports are performed.
      \pre The gauge field should be in staggered order (::STAG).
      \pre The Vector fields should be in CANONICAL order.
     */
    void run(int n, Vector **vout, Vector **vin, const int *dir )
    { pt_1vec(n,(IFloat **)vout, (IFloat **)vin, dir);}
    //! Parallel transports fields of staggered fermionic vectors.
    /*!
      \param n The number of fields.
      \param vout Array of the transported fields.
      \param vin Array of the fields to be transported.
      \param dir Array of directions in which the transports are performed.
      \pre The gauge field should be in staggered order (::STAG).
      \pre The Matrix fields should be in CANONICAL order.      
     */
    void run(int n, Matrix **mout, Matrix **min, const int *dir )
    { pt_mat(n,(IFloat **)mout, (IFloat **)min, dir);}

    //! Not implemented
    void run(Vector *vout, Vector *vin, int dir );

    /*! Computes sum[x] = vect[x] vect[x + hop dir]^dagger
      where the sum is over n_vect vectors and the hop is in a forward direction.
    */
    void vvpd(Vector **vect, int n_vect,
	      const int *dir, int n_dir, int hop, Matrix **sum);

    //! u[x] = v[x+dir] for n_dir forward or backward directions dir.
    void shift_field(Matrix **v, const int *dir, int n_dir,
		     int hop, Matrix **u);

    //! u[-/+nu](x) = U_[-/+nu](x)
    void shift_link(Matrix **u, const int *dir, int n_dir);

    ~ParTransAsqtad();

};

//------------------------------------------------------------------

//! A class describing parallel transports for all sorts of staggered fermions.
/*!
  These are operations of the form
  \f$ u(x) = U_\mu(x) v(x+\mu) \f$
  where \e u and \e v are fermionic vectors or 3x3 matrices.

  This class just reimplements ParTrans; the derived class should be used.
*/
//------------------------------------------------------------------
class ParTransWilsonTypes : public ParTrans
{
 private:
  char *cname;    // Class name.

 public:
  ParTransWilsonTypes(Lattice& latt);            // Lattice object.

  virtual ~ParTransWilsonTypes();

};

//------------------------------------------------------------------

//! A class describing parallel transports for all sorts of staggered fermions.
/*!
  These are operations of the form
  \f$ u(x) = U_\mu(x) v(x+\mu) \f$
  where \e u and \e v are fermionic vectors or 3x3 matrices.

  This class just reimplements ParTrans; the derived class should be used.
*/
//------------------------------------------------------------------
class ParTransGauge : public ParTrans
{
 private:
  char *cname;    // Class name.

 public:
  ParTransGauge(Lattice& latt);            // Lattice object.

  virtual ~ParTransGauge();

  void run(int n, Matrix **mout, Matrix **min, const int *dir )
    { pt_mat(n,(IFloat **)mout, (IFloat **)min, dir);}

};


#endif

CPS_END_NAMESPACE
