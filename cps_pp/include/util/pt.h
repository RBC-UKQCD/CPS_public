#include<config.h>
#ifdef PARALLEL
#include<comms/sysfunc.h>
#endif
CPS_START_NAMESPACE
/*!\file
  \brief  Definition of the parallel transport classes.

  $Id: pt.h,v 1.7 2004-05-12 17:23:57 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-05-12 17:23:57 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/pt.h,v 1.7 2004-05-12 17:23:57 zs Exp $
//  $Id: pt.h,v 1.7 2004-05-12 17:23:57 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: pt.h,v $
//  $Revision: 1.7 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/pt.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------


#ifndef INCLUDED_PT_H
#define INCLUDED_PT_H

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/vector.h>
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


 public:

  static int scope_lock;           // lock that forbids more than
                                   // one ParTrans object to be on
                                   // scope at any time.
  ParTrans(Lattice& latt);         

  virtual ~ParTrans();

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

};

struct gauge_agg{
  int src;
  int dest;
  IFloat mat[18];
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

    IFloat * rcv_buf[2*6];
    IFloat * tmp_buf[2*6];
    IFloat *gauge_field_addr;

    void pt_init(const void *);
    void pt_init_g();
    void pt_delete();
    void pt_delete_g();

    CnvFrmType converted;
    
  public:

    /*!
      \param latt The lattice containing the gauge field on which this
      operation is defined.
      \pre The gauge field should be in staggered order (::STAG).
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
    void run(int n, Vector **vout, Vector **vin, const int *dir );

    //! Parallel transports fields of staggered fermionic vectors.
    /*!
      \param n The number of fields.
      \param vout Array of the transported fields.
      \param vin Array of the fields to be transported.
      \param dir Array of directions in which the transports are performed.
      \pre The gauge field should be in staggered order (::STAG).
      \pre The Matrix fields should be in CANONICAL order.      
     */
    void run(int n, Matrix **mout, Matrix **min, const int *dir );

    //! Not implemented
    void run(Vector *vout, Vector *vin, int dir );

    ~ParTransAsqtad();

};

#endif

CPS_END_NAMESPACE
