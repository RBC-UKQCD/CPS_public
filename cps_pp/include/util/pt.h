#include<config.h>
#ifdef PARALLEL
#include<comms/sysfunc.h>
#endif
CPS_START_NAMESPACE
/*!\file
  \brief  Definition of the Dirac operator classes: DiracOp, DiracOpStagTypes.

  $Id: pt.h,v 1.5 2004-04-27 03:51:16 cwj Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: cwj $
//  $Date: 2004-04-27 03:51:16 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/pt.h,v 1.5 2004-04-27 03:51:16 cwj Exp $
//  $Id: pt.h,v 1.5 2004-04-27 03:51:16 cwj Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: pt.h,v $
//  $Revision: 1.5 $
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
//! A class representing operations on the Dirac operator.
/*!
  This is an abstract base class, so the details specific to the various
  types of fermion action are defined in the derived classes.
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
  ParTrans(Lattice& latt);           // Lattice object.

  virtual ~ParTrans();

};


//------------------------------------------------------------------
//
// DiracOpStagTypes class.
//
//! A class describing the Dirac operator for all sorts of staggered fermions.
/*!
  This is an abstract base class from which the staggered fermion Dirac
  operator classes are derived.

  The staggered fermion is
  \f[
  M_{xy} =
  m_0 - \sum_mu e^{\sum_{i=0}^{\mu-1}x_i} 
  [ U^\dagger_\mu(x) \delta_{x\,y+\mu} - U_\mu(x-\mu) \delta_{x\,y-\mu} ]
  \f]

  \e N.B  The phases are implemented in the gauge field when it is
  converted to staggered (::STAG) storage order.
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

typedef struct gaguge_agg{
  int src;
  int dest;
  IFloat mat[18];
} gauge_agg;

//------------------------------------------------------------------
//
// DiracOpStag is derived from DiracOpStagTypes and is the front
// end for all Dirac operators associated with Staggered fermions.

//! A class describing the Parallel Transport operator for staggered fermions.
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
 public:
  ParTransAsqtad(Lattice& latt);            // Lattice object.
  void run(int n, Vector **vout, Vector **vin, const int *dir );
  void run(int n, Matrix **mout, Matrix **min, const int *dir );
  void run(Vector *vout, Vector *vin, const int dir );

  ~ParTransAsqtad();

};

#endif

CPS_END_NAMESPACE
