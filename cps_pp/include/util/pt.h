#include<config.h>
#ifdef PARALLEL
#include<sysfunc.h>
#endif
CPS_START_NAMESPACE
/*!\file
  \brief  Definition of the Dirac operator classes: DiracOp, DiracOpStagTypes.

  $Id: pt.h,v 1.3 2004-01-14 07:42:59 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-01-14 07:42:59 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/pt.h,v 1.3 2004-01-14 07:42:59 chulwoo Exp $
//  $Id: pt.h,v 1.3 2004-01-14 07:42:59 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.1.2.1.2.1  2003/12/15 18:52:21  cwj
//  *** empty log message ***
//
//  Revision 1.1.2.1  2003/11/05 16:13:21  mike
//  Initial attempt at producing working branch
//
//  Revision 1.2  2003/10/21 17:53:08  chulwoo
//  added asqtad_KS
//  changes for stagerred and asqtad action
//
//  Revision 1.1.1.1  2003/09/18 22:30:58  chulwoo
//  Mike's files for single node QCDOC + Parallel transport
//  I added some hacks for PARALLEL without MPI_SCU
//  PARALLEL=2 set PARALLEL without MPI_SCU
//
//
//  Revision 1.4  2003/08/29 21:02:56  mike
//  Removed MatMInv function as was unnecessary.
//
//  Revision 1.3  2003/08/29 20:28:55  mike
//  Added MInvCG, the multishift CG invertor, used by AlgHmcRHMC.
//
//  Revision 1.2  2003/07/24 16:53:53  zs
//  Addition of documentation via doxygen:
//  doxygen-parsable comment blocks added to many source files;
//  New target in makefile and consequent alterations to configure.in;
//  New directories and files under the doc directory.
//
//  Revision 1.5  2001/08/16 12:54:30  anj
//  Some fixes follosin the float-> IFloat change, mostly of the (variable
//  anme) IFloat_p -> float_p type.  A few fixes to ensure the test
//  programs use the same level of verbosity throughout, and an update of
//  the regression.pl script to make it more useful. Anj
//
//  Revision 1.4  2001/08/16 10:50:29  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.2  2001/06/19 18:13:16  anj
//  Serious ANSIfication.  Plus, degenerate double64.h files removed.
//  Next version will contain the new nga/include/double64.h.  Also,
//  Makefile.gnutests has been modified to work properly, propagating the
//  choice of C++ compiler and flags all the way down the directory tree.
//  The mpi_scu code has been added under phys/nga, and partially
//  plumbed in.
//
//  Everything has newer dates, due to the way in which this first alteration was handled.
//
//  Anj.
//
//  Revision 1.2  2001/05/25 06:16:09  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: pt.h,v $
//  $Revision: 1.3 $
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

  static int scope_lock;           // lock that forbids more than
                                   // one DiracOp object to be on
                                   // scope at any time.

 protected:
  Lattice& lat;                    //!< The lattice..
  Matrix *gauge_field;             //!< Pointer to the gauge field.


 public:
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
