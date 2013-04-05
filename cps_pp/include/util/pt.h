#include<config.h>
#ifdef PARALLEL
#include<comms/sysfunc_cps.h>
#endif
#if TARGET == QCDOC 
#include<util/pt_int.h>
#endif
#ifdef USE_QMP
#include<util/pt_int.h>
#endif
CPS_START_NAMESPACE
/*!\file
  \brief  Definition of the parallel transport classes.

  $Id: pt.h,v 1.24 2013-04-05 17:46:30 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2013-04-05 17:46:30 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/pt.h,v 1.24 2013-04-05 17:46:30 chulwoo Exp $
//  $Id: pt.h,v 1.24 2013-04-05 17:46:30 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: pt.h,v $
//  $Revision: 1.24 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/pt.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------


#ifndef INCLUDED_PT_H
#define INCLUDED_PT_H

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/vector.h>
//#include <comms/scu.h>
CPS_START_NAMESPACE

void pt_init(Lattice &lat);  //!< Initialization for parallel transporters
void pt_init_g();
void pt_delete();
void pt_delete_g();
void pt_mat(int n, Float **mout, Float **min, const int *dir);
void pt_1vec(int n, Float **vout, Float **vin, int const *dir);
void pt_2vec(int n, Float **vout, Float **vin, const int *dir);
//void pt_set_hop_pointer();
int pt_offset(int dir, int hop);
void pt_vvpd(Float **vect, int n_vect, const int *dir,
	     int n_dir, int hop, Float **sum);
void pt_vvpd(Float **vect2, Float ***vect, int n_vect, const int *dir,
	     int n_dir, int hop, Float **sum, int overwrite);
void pt_shift_field(Float **v, const int *dir, int n_dir,
		    int hop, Float **u);
void pt_shift_field_vec(Float **v, const int *dir, int n_dir,
		    int hop, Float **u);
void pt_shift_link(Float **u, const int *dir, int n_dir);

//---------------------------------------------------------------
//Checkerboarding methods
void pt_mat_cb(int n, Float **mout, Float **min, const int *dir, ChkbType cb);  //!<Parallel transport for checkerboarded Matrix fields
void pt_mat_cb(int n, Float **mout, Float **min, const int *dir, ChkbType cb, Float * new_gauge_field);  //!<Parallel transport for checkerboarded Matrix fields
void pt_mat_norm(int n, Float **mout, Float **min, const int *dir, ChkbType cb, Float *gauge);

void pt_1vec_cb(int n, Float **vout, Float **vin, const int *dir, ChkbType cb); //!<Parallel transport for checkerboarded Vector fields
void pt_1vec_cb(int n, Float **vout, Float **vin, const int *dir, ChkbType cb, Float * new_gauge_field); //!<Parallel transport for checkerboarded Vector fields
void pt_1vec_cb(int n, Float *vout, Float **vin, const int *dir, ChkbType cb, int pad); //!<Parallel transport for padded checkerboarded Vector fields
void pt_1vec_cb(int n, Float *vout, Float **vin, const int *dir, ChkbType cb, int pad, Float * new_gauge_field); //!<Parallel transport for padded checkerboarded Vector fields
void pt_1vec_cb_norm(int n, Float **vout, Float **vin, const int *dir,ChkbType cb, Float * gauge);
void pt_1vec_cb_pad(int n, Float *vout, Float **vin, const int *dir,ChkbType cb,int pad, Float * gauge);
//---------------------------------------------------------------


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
      VECT_LEN=6,          //!< Number of Floats in  a Vector
      MATRIX_SIZE=18,  //!< Number of Floats in  a Matrix
      SITE_LEN=72      //!< Number of Floats in four Matrix's
    };

};

#if 0
struct gauge_agg{
  int src;
  int dest;
  IFloat mat[18];
};

//-----------------------------------------------------------------------
//Checkerboarded parallel transport

struct gauge_agg_cb{
  //Index for initial position of field
  int src;
  //Index for "transported field"
  int dest;
  //Index for the padded field.  This assumes that the field will be
  //transported in all 8 possible directions.  Instead of storing
  //each direction in a single block, fields transported to each
  //lattice site from different directions are stored together
  //This allows for faster linear combination in the p4 action.
  int dest2;
  //Index for the gauge link
  int gauge_index;
  //Determines if the gauge link needs to be complex conjugated
  int dagger;
};

//-----------------------------------------------------------------------

struct hop_pointer {
  int src;
  int dest;
};
#endif

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
      \pre The Vector fields should be in ::CANONICAL order.
     */
    void run(int n, Vector **vout, Vector **vin, const int *dir)
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
    void run(int n, Matrix **mout, Matrix **min, const int *dir)
    { pt_mat(n,(IFloat **)mout, (IFloat **)min, dir);}

    //! Not implemented
    void run(Vector *vout, Vector *vin, int dir );

    /*! Computes sum[x] = vect[x] vect[x + hop dir]^dagger
      where the sum is over n_vect vectors and the hop is in a forward direction.
    */
    void vvpd(Vector **vect, int n_vect,
	      const int *dir, int n_dir, int hop, Matrix **sum){
	pt_vvpd((IFloat **)vect,n_vect,dir,n_dir,hop,(IFloat **)sum);
     }

    /*! Computes sum[x] = vect2[x] vect[x + hop dir]^dagger
      where the sum is over n_vect vectors and the hop is in a forward direction.
    */
    void vvpd(Vector **vect2, Vector ***vect, int n_vect,
	      const int *dir, int n_dir, int hop, Matrix **sum, int overwrite){
	pt_vvpd((IFloat **)vect2, (IFloat ***)vect,n_vect,dir,n_dir,hop,
		(IFloat **)sum, overwrite);
     }

    //! u[x] = v[x+dir] for n_dir forward or backward directions dir.
    void shift_field(Matrix **v, const int *dir, int n_dir,
		     int hop, Matrix **u){
        pt_shift_field((IFloat **)v,dir,n_dir,hop,(IFloat **)u);
    }

    //! u[x] = v[x+dir] for n_dir forward or backward directions dir.
    void shift_field_vec(Vector **v, const int *dir, int n_dir,
		     int hop, Vector **u){
        pt_shift_field_vec((IFloat **)v,dir,n_dir,hop,(IFloat **)u);
    }
    

    //! u[-/+nu](x) = U_[-/+nu](x)
    void shift_link(Matrix **u, const int *dir, int n_dir){
        pt_shift_link((IFloat **)u,dir,n_dir);
    }

    /*
    void shift_field_m1m2(Matrix **v, const int *dir, int n_dir,
			  int hop, Matrix **u, Matrix **x, Matrix **y){
      pt_shift_field_m1m2(v,dir,n_dir,hop,u,x,y);
    }
    

    void m1m2(Matrix **v, Matrix **x, const int *dir, int n_dir,
		     int hop, Matrix **u){
      pt_m1m2(v,x,dir,n_dir,hop,u);
    }
    */

    ~ParTransAsqtad();

};

//------------------------------------------------------------------
//! A class describing the Parallel Transport operator for staggered fermions.
/*!
  These are operations of the form
  \f$ u(x) = U_\mu(x) v(x+\mu) \f$
  where \e u is a field defined on sites of a given
  parity and \e v is a field defined on sites of
  the opposite parity.

  The fermions fields must be in staggered  ::STAG order

  The gauge field must be in staggered ::STAG order.
*/
//------------------------------------------------------------------
class ParTransStaggered_cb : public ParTransStagTypes
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
    ParTransStaggered_cb(Lattice& latt);            // Lattice object.

    //! Parallel transports Vector fields defined on single parity sites
    /*!
      \param n The number of fields.
      \param vout Array of the transported fields.
      \param vin Array of the fields to be transported.
      \param dir Array of directions in which the transports are performed.
      \param cb The parity on which the field vin is defined
      \param pad Parameter that determines whether the transported fields should be padded for QCDOC optimization
      \pre The gauge field should be in staggered order (::STAG).
      \pre The Vector fields should be in ::STAG order.
     */
    void run(int n, Vector **vout, Vector **vin, const int *dir,ChkbType cb)
      { pt_1vec_cb(n,(IFloat **)vout, (IFloat **)vin, dir, cb);}

    void run(int n, Vector **vout, Vector **vin, const int *dir,ChkbType cb, IFloat * new_gauge_field)
      { pt_1vec_cb(n,(IFloat **)vout, (IFloat **)vin, dir, cb, new_gauge_field);}

    void run(int n, IFloat * vout, Vector **vin, const int *dir, ChkbType cb, int pad)
      { pt_1vec_cb(n,vout, (IFloat **)vin, dir, cb, pad);}

    void run(int n, IFloat * vout, Vector **vin, const int *dir, ChkbType cb, int pad, IFloat * new_gauge_field)
      { pt_1vec_cb(n,vout, (IFloat **)vin, dir, cb, pad, new_gauge_field);}

    //! Parallel transports Matrix fields defined on single parity sites
    /*!
      \param n The number of fields.
      \param vout Array of the transported fields.
      \param vin Array of the fields to be transported.
      \param dir Array of directions in which the transports are performed.
      \param cb The parity on which the field min is defined
      \param pad Parameter that determines whether the transported fields should be padded for QCDOC optimization
      \pre The gauge field should be in staggered order (::STAG).
      \pre The Matrix fields should be in ::STAG order.      
     */
    void run(int n, Matrix **mout, Matrix **min, const int *dir,ChkbType cb)
    { pt_mat_cb(n,(IFloat **)mout, (IFloat **)min, dir, cb);}

    void run(int n, Matrix **mout, Matrix **min, const int *dir,ChkbType cb, IFloat * new_gauge_field)
    { pt_mat_cb(n,(IFloat **)mout, (IFloat **)min, dir, cb,new_gauge_field);}


    //! Not implemented
    void run(Vector *vout, Vector *vin, int dir );

    ~ParTransStaggered_cb();

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

    //! u[x] = v[x+dir] for n_dir forward or backward directions dir.
    void shift_field(Matrix **v, const int *dir, int n_dir,
		     int hop, Matrix **u){
        pt_shift_field((IFloat **)v,dir,n_dir,hop,(IFloat **)u);
    }


    //! u[-/+nu](x) = U_[-/+nu](x)
    void shift_link(Matrix **u, const int *dir, int n_dir){
        pt_shift_link((IFloat **)u,dir,n_dir);
    }

};


#endif

CPS_END_NAMESPACE
