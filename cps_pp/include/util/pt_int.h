#include <config.h>
/*!\file
  \brief Declaration of functions used by the parallel transport classes.

  $Id: pt_int.h,v 1.5 2004-12-07 05:23:18 chulwoo Exp $
  Why are (at least some of) these not class methods?
*/
#include <util/lattice.h>
CPS_START_NAMESPACE
void pt_init(Lattice &lat);  //!< Initialization for parallel transporters
void pt_init_g();
void pt_delete();
void pt_delete_g();
void pt_mat(int n, IFloat **mout, IFloat **min, const int *dir);
void pt_1vec(int n, IFloat **vout, IFloat **vin, const int *dir);
void pt_2vec(int n, IFloat **vout, IFloat **vin, const int *dir);
void pt_set_hop_pointer();
int pt_offset(int dir, int hop);
void pt_vvpd(Vector **vect, int n_vect, const int *dir,
	     int n_dir, int hop, Matrix **sum);
void pt_shift_field(Matrix **v, const int *dir, int n_dir,
		    int hop, Matrix **u);
void pt_shift_link(Matrix **u, const int *dir, int n_dir);

//---------------------------------------------------------------
//Checkerboarding methods
void pt_mat_cb(int n, IFloat **mout, IFloat **min, const int *dir, ChkbType cb);  //!<Parallel transport for checkerboarded Matrix fields
void pt_1vec_cb(int n, IFloat **vout, IFloat **vin, const int *dir, ChkbType cb); //!<Parallel transport for checkerboarded Vector fields
void pt_1vec_cb(int n, IFloat **vout, IFloat **vin, const int *dir, ChkbType cb, IFloat * new_gauge_field); //!<Parallel transport for checkerboarded Vector fields
void pt_1vec_cb(int n, IFloat *vout, IFloat **vin, const int *dir, ChkbType cb, int pad); //!<Parallel transport for padded checkerboarded Vector fields
void pt_1vec_cb(int n, IFloat *vout, IFloat **vin, const int *dir, ChkbType cb, int pad, IFloat * new_gauge_field); //!<Parallel transport for padded checkerboarded Vector fields
void pt_1vec_cb_norm(int n, IFloat **vout, IFloat **vin, const int *dir,ChkbType cb, IFloat * gauge);
void pt_1vec_cb_pad(int n, IFloat *vout, IFloat **vin, const int *dir,ChkbType cb,int pad, IFloat * gauge);
//---------------------------------------------------------------
CPS_END_NAMESPACE
