#include <config.h>
/*!\file
  \brief Declaration of functions used by the parallel transport classes.

  $Id: pt_int.h,v 1.4 2004-09-02 16:59:06 zs Exp $
  Why are (at least some of) these not class methods?
*/
#include <util/lattice.h>
CPS_START_NAMESPACE
    void pt_init(Lattice &lat);
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
CPS_END_NAMESPACE
