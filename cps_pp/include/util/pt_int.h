#include <config.h>
#include <util/lattice.h>
CPS_START_NAMESPACE
    void pt_init(Lattice &lat);
    void pt_init_g();
    void pt_delete();
    void pt_delete_g();
    void pt_mat(int n, IFloat **mout, IFloat **min, const int *dir);
    void pt_1vec(int n, IFloat **vout, IFloat **vin, const int *dir);
    void pt_2vec(int n, IFloat **vout, IFloat **vin, const int *dir);
CPS_END_NAMESPACE
