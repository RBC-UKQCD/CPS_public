#include <util/lattice/eigcg_controller.h>
#include <alg/eigcg_arg.h>
#include <util/gjp.h>
#include <vector>

#include "eigcg.h"

using namespace std;
USING_NAMESPACE_CPS

EigCG::EigCG(EigCGArg *eigcg_arg, bool use_float)
    :created(eigcg_arg != NULL),
     is_float(use_float)
{
    if(eigcg_arg == NULL) return;
    
    //set eigcg instance
    vector<double> restart(eigcg_arg->restart,
                           eigcg_arg->restart + eigcg_arg->restart_len);

    int vec_len = GJP.VolNodeSites() * GJP.SnodeSites() * 12;
    
    if(use_float) {
        EigCGController<float>::setInstance(eigcg_arg->nev, 
                                            eigcg_arg->m,
                                            eigcg_arg->max_def_len,
                                            eigcg_arg->max_eig_cut,
                                            restart,
                                            eigcg_arg->always_restart,
                                            vec_len);
    } else {
        EigCGController<double>::setInstance(eigcg_arg->nev, 
                                             eigcg_arg->m,
                                             eigcg_arg->max_def_len,
                                             eigcg_arg->max_eig_cut,
                                             restart,
                                             eigcg_arg->always_restart,
                                             vec_len);
    }
}

EigCG::~EigCG()
{
    if(created) {
        if(is_float) {
            EigCGController<float> *eigcg = EigCGController<float>::getInstance();
            eigcg->free();
        } else {
            EigCGController<double> *eigcg = EigCGController<double>::getInstance();
            eigcg->free();
        }
    }
}
