#include <util/lattice/eigcg_controller.h>
#include <alg/eigcg_arg.h>
#include <util/verbose.h>
#include <util/gjp.h>
#include <util/qcdio.h>
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

void EigCG::printH(const string &fn)const
{
    if(!created) return;

    FILE *fp = Fopen(fn.c_str(), "w");

    if(is_float) {
        EigCGController<float> *eigcg = EigCGController<float>::getInstance();
        const int size = eigcg->def_len;

        VRB.Result("EigCG", "printH", "def_len = %d\n", eigcg->def_len);

        for(int i = 0; i < size; ++i) {
            for(int j = 0; j < size; ++j) {
                Fprintf(fp, "%17.10e %17.10e ",
                        real(eigcg->H[i * size + j]),
                        imag(eigcg->H[i * size + j]));
            }
            Fprintf(fp, "\n");
        }
    } else {
        EigCGController<double> *eigcg = EigCGController<double>::getInstance();
        const int size = eigcg->def_len;

        VRB.Result("EigCG", "printH", "def_len = %d\n", eigcg->def_len);

        for(int i = 0; i < size; ++i) {
            for(int j = 0; j < size; ++j) {
                Fprintf(fp, "%17.10e %17.10e ",
                        real(eigcg->H[i * size + j]),
                        imag(eigcg->H[i * size + j]));
            }
            Fprintf(fp, "\n");
        }
    }

    Fclose(fp);
}
