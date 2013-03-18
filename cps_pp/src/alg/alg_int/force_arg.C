#include <util/qcdio.h>
#include <alg/force_arg.h>
#include <comms/glb.h>
#include <util/gjp.h>

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include <omp.h>

CPS_START_NAMESPACE
using namespace std;

void ForceArg::print(Float dt, char *label)const
{
    const char *fname = "print()";
    FILE *fp = Fopen("force.dat", "a");
    
    if(!fp) {
        ERR.FileA(cname, fname, "force.dat");
    }
    Fprintf(fp,"%s L1 = %e L2 = %e Linf = %e dt = %f\n",
            label, L1, L2, Linf, dt);
    Fclose(fp);
}

void ForceArg::glb_reduce()
{
    glb_sum(&L1);
    glb_sum(&L2);
    glb_max(&Linf);

    double links = 4.0 * GJP.VolSites();
    L1 /= links;
    L2 = std::sqrt(L2 / links);
}

// measure the L1/L2/Linf norm, mom is an array of momenta defined
// on each link.
void ForceArg::measure(const Matrix mom[])
{
    L1 = 0;
    L2 = 0;
    Linf = 0;

    int links = GJP.VolNodeSites() * 4;

    // The result of the following will be slightly different when
    // running multiple times on the same data, due to OpenMP.
#pragma omp parallel
    {
        ForceArg f_arg; // threaded

#pragma omp for nowait
        for (int i = 0; i < links; ++i) {
            Float a2 = mom[i].norm();
            Float a = sqrt(a2);

            f_arg.L1 += a;
            f_arg.L2 += a2;
            f_arg.Linf = a > f_arg.Linf ? a : f_arg.Linf;
        }

#pragma omp critical
        {
            this->combine(f_arg);
        }
    } // omp
}

CPS_END_NAMESPACE
