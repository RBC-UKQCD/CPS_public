#include <util/qcdio.h>
#include <alg/force_arg.h>
#include <cstdlib>
#include <cstdio>

CPS_START_NAMESPACE
using namespace std;

void ForceArg::print(Float dt, char *label)const {
    const char *fname = "print()";
    FILE *fp = Fopen("force.dat", "a");
    
    if(!fp) {
        ERR.FileA(cname, fname, "force.dat");
    }
    Fprintf(fp,"%s L1 = %e L2 = %e Linf = %e dt = %f\n",
            label, L1, L2, Linf, dt);
    Fclose(fp);
}

CPS_END_NAMESPACE
