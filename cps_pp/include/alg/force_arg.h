/*
  Class used to return information regarding the MD force.
  L1 - The L1 norm of the force
  L2 - The L2 norm of the force
  Linf - The Linfinity norm of the force

  Placed print routine in this class to minimise rewriting if code evolves.
*/

#ifndef FORCE_ARG_H
#define FORCE_ARG_H

#include <config.h>
#include <alg/enum.h>
#include <cmath>

CPS_START_NAMESPACE

class ForceArg {
public:
    const char *cname;
    Float L1;
    Float L2;
    Float Linf;

public:
    ForceArg()
        :cname("ForceArg"),L1(0.),L2(0.),Linf(0.) {}

    ForceArg(Float _L1, Float _L2, Float _Linf)
        :cname("ForceArg"),L1(fabs(_L1)),L2(fabs(_L2)),Linf(fabs(_Linf)) {}

    ~ForceArg() {}

    void combine(const ForceArg &a) {
        L1 += a.L1;
        L2 += a.L2;
        Linf = Linf > a.Linf ? Linf : a.Linf;
    }

    void glb_reduce();

    // measure the L1/L2/Linf norm, mom is an array of momenta defined
    // on each link.
    void measure(const Matrix mom[]);
    void print(Float dt, char *label)const;
};

CPS_END_NAMESPACE

#endif /* !FORCE_ARG_H */
