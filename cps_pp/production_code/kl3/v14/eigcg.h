// -*- mode:c++; c-basic-offset:4 -*-
#ifndef EIGCG_H_KL3__
#define EIGCG_H_KL3__

#include <alg/eigcg_arg.h>

class EigCG {
public:
    EigCG(cps::EigCGArg *eigcg_arg = NULL, bool use_float = false);
    ~EigCG();
private:
    bool created;
    bool is_float;
};

#endif
