// -*- mode:c++;c-basic-offset:4 -*-
#ifndef INCLUDED_FFORCE_WILSON_TYPE_H__
#define INCLUDED_FFORCE_WILSON_TYPE_H__

#include <config.h>
#include <util/vector.h>
#include <alg/force_arg.h>

CPS_START_NAMESPACE

// FforceWilsonType: computing fermion force from the ''v1/v2''
// vectors. See e.g. Fbfm::EvolveMomFforceBase() for its usage.
//
// 1. Make sure the gauge field has the boundry condition turned
// on. This class doesn't in any way deal with boundry settings.
//
// 2. vec1 and vec2 are 2 auxiliary vectors.  They must be stored in
// (color, spin, s, x, y, z, t) order.
class FforceWilsonType {
public:
    FforceWilsonType(Matrix *momentum, Matrix *gauge_field,
                     Float *vec1, Float *vec2, int Ls, Float dt);
    ~FforceWilsonType();
    ForceArg run();

private:
    void collect_surface(int mu);
    void comm(int mu);
    // these 2 functions are called from within an OpenMP construction
    ForceArg do_internal(int mu, int nthreads);
    ForceArg do_surface(int mu, int nthreads);

private:
    ForceArg f_arg;
    long lcl[5]; // local size
    Matrix *mom;
    Matrix *gauge;
    Float *v1;
    Float *v2;

    // stuff for communication
    Float *sndbuf;
    Float *rcvbuf;

    Float coef;

    // All of the following are computed as number of Floats
    long surf_size[4];
    long v1so[4]; // v1 surface offset in sndbuf/rcvbuf
    long v2so[4]; // v2 surface offset in sndbuf/rcvbuf
    long bufsize; // total size of sndbuf/rcvbuf
};

CPS_END_NAMESPACE

#endif
