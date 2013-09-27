// -*- mode:c++; c-basic-offset:4 -*-
#ifndef PROP_CONTAINER_H_KL3
#define PROP_CONTAINER_H_KL3

#include <vector>
#include <string>
#include <util/gjp.h>
#include <util/verbose.h>
#include <alg/qpropw.h>
#include <omp.h>
#include <cassert>

#include "my_util.h"

// Propagators from all time slices (including the case where the
// source is on the mirrored lattice).
class AllProp {
public:
    enum PREC { SINGLE, DOUBLE };
    AllProp(PREC _p)
        :lcl_vol(cps::GJP.VolNodeSites()),
         t_size_glb(cps::GJP.Sites(3)),
         prec(_p)
    {
        if(prec == SINGLE) {
            ps.assign(t_size_glb, std::vector<cps::WilsonMatrixS>());
            as.assign(t_size_glb, std::vector<cps::WilsonMatrixS>());
        } else {
            pd.assign(t_size_glb, std::vector<cps::WilsonMatrix>());
            ad.assign(t_size_glb, std::vector<cps::WilsonMatrix>());
        }
    }

    // return a WilsonMatrix according to the type of propagators
    //
    // t: source location (time slice), for P or A propagators this
    // can be any value in [0, t_size_glb), for P+A/P-A propagators
    // the size is doubled: t in [0, 2 * t_size_glb).
    //
    // i: sink location (4D corrdinate), for P or A propagators this
    // can be anything in [0, lcl_vol), for P+A/P-A propagators
    // this can be anything in [0, 2 * lcl_vol).
    const cps::WilsonMatrix operator()(size_t i, size_t t, PROP_TYPE ptype)const {
        switch(ptype) {

        case PROP_P:
            assert(t < t_size_glb);
            assert(i < lcl_vol);
            return prec == SINGLE
                ? cps::WilsonMatrix(ps[t][i])
                : pd[t][i];

        case PROP_A:
            assert(t < t_size_glb);
            assert(i < lcl_vol);
            return prec == SINGLE
                ? cps::WilsonMatrix(as[t][i])
                : ad[t][i];

        case PROP_PA:
            {
                assert(t < 2 * t_size_glb);
                assert(i < 2 * lcl_vol);

                bool add = t < t_size_glb == i < lcl_vol;
                if(t >= t_size_glb) t -= t_size_glb;
                if(i >= lcl_vol) i -= lcl_vol;

                cps::WilsonMatrix pi, ai;
                if(prec == SINGLE) {
                    pi = ps[t][i];
                    ai = as[t][i];
                } else {
                    pi = pd[t][i];
                    ai = ad[t][i];
                }

                return 0.5 * (add ? pi + ai : pi - ai);
            }

        default:
            assert(false);
        }
    }

    // Test if a certain type of propagator is NOT calculated on a
    // given time slice.
    //
    // P+A/P-A propagator on a given time slice requires both periodic
    // and antiperiodic propagators.
    bool empty(size_t t, PROP_TYPE ptype)const {
        switch(ptype) {
        case PROP_P:

            assert(t < t_size_glb);
            return prec == SINGLE
                ? ps[t].empty()
                : pd[t].empty();

        case PROP_A:
            assert(t < t_size_glb);
            return prec == SINGLE
                ? as[t].empty()
                : ad[t].empty();

        case PROP_PA:
            assert(t < 2 * t_size_glb);
            if(t >= t_size_glb) t -= t_size_glb;

            if(prec == SINGLE) {
                return ps[t].empty() || as[t].empty();
            } else {
                return pd[t].empty() || ad[t].empty();
            }

        default:
            assert(false);
        }
    }

    void clearP(void) {
        for(unsigned i = 0; i < ps.size(); ++i) {
            ps[i].clear();
        }
        for(unsigned i = 0; i < pd.size(); ++i) {
            pd[i].clear();
        }
    }
    void clearA(void) {
        for(unsigned i = 0; i < as.size(); ++i) {
            as[i].clear();
        }
        for(unsigned i = 0; i < ad.size(); ++i) {
            ad[i].clear();
        }
    }
    void clear(void) {
        clearP();
        clearA();
    }

    // Add a propagator where the source is located at time slice t.
    // If periodic == true then it will be treated as a P-boundary
    // condition propagator, otherwise it will be treated as an
    // A-boundary condition propagator.
    //
    // Potentially transforms the propagator to single precision to
    // save some memory.
    void add(cps::QPropW &qp, size_t t, bool periodic) {
        if(prec == SINGLE) {
            std::vector<cps::WilsonMatrixS> &wm = periodic ? ps[t] : as[t];
            assert(wm.empty());

            wm.resize(lcl_vol);
#pragma omp parallel for
            for(size_t i = 0; i < lcl_vol; ++i) {
                wm[i] = qp[i];
            }
        } else {
            std::vector<cps::WilsonMatrix> &wm = periodic ? pd[t] : ad[t];
            assert(wm.empty());

            wm.resize(lcl_vol);
#pragma omp parallel for
            for(size_t i = 0; i < lcl_vol; ++i) {
                wm[i] = qp[i];
            }
        }
    }

    // store all propagators I have.
    //
    //! IMPORTANT: This function assumes gauge fixed wall source!!!
    void store_all(const std::string &fn_stem, double mass, int traj)const;

private:
    void store(const std::string &fn,
               const std::vector<cps::WilsonMatrix> &prop,
               int t)const;
    void store_qio(const std::string &fn,
               const std::vector<cps::WilsonMatrix> &prop,
               int t)const;
private:
    const size_t lcl_vol;
    const size_t t_size_glb;

    const PREC prec;
    std::vector<std::vector<cps::WilsonMatrixS> > ps; // P prop (single)
    std::vector<std::vector<cps::WilsonMatrixS> > as; // A prop (single)
    std::vector<std::vector<cps::WilsonMatrix> >  pd; // P prop (double)
    std::vector<std::vector<cps::WilsonMatrix> >  ad; // A prop (double)
};

#endif
