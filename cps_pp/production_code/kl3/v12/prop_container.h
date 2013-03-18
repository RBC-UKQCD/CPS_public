// -*- mode:c++; c-basic-offset:4 -*-
#ifndef PROP_CONTAINER_H_KL3
#define PROP_CONTAINER_H_KL3

#include <util/gjp.h>
#include <util/verbose.h>
#include <alg/qpropw.h>
#include <omp.h>
#include <cassert>

#include "my_util.h"

class PropAP {
public:
    PropAP(const std::vector<cps::WilsonMatrix> &f,
           const std::vector<cps::WilsonMatrix> &s)
        :fst(&f),snd(&s),
         lcl_vol(cps::GJP.VolNodeSites()) {
    }
    const cps::WilsonMatrix &operator[](size_t i)const {
        if(i < lcl_vol) return (*fst)[i];
        else return (*snd)[i - lcl_vol];
    }
private:
    const std::vector<cps::WilsonMatrix> *fst;
    const std::vector<cps::WilsonMatrix> *snd;
    // If the source is on the local lattice, then the 1st half is
    // (P+A)/2 and the 2nd half is (P-A)/2. If the source is on the
    // mirrored lattice, then the 1st half is (P-A)/2 and the 2nd half
    // is (P+A)/2.
    size_t lcl_vol;
};

// Propagators from all time slices (including the case where the
// source is on the mirrored lattice).
class AllProp {
public:
    AllProp()
        :lcl_vol(cps::GJP.VolNodeSites()),
         t_size_glb(cps::GJP.TnodeSites() * cps::GJP.Tnodes()),
         p(t_size_glb),
         m(t_size_glb)
    {
        for(unsigned i = 0; i < t_size_glb; ++i) {
            prop.push_back(PropAP(p[i], m[i]));
        }
        for(unsigned i = 0; i < t_size_glb; ++i) {
            prop.push_back(PropAP(m[i], p[i]));
        }
    }

    const PropAP &operator[](size_t i)const {
        return prop[i];
    }

    bool empty(size_t t)const {
        assert(t < 2 * t_size_glb);
        if(t >= t_size_glb) t -= t_size_glb;
        return p[t].empty();
    }

    // Add a propagator where the source is located at time slice t.
    // If periodic == true then it will be treated as a P-boundary
    // condition propagator, otherwise it will be treated as an
    // A-boundary condition propagator.
    void add(cps::QPropW &qp, size_t t, bool periodic) {
        std::vector<cps::WilsonMatrix> &pt = p[t];
        std::vector<cps::WilsonMatrix> &mt = m[t];
        
        if(pt.size() == 0) {
            pt.assign(lcl_vol, cps::WilsonMatrix(0));
        }
        if(mt.size() == 0) {
            mt.assign(lcl_vol, cps::WilsonMatrix(0));
        }
#pragma omp parallel for
        for(size_t i = 0; i < lcl_vol; ++i) {
            pt[i] += qp[i];
            if(periodic) {
                mt[i] += qp[i];
            } else {
                mt[i] -= qp[i];
            }
        }
    }

    void avg(void) {
        for(size_t t = 0; t < t_size_glb; ++t) {
            std::vector<cps::WilsonMatrix> &pt = p[t];
            std::vector<cps::WilsonMatrix> &mt = m[t];

            assert(pt.size() == mt.size());
            if(pt.size() == 0) continue;
#pragma omp parallel for
            for(size_t i = 0; i < lcl_vol; ++i) {
                pt[i] *= 0.5;
                mt[i] *= 0.5;
            }
        }
    }
private:
    const size_t lcl_vol;
    const size_t t_size_glb;

    std::vector<std::vector<cps::WilsonMatrix> > p; // (P+A)/2
    std::vector<std::vector<cps::WilsonMatrix> > m; // (P-A)/2

    std::vector<PropAP> prop;
};


    // This function is not supported and must be checked again when
    // use.
//     void apply_mom(const double mom[4]) {
//         if(mom[3] != 0) {
//             fprintf(stderr, "Adding momentum in t direction is not supported.\n");
//             exit(-1);
//         }

//         const int lcl[4] = {
//             cps::GJP.XnodeSites(), cps::GJP.YnodeSites(),
//             cps::GJP.ZnodeSites(), cps::GJP.TnodeSites(),
//         };

//         const int shift[4] = {
//             cps::GJP.XnodeSites() * cps::GJP.XnodeCoor(),
//             cps::GJP.YnodeSites() * cps::GJP.YnodeCoor(),
//             cps::GJP.ZnodeSites() * cps::GJP.ZnodeCoor(),
//             cps::GJP.TnodeSites() * cps::GJP.TnodeCoor(),
//         };

//         const int glb[4] = {
//             cps::GJP.XnodeSites() * cps::GJP.Xnodes(),
//             cps::GJP.YnodeSites() * cps::GJP.Ynodes(),
//             cps::GJP.ZnodeSites() * cps::GJP.Znodes(),
//             cps::GJP.TnodeSites() * cps::GJP.Tnodes(),
//         };

//         cps::VRB.Result("PropAP", "apply_mom", "mom = %.3e %.3e %.3e %.3e\n",
//                         mom[0], mom[1], mom[2], mom[3]);

//         const size_t lcl_vol = cps::GJP.VolNodeSites();

// #pragma omp parallel for
//         for(int i = 0; i < lcl_vol; ++i) {
//             int glb_x[4];
//             compute_coord(glb_x, lcl, shift, i);

//             const double PI = 3.1415926535897932384626433832795028842;
//             double alpha = 0.;
//             for(int mu = 0; mu < 4; ++mu) {
//                 alpha -= mom[mu] * 2.0 * PI * glb_x[mu] / glb[mu];
//             }

//             wm[i]           *= cps::Rcomplex(std::cos(alpha), std::sin(alpha));
//             wm[i + lcl_vol] *= cps::Rcomplex(std::cos(alpha), std::sin(alpha));
//         }
//     }

#endif
