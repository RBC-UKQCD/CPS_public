// -*- mode:c++;c-basic-offset:4 -*-
#ifndef INCLUDED_PT_MAT_H__
#define INCLUDED_PT_MAT_H__

#include <config.h>
/*! \file
  \brief  Generic parallel transport.
  
  $Id: pt_mat.h,v 1.5 2013-03-18 19:33:13 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2013-03-18 19:33:13 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/pt_mat.h,v 1.5 2013-03-18 19:33:13 chulwoo Exp $
//  $Id: pt_mat.h,v 1.5 2013-03-18 19:33:13 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: pt_mat.h,v $
//  $Revision: 1.5 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/pt_mat.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
#include <util/gjp.h>
#include <comms/scu.h>
#include <util/omp_wrapper.h>
#include <vector>
#include <cstddef>

CPS_START_NAMESPACE


// namespace pt_generic: contains a single useful routine pt() for
// generic parallel transport. Roughly this is a rewrite of pt_mat()
// (you can find equivalent code in
// src/util/parallel_transport/pt_base/noarch/pt.C) but use threaded
// calculation and parallel communication.
//
// pt() does the following:
// 1. For forward direction (dir = 0, 2, 4, 6 correspond to +X, +Y, +Z, +T),
// out(x) = U_\mu(x) in(x+\mu) or out(x-\mu) = U_\mu(x-\mu) in(x)
//
// 2. For backward direction (dir = 1, 3, 5, 7 correspond to -X, -Y, -Z, -T),
// out(x) = U_\mu(x-\mu)^\dag in(x-\mu) or out(x+\mu) = U_\mu(x)^\dag in(x)
namespace pt_generic {
    enum { NUM_DIR = 8 };

    static inline int idx_4d(const int x[4], const int lx[4]) {
        int ret = 0;
        for(int i = 3; i >= 0; --i) {
            ret = ret * lx[i] + x[i];
        }
        return ret;
    }

    static inline int idx_4d_surf(const int x[4], const int lx[4], int mu) {
        int ret = 0;
        for(int i = 3; i >= 0; --i) {
            if(i == mu) continue;
            ret = ret * lx[i] + x[i];
        }
        return ret;
    }

    // Compute all size and offset data need for this batch of parallel transport.
    // Everything is measured in # of elements of the Site object.
    static void compute_comm(std::vector<int> &comm_size,
                             std::vector<int> &comm_off,
                             int total[],
                             int N, const int dir[])
    {
        int lcl[4] = { GJP.XnodeSites(), GJP.YnodeSites(), GJP.ZnodeSites(), GJP.TnodeSites() };
        int vol4d = GJP.VolNodeSites();
        int surf[4] = {vol4d / lcl[0], vol4d / lcl[1], vol4d / lcl[2], vol4d / lcl[3] };

        for(int i = 0; i < NUM_DIR; ++i) total[i] = 0;

        comm_size.resize(N);
        comm_off.resize(N);

        for(int i = 0; i < N; ++i) {
            int mu = dir[i];
            comm_size[i] = surf[mu / 2];
            comm_off[i] = total[mu];
            total[mu] += comm_size[i];
        }
    }

    static void thread_work_partial(int nwork, int me, int nthreads,
                                    int &mywork, int &myoff)
    {
        int basework = nwork / nthreads;
        int backfill = nthreads - (nwork % nthreads);
        mywork = (nwork + me) / nthreads;
        myoff  = basework * me;
        if ( me > backfill ) 
            myoff += (me-backfill);
    }

    // collect surface data
    template<typename Link, typename Site>
    void collect(Site buf[], Link gauge[], Site vol[], int mu, bool do_send)
    {
        int lcl[4] = { GJP.XnodeSites(), GJP.YnodeSites(), GJP.ZnodeSites(), GJP.TnodeSites() };
        int low[4] = { 0, 0, 0, 0 };
        int high[4] = { lcl[0], lcl[1], lcl[2], lcl[3] };

        bool backwards = mu % 2 ? true : false;
        mu /= 2;

        low[mu] = do_send == backwards ? lcl[mu] - 1 : 0;
        high[mu] = low[mu] + 1;

        const int hl[4] = {high[0] - low[0],
                           high[1] - low[1],
                           high[2] - low[2],
                           high[3] - low[3] };
        const int hl_sites = hl[0] * hl[1] * hl[2] * hl[3];

        int me = omp_get_thread_num();
        int nthreads = omp_get_num_threads();
        int mywork, myoff;

        thread_work_partial(hl_sites, me, nthreads, mywork, myoff);
        for(int i = 0; i < mywork; ++i) {
            int x[4];
            int tmp = i + myoff;
            x[0] = tmp % hl[0] + low[0]; tmp /= hl[0];
            x[1] = tmp % hl[1] + low[1]; tmp /= hl[1];
            x[2] = tmp % hl[2] + low[2]; tmp /= hl[2];
            x[3] = tmp % hl[3] + low[3];

            int off_4d = idx_4d(x, lcl);
            int off_3d = idx_4d_surf(x, lcl, mu);
        
            if(do_send) {
                if(backwards) {
                    Link tmp;
                    tmp.Dagger(gauge[mu + 4 * off_4d]);
                    buf[off_3d] = tmp * vol[off_4d];
                } else {
                    buf[off_3d] = vol[off_4d];
                }
            } else {
                if(backwards) {
                    vol[off_4d] = buf[off_3d];
                } else {
                    vol[off_4d] = gauge[mu + 4 * off_4d] * buf[off_3d];
                }
            }
        }
    }

    template<typename Link, typename Site>
    void compute_internal(Site out[], Link gauge[], Site in[], int mu, int nthreads)
    {
        int lcl[4] = { GJP.XnodeSites(), GJP.YnodeSites(), GJP.ZnodeSites(), GJP.TnodeSites() };
        int low[4] = { 0, 0, 0, 0 };
        int high[4] = { lcl[0], lcl[1], lcl[2], lcl[3] };

        bool backwards = mu % 2 ? true : false;
        mu /= 2;

        low[mu] = backwards ? 0 : 1;
        high[mu] = low[mu] + lcl[mu] - 1;

        const int hl[4] = {high[0] - low[0],
                           high[1] - low[1],
                           high[2] - low[2],
                           high[3] - low[3] };
        const int hl_sites = hl[0] * hl[1] * hl[2] * hl[3];

        int me = omp_get_thread_num();
        int mywork, myoff;

        thread_work_partial(hl_sites, me, nthreads, mywork, myoff);
        for(int i = 0; i < mywork; ++i) {
            int x[4];
            int tmp = i + myoff;
            x[0] = tmp % hl[0] + low[0]; tmp /= hl[0];
            x[1] = tmp % hl[1] + low[1]; tmp /= hl[1];
            x[2] = tmp % hl[2] + low[2]; tmp /= hl[2];
            x[3] = tmp % hl[3] + low[3];

            int y[4] = {x[0], x[1], x[2], x[3]};
            y[mu] = backwards ? x[mu] + 1 : x[mu] - 1;
            int off_x = idx_4d(x, lcl);
            int off_y = idx_4d(y, lcl);

            if(backwards) {
                Link tmp;
                tmp.Dagger(gauge[4 * off_x + mu]);
                out[off_y] = tmp * in[off_x];
            } else {
                out[off_y] = gauge[4 * off_y + mu] * in[off_x];
            }
        }
    }

    template<typename Link, typename Site>
    void pt(int N, Site *out[], Link *gauge[], Site *in[], const int dir[])
    {
        // 1. compute comm space requirement and pointer shift
        std::vector<int> comm_size;
        std::vector<int> comm_off;
        int total[NUM_DIR];

        compute_comm(comm_size, comm_off, total, N, dir);

        // 2. allocate comm space
        Site *sndbuf[NUM_DIR];
        Site *rcvbuf[NUM_DIR];

        for(int i = 0; i < NUM_DIR; ++i) {
            if(total[i] == 0) { // no communication in this direction
                sndbuf[i] = rcvbuf[i] = NULL;
            } else {
                sndbuf[i] = new Site[total[i]];
                rcvbuf[i] = new Site[total[i]];
            }
        }

        // 3. compute/collect surface sites for communication
#pragma omp parallel
        {
            for(int d = 0; d < N; ++d) {
                collect(sndbuf[dir[d]] + comm_off[d], gauge[d], in[d], dir[d], true);
            }
        }

        // 4. single threaded comm
        // If you fix the multi thread MPI problem on BG/Q you can do multithreaded comm here.
        for(int mu = 0; mu < 4; ++mu) {
            // forwards                
            if(total[2 * mu]) {
                getPlusData((Float *)rcvbuf[2 * mu], (Float *)sndbuf[2 * mu],
                            sizeof(Site) * total[2 * mu] / sizeof(Float), mu);
            }
            // backwards
            if(total[2 * mu + 1]) {
                getMinusData((Float *)rcvbuf[2 * mu + 1], (Float *)sndbuf[2 * mu + 1],
                             sizeof(Site) * total[2 * mu + 1] / sizeof(Float), mu);
            }
        }

#pragma omp parallel
        {
            int nthreads = omp_get_num_threads();
            // int me = omp_get_thread_num();

            // 5. compute internal parallel transport
            for(int d = 0; d < N; ++d) {
                compute_internal(out[d], gauge[d], in[d], dir[d], nthreads);
            }
            //#pragma omp barrier

            // 6. compute/collect surface sites 
            for(int d = 0; d < N; ++d) {
                collect(rcvbuf[dir[d]] + comm_off[d], gauge[d], out[d], dir[d], false);
            }
        }

        for(int i = 0; i < NUM_DIR; ++i) {
            delete[] sndbuf[i];
            delete[] rcvbuf[i];
        }
    }
}

CPS_END_NAMESPACE

#endif
