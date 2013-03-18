// -*- mode:c++; c-basic-offset: 4 -*-
#ifndef INCLUDED_BFM_MIXED_SOLVER_HT_H
#define INCLUDED_BFM_MIXED_SOLVER_HT_H

#include <util/lattice/bfm_evo.h>
#include <alg/enum_int.h>

#include <chroma.h>
#include <omp.h>
#include <pthread.h>

#include <cstdio>
#include <cstdlib>
#include <unistd.h>

#include <sys/stat.h>
#include <sys/time.h>

namespace mixed_cg {

    // check if 2 instances of bfm agree on what they are going to do.
    template<typename Float_out, typename Float_in>
    inline bool check(bfm_evo<Float_out> &bfm_out,
               bfm_evo<Float_in> &bfm_in)
    {
        if(bfm_out.node_latt[0] != bfm_in.node_latt[0]) return false;
        if(bfm_out.node_latt[1] != bfm_in.node_latt[1]) return false;
        if(bfm_out.node_latt[2] != bfm_in.node_latt[2]) return false;
        if(bfm_out.node_latt[3] != bfm_in.node_latt[3]) return false;
        if(bfm_out.Ls != bfm_in.Ls) return false;
        if(bfm_out.precon_5d != bfm_in.precon_5d) return false;
        return true;
    }

    // Convert between single/double precision bfm fermions
    template<typename Float_out, typename Float_in>
    inline void threaded_convFermion(Fermion_t out, Fermion_t in,
                                     bfm_evo<Float_out> &bfm_out,
                                     bfm_evo<Float_in> &bfm_in)
    {
        int me = bfm_out.thread_barrier();

        if(!check(bfm_out, bfm_in)) {
            if(bfm_out.isBoss() && !me) {
                printf("Output/Input fermions don't match.\n");
            }
            exit(-1);
        }

        // Simple copy, this shouldn't be called, right?
        if(sizeof(Float_out) == sizeof(Float_in)) {
            bfm_out.copy(out, in);
            return;
        }

        // otherwise, we do the conversion
        //
        // Note: this function is running under threaded environment.
        int Nspinco=12;
        int out_i_inc = bfm_out.simd() * 2;
        int  in_i_inc =  bfm_in.simd() * 2;

        int out_lat[5] = { bfm_out.node_latt[0],
                           bfm_out.node_latt[1],
                           bfm_out.node_latt[2],
                           bfm_out.node_latt[3],
                           bfm_out.Ls };

        int vol5d_out = out_lat[0] * out_lat[1] * out_lat[2] * out_lat[3] * out_lat[4];


        int thrlen, throff;
        bfm_out.thread_work_nobarrier(vol5d_out, me, thrlen, throff);

        Float_out *outf = (Float_out *)out;
        Float_in  * inf = (Float_in  *)in;

        for(int site = throff; site < throff + thrlen; ++site) {
            int x[4], s;
            int si = site;
            x[0] = si % out_lat[0];    si = si / out_lat[0];
            x[1] = si % out_lat[1];    si = si / out_lat[1];
            x[2] = si % out_lat[2];    si = si / out_lat[2];
            x[3] = si % out_lat[3];
            s    = si / out_lat[3];
        
            // both in and out must have the same preconditioning scheme.
            int sp = bfm_out.precon_5d ? s : 0;
            if ( (x[0]+x[1]+x[2]+x[3] + sp &0x1) == 1) {
                int out_base = bfm_out.bagel_idx5d(x, s, 0, 0, Nspinco, 1);
                int  in_base =  bfm_in.bagel_idx5d(x, s, 0, 0, Nspinco, 1);

                for ( int co=0;co<Nspinco;co++ ) { 
                    for ( int reim=0;reim<2;reim++ ) {
                        int out_id = out_base + reim + co * out_i_inc;
                        int  in_id =  in_base + reim + co *  in_i_inc;

                        outf[out_id] = inf[in_id];
                    }}//co,reim
            }//cb
        }//xyzts
    }

    // Reinitialize communication subsystem.
    // Check bfmcommspi.C in the bfm package to see if this can be avoided.
    template<typename Float_new, typename Float_old>
    inline void switch_comm(bfm_evo<Float_new> &bfm_new, bfm_evo<Float_old> &bfm_old)
    {
        if(static_cast<void *>(&bfm_new)
           == static_cast<void *>(&bfm_old)) return;
        int me = bfm_old.thread_barrier();

        if(!me) {
            bfm_old.comm_end();
            bfm_new.comm_init();
            // Question: how do we propagate the reinitialized information
            // to other threads?
            // Answer: thread barrier does this.
        }
        bfm_new.thread_barrier();
    }

    // Both sol and src are double precision fermions. Single precision
    // solver is only used internally.
    //
    // Things to be set before using this function:
    //
    // double precision solver mass, stopping condition, max iteration
    // number.
    //
    // single precision solver mass, stopping condition, max iteration
    // number.
    //
    // Gauge field must be initialized for both double and single prec
    // solvers.
    //
    // the communication subsystem must be ready for bfm_d to use (due to
    // the way bfmcommspi is written, one must reinitialize the
    // communication object when switching between single and double
    // precisions).
    //
    // max_cycle: the maximum number of restarts will be performed.
    inline int threaded_cg_mixed_MdagM(Fermion_t sol, Fermion_t src,
                                       bfm_evo<double> &bfm_d, bfm_evo<float> &bfm_f,
                                       int max_cycle, cps::InverterType itype = cps::CG,
                                       // the following parameters are for deflation
                                       multi1d<Fermion_t[2]> *evec = NULL,
                                       multi1d<float> *eval = NULL,
                                       int N = 0)
    {
        int me = bfm_d.thread_barrier();

        double frsd = bfm_f.residual;

        Fermion_t src_d = bfm_d.threadedAllocFermion();
        Fermion_t tv1_d = bfm_d.threadedAllocFermion();
        Fermion_t tv2_d = bfm_d.threadedAllocFermion();
        Fermion_t sol_f = bfm_f.threadedAllocFermion();
        Fermion_t src_f = bfm_f.threadedAllocFermion();

        double src_norm = bfm_d.norm(src);
        double stop = src_norm * bfm_d.residual * bfm_d.residual;

        int iter = 0;
        for(int i = 0; i < max_cycle; ++i) {
            // compute double precision rsd and also new RHS vector.
            bfm_d.Mprec(sol,   tv1_d, src_d, 0, 0);
            bfm_d.Mprec(tv1_d, tv2_d, src_d, 1, 0); // tv2_d = MdagM * sol
            double norm = bfm_d.axpy_norm(src_d, tv2_d, src, -1.);

            if(bfm_f.isBoss() && !me) {
                printf("cg_mixed_MdagM: iter = %d rsd = %17.10e(d) stop = %17.10e(d)\n",
                       i, norm, stop);
            }

            // my ad hoc stopping condition
            if(norm < 100. * stop) break;

            // will this cause a deadlock when combined with the
            // condition above?  i.e., will we lose a factor of huge
            // factor in the accuracy of rsd when converting from
            // single to double?
            while(norm * bfm_f.residual * bfm_f.residual < stop) bfm_f.residual *= 2;

            threaded_convFermion(src_f, src_d, bfm_f, bfm_d);
            switch_comm(bfm_f, bfm_d);

            bfm_f.set_zero(sol_f);
            switch(itype) {
            case cps::CG:
                if(evec && eval && N) {
                    bfm_f.deflate(sol_f, src_f, evec, eval, N);
                }
                iter += bfm_f.CGNE_prec_MdagM(sol_f, src_f);
                break;
            case cps::EIGCG:
                iter += bfm_f.Eig_CGNE_prec(sol_f, src_f);
                break;
            default:
                if(bfm_f.isBoss() && !me) {
                    printf("cg_mixed_MdagM: unsupported inverter type.\n");
                }
                exit(-1);
            }

            switch_comm(bfm_d, bfm_f);
            threaded_convFermion(tv1_d, sol_f, bfm_d, bfm_f);

            bfm_d.axpy(sol, tv1_d, sol, 1.);
        }

        bfm_d.threadedFreeFermion(src_d);
        bfm_d.threadedFreeFermion(tv1_d);
        bfm_d.threadedFreeFermion(tv2_d);

        bfm_f.threadedFreeFermion(sol_f);
        bfm_f.threadedFreeFermion(src_f);

        iter += bfm_d.CGNE_prec_MdagM(sol, src);

        double sol_norm = bfm_d.norm(sol);
        if(bfm_d.isBoss() && !me) {
            printf("cg_mixed_MdagM: final sol norm = %17.10e\n", sol_norm);
        }

        bfm_f.residual = frsd;
        return iter;
    }

    // apply single precision solver to double precision vectors. Both
    // sol_d and src_d are in double precision.
    //
    // sol_f and src_f are auxiliary fermions in single
    // precision. Their content will be overriden after calling this
    // function.
    //
    // If import_guess == true, then we import sol_d as an initial
    // guess, other we do a zero start CG.
    inline int cg_single_prec(Fermion_t sol_d, Fermion_t src_d,
                              Fermion_t sol_f, Fermion_t src_f,
                              bfm_evo<double> &bfm_d, bfm_evo<float> &bfm_f)
    {
        threaded_convFermion(src_f, src_d, bfm_f, bfm_d);
        switch_comm(bfm_f, bfm_d);
        
        bfm_f.set_zero(sol_f);
        int iter = bfm_f.CGNE_prec_MdagM(sol_f, src_f);
        
        switch_comm(bfm_d, bfm_f);
        threaded_convFermion(sol_d, sol_f, bfm_d, bfm_f);
        return iter;
    }

    // cg_MdagM_single_precnd: Nested CG, single precision solver is
    // used as a preconditioner.
    //
    // Calling interface is the same as threaded_cg_mixed_MdagM().
    inline int cg_MdagM_single_precnd(Fermion_t sol, Fermion_t src,
                                      bfm_evo<double> &bfm_d, bfm_evo<float> &bfm_f)
    {
        int me = bfm_d.thread_barrier();

        double frsd = bfm_f.residual;

        Fermion_t r = bfm_d.threadedAllocFermion();
        Fermion_t minvr = bfm_d.threadedAllocFermion();
        Fermion_t d = bfm_d.threadedAllocFermion();
        Fermion_t ad = bfm_d.threadedAllocFermion();
        Fermion_t aad = bfm_d.threadedAllocFermion();
        Fermion_t tv1 = bfm_d.threadedAllocFermion();
        Fermion_t sol_f = bfm_f.threadedAllocFermion();
        Fermion_t src_f = bfm_f.threadedAllocFermion();

        const double src_norm = bfm_d.norm(src);
        const double stop = src_norm * bfm_d.residual * bfm_d.residual;

        int iter_s = 0;

        bfm_d.Mprec(sol,  ad, tv1, 0, 0);
        bfm_d.Mprec(ad , aad, tv1, 1, 0); // aad = MdagM * x0 (double prec)
        bfm_d.axpy(r, aad, src, -1.0); // r0 = b - MdagM * x0

        iter_s += cg_single_prec(minvr, r, sol_f, src_f, bfm_d, bfm_f);
        bfm_d.copy(d, minvr); // d0 = (M'dagM')^(-1) * r0 
        double rtminvr = bfm_d.inner_real(r, minvr);

        int k = 1;
        for(; k <= bfm_d.max_iter; ++k) {
            double dtad = bfm_d.Mprec(d, ad, tv1, 0, 1);
            bfm_d.Mprec(ad, aad, tv1, 1, 0);  // aad = MdagM * d[k] (double prec)

            double alpha = rtminvr / dtad;
            bfm_d.axpy(sol, d, sol, alpha);
            double rsd = bfm_d.axpy_norm(r, aad, r, -alpha);
            
            // check watch file
            FILE *fp = fopen("stop.file", "r");
            if(fp) {
                fclose(fp);
                printf("Found watchfile stop.file\n");
                return 1;
            }

            // check stopping condition
            if(rsd < stop) {
                // compute true residual
                bfm_d.Mprec(sol,  ad, tv1, 0, 0);
                bfm_d.Mprec(ad , aad, tv1, 1, 0);
                double true_rsd = bfm_d.axpy_norm(tv1, aad, src, -1.0);

                if(bfm_d.isBoss() && !me) {
                    printf("cg_MdagM_single_precnd: converged in %d(d)+%d(s) iterations.\n", k, iter_s);
                    printf("cg_MdagM_single_precnd: true residual = %17.10e.\n", sqrt(true_rsd/src_norm));
                }
                break;
            }

            iter_s += cg_single_prec(minvr, r, sol_f, src_f, bfm_d, bfm_f);
            double tmp = bfm_d.inner_real(r, minvr);
            double beta = tmp / rtminvr;
            rtminvr = tmp;

            bfm_d.axpby(d, minvr, d, 1.0, beta);
        }

        if(k > bfm_d.max_iter) {
            if(bfm_d.isBoss() && !me) {
                printf("cg_MdagM_single_precnd: CG not converged in %d(d)+%d(s) iterations.\n", k, iter_s);
            }
        }

        bfm_d.threadedFreeFermion(r);
        bfm_d.threadedFreeFermion(minvr);
        bfm_d.threadedFreeFermion(d);
        bfm_d.threadedFreeFermion(ad);
        bfm_d.threadedFreeFermion(aad);
        bfm_d.threadedFreeFermion(tv1);
        bfm_f.threadedFreeFermion(sol_f);
        bfm_f.threadedFreeFermion(src_f);

        bfm_d.CGNE_prec_MdagM(sol, src);

        bfm_f.residual = frsd;
        return k + iter_s;
    }
    
    inline int threaded_cg_mixed_M(Fermion_t sol[2], Fermion_t src[2],
                                   bfm_evo<double> &bfm_d, bfm_evo<float> &bfm_f,
                                   int max_cycle, cps::InverterType itype = cps::CG,
                                   // the following parameters are for deflation
                                   multi1d<Fermion_t[2]> *evec = NULL,
                                   multi1d<float> *eval = NULL,
                                   int N = 0)
    {
        int me = bfm_d.thread_barrier();
        Fermion_t be = bfm_d.threadedAllocFermion();
        Fermion_t bo = bfm_d.threadedAllocFermion();
        Fermion_t ta = bfm_d.threadedAllocFermion(); 
        Fermion_t tb = bfm_d.threadedAllocFermion(); 

        double nsrc = bfm_d.norm(src[0]) + bfm_d.norm(src[1]);
        if (bfm_d.isBoss() && !me ) {
            printf("threaded_cg_mixed_M: source norm is %17.10e\n", nsrc);
        }

        // eo preconditioning
        bfm_d.MooeeInv(src[Even], ta, DaggerNo);
        bfm_d.Meo     (ta, tb, Odd, DaggerNo); // tb == Moe Mee^{-1} src[e]
        bfm_d.axpy    (ta, tb, src[Odd], -1.0);
        bfm_d.Mprec   (ta, bo, tb, DaggerYes); // bo = Mprec^dag (src[o] - Moe Mee^{-1} src[e])

        int iter = threaded_cg_mixed_MdagM(sol[Odd], bo, bfm_d, bfm_f, max_cycle, itype, evec, eval, N);
  
        bfm_d.Meo(sol[Odd], ta, Even, DaggerNo);
        bfm_d.axpy(tb, ta, src[Even], -1.0);
        bfm_d.MooeeInv(tb, sol[Even], DaggerNo);
  
        double nsol = bfm_d.norm(sol[0]) + bfm_d.norm(sol[1]);

        // compute final residual
        Fermion_t tmp[2] = {be, bo};
        bfm_d.Munprec(sol, tmp, ta, DaggerNo);

        double ndiff = 0.;
        for(int i = 0; i < 2; ++i) {
            bfm_d.axpy(tb, tmp[i], src[i], -1.0);
            ndiff += bfm_d.norm(tb);
        }

        if (bfm_d.isBoss() && !me ) {
            printf("threaded_cg_mixed_M: unprec sol norm = %17.10e, residual = %17.10e\n",
                   nsol, sqrt(ndiff / nsrc));
        }
  
        bfm_d.threadedFreeFermion(be);
        bfm_d.threadedFreeFermion(bo);
        bfm_d.threadedFreeFermion(ta);
        bfm_d.threadedFreeFermion(tb);

        return iter;
    }
}

#endif
