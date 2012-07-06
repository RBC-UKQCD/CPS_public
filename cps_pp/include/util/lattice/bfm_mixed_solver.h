// -*- mode:c++; c-basic-offset: 4 -*-
#ifndef INCLUDED_BFM_MIXED_SOLVER_HT_H
#define INCLUDED_BFM_MIXED_SOLVER_HT_H

#include <util/lattice/bfm_evo.h>

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
                                       int max_cycle)
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
            // double sol_norm = bfm_d.norm(sol);

            // compute double precision rsd and also new RHS vector.
            bfm_d.Mprec(sol,   tv1_d, src_d, 0, 0);
            bfm_d.Mprec(tv1_d, tv2_d, src_d, 1, 0); // tv2_d = MdagM * sol
            double norm = bfm_d.axpy_norm(src_d, tv2_d, src, -1.);

            if(bfm_f.isBoss() && !me) {
                printf("cg_mixed_MdagM: iter = %d rsd = %17.10e(d) stop = %17.10e(d)\n",
                       i, norm, stop);
            }

            // my ad hoc stopping condition
            // if(norm < stop) break;
            if(norm < 10000. * stop) break;

            // will this cause a deadlock when combined with the
            // condition above?  i.e., will we lose a factor of
            // 10000/4. in the accuracy of rsd when converting from
            // single to double?
            while(norm * bfm_f.residual * bfm_f.residual < stop) bfm_f.residual *= 2;

            threaded_convFermion(src_f, src_d, bfm_f, bfm_d);
            switch_comm(bfm_f, bfm_d);

            bfm_f.set_zero(sol_f);
            iter += bfm_f.CGNE_prec_MdagM(sol_f, src_f);

            // double norm_sol_f_1 = bfm_f.norm(sol_f);
            switch_comm(bfm_d, bfm_f);
            threaded_convFermion(tv1_d, sol_f, bfm_d, bfm_f);
            // double norm_sol_f_2 = bfm_d.norm(tv1_d);
            // if(bfm_d.isBoss() && !me) {
            //     printf("cg_mixed_MdagM: iter = %d sol norm = %17.10e(f) %17.10e(d)\n",
            //            i, norm_sol_f_1, norm_sol_f_2);
            // }

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

    struct mixed_cg_thread_arg {
        bfm_evo<double> &bfm_d;
        bfm_evo<float>  &bfm_f;

        Fermion_t src; // double precision
        Fermion_t sol; // double precision

        int max_cycle;
    };

    inline void *mixed_cg_thread_function(void *myarg)
    {
        mixed_cg_thread_arg *args = (mixed_cg_thread_arg *)myarg;
    
        bfm_evo<double> &bfm_d = args->bfm_d;
        bfm_evo<float>  &bfm_f = args->bfm_f;
        Fermion_t src = args->src;
        Fermion_t sol = args->sol;
        int max_cycle = args->max_cycle;

        args->max_cycle = threaded_cg_mixed_MdagM(sol, src, bfm_d, bfm_f, max_cycle);

        return myarg;
    }

    inline void mixed_cg_spawn_and_run(void *(*thread_func)(void *), void *arg)
    {
        int THREADS=64;
        pthread_t *handles = (pthread_t *) malloc(THREADS*sizeof(pthread_t));

        for(int t =0; t<THREADS; t++ ) {
            if ( pthread_create(&handles[t],(pthread_attr_t *)NULL, thread_func, arg))
                perror("pthread create");
        }

        for(int t =0; t<THREADS; t++ ) {
            if ( pthread_join(handles[t],NULL)) {
                perror("pthread_join");
            }
        }

        free (handles);
    }

    // non-threaded version, just a wrapper
    inline int cg_mixed_MdagM(Fermion_t sol, Fermion_t src,
                              bfm_evo<double> &bfm_d, bfm_evo<float> &bfm_f,
                              int max_cycle)
    {
        mixed_cg_thread_arg arg = { bfm_d, bfm_f, src, sol, max_cycle };
        mixed_cg_spawn_and_run(mixed_cg_thread_function, (void *)&arg);
        return arg.max_cycle;
    }

};

#endif
