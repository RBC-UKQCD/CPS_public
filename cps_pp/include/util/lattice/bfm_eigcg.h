/* -*- mode:c++; c-basic-offset:4 -*- */
#ifndef BFM_EIGCG_H
#define BFM_EIGCG_H

/* EigCG code for Mpc^dagger Mpc solver 2012 by Qi
 * copied from the CG code (2012 ?version in bfm) and made modifications and 
 * adding the eigCG part
 */

#include <bfm.h>
#include <bfm_qdp.h>

#include <util/lattice/bfm_evo.h>
#include <util/lattice/eigcg_controller.h>

#include <stdio.h>
#include <cmath>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <cstring>
#include <cassert>

#ifdef LIST_ENGINE
#include <pprefetch.h>
#endif

#include<qmp.h>

void eigen_solver(double *A, double *EV, double *E, int n);
template<class Float> void matrix_dgemm (const int M,const int N, const int K, Float **A, const double *B, double *C);
void min_eig_index(int *INDEX, int nev,double *EIG, int n);
void invert_H_matrix(complex<double> *data, int n); //if n is large enough, must be parallerized!!!
template<class Float> void eigcg_vec_mult(Float* V, const int m, double *QZ, const int n, const int f_size_cb, const int nthread, const int me);

template <class Float>
Fermion_t bfm_evo<Float>::threadedAllocCompactFermion   (int mem_type)
{
    Fermion_t ret;
    int me = this->thread_barrier();
    if ( me == 0 ) {
        ret = this->allocCompactFermion(mem_type);
    }
    ret = this->thread_bcast(me,ret);
    this->thread_barrier();
    return ret;
}

template <class Float>
Fermion_t bfm_evo<Float>::allocCompactFermion   (int mem_type)
{
    int words = 24 * this->node_cbvol;
    Float *ferm = (Float *)bfm_alloc(words*this->cbLs*sizeof(Float),mem_type);

    if(ferm == 0){
        if ( this->isBoss() ) printf("bfmbase::allocFermion\n");
        fflush(stdout);
        exit(-1);
    }
    return (Fermion_t)ferm;
}

template <class Float>
void* bfm_evo<Float>::threaded_alloc(int length, int mem_type)
{
    void *v;
    int me = this->thread_barrier();
    if(me==0)
        {
            v = bfm_alloc(length, mem_type);
            if(v== NULL){
                if ( this->isBoss() ) printf("bfmbase::allocFermion failed\n");
                fflush(stdout);
                exit(-1);
            }
        }
    v = this->thread_bcast(me, v);
    this->thread_barrier();
    return v;
}

template <class Float>
void bfm_evo<Float>::threaded_free(void* handle)
{
    int me = this->thread_barrier();
    if ( me == 0 ) { 
        Float * ferm = ( Float * ) handle;
        bfm_free(ferm);
    }
    this->thread_barrier();
}

template<class Float>
int bfm_evo<Float>::EIG_CGNE_M(Fermion_t solution[2], Fermion_t source[2])
{
    int me = this->thread_barrier();
    Fermion_t src = this->threadedAllocFermion(); 
    Fermion_t tmp = this->threadedAllocFermion(); 
    Fermion_t Mtmp= this->threadedAllocFermion(); 

    // src_o = Mdag * (source_o - Moe MeeInv source_e)
    this->MooeeInv(source[Even],tmp,DaggerNo);
    this->Meo(tmp,src,Odd,DaggerNo);
    this->axpy(tmp,src,source[Odd],-1.0);
    this->Mprec(tmp,src,Mtmp,DaggerYes);  
  
    int iter = this->Eig_CGNE_prec(solution[Odd], src);

    // sol_e = M_ee^-1 * ( src_e - Meo sol_o )...
    this->Meo(solution[Odd],tmp,Even,DaggerNo);
    this->axpy(src,tmp,source[Even],-1.0);
    this->MooeeInv(src,solution[Even],DaggerNo);
  
    this->threadedFreeFermion(tmp);
    this->threadedFreeFermion(src);
    this->threadedFreeFermion(Mtmp);

    return iter;
}

template<class Float>
int bfm_evo<Float>::Eig_CGNE_prec(Fermion_t psi, Fermion_t src)
{
    double f;
    double cp,c,a,d,b;
    int me = this->thread_barrier();

    Fermion_t p   = this->threadedAllocFermion(mem_fast); 
    Fermion_t tmp = this->threadedAllocFermion(mem_fast); 
    Fermion_t mp  = this->threadedAllocFermion(mem_fast); 
    Fermion_t mmp = this->threadedAllocFermion(mem_fast); 
    Fermion_t r   = this->threadedAllocCompactFermion(mem_fast); 
    Fermion_t mmp_prev = this->threadedAllocCompactFermion(mem_fast); 


    //NOTICE!!!!
    //To save memory, here we should use float!! no matter what Float is.
    //But then we need to modify a few functions like this->axpy, this->inner
    //to take float in mix with double in bfmBase
    EigCGController<Float>* eigcg = EigCGController<Float>::getInstance();

    assert(eigcg->get_vec_len() == this->cbLs*24*this->node_cbvol);

    complex<double> *invH=NULL;
    if(eigcg->def_len > 0) {
        struct timeval proj_start,proj_end,proj_diff;
        gettimeofday(&proj_start,NULL);

        if(eigcg->def_len<eigcg->get_max_def_len()) {
            invH = (complex<double>*)threaded_alloc(eigcg->def_len*eigcg->def_len*sizeof(complex<double>));

            //calculate invH from eigcg->H
            //---possiblly need to improve--
            if(me==0) {
                for(int i=0;i<eigcg->def_len;i++)
                    for(int j=0;j<eigcg->def_len;j++)
                        invH[i*eigcg->def_len+j]=eigcg->H[i*eigcg->get_max_def_len()+j];
                invert_H_matrix(invH, eigcg->def_len);
            }
            this->thread_barrier();
            //------
        } else {
            invH = eigcg->H;
        }

        complex<double> *Ub = NULL;
        complex<double> *invHUb = NULL;
        Ub = (complex<double>*)threaded_alloc(eigcg->def_len*sizeof(complex<double>));
        invHUb = (complex<double>*)threaded_alloc(eigcg->def_len*sizeof(complex<double>));
        //Do initial projection
        //low modes guess = U*invH*U^dag*src;
        for(int n=0; n<eigcg->def_len; n++) {
            complex<double> dt = this->inner((Fermion_t)(eigcg->getU(n)), src);
            if(me==0)Ub[n] = dt;
        }
        if(me==0) {
            for(int i=0;i<eigcg->def_len;i++) {
                invHUb[i]=0.0;
                for(int j=0;j<eigcg->def_len;j++)
                    invHUb[i] += invH[i*eigcg->def_len+j]*Ub[j]; //notice the index ofr invH.
            }
        }

        // FIXME: I don't know why copying stuff directly to psi doesn't work.
        this->set_zero(tmp);
        for(int i=0;i<eigcg->def_len;i++) {
            this->zaxpy(tmp, (Fermion_t)(eigcg->getU(i)), tmp, invHUb[i]); 
            // this->copy(tmp,(Fermion_t)(eigcg->getU(i)));
            // this->scale(tmp,invHUb[i].real(), invHUb[i].imag());
            // this->axpy(psi,tmp,psi,1.0);
        }
        this->copy(psi, tmp);

        threaded_free(Ub);
        threaded_free(invHUb);

        gettimeofday(&proj_end,NULL);
        timersub(&proj_end,&proj_start,&proj_diff);

        if ( this->isBoss() && this->verbose && !me ) {
            printf("bfmbase::eig_CGNE_prec projection of lowmodes in %d.%6.6d s\n",proj_diff.tv_sec,proj_diff.tv_usec);
        }
    }

    //Initial residual computation & set up
    d=this->Mprec(psi,mp,tmp,DaggerNo);
    b=this->Mprec(mp,mmp,tmp,DaggerYes);

    this->axpy (r, mmp, src,-1.0);
    this->copy(p,r);
    cp= this->norm(r);

    Float ssq =  this->norm(src);
    Float rsq =  this->residual* this->residual*ssq;

    //Check if guess with lowmodes is really REALLY good :)
    if ( cp <= rsq ) {
        if ( this->verbose && this->isBoss() && !me ) {
            printf("bfmbase::CGNE_prec k=0 converged - suspiciously nice guess %le %le\n",cp,rsq);
        }
        this->threadedFreeFermion(tmp);
        this->threadedFreeFermion(p);
        this->threadedFreeFermion(mp);
        this->threadedFreeFermion(mmp);
        this->threadedFreeFermion(mmp_prev);
        this->threadedFreeFermion(r);

        if(eigcg->def_len>0 && eigcg->def_len<eigcg->get_max_def_len()) {
            threaded_free(invH);
        }
        return 0;
    }

    if ( this->verbose && this->isBoss() && !me ) {
        printf("bfmbase::Eig_CGNE_prec k=0 residual %le rsq %le\n",cp,rsq);
    }

    double restart_rsq = 0.0;
    int restart_pos = 0;
    if(eigcg->restart.size()>0) {
        restart_rsq = eigcg->restart[restart_pos]*eigcg->restart[restart_pos]*ssq;
    }

    const bool do_eigcg = (eigcg->get_nev()>0 && eigcg->def_len<eigcg->get_max_def_len());
    int i_eig = 0;
    a = 1.0;
    b = 0.0;
    double a_old = a;
    double b_old = b;

    double *T=NULL;
    double *Tt=NULL;
    double *Q=NULL;
    double *QZ=NULL;
    double *H=NULL;
    double *Y=NULL;
    double *EIG=NULL;
    int *INDEX=NULL;
    double *M=NULL;
    if(do_eigcg) {
        int m=eigcg->get_m();
        int nev=eigcg->get_nev();
        T = (double*)threaded_alloc(m*(m+1)/2*sizeof(double));
        // mxm real symetric matrix, lower triangular only
        Tt = (double*)threaded_alloc(m*(m+1)/2*sizeof(double));
        // (m)x(m) real symetric matrix, lower triangular only
        //T(i,j) = T[ i(i+1)/2+j]
        //NOTICE!! T is actually very sparse matrix, so it can be speed up by a very large factor if wanted
        //		Y = new double[m*m];//m*m real matrix, to store eigen vectors when needed
        Y = (double*)threaded_alloc(m*m*sizeof(double));
        //m*m real matrix, to store eigen vectors when needed
        Q = (double*)threaded_alloc(m*2*nev*sizeof(double));
        //m*(2nev) real matrix
        QZ = (double*)threaded_alloc(m*2*nev*sizeof(double));
        //m*(2nev) real matrix
        H = (double*)threaded_alloc(2*nev*(2*nev+1)/2*sizeof(double));
        //(2nev)*(2nev) real matrix, symetric, lower part
        EIG = (double*)threaded_alloc(m*sizeof(double));
        M = (double*)threaded_alloc(2*nev*sizeof(double));
        //eigen value
        INDEX = (int*)threaded_alloc(nev*sizeof(int));
        //this vector is not necessay if we have a good eig solver for nev low
        if(me==0) {
            for(int i=0;i<m*(m+1)/2;i++) {
                T[i] = Tt[i] = 0.0;
            }
            for(int i = 0; i < m * m; ++i) {
                Y[i] = 0.0;
            }
            for(int i = 0; i < m * 2 * nev; ++i) {
                Q[i] = QZ[i] = 0.0;
            }
            for(int i = 0; i < 2*nev*(2*nev+1)/2; ++i) {
                H[i] = 0.0;
            }
            for(int i = 0; i < m; ++i) {
                EIG[i] = 0.0;
            }
            for(int i = 0; i < 2 * nev; ++i) {
                M[i] = 0.0;
            }
            for(int i = 0; i < nev; ++i) {
                INDEX[i] = 0;
            }
        }
    }

    struct timeval start,stop;
    gettimeofday(&start,NULL);

    double eigProj_flops = 0.0;
    double t_lancos_proj = 0.0; 
    double t_lancos = 0.0;
    double t = 0.0;//CG time

    for (int k=1;k<=this->max_iter;k++) {
        this->iter=k;

        if(do_eigcg && i_eig==eigcg->get_m()) {
            this->copy(mmp_prev,mmp);
        }

        c=cp; //this->norm_r
        d = this->Mprec(p,mp,tmp,0,1);
        double qq=this->Mprec(mp,mmp,tmp,1); 

        if(do_eigcg) {
            struct timeval proj_start,proj_end,proj_diff;
            gettimeofday(&proj_start,NULL);

            const int m = eigcg->get_m();
            const int nev=eigcg->get_nev();
            const int nev2 = 2*nev;
            int rank;
            //T(i_eig-1,i_eig-1)
            if(me==0) {
                if(i_eig!=0)T[(i_eig-1)*i_eig/2+i_eig-1]=1.0/a+b_old/a_old;
            }
            if(i_eig==m) {
                //Yb need lowest nev eigen vector of T(m-1)
                //Y need lowest nev eigen vector of T(m)
                if(me==0) {
                    for(int i=0;i<m*(m-1)/2;i++)Tt[i]=T[i];
                    eigen_solver(Tt,Y,EIG,m-1);//NOTICE, this is NOT needed, only need to calculate the lowest nev, not all m >2*nev
                    //NOTICE, Y transpose is the eigen vectors.
                    min_eig_index(INDEX,nev,EIG,m-1);
                    //Q(nev:2nev-1)=Yb; with Yb last row zero
                    for(int i=nev;i<2*nev;i++) {
                        //Y transpose is eigen vector
                        //for(int j=0;j<m-1;j++)Q[j*2*nev+i]=Y[j+INDEX[i-nev]*(m-1)];
                        //Y is eigen vector
                        for(int j=0;j<m-1;j++)Q[j*2*nev+i]=Y[j*(m-1)+INDEX[i-nev]];
                        Q[(m-1)*2*nev+i]=0.0;
                    }
                    for(int i=0;i<m*(m+1)/2;i++)Tt[i]=T[i];
                    eigen_solver(Tt,Y,EIG,m);//NOTICE, this is NOT needed, only need to calculate the lowest nev, not all m >2*nev
                    min_eig_index(INDEX,nev,EIG,m);
                    //Q(0:nev-1)=Yb; 
                    for(int i=0;i<nev;i++) {
                        //for(int j=0;j<m;j++)Q[j*2*nev+i]=Y[j+INDEX[i]*m];
                        //Y is eigen vector
                        for(int j=0;j<m;j++)Q[j*2*nev+i]=Y[j*m+INDEX[i]];
                    }

                    //Q=orth([Y,Yb]); with YB last row zero
                    //rank(Q) may be smaller than 2*nev. remove these
                    //should be optimized. maybe save row first for Q
                    rank=nev;
                    for(int i=nev;i<2*nev;i++) {
                        for(int j=0;j<rank;j++) {
                            double xy=0.0;
                            for(int k=0;k<m;k++)xy+=Q[k*2*nev+i]*Q[k*2*nev+j];
                            for(int k=0;k<m;k++)Q[k*2*nev+i]-=xy*Q[k*2*nev+j];
                        }
                        //this->normalize
                        double xx=0.0;
                        for(int k=0;k<m;k++)xx+=Q[k*2*nev+i]*Q[k*2*nev+i];
                        if(xx>1e-16) {
                            xx=1.0/std::sqrt(xx);
                            for(int k=0;k<m;k++)Q[k*2*nev+rank]=xx*Q[k*2*nev+i];
                            rank++; 
                        }
                    }
                    //printf("Rank of Q = %d\n",rank);
                    //H=Q' * T * Q
                    for(int i=0;i<rank;i++) {
                        for(int j=0;j<=i;j++) {
                            H[i*(i+1)/2+j]=0.0;
                            for(int l=0;l<m;l++)
                                for(int k=0;k<m;k++) {
                                    if(k<=l)H[i*(i+1)/2+j]+=Q[l*nev2+i]*T[l*(l+1)/2+k]*Q[k*nev2+j];
                                    else H[i*(i+1)/2+j]+=Q[l*nev2+i]*T[k*(k+1)/2+l]*Q[k*nev2+j];
                                }
                        }
                    }
                    //[Z,M]=eig(H)
                    eigen_solver(H,Y,M,rank);
                    for(int i=rank;i<2*nev;i++)M[i]=0.0;//set M[i>=rank] to zero
                    //Notice, Y transpose is eigenvectos.
                    if(0)for(int i=0;i<rank;i++)printf("eig %d : %e \n",i,M[i]);
				
                    //V=V*(Q*Z)
                    //1.QZ=Q*Z
                    //transpoze QZ here to speed up the later calculation
                    for(int j=0;j<rank;j++) {
                        for(int i=0;i<m;i++) {
                            QZ[i+m*j]=0.0;
                            for(int k=0;k<rank;k++)QZ[i+m*j]+=Q[i*nev2+k]*Y[j+k*rank];
                        }
                    }
                }
                void* ptmp = &rank;
                ptmp = this->thread_bcast(me, ptmp);
                rank = *((int*)(ptmp)); //broad cast rank that is only calculated on me==0
                this->thread_barrier();

                i_eig=rank;
                //2.V=V*QZ need to be implement very efficiently
                //. The way we save QZ is transpozed to column first 
                struct timeval proj_start_1,proj_end_1,proj_diff_1;
                gettimeofday(&proj_start_1,NULL);

                // eigcg_vec_mult(eigcg->getV(0),m,QZ,rank,eigcg->get_vec_len(), this->nthread, me);
                eigcg_vec_mult2(eigcg->getV(0),m,QZ,rank,eigcg->get_vec_len(), this->nthread, me, *this);
                // eigcg_vec_mult3(eigcg->getV(0),m,QZ,rank,eigcg->get_vec_len(), this->nthread, me, *this);

                eigProj_flops += 2*eigcg->get_vec_len()*m*rank;
                gettimeofday(&proj_end_1,NULL);
                timersub(&proj_end_1,&proj_start_1,&proj_diff_1);
                t_lancos_proj += proj_diff_1.tv_sec*1.0E6 + proj_diff_1.tv_usec;

                if(me==0) {
                    //T=M
                    for(int i=0;i<m*(m+1)/2;i++)T[i]=0.0;
                    for(int i=0;i<rank;i++)T[i*(i+1)/2+i]=M[i];
                }
                //w=mmp-beta*mmp_prev
                this->axpy(mmp_prev,mmp_prev,mmp,-b);
                //T(i_eig+1,1:i_eig)=w^dag * V/sqrt(rsq) 
                //T is symmetric and REAL !! TESTED
                for(int ii=0;ii<i_eig;ii++) {
                    //T(i_eig,ii)
                    double dt=this->inner_real(mmp_prev,eigcg->getV(ii));
                    if(me==0) {
                        T[i_eig*(i_eig+1)/2+ii] = dt;
                        T[i_eig*(i_eig+1)/2+ii]/=std::sqrt(cp);
                    }
                }
            } else {
                if(me==0) {
                    if(i_eig!=0) {
                        //T(i_eig,i_eig-1)
                        T[i_eig*(i_eig+1)/2+i_eig-1]=-sqrt(b)/a;
                    }
                }
            }
            //V[i_eig]=r/sqrt(rsq);
            this->copy(eigcg->getV(i_eig),r);
            this->scale(eigcg->getV(i_eig),1.0/std::sqrt(cp));
            i_eig++;

            gettimeofday(&proj_end,NULL);
            timersub(&proj_end,&proj_start,&proj_diff);
            t_lancos += proj_diff.tv_sec*1e6+proj_diff.tv_usec;
        }

        a_old = a;
        a = c/d; //alpha
        b_old = b;
        b  = a*(a*qq-d)/c; //beta
        cp = this->cg_psi_p(psi,p,r,mmp,a,b);
        //r=r-a*mmp;
        //psi=psi+a*p;
        //p=r+b*p;
        //cp = this->norm_r

        //restart condition
        if(cp <= restart_rsq && eigcg->def_len >0 && (eigcg->def_len==eigcg->get_max_def_len() || eigcg->is_always_restart())) {
            struct timeval proj_start,proj_end,proj_diff;
            gettimeofday(&proj_start,NULL);
            restart_pos++;
            if(restart_pos >= eigcg->restart.size())restart_rsq=0.0;
            else restart_rsq = eigcg->restart[restart_pos]*eigcg->restart[restart_pos]*ssq;
            complex<double> *Ub = NULL;
            complex<double> *invHUb = NULL;
            Ub = (complex<double>*)threaded_alloc(eigcg->def_len*sizeof(complex<double>));
            invHUb = (complex<double>*)threaded_alloc(eigcg->def_len*sizeof(complex<double>));
            //Do projection on A e = r = b-Ax
            //low modes guess = U*invH*U^dag*r;
            for(int n=0; n<eigcg->def_len; n++) {
                complex<double> dt = this->inner((Fermion_t)(eigcg->getU(n)), r);
                if(me==0)Ub[n] = dt;
            }
            if(me==0) {
                for(int i=0;i<eigcg->def_len;i++) {
                    invHUb[i]=0.0;
                    for(int j=0;j<eigcg->def_len;j++)
                        invHUb[i] += invH[i*eigcg->def_len+j]*Ub[j]; //notice the index ofr invH.
                }
            }
            //set tmp=0.0;
            this->set_zero(tmp);
            for(int i=0;i<eigcg->def_len;i++) {
                this->zaxpy(tmp, (Fermion_t)(eigcg->getU(i)), tmp, invHUb[i]);
                // this->copy(mp,(Fermion_t)(eigcg->getU(i)));
                // this->scale(mp,invHUb[i].real(), invHUb[i].imag());
                // this->axpy(tmp,mp,tmp,1.0);
            }
            threaded_free(Ub);
            threaded_free(invHUb);

            this->axpy(psi,tmp,psi,1.0);//update solution

            d=this->Mprec(psi,mp,tmp,DaggerNo);
            b=this->Mprec(mp,mmp,tmp,DaggerYes);

            this->axpy (r, mmp, src,-1.0);
            this->copy(p, r);
            cp= this->norm(r);

            i_eig = 0;
            a = 1.0;
            b = 0.0;
            a_old = a;
            b_old = b;

            gettimeofday(&proj_end,NULL);
            timersub(&proj_end,&proj_start,&proj_diff);

            if ( this->isBoss() && !me && this->verbose) {
                printf("bfmbase::eig_CGNE_prec restart with projection of lowmodes in %d.%6.6d s\n",proj_diff.tv_sec,proj_diff.tv_usec);
            }
            t -= proj_diff.tv_sec*1e6+proj_diff.tv_usec;
        }

        if ( k%50==0 && this->verbose && this->isBoss() && !me ) {
            printf("residual[%d] = %le \n", k, cp);
        }

        // Stopping condition
        if ( cp <= rsq ) { 

            gettimeofday(&stop,NULL);
            struct timeval diff;
            timersub(&stop,&start,&diff);

            if ( this->isBoss() && !me ) printf("bfmbase::CGNE_prec converged in %d iterations\n",k);
            if ( this->isBoss() && !me ) printf("bfmbase::CGNE_prec converged in total %d.%6.6d s\n",diff.tv_sec,diff.tv_usec);


            double flops = this->mprecFlops()*2.0 + 2.0*this->axpyNormFlops() + this->axpyFlops()*2.0;
            flops = flops * k;

            t += diff.tv_sec*1.0E6 + diff.tv_usec;
            t -= t_lancos;
            if ( this->isBoss()&& !me ) {
                //printf("bfmbase::Eig_CGNE_prec: %d mprec flops/site\n",mprecFlopsPerSite());
                printf("bfmbase::Eig_CGNE_prec: CG part: %le s for %le flops\n",t/1e6, flops);
                printf("bfmbase::Eig_CGNE_prec: CG part: %le Mflops per node\n",flops/t);
                if(do_eigcg) {
                    printf("bfmbase::eig_CGNE_prec lancos total in %f s\n", t_lancos/1e6);
                    printf("bfmbase::eig_CGNE_prec lancos project in %f s and %f Mflops \n", t_lancos_proj/1e6, eigProj_flops/t_lancos_proj);
                }
            }

            this->Mprec(psi,mp,tmp,0);
            this->Mprec(mp,mmp,tmp,1); 
            this->axpy(tmp,src,mmp,-1.0);
            double true_residual = sqrt(this->norm(tmp)/this->norm(src));
            if ( this->isBoss() && !me ) {
                printf("bfmbase::CGNE_prec: true residual is %le \n",true_residual);
            }
            break;
        }
    }
    if(eigcg->def_len>0 && eigcg->def_len<eigcg->get_max_def_len()) {
        threaded_free(invH);
    }

    //add vectors V to U (deflation space) when necessary (have more low modes coming)
    if(do_eigcg) {
        struct timeval addV_start,addV_end,addV_diff;
        gettimeofday(&addV_start,NULL);
        int *low = new int[2*eigcg->get_nev()];
        //add low modes from the lowest first
        //sort
        int s;
        for(s=0;s<2*eigcg->get_nev();s++)low[s]=s;
        for(int i=0;i<2*eigcg->get_nev()-1;i++) {
            for(int j=2*eigcg->get_nev()-1;j>i;j--) {
                if(M[low[j]]<M[low[j-1]]){s=low[j];low[j]=low[j-1];low[j-1]=s;}
            }
        }

        //add
        for(int i=0;i<2*eigcg->get_nev();i++) {
            if(this->isBoss() && this->verbose && !me ) {
                printf("eigen value %d is M[%d] = %e \n",i, low[i], M[low[i]]);
            }
            //update deflation space U from V
            if(M[low[i]]<eigcg->get_max_eig_cut() && M[low[i]]>1e-20 && eigcg->def_len<eigcg->get_max_def_len()) {
                //remember to set M[i>=rank] to zero

                //update H 
                //U[def_len]->CopyVec(V[low[i]],vec_len);
                //double *fvptr = (Float *)V[low[i]];
                //for(int ii=0;ii<vec_len;ii++){U[def_len][ii]=(float)(fvptr[ii]);}//fvptr[ii]=(Float)(U[def_len][ii]);//Make low accuracy
                this->copy((Fermion_t)eigcg->getU(eigcg->def_len), (Fermion_t)eigcg->getV(low[i]));
                this->copy(p,(Fermion_t)eigcg->getV(low[i]));
                this->Mprec(p,mp,tmp,DaggerNo);
                this->Mprec(mp,mmp,tmp,DaggerYes);
                //for(int j=0;j<def_len;j++)H[j*max_def_len+def_len]=U[j]->CompDotProductGlbSum(AU,vec_len);
                for(int j=0;j<eigcg->def_len;j++) {
                    complex<double> dt = this->inner((Fermion_t)(eigcg->getU(j)),mmp);
                    if(me==0)eigcg->H[j*eigcg->get_max_def_len()+eigcg->def_len]=dt;
                }
                if(me==0) {
                    for(int j=0;j<eigcg->def_len;j++) {
                        eigcg->H[eigcg->def_len*eigcg->get_max_def_len()+j]
                            =std::conj(eigcg->H[j*eigcg->get_max_def_len()+eigcg->def_len]);
                    }
                }
                complex<double> dt = this->inner((Fermion_t)(eigcg->getV(low[i])),mmp);
                if(me==0) {
                    eigcg->H[eigcg->def_len*eigcg->get_max_def_len()+eigcg->def_len] =  dt;
                    eigcg->def_len++;	
                }
                if(this->isBoss() && this->verbose && !me) {
                    printf("def_len = %d\n",eigcg->def_len);
                }
            }
        }
        delete [] low;
        if(me==0) {
            if(eigcg->def_len == eigcg->get_max_def_len())
                invert_H_matrix(eigcg->H,eigcg->get_max_def_len());
        }
        gettimeofday(&addV_end,NULL);
        timersub(&addV_end,&addV_start,&addV_diff);
        if(this->isBoss() && !me) {
            printf("bfmbase::Eig_CGNE_prec: time to add V to deflation space %d.%6.6d \n",
                   addV_diff.tv_sec, addV_diff.tv_usec);
        }
    }

    if(do_eigcg) {
        threaded_free(T);
        threaded_free(Tt);
        threaded_free(Y);
        threaded_free(Q);
        threaded_free(QZ);
        threaded_free(H);
        threaded_free(EIG);
        threaded_free(M);
        threaded_free(INDEX);
    }

    this->threadedFreeFermion(tmp);
    this->threadedFreeFermion(p);
    this->threadedFreeFermion(mp);
    this->threadedFreeFermion(mmp);
    this->threadedFreeFermion(mmp_prev);
    this->threadedFreeFermion(r);
#ifdef LIST_ENGINE
    if (this->list_engine ) L1P_PatternUnconfigure();
#endif

    if ( this->iter>=this->max_iter && this->isBoss() && !me ) {
        printf("bfmbase::CGNE_prec: CG not converged \n");
        return -1;
    }
    else {
        return this->iter;
    }
}

#ifndef BLOCK_SIZE
#define BLOCK_SIZE ((int) 4)
#endif
template<class Float> void basic_dgemm (const int KB, const int M, const int N, const int K, Float **A, const double *B, double *C, const int NB)
{
    int i, j, k;
    register const double *Bp;
    register Float **pA=A;
    switch(K)
        {
        case 1:
            for (i = 0; i < M; ++i) 
                {
                    Bp=B;
                    for (j = 0; j < N; ++j) 
                        {
                            *(C + i*NB + j) += pA[0][i]*(*(Bp++));
                        }
                }
            break;
        case 2:
            for (i = 0; i < M; ++i) 
                {
                    Bp=B;
                    for (j = 0; j < N; ++j) 
                        {
                            register double cij = *(C + i*NB + j);
                            cij += pA[0][i]*(*(Bp++));
                            cij += pA[1][i]*(*(Bp++));
                            *(C + i*NB + j) = cij;
                        }
                }
            break;
        case 3:
            for (i = 0; i < M; ++i) 
                {
                    Bp=B;
                    for (j = 0; j < N; ++j) 
                        {
                            register double cij = *(C + i*NB + j);
                            cij += pA[0][i]*(*(Bp++));
                            cij += pA[1][i]*(*(Bp++));
                            cij += pA[3][i]*(*(Bp++));
                            *(C + i*NB + j) = cij;
                        }
                }
            break;
        case 4:
            for (i = 0; i < M; ++i) 
                {
                    Bp=B;
                    for (j = 0; j < N; ++j) 
                        {
                            register double cij = *(C + i*NB + j);
                            cij += pA[0][i]*(*(Bp++));
                            cij += pA[1][i]*(*(Bp++));
                            cij += pA[2][i]*(*(Bp++));
                            cij += pA[3][i]*(*(Bp++));
                            *(C + i*NB + j) = cij;
                        }
                }
            break;
        default:
            break;
        }
}

template<class Float> void unrolled_dgemm (const int lda, Float **A, const double *B, double *C, const int N)
{
    const double * pB=B;
    register const Float A0_0 = A[0][0];
    register const double B0_0 = *(pB++);
    register const Float A0_1 = A[1][0];
    register const double B1_0 = *(pB++);
    register const Float A0_2 = A[2][0];
    register const double B2_0 = *(pB++);
    register const Float A0_3 = A[3][0];
    register const double B3_0 = *(pB++);
    register const Float A1_0 = A[0][1];
    register const double B0_1 = *(pB++);
    register const Float A1_1 = A[1][1];
    register const double B1_1 = *(pB++);
    register const Float A1_2 = A[2][1];
    register const double B2_1 = *(pB++);
    register const Float A1_3 = A[3][1];
    register const double B3_1 = *(pB++);
    register const Float A2_0 = A[0][2];
    register const double B0_2 = *(pB++);
    register const Float A2_1 = A[1][2];
    register const double B1_2 = *(pB++);
    register const Float A2_2 = A[2][2];
    register const double B2_2 = *(pB++);
    register const Float A2_3 = A[3][2];
    register const double B3_2 = *(pB++);
    register const Float A3_0 = A[0][3];
    register const double B0_3 = *(pB++);
    register const Float A3_1 = A[1][3];
    register const double B1_3 = *(pB++);
    register const Float A3_2 = A[2][3];
    register const double B2_3 = *(pB++);
    register const Float A3_3 = A[3][3];
    register const double B3_3 = *(pB++);

    register double c0_0 = *(C + 0*N + 0);
    register double c0_1 = *(C + 0*N + 1);
    register double c0_2 = *(C + 0*N + 2);
    register double c0_3 = *(C + 0*N + 3);
    register double c1_0 = *(C + 1*N + 0);
    register double c1_1 = *(C + 1*N + 1);
    register double c1_2 = *(C + 1*N + 2);
    register double c1_3 = *(C + 1*N + 3);
    register double c2_0 = *(C + 2*N + 0);
    register double c2_1 = *(C + 2*N + 1);
    register double c2_2 = *(C + 2*N + 2);
    register double c2_3 = *(C + 2*N + 3);
    register double c3_0 = *(C + 3*N + 0);
    register double c3_1 = *(C + 3*N + 1);
    register double c3_2 = *(C + 3*N + 2);
    register double c3_3 = *(C + 3*N + 3);

    c0_0 += A0_0 * B0_0;
    c0_0 += A0_1 * B1_0;
    c0_0 += A0_2 * B2_0;
    c0_0 += A0_3 * B3_0;
    *(C + 0*N + 0) = c0_0;

    c0_1 += A0_0 * B0_1;
    c0_1 += A0_1 * B1_1;
    c0_1 += A0_2 * B2_1;
    c0_1 += A0_3 * B3_1;
    *(C + 0*N + 1) = c0_1;

    c0_2 += A0_0 * B0_2;
    c0_2 += A0_1 * B1_2;
    c0_2 += A0_2 * B2_2;
    c0_2 += A0_3 * B3_2;
    *(C + 0*N + 2) = c0_2;

    c0_3 += A0_0 * B0_3;
    c0_3 += A0_1 * B1_3;
    c0_3 += A0_2 * B2_3;
    c0_3 += A0_3 * B3_3;
    *(C + 0*N + 3) = c0_3;

    c1_0 += A1_0 * B0_0;
    c1_0 += A1_1 * B1_0;
    c1_0 += A1_2 * B2_0;
    c1_0 += A1_3 * B3_0;
    *(C + 1*N + 0) = c1_0;

    c1_1 += A1_0 * B0_1;
    c1_1 += A1_1 * B1_1;
    c1_1 += A1_2 * B2_1;
    c1_1 += A1_3 * B3_1;
    *(C + 1*N + 1) = c1_1;

    c1_2 += A1_0 * B0_2;
    c1_2 += A1_1 * B1_2;
    c1_2 += A1_2 * B2_2;
    c1_2 += A1_3 * B3_2;
    *(C + 1*N + 2) = c1_2;

    c1_3 += A1_0 * B0_3;
    c1_3 += A1_1 * B1_3;
    c1_3 += A1_2 * B2_3;
    c1_3 += A1_3 * B3_3;
    *(C + 1*N + 3) = c1_3;

    c2_0 += A2_0 * B0_0;
    c2_0 += A2_1 * B1_0;
    c2_0 += A2_2 * B2_0;
    c2_0 += A2_3 * B3_0;
    *(C + 2*N + 0) = c2_0;

    c2_1 += A2_0 * B0_1;
    c2_1 += A2_1 * B1_1;
    c2_1 += A2_2 * B2_1;
    c2_1 += A2_3 * B3_1;
    *(C + 2*N + 1) = c2_1;

    c2_2 += A2_0 * B0_2;
    c2_2 += A2_1 * B1_2;
    c2_2 += A2_2 * B2_2;
    c2_2 += A2_3 * B3_2;
    *(C + 2*N + 2) = c2_2;

    c2_3 += A2_0 * B0_3;
    c2_3 += A2_1 * B1_3;
    c2_3 += A2_2 * B2_3;
    c2_3 += A2_3 * B3_3;
    *(C + 2*N + 3) = c2_3;

    c3_0 += A3_0 * B0_0;
    c3_0 += A3_1 * B1_0;
    c3_0 += A3_2 * B2_0;
    c3_0 += A3_3 * B3_0;
    *(C + 3*N + 0) = c3_0;

    c3_1 += A3_0 * B0_1;
    c3_1 += A3_1 * B1_1;
    c3_1 += A3_2 * B2_1;
    c3_1 += A3_3 * B3_1;
    *(C + 3*N + 1) = c3_1;

    c3_2 += A3_0 * B0_2;
    c3_2 += A3_1 * B1_2;
    c3_2 += A3_2 * B2_2;
    c3_2 += A3_3 * B3_2;
    *(C + 3*N + 2) = c3_2;

    c3_3 += A3_0 * B0_3;
    c3_3 += A3_1 * B1_3;
    c3_3 += A3_2 * B2_3;
    c3_3 += A3_3 * B3_3;
    *(C + 3*N + 3) = c3_3;
}

//A(M*K) B(K*N)= C(M*N)
template<class Float> void matrix_dgemm (const int M,const int N, const int K, Float **A, const double *B, double *C)
{
    const int m_cblocks = M / BLOCK_SIZE + (M%BLOCK_SIZE ? 1 : 0);
    const int k_cblocks = K / BLOCK_SIZE + (K%BLOCK_SIZE ? 1 : 0);
    const int n_cblocks = N / BLOCK_SIZE + (N%BLOCK_SIZE ? 1 : 0);
    double *pt;
    int ci, cj, ck;
    for (ci = 0; ci < m_cblocks; ++ci) {
        const int i = ci * BLOCK_SIZE;
        const int M_cblock = (i+BLOCK_SIZE>M? M-i:BLOCK_SIZE);
        for (cj = 0; cj < n_cblocks; ++cj) {
            const int j = cj * BLOCK_SIZE;
            const int N_cblock = (j+BLOCK_SIZE>N? N-j:BLOCK_SIZE);
            for (ck = 0; ck < k_cblocks; ++ck) {
                const int k = ck * BLOCK_SIZE;
                const int K_cblock = (k+BLOCK_SIZE>K? K-k:BLOCK_SIZE);
                if((M_cblock==BLOCK_SIZE)&&(N_cblock==BLOCK_SIZE)&&(K_cblock==BLOCK_SIZE))
                    {
                        //unrolled_dgemm (K, A + i*K + k, B + k + j*K, C + i*N + j, N);
                        unrolled_dgemm (K, A+k, B + k*BLOCK_SIZE+ j*K, C + i*N + j, N);
                    }
                else
                    {
                        //basic_dgemm (K,M_cblock, N_cblock, K_cblock, A + i*K + k, B + k + j*K, C + i*N + j,N);
                        basic_dgemm (K,M_cblock, N_cblock, K_cblock, A + k, B + k*N_cblock+ j*K, C + i*N + j,N);
                    }
            }
        }
        //replace A with C
        pt=C+i*N;
        for(ck=0;ck<M_cblock;ck++)
            for(cj=0;cj<N;cj++)
                A[cj][ck]=(*(pt++));

        for(cj=0;cj<K;cj++)A[cj]+=M_cblock;
    }
}

template<class Float>
void eigcg_vec_mult(Float* V, const int m, double *QZ, const int n, const int f_size_cb, const int nthread, const int me)
//QZ is saved in column major format. 
//perform V = V*QZ;
{
    //implementation 4
    //reorder QZ to 4x4 blocks
    Float **Vptr = new Float*[m];
    assert(f_size_cb%nthread==0);
    int each = f_size_cb/nthread;
    for(int i=0;i<m;i++)Vptr[i] = V+i*f_size_cb+me*each;

    double *aux = new double[m*n];
    double *pt=aux;
    const int BS = 4;
    const int m_cblocks = m / BS + (m%BS ? 1 : 0);
    const int n_cblocks = n / BS + (n%BS ? 1 : 0);
    for(int c=0;c<n_cblocks;++c)// row direction
        {
            const int i=c*BS;
            for(int r=0;r<m_cblocks;++r) //column direction
                {
                    const int j=r*BS;
                    for(int ci=i;(ci<i+BS) && (ci<n); ci++)
                        for(int rj=j;(rj<j+BS)&&(rj<m);rj++)
                            *(pt++)=QZ[ci*m+rj];
                }
        }
    const int BBSS=24;
    int xlen=BBSS*n;
    double *x = new double[BBSS*n];
    for(int block=0;block<each/BBSS;block++)
        {
            memset(x,0,xlen*sizeof(double));
            matrix_dgemm(BBSS,n,m,Vptr,aux,x);
        }
    delete [] aux;
    delete [] x;
    delete [] Vptr;
}

template<class Float>
void eigcg_vec_mult2(Float* V, const int m, double *QZ, const int n, const int f_size_cb,
                     const int nthread, const int me,
                     bfm_evo<Float> &bfmobj)
//QZ is saved in column major format. 
//perform V = V*QZ;
{
    std::vector<Fermion_t> ret(n, NULL);
    for(int i = 0; i < n; ++i) {
        ret[i] = bfmobj.threadedAllocFermion();
    }

    for(int i = 0; i < n; ++i) {
        bfmobj.set_zero(ret[i]);
        for(int j = 0; j < m; ++j) {
            bfmobj.axpy(ret[i],
                        (Fermion_t)(V + j * f_size_cb), ret[i],
                        QZ[i * m + j]);
        }
    }

    for(int i = 0; i < n; ++i) {
        bfmobj.copy((Fermion_t)(V + i * f_size_cb), ret[i]);
        bfmobj.threadedFreeFermion(ret[i]);
    }
}

#ifndef ALIGNIT
#define ALIGNIT(A) (double *)( (((uint64_t)A) + 31)& (~0x1FUL) );
#endif

template <class Float>
void myaxpy(Float *r, Float *x, Float *y, double a, int len)
{
  double aa_b[4+4]; // 64 bytes
  double *aa = ALIGNIT(aa_b);
  for(int i=0;i<4;i++) aa[i] = a;

  if ( sizeof(Float) == sizeof(double) ) {
    vmx_vaxpy((double *)r,(double *)aa,(double *)x,(double *)y,len);
  } else { 
    vmx_vaxpy_s((float *)r,(double *)aa,(float *)x,(float *)y,len);
  }

  // thread_barrier();

  return;
}

template<class Float>
void eigcg_vec_mult3(Float* V, const int m, double *QZ, const int n, const int f_size_cb,
                     const int nthread, const int me,
                     bfm_evo<Float> &bfmobj)
//QZ is saved in column major format. 
//perform V = V*QZ;
{
    std::vector<Float *> ret(n, NULL);
    for(int i = 0; i < n; ++i) {
        ret[i] = (Float *)bfmobj.threadedAllocFermion();
        bfmobj.set_zero(ret[i]);
    }

    assert(nthread >= n);
    int mywork, myoff;
    bfmobj.thread_work_partial_nobarrier(f_size_cb / 24, me / n, nthread / n, mywork, myoff);
    myoff *= 24;

    if(me < nthread / n * n) {
        int i = me % n;
        for(int j = 0; j < m; ++j) {
            myaxpy(ret[i] + myoff, // &r
                   V + j * f_size_cb + myoff, // &x
                   ret[i] + myoff, // &y
                   QZ[i * m + j], // a
                   mywork); // length
        }
    }

    for(int i = 0; i < n; ++i) {
        bfmobj.copy((Fermion_t)(V + i * f_size_cb), ret[i]);
        bfmobj.threadedFreeFermion(ret[i]);
    }
}

#endif
