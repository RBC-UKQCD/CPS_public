//#define TEST

#ifdef TEST
#include<iostream>
#endif

#include <unistd.h>

#include <config.h>
#include <util/time_cps.h>
#include <util/dirac_op.h>
#include <util/lattice.h>
#include <util/smalloc.h>
#include <util/vector.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/error.h>
#include <util/checksum.h>
#include <comms/glb.h>
#include <math.h>
#include <stdio.h>
#include <comms/sysfunc_cps.h>
#include <cstring>
#include <util/vector_template_func.h>
CPS_START_NAMESPACE

extern "C" { 
	void invcg_r_norm(IFloat *resa, IFloat *scale, IFloat *mult, IFloat *add, int ncvec, IFloat *norm);
	void invcg_xp_update(IFloat *out1, IFloat *out2, IFloat *A, IFloat *B, IFloat *mult, IFloat *add, int size);
	void eigen_solver(IFloat *A, IFloat *EV, IFloat *E, int n);
	void matrix_dgemm (const int M,const int N, const int K, double **A, const double *B, double *C);
//	void matrix_dgemm (const int M,const int N, const int K, double *A, const double *B, double **C);
}

// The granularity used in the interleaving
#define GRAN 12
void min_eig_index(int *INDEX, int nev,Float *EIG, int n)
{
	//complexity=n*nev; can be reduced to n
	int max;
	for(int i=0;i<nev;i++)INDEX[i]=i;
	for(int i=nev;i<n;i++)
	{
		//find max in INDEX
		max=0;
		for(int j=1;j<nev;j++)if(EIG[INDEX[j]]>EIG[INDEX[max]])max=j;
		//change the max one
		if(EIG[i]<EIG[INDEX[max]])INDEX[max]=i;
	}
	
}

void invcg_copy_rnorm(Float *v, Float rsq, Float *X, int len)
{
	register Float a=1.0/sqrt(rsq);
	Float *vp=v;
	Float *rp=X;
	for(int i=0;i<len;i++)
	{
		*(vp++)=a*(*(rp++));
		*(vp++)=a*(*(rp++));
		*(vp++)=a*(*(rp++));
		*(vp++)=a*(*(rp++));
		*(vp++)=a*(*(rp++));
		*(vp++)=a*(*(rp++));
		*(vp++)=a*(*(rp++));
		*(vp++)=a*(*(rp++));
		*(vp++)=a*(*(rp++));
		*(vp++)=a*(*(rp++));
		*(vp++)=a*(*(rp++));
		*(vp++)=a*(*(rp++));
		rp+=12;
	}
}

void eigcg_vec_mult(Vector **V,int m, Float *QZ, int n, int f_size_cb)
//QZ is saved in column major format. 
{
	Float **Vptr=(Float* *)smalloc(m*sizeof(Float *));
	for(int i=0;i<m;i++)Vptr[i]=(Float *)V[i];

	//implementation 3. unroll. 96.1 seconds ~75Mflops
//	int BLOCK_SIZE=12;
//	int xlen=BLOCK_SIZE*n;
//	Float *x=(Float *)fmalloc(BLOCK_SIZE*n*sizeof(Float));
//	for(int block=0;block<f_size_cb/BLOCK_SIZE;block++)
//	{
//		std::memset(x,0,xlen*sizeof(Float));
//		matrix_dgemm(BLOCK_SIZE,n,m,Vptr,QZ,x);
//	}
//

	//implementation 4
	//reorder QZ to 4x4 blocks
	Float *aux=(Float*)fmalloc(m*n*sizeof(Float));
	Float *pt=aux;
   	int BLOCK_SIZE=4;
	const int m_cblocks = m / BLOCK_SIZE + (m%BLOCK_SIZE ? 1 : 0);
	const int n_cblocks = n / BLOCK_SIZE + (n%BLOCK_SIZE ? 1 : 0);
	for(int c=0;c<n_cblocks;++c)// row direction
	{
		const int i=c*BLOCK_SIZE;
		for(int r=0;r<m_cblocks;++r) //column direction
		{
			const int j=r*BLOCK_SIZE;
			for(int ci=i;(ci<i+BLOCK_SIZE) && (ci<n); ci++)
				for(int rj=j;(rj<j+BLOCK_SIZE)&&(rj<m);rj++)
					*(pt++)=QZ[ci*m+rj];
		}
	}
	BLOCK_SIZE=12;
	int xlen=BLOCK_SIZE*n;
	Float *x=(Float *)fmalloc(BLOCK_SIZE*n*sizeof(Float));
	for(int block=0;block<f_size_cb/BLOCK_SIZE;block++)
	{
		std::memset(x,0,xlen*sizeof(Float));
		matrix_dgemm(BLOCK_SIZE,n,m,Vptr,aux,x);
	}
	sfree(aux);

	sfree(x);
	sfree(Vptr);
}

//U saves the deflation space vectors. H=U^dag*A*U, invH=inv(H), def_len is the number of deflation space vectors
//set restart=0 will never restart, or else it will restart once when relative residule < restart
//return number of iterations to converge
//V returns the calculated eig vectors(length m, only use the first 2*nev), and M are the eigen values(length 2*nev). When we do deflation, we should only pick those has small eigen values to U
int DiracOp::InvEigCg(Vector *sol, Vector *src, Float *true_res, const int nev, const int m, Vector **V, const int vec_len, Float *M, float **U, Rcomplex *invH, const int def_len, const Float *restart, const int restart_len)
{
	
	char *fname = "InvEigCg(V*,V*,F,F*)";
	VRB.Func(cname,fname);

	if(nev>0 && m<=2*nev)ERR.General(cname,fname,"m should larger than 2*nev\n");

	int f_size_cb;     // Node checkerboard size of the fermion field
	// Set the node checkerboard size of the fermion field
	//------------------------------------------------------------------
	if(lat.Fclass() == F_CLASS_CLOVER) {
		f_size_cb = GJP.VolNodeSites() * lat.FsiteSize() / 2;
	 } else {
		f_size_cb = GJP.VolNodeSites() * lat.FsiteSize() / (lat.FchkbEvl()+1);
	 }

	if(vec_len!=f_size_cb)ERR.General(cname,fname,"vector length V does not match the length of solution and src vectors!\n");


	int iter=0; //Current number of CG iterations
	int max_iter=dirac_arg->max_num_iter; //max iteration number

	if (f_size_cb % GRAN != 0)ERR.General(cname,fname,"Field length %d is not a multiple of granularity %d\n", GRAN, f_size_cb);

	//calculate source norm
	Float src_norm_sq = src->NormSqNode(f_size_cb);
    DiracOpGlbSum(&src_norm_sq);

	VRB.Flow(cname,fname, "nev = %d, m= %d\n", nev, m);
	VRB.Flow(cname,fname, "Deflation length = %d \n", def_len);
	for(int i=0;i<restart_len;i++)VRB.Flow(cname,fname, "restart condition = %e\n", restart[i]);

	Float stp_cnd = src_norm_sq * dirac_arg->stop_rsd * dirac_arg->stop_rsd;
	VRB.Flow(cname,fname, "stp_cnd =%e\n", IFloat(stp_cnd));

	// Allocate memory for the solution/residual field.
	//------------------------------------------------------------------
	IFloat *X = (IFloat *) fmalloc(cname,fname,"X",2*f_size_cb * sizeof(Float));

	// Allocate memory for the direction vector dir.
	//------------------------------------------------------------------
	Vector *dir;
	if(GJP.VolNodeSites() >4096) 
		dir = (Vector *) smalloc(cname,fname,"dir",f_size_cb * sizeof(Float));
	else dir = (Vector *) fmalloc(cname,fname,"dir",f_size_cb * sizeof(Float));

	// Allocate mem. for the result vector of matrix multiplication mmp.
	//------------------------------------------------------------------
	Vector *mmp;
	if(GJP.VolNodeSites() >4096) 
		mmp = (Vector *) smalloc(cname,fname,"mmp",f_size_cb * sizeof(Float));
	else mmp = (Vector *) fmalloc(cname,fname,"mmp",f_size_cb * sizeof(Float));

	Vector *mmp_prev=NULL;
	//eigCG part
	if(nev>0)
	{
		mmp_prev = (Vector *) smalloc(cname,fname,"mmp",f_size_cb * sizeof(Float));
	}
	//eigCG end

	Rcomplex *Ub=NULL,*invHUb=NULL;
	if(def_len>0)
	{
		Ub=(Rcomplex *) fmalloc(cname,fname,"Ub",2*def_len*sizeof(Float));
		invHUb=(Rcomplex *) fmalloc(cname,fname,"invHUb",2*def_len*sizeof(Float));
	}


	//initial solution guess
	sol->VecZero(f_size_cb);
	if(def_len>0)
	{
		//sol = U*invH*U^dag*src;
		for(int ii=0;ii<def_len;ii++)
		{
			//in CPS, dot(a,b)=a^* * b
			//Ub[ii]=U[ii]->CompDotProductGlbSum(src,f_size_cb);
			//NOTICE!!! Should be improved to do a single glb_sum for all Ub
			//use function CompDotProductNode, after loop, do glb_sum(Ub,2*def_len)
			Float c_r, c_i;
			compDotProduct<float, Float>(&c_r, &c_i, U[ii], (Float *)src,f_size_cb);
			glb_sum_five(&c_r);
			glb_sum_five(&c_i);
			Ub[ii]=Complex(c_r,c_i);
		}
		int index=0;
		for(int ii=0;ii<def_len;ii++)
		{
			invHUb[ii]=0.0;
			for(int jj=0;jj<def_len;jj++)
				invHUb[ii]+=invH[index++]*Ub[jj];
		}
		for(int ii=0;ii<def_len;ii++)
		{
			//sol->CTimesV1PlusV2(invHUb[ii],U[ii],sol,f_size_cb);
			cTimesV1PlusV2<Float,float,Float>((Float *)sol, real(invHUb[ii]), imag(invHUb[ii]), U[ii],(Float *)sol,f_size_cb);

		}
#ifdef TEST
		for(int ii=0;ii<def_len;ii++)
		{
			Float xx=U[ii]->NormSqNode(f_size_cb);
			DiracOpGlbSum(&xx);
			std::cout<<"U[ii] norm = "<<xx<<std::endl;
		}
		std::cout<<"Ub vector"<<std::endl;
		for(int i=0;i<def_len;i++)
			std::cout<<Ub[i].real()<<'\t';
		std::cout<<std::endl;
		std::cout<<"inv(H) matrix:"<<std::endl;
		for(int i=0;i<def_len;i++)
		{
			for(int j=0;j<def_len;j++)
				std::cout<<invH[i*def_len+j].real()<<'\t';
			std::cout<<std::endl;
		}
		Float xx=sol->NormSqNode(f_size_cb);
		DiracOpGlbSum(&xx);
		std::cout<<"inital guess norm = "<<xx<<std::endl;
#endif
	}

	//dir = res(part of X vector) = src - MatPcDagMatPc * sol
    MatPcDagMatPc(mmp, sol);
    dir->CopyVec(src, f_size_cb);
    dir->VecMinusEquVec(mmp, f_size_cb);

	//aux pointers
    IFloat *Fsol = (IFloat*)sol;
    IFloat *Fdir = (IFloat*)dir;
    IFloat *Fmmp = (IFloat*)mmp;

    // Interleave solution and residual
    IFloat *Xptr = X;
    for (int j=0; j<f_size_cb/GRAN;j++) {
      for (int i=0; i<GRAN; i++) *Xptr++ = *(Fsol+j*GRAN+i); //initial solution
      for (int i=0; i<GRAN; i++) *Xptr++ = *(Fdir+j*GRAN+i); //residule
    }

	Float res_norm_sq_prv,res_norm_sq_cur;
	Float alpha,beta,pAp;

    res_norm_sq_cur = dir->NormSqNode(f_size_cb);
    DiracOpGlbSum(&res_norm_sq_cur);
	VRB.Flow(cname,fname, "|res[initial]|^2 = %e\n", IFloat(res_norm_sq_cur));

	sync();

	//eigCG part 
	int i_eig=0;
	alpha=1.0; //avoid 0/0 at the first eig iteration 
	beta=0.0;
	Float alpha_old;
	Float beta_old;
	Float *T=NULL;
	Float *Tt=NULL;
	Float *Q=NULL;
	Float *QZ=NULL;
	Float *H=NULL;
	Float *Y=NULL;
	Float *EIG=NULL;
	int *INDEX=NULL;
	if(nev>0)
	{
		T=(Float *) fmalloc(cname,fname,"T", m*(m+1)/2*sizeof(Float));// mxm real symetric matrix, lower triangular only
		Tt=(Float *) fmalloc(cname,fname,"T", m*(m+1)/2*sizeof(Float));// (m)x(m) real symetric matrix, lower triangular only
		//T(i,j) = T[ i(i+1)/2+j]
		//NOTICE!! T is actually very sparse matrix, so it can be speed up by a very large factor
		for(int i=0;i<m*(m+1)/2;i++)T[i]=0.0;
		Y=(Float *) fmalloc(cname,fname,"Y", m*m*sizeof(Float));//m*m real matrix, to store eigen vectors when needed
		Q=(Float *) fmalloc(cname,fname,"Q", m*2*nev*sizeof(Float));//m*(2nev) real matrix
		QZ=(Float *) fmalloc(cname,fname,"QZ", m*2*nev*sizeof(Float));//m*(2nev) real matrix
		H=(Float *) fmalloc(cname,fname,"H", 2*nev*(2*nev+1)/2*sizeof(Float));//(2nev)*(2nev) real matrix, symetric, lower part
		EIG=(Float *)fmalloc(cname,fname,"EIG",m*sizeof(Float));

		INDEX=(int *)fmalloc(cname,fname,"INDEX",nev*sizeof(int));//this vector is not necessay if we have a good eig solver for nev low
	}
	//
	//eigCG part end

	int restarted=0;
	Float *restartcond;
	if(restart_len>0)
	{
		restartcond=(Float *)fmalloc(cname,fname,"restartcond",restart_len*sizeof(Float));
		for(int i=0;i<restart_len;i++)restartcond[i]=restart[i]*restart[i]*src_norm_sq;
	}
	Float eigTotal=0.0;
	Float eigProj=0.0;
	Float defTime=0.0;
	Float total_time=0.0;
	total_time-=dclock();
    Float linalg_flops = 0;
    Float eigProj_flops = 0;
	Float linalg_time=0;
	CGflops=0;
	for(iter=0;iter<max_iter;iter++)
	{
		//eigCG part
		if(nev>0 && i_eig==m)mmp_prev->CopyVec(mmp,f_size_cb);
		//eig CG end
		
		res_norm_sq_prv = res_norm_sq_cur;
		MatPcDagMatPc(mmp,dir,&pAp);
		DiracOpGlbSum(&pAp);
		if(pAp==0)break;

		if(nev>0)
		{
			eigTotal-=dclock();
			int nev2=2*nev;
			//T(i_eig-1,i_eig-1)
			if(i_eig!=0)T[(i_eig-1)*i_eig/2+i_eig-1]=1.0/alpha+beta_old/alpha_old;
			if(i_eig==m)
			{
				//Yb need lowest nev eigen vector of T(m-1)
				//Y need lowest nev eigen vector of T(m)
#ifdef TEST
				std::cout<<" T matrix: "<<std::endl;
				for(int i=0;i<m;i++)
				{
					for(int j=0;j<m;j++)
					{
						if(i>=j)std::cout<<T[i*(i+1)/2+j]<<'\t';
						else std::cout<<T[j*(j+1)/2+i]<<'\t';
					}
					std::cout<<std::endl;
				}
#endif
				for(int i=0;i<m*(m-1)/2;i++)Tt[i]=T[i];
				eigen_solver(Tt,Y,EIG,m-1);//NOTICE, this is NOT needed, only need to calculate the lowest nev, not all m >2*nev
				//NOTICE, Y transpose is the eigen vectors.
				min_eig_index(INDEX,nev,EIG,m-1);
				//Q(nev:2nev-1)=Yb; with Yb last row zero
				for(int i=nev;i<2*nev;i++)
				{
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
				for(int i=0;i<nev;i++)
				{
					//for(int j=0;j<m;j++)Q[j*2*nev+i]=Y[j+INDEX[i]*m];
					//Y is eigen vector
					for(int j=0;j<m;j++)Q[j*2*nev+i]=Y[j*m+INDEX[i]];
				}
				
				//Q=orth([Y,Yb]); with YB last row zero
				//rank(Q) may be smaller than 2*nev. remove these
				//should be optimized. maybe save row first for Q
				int rank=nev;
				for(int i=nev;i<2*nev;i++)
				{
					for(int j=0;j<rank;j++)
					{
						Float xy=0.0;
						for(int k=0;k<m;k++)xy+=Q[k*2*nev+i]*Q[k*2*nev+j];
						for(int k=0;k<m;k++)Q[k*2*nev+i]-=xy*Q[k*2*nev+j];
					}
					//normalize
					Float xx=0.0;
					for(int k=0;k<m;k++)xx+=Q[k*2*nev+i]*Q[k*2*nev+i];
					if(xx>1e-16)
					{
						xx=1.0/sqrt(xx);
						for(int k=0;k<m;k++)Q[k*2*nev+rank]=xx*Q[k*2*nev+i];
						rank++;
					}
				}
				VRB.Flow(cname,fname,"Rank of Q = %d\n",rank);
				//H=Q' * T * Q
				for(int i=0;i<rank;i++)
					for(int j=0;j<=i;j++)
					{
						H[i*(i+1)/2+j]=0.0;
						for(int l=0;l<m;l++)
							for(int k=0;k<m;k++)
							{
								if(k<=l)H[i*(i+1)/2+j]+=Q[l*nev2+i]*T[l*(l+1)/2+k]*Q[k*nev2+j];
								else H[i*(i+1)/2+j]+=Q[l*nev2+i]*T[k*(k+1)/2+l]*Q[k*nev2+j];
							}
					}

#ifdef TEST
				std::cout<<"H matrix:"<<std::endl;
				for(int i=0;i<rank;i++)
				{
					for(int j=0;j<rank;j++)
					{
						if(i>=j)std::cout<<H[i*(i+1)/2+j]<<'\t';
						else std::cout<<H[j*(j+1)/2+i]<<'\t';
					}
					std::cout<<std::endl;
				}
#endif
				//[Z,M]=eig(H)
				eigen_solver(H,Y,M,rank);
				for(int i=rank;i<2*nev;i++)M[i]=0.0;//set M[i>=rank] to zero
				//NOtice, Y transpose is eigenvectos.
				for(int i=0;i<rank;i++)VRB.Flow(cname,fname,"eig %d : %e \n",i,M[i]);
				
				//V=V*(Q*Z)
				//1.QZ=Q*Z
				//transpoze QZ here to speed up the later calculation
				for(int j=0;j<rank;j++)
				for(int i=0;i<m;i++)
				{
					QZ[i+m*j]=0.0;
					//for(int k=0;k<rank;k++)QZ[i+m*j]+=Q[i*nev2+k]*Y[j*rank+k];
					for(int k=0;k<rank;k++)QZ[i+m*j]+=Q[i*nev2+k]*Y[j+k*rank];
				}

#ifdef TEST
				std::cout<<"QZ matrix:"<<std::endl;
				for(int i=0;i<m;i++)
				{
					for(int j=0;j<rank;j++)
					{
						std::cout<<QZ[i+j*m]<<'\t';
					}
					std::cout<<std::endl;
				}
#endif
				eigProj-=dclock();
				//2.V=V*QZ need to be implement very efficiently
				//. The way we save QZ is transpozed to column first 
				eigcg_vec_mult(V,m,QZ,rank,f_size_cb);
				eigProj+=dclock();
				eigProj_flops += 2*f_size_cb*m*rank;

				//T=M
				for(int i=0;i<m*(m+1)/2;i++)T[i]=0.0;
				for(int i=0;i<rank;i++)T[i*(i+1)/2+i]=M[i];
				i_eig=rank;
				//w=mmp-beta*mmp_prev
				mmp_prev->FTimesV1PlusV2(-beta,mmp_prev,mmp,f_size_cb);
				//T(i_eig+1,1:i_eig)=w^dag * V/sqrt(rsq) 
				//T is symmetric and REAL !! TESTED
				for(int ii=0;ii<i_eig;ii++)
				{
					//T(i_eig,ii)
					T[i_eig*(i_eig+1)/2+ii]=mmp_prev->ReDotProductGlbSum(V[ii],f_size_cb); //again these global sum can be donce by once to speed up
					T[i_eig*(i_eig+1)/2+ii]/=sqrt(res_norm_sq_cur);
				}
			}
			else
			{
				if(i_eig!=0)
				{
					//T(i_eig,i_eig-1)
					T[i_eig*(i_eig+1)/2+i_eig-1]=-sqrt(beta)/alpha;
				}
			}
			//V[i_eig]=r/sqrt(rsq);
			Float *vptr = (Float *)V[i_eig];
			invcg_copy_rnorm(vptr, res_norm_sq_cur, X+GRAN, f_size_cb/GRAN);
			i_eig++;
			eigTotal+=dclock();
		}
		
		alpha_old=alpha; //eigCG part
		alpha = -res_norm_sq_prv/pAp;
		// res = - alpha * (MatPcDagMatPc * dir) + res;
		// res_norm_sq_cur = res * res

		//test
		//Float test=((Vector *)X)->NormSqGlbSum(f_size_cb*2);
		//VRB.Flow(cname,fname,"X norm=%e \n",test);
		//test=mmp->NormSqGlbSum(f_size_cb);
		//VRB.Flow(cname,fname,"mmp norm=%e \n",test);

		linalg_time-=dclock();
		invcg_r_norm(X+GRAN, &alpha, Fmmp, X+GRAN, f_size_cb/GRAN, &res_norm_sq_cur);
		linalg_time+=dclock();
		DiracOpGlbSum(&res_norm_sq_cur);
		linalg_flops+=f_size_cb*4;
		alpha = -alpha;
		beta_old=beta; //eigCG part
		beta = res_norm_sq_cur / res_norm_sq_prv;
		//VRB.Flow(cname,fname,"a=%e, b=%e, pAp=%e \n",alpha,beta,pAp);
		// sol = alpha * dir + sol;
		// dir = beta * dir + res;
		linalg_time-=dclock();
		invcg_xp_update(X, Fdir, &alpha, &beta, Fdir, X, f_size_cb/GRAN);
		linalg_time+=dclock();
		linalg_flops+=f_size_cb*4;
		
		//consider restarting the init-CG once
		if(restarted<restart_len && def_len>0 && res_norm_sq_cur<restartcond[restarted])
		{
			defTime-=dclock();
			restarted++;
			VRB.Flow(cname,fname,"eigCG restarted at res_norm_sq_cur= %e\n",res_norm_sq_cur);
			//reuse dir vector as the initial guess vector of A*e=r ,with dir=e=U*invH*U^dag*r
			//reuse mmp for r in X
			//use sol for x in X
			Xptr = X;
			for (int j=0; j<f_size_cb/GRAN; j++) {
			  for (int i=0; i<GRAN; i++) *(Fsol+j*GRAN+i)=*Xptr++; // solution
			  for (int i=0; i<GRAN; i++) *(Fmmp+j*GRAN+i)=*Xptr++; // residule
			}

			for(int ii=0;ii<def_len;ii++)
			{
				//Ub[ii]=U[ii]->CompDotProductGlbSum(mmp,f_size_cb);
				//NOTICE!!! Should be improved to a single glb_sum for all Ub
				Float c_r, c_i;
				compDotProduct<float, Float>(&c_r, &c_i, U[ii], (Float *)mmp,f_size_cb);
				glb_sum_five(&c_r);
				glb_sum_five(&c_i);
				Ub[ii]=Complex(c_r,c_i);
			}
			int index=0;
			for(int ii=0;ii<def_len;ii++)
			{
				invHUb[ii]=0.0;
				for(int jj=0;jj<def_len;jj++)
					invHUb[ii]+=invH[index++]*Ub[jj];
			}
			dir->VecZero(f_size_cb);
			for(int ii=0;ii<def_len;ii++)
			{
				//dir->CTimesV1PlusV2(invHUb[ii],U[ii],dir,f_size_cb);
				cTimesV1PlusV2<Float,float,Float>((Float *)dir, real(invHUb[ii]), imag(invHUb[ii]), U[ii],(Float *)dir,f_size_cb);
			}

			sol->VecAddEquVec(dir,f_size_cb); //get new sol
			//reuse the first half of X to save M^dag*M*e
			Vector *PAe=(Vector *)X;

			MatPcDagMatPc(PAe,dir);
			mmp->VecMinusEquVec(PAe,f_size_cb); //get new res

			res_norm_sq_cur = mmp->NormSqNode(f_size_cb);
			DiracOpGlbSum(&res_norm_sq_cur);

			dir->CopyVec(mmp,f_size_cb);

			Xptr = X;
			for (int j=0; j<f_size_cb/GRAN;j++) {
			  for (int i=0; i<GRAN; i++) *Xptr++ = *(Fsol+j*GRAN+i); //new initial solution
			  for (int i=0; i<GRAN; i++) *Xptr++ = *(Fmmp+j*GRAN+i); //new residule
			}

			defTime+=dclock();
		}
		
		VRB.Flow(cname,fname, "|res[%d]|^2 = %e\n", iter, IFloat(res_norm_sq_cur));
		if(res_norm_sq_cur <= stp_cnd) break;
	}
	total_time+=dclock();
	
	VRB.Result(cname,fname,"1. Time on CG : %e  seconds in %e flops\n", total_time-eigTotal-defTime,((Float)CGflops+linalg_flops)/(total_time-eigTotal-defTime));
	VRB.Result(cname,fname,"1.x CG linear algebra : %e flops / %e seconds = %e flops\n", linalg_flops, linalg_time, linalg_flops/linalg_time);
	VRB.Result(cname,fname,"2. Total Time on eig part: %e seconds \n", eigTotal);
	if(nev>0)VRB.Result(cname,fname,"2.x projection part of eig part : %e flops / %e seconds = %e flops\n", eigProj_flops, eigProj, eigProj_flops/eigProj);
	if(def_len>0)VRB.Result(cname,fname,"3. deflation(restart) time : %e seconds\n", defTime);


	if(iter == max_iter)
		VRB.Warn(cname,fname, "CG reached max iterations = %d. |res|^2 = %e\n", iter+1, IFloat(res_norm_sq_cur) );

    Xptr = X-GRAN;
    for (int j=0; j<f_size_cb; j++) {
      if (j%GRAN==0) Xptr += GRAN;
      *(Fsol++) = *(Xptr++);
    }
    MatPcDagMatPc(mmp, sol);
    dir->CopyVec(src, f_size_cb);
    dir->VecMinusEquVec(mmp, f_size_cb);
    res_norm_sq_cur = dir->NormSqNode(f_size_cb);
    DiracOpGlbSum(&res_norm_sq_cur);
	*true_res = res_norm_sq_cur/src_norm_sq;
	*true_res = sqrt(*true_res);
    VRB.Result(cname,fname, "True |res| / |src| = %e, iter = %d\n", IFloat(*true_res), iter);

	// Free memory
	sfree(cname,fname, "mmp", mmp);
	sfree(cname,fname, "dir", dir);
	sfree(cname,fname, "X", X);
	if(def_len>0)
	{
		sfree(cname,fname, "Ub", Ub);
		sfree(cname,fname, "invHUb", invHUb);
	}
	//eigCG part
	if(nev>0)
	{
		sfree(cname,fname,"mmp_prev",mmp_prev);
	}
	//eigCG part end
	if(restart_len>0)sfree(cname,fname,"restartcond",restartcond);

	sync();

	return iter;
}

CPS_END_NAMESPACE
