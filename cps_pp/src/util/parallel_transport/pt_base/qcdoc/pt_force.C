
/*!\file
  \brief  Implementation of Fasqtad::EvolveMomFforce.

  $Id: pt_force.C,v 1.7 2012-05-10 05:51:23 chulwoo Exp $
*/
//--------------------------------------------------------------------

#include <string.h>
#include "asq_data_types.h"
#include "pt_int.h"

extern "C"{
void asq_force_cross2dag(PTvector *chi, PTvector *phi, matrix *result,
                      int counter, double *fac);
void asq_force_cross2dag_s(PTvector *chi, PTvector *phi, matrix *result,
                      int counter, float *fac);
void asq_vaxpy3(PTvector *res,Float *scale,PTvector *mult,PTvector *add, int ncvec);
void asq_vaxpy3_s(PTvector *res,Float *scale,PTvector *mult,PTvector *add, int ncvec);
}

void *pt_alloc(size_t request, const char *cname, const
char *fname, const char *vname){
    return qalloc(QCOMMS,request);
}

#ifdef ASQD_SINGLE
#define asq_force_cross2dag(A,B,C,D,E) asq_force_cross2dag_s(A,B,C,D,E)
#define asq_vaxpy3(A,B,C,D,E) asq_vaxpy3_s(A,B,C,D,E)
#endif



#undef PROFILE

// N.B. No optimising provision is made if any of the asqtad coefficients
// are zero.

void PT::asqtad_force(AsqDArg *asq_arg, matrix *mom, Float *X, Float dt){

    char *fname = "asqtad_force()";
//    VRB.Func(cname,fname);

#ifdef PROFILE    
    Flops = 0;
    Float dtime = -dclock();
//    struct timeval start,end;
//    gettimeofday(&start, NULL);
#endif

    const int VAXPY_UNROLL = 6; // For optimised linear algebra functions.
    const int vax_len = vol*6/VAXPY_UNROLL;
    
    // The number of directions to do simultaneously; N can be 1, 2 or 4.
    int N = 4;  
    if (vol>512) N = 1;       // This is to fit the vectors 
//    else if (vol>256) N = 2;  // into EDRAM on QCDOC.

//    VRB.Flow(cname,fname,"N=%d\n",N);

    enum{plus=0, minus=1, n_sign=2};

    // Array in which to accumulate the force term:
    // this must be initialised to zero 

#if 0
    matrix **force = (matrix**)pt_amalloc(pt_alloc, sizeof(matrix), 2, 4, vol);
    memset( (char *)force,0,sizeof(matrix)*4*vol);
#endif
    matrix *force[4];
    for(int i =0;i<4;i++){
      force[i] = (matrix *)FastAlloc(sizeof(matrix)*vol);
      memset( (char *)force[i],0,sizeof(matrix)*vol);
    }


#if 0
    // vector arrays for which we must allocate memory

    PTvector ***Pnu = (PTvector***)pt_amalloc(pt_alloc, sizeof(PTvector), 3, n_sign, N, vol);
    
    PTvector ****P3 = (PTvector****)pt_amalloc(pt_alloc, sizeof(PTvector), 4, n_sign, n_sign, N, vol);

    PTvector ****Prhonu = (PTvector****)pt_amalloc(pt_alloc, sizeof(PTvector), 4, n_sign, n_sign, N, vol);

    PTvector *****P5 = (PTvector*****)pt_amalloc(pt_alloc, sizeof(PTvector), 5, n_sign, n_sign, n_sign, N, vol);

    PTvector ******P7 = (PTvector******)pt_amalloc(pt_alloc, sizeof(PTvector), 6, n_sign, n_sign, n_sign, n_sign, N, vol);

    PTvector ******Psigma7 = (PTvector******)pt_amalloc(pt_alloc, sizeof(PTvector), 6, n_sign, n_sign, n_sign, n_sign, N, vol);
#else
    PTvector *Pnu[2][4];
    PTvector *P3[2][2][4];
    PTvector *Prhonu[2][2][4];
    PTvector *P5[2][2][2][4];
    PTvector *P7[2][2][2][2][4];
    PTvector *Psigma7[2][2][2][2][4];
    PTvector *Pnununu[4];
    PTvector *Pnunu[2][4];
    PTvector *Pnu5[2][2][4];
    PTvector *Pnu3[2][2][4];
    PTvector *Prho5[2][2][2][4];
    PTvector *Psigmarhonu[2][2][2][4];

    size_t vec_size = 6*vol*sizeof(Float);
    for(int i = 0;i<N;i++){
      for(int j = 0;j<n_sign;j++){
        Pnu[j][i] = (PTvector *)Alloc(cname,fname,"Pnu",vec_size,QCOMMS);
//        printf("Pnu[%d][%d]=%p\n",j,i,Pnu[j][i]);
        for(int k = 0;k<n_sign;k++){
          P3[k][j][i] = (PTvector *)Alloc(cname,fname,"Pnu",vec_size,QCOMMS);
          Prhonu[k][j][i] = (PTvector *)Alloc(cname,fname,"Pnu",vec_size,QCOMMS);
          for(int l = 0;l<n_sign;l++){
            P5[l][k][j][i] = (PTvector *)FastAlloc(vec_size);
            for(int m = 0;m<n_sign;m++){
              P7[m][l][k][j][i] = (PTvector *)FastAlloc(vec_size);
              Psigma7[m][l][k][j][i] = (PTvector *)FastAlloc(vec_size);
            }
            Prho5[l][k][j][i] = Psigma7[0][l][k][j][i];
            Psigmarhonu[l][k][j][i] = Psigma7[0][l][k][j][i];
          }
          Pnu3[k][j][i] = P7[0][0][k][j][i];
          Pnu5[k][j][i] = P7[0][0][k][j][i];
        }
        Pnunu[j][i] = Psigma7[0][0][0][j][i];
      }
      Pnununu[i] = Prhonu[0][0][i];
    }
#endif


    
#if 0
    // These vectors can be overlapped with previously allocated memory
    
    PTvector *Pnununu[4] = Prhonu[0][0];
    PTvector *Pnunu[2][4] = Psigma7[0][0][0];;
    PTvector ****Pnu5 = P7[0][0];
    PTvector ****Pnu3 = P7[0][0];
    PTvector *****Prho5 = Psigma7[0];
    PTvector *****Psigmarhonu = Psigma7[0];
#endif


    // input/output arrays for the parallel transport routines
    Float *vin[n_sign*4], *vout[n_sign*4];
    int dir[n_sign*4];
	
    int mu[4], nu[4], rho[4], sigma[4];   // Sets of directions
    int w;                                // The direction index 0...N-1
    int ms, ns, rs, ss;                   // Sign of direction
    bool done[4] = {false,false,false,false};  // Flags to tell us which 
                                               // nu directions we have done.
//    ParTransAsqtad parallel_transport(*this);
	    
#if 0
   dtime = -dclock();
#endif
    for (int m=0; m<4; m+=N){                     	    // Loop over mu
	for(w=0; w<N; w++) mu[w] = (m+w)%4; 

	for (int n=m+1; n<m+4; n++){                        // Loop over nu
	    for(w=0; w<N; w++) nu[w] = (n+w)%4;

	    // Pnu = U_nu X

	    for(int i=0; i<N; i++){
		vin[i] = vin[i+N] = X;
		dir[n_sign*i] = n_sign*nu[i]+plus;        // nu_i
		dir[n_sign*i+1] = n_sign*nu[i]+minus;    // -nu_i
		vout[n_sign*i] = Pnu[minus][i];
		vout[n_sign*i+1] = Pnu[plus][i];
	    }
	    vec(n_sign*N, vout, vin, dir);

	    // P3 = U_mu Pnu
	    // ms is the nu sign index, ms is the mu sign index,
	    // w is the direction index
	    for(int i=0; i<N; i++){
		dir[n_sign*i] = n_sign*mu[i]+plus;        // mu_i
		dir[n_sign*i+1] = n_sign*mu[i]+minus;    // -mu_i
	    }
	    for(ns=0; ns<n_sign; ns++){               // ns is the sign of nu
		for(int i=0; i<N; i++){
		    vin[n_sign*i] = vin[n_sign*i+1] = Pnu[ns][i];
		    vout[n_sign*i] = P3[plus][ns][i];
		    vout[n_sign*i+1] = P3[minus][ns][i];
		}
		vec(n_sign*N, vout, vin, dir);
	    }
	    
	    for(w=0; w<N; w++)
		for(ns=0; ns<n_sign; ns++){
		    force_product_sum(P3[plus][ns][w], Pnu[ns][w],
				      -asq_arg->c3,
				      force[mu[w]]);
		}

	    for(int r=n+1; r<n+4; r++){                     // Loop over rho
		bool nextr = false;
		for(w=0; w<N; w++){
		    rho[w] = (r+w)%4;		
		    if(rho[w]==mu[w]){
			nextr = true;
			break;
		    }
		}
		if(nextr) continue;

		for(w=0; w<N; w++){                         // sigma
		    for(int s=rho[w]+1; s<rho[w]+4; s++){
			sigma[w] = s%4;
			if(sigma[w]!=mu[w] && sigma[w]!=nu[w]) break;
		    }
		}

		// Prhonu = U_rho Pnu 

		for(int i=0; i<N; i++){
		    dir[n_sign*i] = n_sign*rho[i]+plus;        
		    dir[n_sign*i+1] = n_sign*rho[i]+minus;    
		}
		for(ns=0; ns<n_sign; ns++){
		    for(int i=0; i<N; i++){
			vin[n_sign*i] = vin[n_sign*i+1] = Pnu[ns][i];
			vout[n_sign*i] = Prhonu[ns][minus][i];
			vout[n_sign*i+1] = Prhonu[ns][plus][i];
		    }
		    vec(n_sign*N, vout, vin, dir);
		}

		// P5 = U_mu Prhonu

		for(int i=0; i<N; i++){
		    dir[n_sign*i] = n_sign*mu[i]+plus;        
		    dir[n_sign*i+1] = n_sign*mu[i]+minus;    
		}
		for(ns=0; ns<n_sign; ns++) for(rs=0; rs<n_sign; rs++) {
		    for(int i=0; i<N; i++){
			vin[n_sign*i] = vin[n_sign*i+1] = Prhonu[ns][rs][i];
			vout[n_sign*i] = P5[plus][ns][rs][i];
			vout[n_sign*i+1] = P5[minus][ns][rs][i];
		    }
		    vec(n_sign*N, vout, vin, dir);
		}

		// F_mu += P5 Prhonu^dagger
		      
	        for(w=0; w<N; w++)
		    for(ns=0; ns<n_sign; ns++) for(rs=0; rs<n_sign; rs++)
			force_product_sum(P5[plus][ns][rs][w],
					  Prhonu[ns][rs][w],
					  asq_arg->c5,
					  force[mu[w]]);

		// Psigmarhonu = U_sigma P_rhonu
		
		for(int i=0; i<N; i++){
		    dir[n_sign*i] = (n_sign*sigma[i]);        
		    dir[n_sign*i+1] = (n_sign*sigma[i]+1);    
		}
		for(ns=0; ns<n_sign; ns++) for(rs=0; rs<n_sign; rs++){
		    for(int i=0; i<N; i++){
			vin[n_sign*i] = vin[n_sign*i+1] = Prhonu[ns][rs][i];
			vout[n_sign*i] = Psigmarhonu[ns][rs][minus][i];
			vout[n_sign*i+1] = Psigmarhonu[ns][rs][plus][i];
		    }
		    vec(n_sign*N, vout, vin, dir);
		}

		// P7 = U_mu P_sigmarhonu
		for(int i=0; i<N; i++){
		    dir[n_sign*i] = n_sign*mu[i]+plus;        
		    dir[n_sign*i+1] = n_sign*mu[i]+minus;    
		}
		for(ns=0; ns<n_sign; ns++) for(rs=0; rs<n_sign; rs++) for(ss=0; ss<n_sign; ss++){
		    for(int i=0; i<N; i++){
			vin[n_sign*i] = vin[n_sign*i+1] = Psigmarhonu[ns][rs][ss][i];
			vout[n_sign*i] = P7[plus][ns][rs][ss][i];
			vout[n_sign*i+1] = P7[minus][ns][rs][ss][i];
		    }
		    vec(n_sign*N, vout, vin, dir);
		}

		// F_mu -= P7 Psigmarhonu^\dagger
		
	        for(w=0; w<N; w++)
		    for(ns=0; ns<n_sign; ns++) for(rs=0; rs<n_sign; rs++) for(ss=0; ss<n_sign; ss++)
			force_product_sum(P7[plus][ns][rs][ss][w],
					  Psigmarhonu[ns][rs][ss][w],
					  -asq_arg->c7,
					  force[mu[w]]);

		// F_sigma += P7 Psigmarhonu^\dagger
		// N.B. this is the same as one of the previous products.
		
	        for(w=0; w<N; w++)
		    for(ns=0; ns<n_sign; ns++) for(rs=0; rs<n_sign; rs++) 
			force_product_sum(P7[plus][ns][rs][minus][w],
					  Psigmarhonu[ns][rs][minus][w],
					  asq_arg->c7,
					  force[sigma[w]]);

		// F_sigma += Psigmarhonu P7^\dagger
		
	        for(w=0; w<N; w++)
		    for(ns=0; ns<n_sign; ns++) for(rs=0; rs<n_sign; rs++) 
			force_product_sum(Psigmarhonu[ns][rs][minus][w],
					  P7[minus][ns][rs][minus][w],
					  asq_arg->c7,
					  force[sigma[w]]);

		// Psigma7 = U_sigma P7 
		for(int i=0; i<N; i++){
		    dir[n_sign*i] = (n_sign*sigma[i]);        
		    dir[n_sign*i+1] = (n_sign*sigma[i]+1);    
		}
		for(ms=0; ms<n_sign; ms++) for(ns=0; ns<n_sign; ns++) for(rs=0; rs<n_sign; rs++){
		    for(int i=0; i<N; i++){
			vin[n_sign*i] = P7[ms][ns][rs][plus][i];
			vin[n_sign*i+1] = P7[ms][ns][rs][minus][i];
			vout[n_sign*i] = Psigma7[ms][ns][rs][plus][i];
			vout[n_sign*i+1] = Psigma7[ms][ns][rs][minus][i];
		    }
		    vec(n_sign*N, vout, vin, dir);
		}

		// F_sigma += Fsigma7 Frhonu^\dagger

	        for(w=0; w<N; w++)
		    for(ns=0; ns<n_sign; ns++) for(rs=0; rs<n_sign; rs++) 
			force_product_sum(Psigma7[plus][ns][rs][plus][w],
					  Prhonu[ns][rs][w],
					  asq_arg->c7,
					  force[sigma[w]]);

		// F_sigma += Frhonu Fsigma7^\dagger

	        for(w=0; w<N; w++)
		    for(ns=0; ns<n_sign; ns++) for(rs=0; rs<n_sign; rs++) 
			force_product_sum(Prhonu[ns][rs][w],
					  Psigma7[minus][ns][rs][plus][w],
					  asq_arg->c7,
					  force[sigma[w]]);

		// P5 += c_7/c_5 Psigma7

		if(asq_arg->c5!=0.0){
		    Float c75 = asq_arg->c7/asq_arg->c5;
		    for(ms=0; ms<n_sign; ms++) for(ns=0; ns<n_sign; ns++) for(rs=0; rs<n_sign; rs++) for(ss=0; ss<n_sign; ss++) for(w=0; w<N; w++)
			
			asq_vaxpy3(P5[ms][ns][rs][w],&c75, Psigma7[ms][ns][rs][ss][w], P5[ms][ns][rs][w], vax_len);
			Flops += 2*vol*VECT_LEN*N*n_sign*n_sign*n_sign*n_sign;
		}
		// F_rho -= P5 Prhonu^\dagger
	        for(w=0; w<N; w++)
		    for(ns=0; ns<n_sign; ns++)
			force_product_sum(P5[plus][ns][minus][w],
					  Prhonu[ns][minus][w],
					  -asq_arg->c5,
					  force[rho[w]]);

		// F_rho -= Prhonu P5^\dagger
		    
	        for(w=0; w<N; w++)
		    for(ns=0; ns<n_sign; ns++)
			force_product_sum(Prhonu[ns][minus][w],
					  P5[minus][ns][minus][w],
					  -asq_arg->c5,
					  force[rho[w]]);

		// Prho5 = U_rho P5

		for(int i=0; i<N; i++){
		    dir[n_sign*i] = n_sign*rho[i]+plus;        
		    dir[n_sign*i+1] = n_sign*rho[i]+minus;    
		}
		for(ms=0; ms<n_sign; ms++) for(ns=0; ns<n_sign; ns++){
		    for(int i=0; i<N; i++){
			vin[n_sign*i] = P5[ms][ns][plus][i];
			vin[n_sign*i+1] = P5[ms][ns][minus][i];
			vout[n_sign*i] = Prho5[ms][ns][plus][i];
			vout[n_sign*i+1] = Prho5[ms][ns][minus][i];
		    }
		    vec(n_sign*N, vout, vin, dir);
		}

		// F_rho -= Prho5 Pnu^\dagger
		
	        for(w=0; w<N; w++)
		    for(ns=0; ns<n_sign; ns++)
			force_product_sum(Prho5[plus][ns][plus][w],
					  Pnu[ns][w],
					  -asq_arg->c5,
					  force[rho[w]]);

		// F_rho -= Pnu Prho5^\dagger
		
	        for(w=0; w<N; w++)
		    for(ns=0; ns<n_sign; ns++)
			force_product_sum(Pnu[ns][w],
					  Prho5[minus][ns][plus][w],
					  -asq_arg->c5,
					  force[rho[w]]);
		
		// P3 += c_5/c_3 Prho5

		if(asq_arg->c3!=0.0){
		    Float c53 = asq_arg->c5/asq_arg->c3;
		    for(ms=0; ms<n_sign; ms++) for(ns=0; ns<n_sign; ns++) for(rs=0; rs<n_sign; rs++) for(w=0; w<N; w++)
			asq_vaxpy3(P3[ms][ns][w],&c53,Prho5[ms][ns][rs][w], P3[ms][ns][w], vax_len);
		Flops += 2*vol*VECT_LEN*N*n_sign*n_sign*n_sign;
		}
		
	    } // rho+sigma loop

	    // Pnunu = U_nu Pnu

	    for(int i=0; i<N; i++){
		dir[n_sign*i] = n_sign*nu[i]+plus;        
		dir[n_sign*i+1] = n_sign*nu[i]+minus;    
	    }
	    for(int i=0; i<N; i++){
		vin[n_sign*i] = Pnu[minus][i];
		vin[n_sign*i+1] = Pnu[plus][i];
		vout[n_sign*i] = Pnunu[minus][i];
		vout[n_sign*i+1] = Pnunu[plus][i];
	    }
	    vec(n_sign*N, vout, vin, dir);

	    // P5 = U_mu Pnunu

	    for(int i=0; i<N; i++){
		dir[n_sign*i] = n_sign*mu[i]+plus;        
		dir[n_sign*i+1] = n_sign*mu[i]+minus;    
	    }
	    for(ns=0; ns<n_sign; ns++){
		for(int i=0; i<N; i++){
		    vin[n_sign*i] = Pnunu[ns][i];
		    vin[n_sign*i+1] = Pnunu[ns][i];
		    vout[n_sign*i] = P5[plus][ns][0][i];
		    vout[n_sign*i+1] = P5[minus][ns][0][i];
		}
		vec(n_sign*N, vout, vin, dir);
	    }

	    // F_mu += P5 Pnunu^\dagger

	    for(w=0; w<N; w++)
		for(ns=0; ns<n_sign; ns++)
		    force_product_sum(P5[plus][ns][0][w],
				      Pnunu[ns][w],
			          asq_arg->c6,
				      force[mu[w]]);

	    // F_nu -= P5 Pnunu^\dagger
	    // N.B. this is the same as one of the previous products
	    
	    for(w=0; w<N; w++)
		force_product_sum(P5[plus][minus][0][w],
				  Pnunu[minus][w],
			 	  -asq_arg->c6,
				  force[nu[w]]);
	    
	    // F_nu -= Pnunu P5^\dagger
	    
	    for(w=0; w<N; w++)
		force_product_sum(Pnunu[minus][w],
				  P5[minus][minus][0][w],
			 	  -asq_arg->c6,
				  force[nu[w]]);

	    // Pnu5 = U_nu P5

	    for(int i=0; i<N; i++){
		dir[n_sign*i] = n_sign*nu[i]+plus;        
		dir[n_sign*i+1] = n_sign*nu[i]+minus;    
	    }
	    for(ms=0; ms<n_sign; ms++){
		for(int i=0; i<N; i++){
		    vin[n_sign*i] =   P5[ms][plus][0][i]; 
		    vin[n_sign*i+1] = P5[ms][minus][0][i];
		    vout[n_sign*i] =   Pnu5[ms][plus][i];
		    vout[n_sign*i+1] = Pnu5[ms][minus][i];
		}
		vec(n_sign*N, vout, vin, dir);
	    }

	    // F_nu -= Pnu5 Pnu^\dagger

	    for(w=0; w<N; w++)
		force_product_sum(Pnu5[plus][plus][w],
				  Pnu[plus][w],
			 	  -asq_arg->c6,
				  force[nu[w]]);

	    // F_nu -= Pnu Pnu5^\dagger

	    for(w=0; w<N; w++)
		force_product_sum(Pnu[plus][w],
				  Pnu5[minus][plus][w],
			 	  -asq_arg->c6,
				  force[nu[w]]);

	    // P3 += c_L/c_3 Pnu5

	    if(asq_arg->c3!=0.0){
		    Float cl3 = asq_arg->c6/asq_arg->c3;
		for(ms=0; ms<n_sign; ms++) for(ns=0; ns<n_sign; ns++) for(w=0; w<N; w++)
		    asq_vaxpy3(P3[ms][ns][w], &cl3, Pnu5[ms][ns][w], P3[ms][ns][w], vax_len);
		Flops += 2*vol*VECT_LEN*N*n_sign*n_sign;
	    }
	    // F_nu += P3 Pnu^\dagger

	    for(w=0; w<N; w++)
		force_product_sum(P3[plus][minus][w],
				  Pnu[minus][w],
			 	  asq_arg->c3,
				  force[nu[w]]);

	    // F_nu +=  Pnu P3^\dagger

	    for(w=0; w<N; w++)
		force_product_sum(Pnu[minus][w],
				  P3[minus][minus][w],
			 	  asq_arg->c3,
				  force[nu[w]]);
	    
	    // Pnu3 = U_nu P3

	    for(int i=0; i<N; i++)
		dir[i] = n_sign*nu[i]+plus;        
	    for(ms=0; ms<n_sign; ms++){
		for(int i=0; i<N; i++){
		    vin[i] = P3[ms][plus][i]; 
		    vout[i] = Pnu3[ms][plus][i];
		}
		vec(N, vout, vin, dir);
	    }

	    // F_nu += Pnu3 X^\dagger

	    for(w=0; w<N; w++)
		force_product_sum(Pnu3[plus][plus][w], X,
			 	  asq_arg->c3,
				  force[nu[w]]);

	    // F_nu += X Pnu3^\dagger

	    for(w=0; w<N; w++)
		force_product_sum(X, Pnu3[minus][plus][w], 
			 	  asq_arg->c3,
				  force[nu[w]]);

	    // This stuff is to be done once only for each value of nu[w].
	    // Look for  N nu's that haven't been done before.

    	    bool nextn = false;
	    for(w=0; w<N; w++)
		if(done[nu[w]]){
		    nextn = true;
		    break;
		}
	    if(nextn) continue;
	    for(w=0; w<N; w++) done[nu[w]] = true;

	    // Got N new nu's, so do some stuff...
	    
	    // F_nu += Pnu X^\dagger

	    for(w=0; w<N; w++)
		force_product_sum(Pnu[minus][w], X,
			 	  asq_arg->c1,
				  force[nu[w]]);

	    // F_nu += Pnunu Pnu^\dagger

	    for(w=0; w<N; w++)
		force_product_sum(Pnunu[minus][w], Pnu[plus][w],
			 	  -asq_arg->c2,
				  force[nu[w]]);

	    // F_nu += Pnu Pnunu^\dagger
	    
	    for(w=0; w<N; w++)
		force_product_sum(Pnu[minus][w], Pnunu[plus][w],
			 	  asq_arg->c2,
				  force[nu[w]]);

	    // Pnununu = U_nu Pnunu

	    for(int i=0; i<N; i++){
		dir[i] = n_sign*nu[i]+plus;        
		vin[i] = Pnunu[minus][i]; 
		vout[i] = Pnununu[i];
	    }
	    vec(N, vout, vin, dir);
	    
	    // F_nu += Pnununu X^\dagger
		
	    for(w=0; w<N; w++)
		force_product_sum(Pnununu[w], X,
			 	  asq_arg->c2,
				  force[nu[w]]);
		
	    

	} // nu loop
    } // mu loop


    // Now that we have computed the force, we can update the momenta

    update_momenta(force, dt, mom);

#if 1
    for(int i = 0;i<N;i++){
      for(int j = 0;j<n_sign;j++){
        Free(Pnu[j][i]);
        for(int k = 0;k<n_sign;k++){
          Free(P3[k][j][i]);
          Free(Prhonu[k][j][i]);
          for(int l = 0;l<n_sign;l++){
            Free(P5[l][k][j][i]);
            for(int m = 0;m<n_sign;m++){
              Free(P7[m][l][k][j][i]);
              Free(Psigma7[m][l][k][j][i]);
            }
          }
        }
      }
    }
#else

    // Tidy up

    Free(Pnu);
    Free(P3);
    Free(Prhonu);
    Free(P5);
    Free(P7);
    Free(Psigma7);
#endif
#if 0
    Free(force);
#else
    for(int i =0;i<4;i++){
      Free(force[i]);
    }
#endif

#ifdef PROFILE
//    gettimeofday(&end,NULL);
    dtime +=dclock();
    print_flops(fname, Flops, dtime);
#endif

}

#undef PROFILE
void PT::force_product_sum(PTvector *v, PTvector *w,
				    Float coeff, matrix *f){

//  char *fname = "force_product_sum(*V,*V,F,*M)";
  Flops +=78*vol;
  unsigned long v2 = (unsigned long)v;
  if( qalloc_is_fast(v) &&
      qalloc_is_fast(w) &&
      qalloc_is_fast(f) )
    v2 = v2 - 0xb0000000 + 0x9c000000;
  
#ifdef PROFILE
  Float dtime = -dclock();
#endif
  IFloat coeff2 = 2.0*coeff;
//  printf("force_product_sum::coeff =%e\n",coeff);
//  Float *f_p = (Float *)f;
//  printf("f = %0.3e %0.3e %0.3e %0.3e %0.3e %0.3e\n",
//        f_p[0],f_p[1],f_p[2],f_p[3],f_p[4],f_p[5]);
  asq_force_cross2dag((PTvector *)v2, w, f, vol/2, &coeff2);
//  printf("f = %0.3e %0.3e %0.3e %0.3e %0.3e %0.3e\n",
 //         f_p[0],f_p[1],f_p[2],f_p[3],f_p[4],f_p[5]);
#ifdef PROFILE
  dtime += dclock();
  print_flops(cname,fname,78*vol,dtime);
#endif
}

inline int parity(int *n){
  return( (PT::evenodd + n[0]+n[1]+n[2]+n[3])%2);
}

void PT::update_momenta(matrix **force, Float dt, matrix *mom) {

    matrix mf, mfd;
    double dt_tmp;

    int s[4];
    for(s[3]=0; s[3]<size[3]; s[3]++)
	for(s[2]=0; s[2]<size[2]; s[2]++)
	    for(s[1]=0; s[1]<size[1]; s[1]++)
		for(s[0]=0; s[0]<size[0]; s[0]++){

//		    matrix *ip = mom+LexGauge(s,0);

		    for (int mu=0; mu<4; mu++){			
			mf = force[mu][LexVector(s)];
			mf.TrLessAntiHermMatrix();
//			mf *= 0.5;	
			if(parity(s)) dt_tmp =-dt;
			else dt_tmp = dt;
		        matrix *ip = mom+LexGauge(s,mu);
			(ip)->fTimesV1Plus(0.5*dt_tmp,mf);
#if 0
                    Float *ip_f =(Float *)ip;
                    printf("update_momenta[%d](%d %d %d %d)(%d)= (%0.3e %0.3e)\n",
                    LexGauge(s,mu),s[0],s[1],s[2],s[3],mu,ip_f[1],ip_f[2]);
#endif
		    }
		}
    Flops += vol*54;	
    
}


