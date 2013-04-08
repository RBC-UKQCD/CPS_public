#include<config.h>

#ifdef USE_FFTW
#include <fftw3.h>

#ifndef INCLUDED_FFT_H_
#define INCLUDED_FFT_H_

#include <util/data_types.h>
#include <util/site.h>
#include <comms/scu.h>

CPS_START_NAMESPACE

//////////////////////////////////////////////////////////
/*! @file

  @brief  a FFT implementation

  @author Taku Izubuchi
  @date   2007-05-30

  TODOs:

  - There is unnecessary local memory copy in `shift()' and `unshift()'.
    - When there is more than one direction, mu, which doesn't satisfy
         nlr % N == 0
      where
         nlr = GJP.VolSites() / GJP.NodeSites(mu) * [num of Complex on a site, such as 144]
	 N   = GJP.Nodes(mu),
     the `FFT_one_dir()' bangs.

     This limitation is curable by padding data, for example, by
     increasing number of site Complex.  More memory/computation
     preserving way is to pad a chunk in `skew' section and adjust
     accordingly in `unskew()' automatically.

     This is not done for now as the day after tomorrow is the
     deadline of the lattice conference and I will have to write
     something. Tom is very mad at me also.
  
  -----------------------------------------------------------

   CVS keywords
 
   $Author: chulwoo $
   $Date: 2013-04-08 20:50:00 $
   $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/fft.h,v 1.2 2013-04-08 20:50:00 chulwoo Exp $
   $Id: fft.h,v 1.2 2013-04-08 20:50:00 chulwoo Exp $
   $Name: not supported by cvs2svn $
   $Locker:  $
   $RCSfile: fft.h,v $
   $Revision: 1.2 $
   $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/fft.h,v $
   $State: Exp $  */

/*----------------------------------------------------------*/


//------------------------------------------------------------------
//!  A class implementing  FFT in one direction.
/*!
  This class perform FFT in the following way:
  1.  gather a whole dimension into a node
  2.  do local FFT
  3.  redistribute results into nodes
*/

class FFT_one_dir
{
private:
  const int mu;
  //!< The direction this FFT is to be done
  const int N;
  //!< The number of node in the direction
  size_t P;
  //!< The Coordinate of this node in the direction
  Float* fp_tmp;
  //!< pointer to temporary vector

  // Input Data is aligned as  D[nl][n][nr]

  const size_t nl;
  //!< The number of  index(es), which is(are) left of the index for the direction, in a node.
  
  const int n;
  //!< The number of sites in the direction in a node. 
  const size_t nr;
  //!< The number of last index(es), which is(are) right of the index for the direction, in a node

  size_t nlr;     //!< nl*nr
  size_t nlr_N;   //!< nl*nr/N
  size_t nnlr_N;  //!< n*nl*nr/N
  
  char* cname;
  

  //Large Data
  void getPlusDataL(IFloat *rcv_buf, IFloat *send_buf, int len, int mu)
  {
    const int MAX_LENGTH = 4096;
    const int blk(MAX_LENGTH);
    const int N (len / MAX_LENGTH);
      
    for(int i=0;i<N; ++i)
      getPlusData(rcv_buf+i*blk, send_buf+i*blk,
		  MAX_LENGTH, mu );
    if(N*blk < len)
      getPlusData(rcv_buf+N*blk, send_buf+N*blk,
		  len-N*blk, mu );
  }

  void getMinusDataL(IFloat *rcv_buf, IFloat *send_buf, int len, int mu)
  {
    const int MAX_LENGTH = 4096;
    const int blk(MAX_LENGTH);
    const int N (len / MAX_LENGTH);
      
    for(int i=0;i<N; ++i)
      getMinusData(rcv_buf+i*blk, send_buf+i*blk,
		  MAX_LENGTH, mu );
    if(N*blk < len)
      getMinusData(rcv_buf+N*blk, send_buf+N*blk,
		  len-N*blk, mu );
  }


  

  //!  "Skew" data 
  /*!
    Reorder index in a node 

    @param[in]   fpin   fpin [nl][n][nr],
    @param[out] fpout   fpout[N][n][nlr/N]

    also shift index I according to the node coordinate `P'
    
  */
#define Fpin( il_, i_,  ir_ )  (fpin+ (ir_)+nr*( (i_ ) +n *(il_) ))
#define Fpout( I_, i_, ilr_N_ )  (fpout+ (ilr_N_)+nlr_N*( (i_) +n*(I_)  ) )  
  void skew(Float *fpin, Float *fpout)
  {
   
    //  convert input data so that `n' is the first index
    for(int il=0;il<nl;++il)
      for(int i=0;i<n;++i) 
	for(int ir=0;ir<nr;++ir) {
	  int ilr= ir+nr*il;
	  int I  = ilr / nlr_N;
	  int ilr_N = ilr % nlr_N;

	  *(Fpout( (N+I-P)%N, i, ilr_N )) = *(Fpin(il, i, ir));
	}
  }
#undef Fpin
#undef Fpout
  
  //!  "Shift" data 
  /*!
    Pass data between nodes

    @param[in/out]  fp1[N][n][nlr/N]
    @param[tmp]     fp2[N][n][nlr/N]  working vector

  */
#define Fp1(I_,iilr_N_)  (fp1 + (iilr_N_)+nnlr_N*(I_) )
#define Fp2(I_,iilr_N_)  (fp2 + (iilr_N_)+nnlr_N*(I_) )
  void shift( Float *fp1, Float *fp2 )
  {
    // just copy the first line
    moveMem( Fp2(0,0),  Fp1(0,0), nnlr_N*sizeof(IFloat) );

#if 0
    /*
      slower version uses getMinusData for all communication
    */
    for(int I=1; I<N;++I){
      getMinusDataL( Fp2(I,0), Fp1(I,0),
		    (N-I)*nnlr_N, mu ) ;
      if(I!=N-1)
	moveMem(  Fp1(I+1, 0), Fp2(I+1, 0),
		  (N-I-1)*nnlr_N*sizeof(IFloat) );      
    }
#else
    /*
      faster version uses get{Minus,Plus}Data for {first,second} half
    */
      const int M=N/2;
    // The first half uses getMinus
    for(int I=1; I<M;++I){
      getMinusDataL( Fp2(I,0), Fp1(I,0),
		    (M-I)*nnlr_N, mu ) ;
      moveMem(  Fp1(I+1, 0), Fp2(I+1, 0),
		  (M-I-1)*nnlr_N*sizeof(IFloat) );      
    }
    // The second half uses getPlus
    for(int I=M; I<N;++I){
      getPlusDataL( Fp2(M,0), Fp1(M,0),
		   (N-I)*nnlr_N, mu ) ;
      if(I!=N-1)
	moveMem(  Fp1(M, 0), Fp2(M, 0),
		  (N-I-1)*nnlr_N*sizeof(IFloat) );      
    }
#endif
    // The next local shift can be skipped by noting
    // the local shift (and flip parity) could be
    // compensate in fourier space,
    // but let's be safe and lazy for now.
    for(int I=0; I<N;++I)
      moveMem(	Fp1(I, 0), Fp2( (P-I+N)%N, 0 ), nnlr_N*sizeof(IFloat) );
  }
#undef Fp1
#undef Fp2


  //!  fft in a node (lazy & slow)
  /*!
    @param[in]  fpin[N][n][nlr/N] 
    @param[out] fpout[N][n][nlr/N] 
    @param[in]  fft_forward : 0 for backward fourier transformation

    call FFTW version 3 
  */


void fft(int n, double theta, double ar[], double ai[], 
        double tmpr[], double tmpi[])
{
    int radix, n_radix, j, m, r;
    double xr, xi, wr, wi;

    if (n <= 1) return;
    /* ---- factorization ---- */
    for (radix = 2; radix * radix <= n; radix++) {
        if (n % radix == 0) break;
    }
    if (n % radix != 0) radix = n;
    n_radix = n / radix;
    /* ---- butterflies ---- */
    for (j = 0; j < n_radix; j++) {
        for (m = 0; m < radix; m++) {
            xr = ar[j];
            xi = ai[j];
            for (r = n_radix; r < n; r += n_radix) {
                wr = cos(theta * m * r);
                wi = sin(theta * m * r);
                xr += wr * ar[r + j] - wi * ai[r + j];
                xi += wr * ai[r + j] + wi * ar[r + j];
            }
            wr = cos(theta * m * j);
            wi = sin(theta * m * j);
            tmpr[m * n_radix + j] = xr * wr - xi * wi;
            tmpi[m * n_radix + j] = xi * wr + xr * wi;
        }
    }
    for (r = 0; r < n; r += n_radix) {
        fft(n_radix, theta * radix, &tmpr[r], &tmpi[r], ar, ai);
    }
    for (j = 0; j < n_radix; j++) {
        for (m = 0; m < radix; m++) {
            ar[radix * j + m] = tmpr[n_radix * m + j];
            ai[radix * j + m] = tmpi[n_radix * m + j];
        }
    }
}
  

  //!  fft in a node
  /*!
    @param[in]  fpin[N][n][nlr/N] 
    @param[out] fpout[N][n][nlr/N] 
    @param[in]  fft_forward : 0 for backward fourier transformation

    call FFTW version 3 
  */
  void local_fft(Float *fpin, Float *fpout, int fft_forward )
  {
    //fftw_complex *in, *out;
    //    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N*Nr);
    //    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N*Nr);

    const int rank=1;
    int nlen[1]={N*n};
    const int howmany=nlr_N/2; // in the unit of complex
    const int *inembed=NULL;
    const int istride=nlr_N/2; // in the unit of complex
    const int idist=1;

    const int *onembed=NULL;
    const int ostride=nlr_N/2; // in the unit of complex
    const int odist=1;
  
#if 0
    printf("nlen=%d %d %d, nl=%d,nr=%d,howmany=%d\n",
	   N,n, nlen[0],nl,nr,howmany);

#define Fpin( I_,iilr_N_)  (fpin + (iilr_N_)+nnlr_N*(I_) )
    for(int I=0;I<N;++I) for(int iilr_N=0;iilr_N<nnlr_N;++iilr_N)
	printf("IN %d %d %f\n",I, iilr_N, *(Fpin(I,iilr_N)) );
#undef Fpin
#endif
    
    int FFTW_F_or_B;
    if( fft_forward ) FFTW_F_or_B = FFTW_FORWARD;
    else  FFTW_F_or_B = FFTW_BACKWARD;
    
    fftw_plan  plan =
      fftw_plan_many_dft(rank, nlen, howmany,
			 (fftw_complex*) fpin,  inembed, istride, idist,
			 (fftw_complex*) fpout, onembed, ostride, odist,
			 FFTW_F_or_B,  FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    //    fftw_free(in); fftw_free(out);

  }

  //!  "Unshift" data
  /*!
    Pass data between nodes

    @param[in/out]  fp1[N][n][nlr/N]
    @param[tmp]     fp2[N][n][nlr/N]  working vector

  */
#define Fp1(I_,iilr_N_)  (fp1 + (iilr_N_)+nnlr_N*(I_) )
#define Fp2(I_,iilr_N_)  (fp2 + (iilr_N_)+nnlr_N*(I_) )
  void unshift( Float *fp1, Float *fp2 )
  {
    // The next local shift can be skipped by noting
    // the local shift (and flip parity) could be
    // compensate in fourier space,
    // but let's be safe and lazy for now.
    for(int I=0; I<N;++I)
      moveMem(	Fp1(I, 0), Fp2( (P-I+N)%N, 0 ), nnlr_N*sizeof(IFloat) );

#if 0
    /*
      slower version uses getMinusData for all communication
    */
    for(int I=N-1; I>=1;--I){
      getPlusDataL( Fp2(I,0), Fp1(I,0),
		   (N-I)*nnlr_N, mu ) ;
      moveMem(  Fp1(I, 0), Fp2(I, 0),
		(N-I)*nnlr_N*sizeof(IFloat) );      
    }
#else
    /*
      faster version uses get{Minus,Plus}Data for {first,second} half
    */
      const int M=N/2;
    // The first half uses getPlusMinus
    for(int I=M-1; I>=1;--I){
      getPlusDataL( Fp2(I,0), Fp1(I,0),
		    (M-I)*nnlr_N, mu ) ;
      moveMem(  Fp1(I, 0), Fp2(I, 0),
		  (M-I)*nnlr_N*sizeof(IFloat) );      
    }
    // The second half uses getMinus
    for(int I=N-1; I>=M;--I){
      getMinusDataL( Fp2(M,0), Fp1(M,0),
		    (N-I)*nnlr_N, mu ) ;
      moveMem(  Fp1(M, 0), Fp2(M, 0),
		  (N-I)*nnlr_N*sizeof(IFloat) );      
    }
#endif
    // just copy the first line
    moveMem( Fp2(0,0),  Fp1(0,0), nnlr_N*sizeof(IFloat) );
  }
#undef Fp1
#undef Fp2


  //!  "Unshift" data
  /*!
    Pass back data between nodes

    @param[in/out]  fp1[N][n][nlr/N]
    @param[tmp]     fp2[N][n][nlr/N]  working vector

  */
#define Fpin(  I_, i_, ilr_N_ ) (fpin+ (ilr_N_)+nlr_N*( (i_) +n*(I_)  ) )  
#define Fpout( il_, i_,  ir_ )  (fpout+ (ir_)+nr*( (i_ ) +n *(il_) ))
  void unskew(Float *fpin, Float *fpout)
  {
    //  convert input data so that `n' is the first index
    for(int il=0;il<nl;++il)
      for(int i=0;i<n;++i) 
	for(int ir=0;ir<nr;++ir) {
	  int ilr= ir+nr*il;
	  int I  = ilr / nlr_N;
	  int ilr_N = ilr % nlr_N;
	  *(Fpout(il, i, ir))=  *(Fpin( (N+I-P)%N, i, ilr_N )) ; 
	}
  }
#undef Fpin
#undef Fpout
  
public:


  //!  FFT in one direction on a node
  /*!

  Y[k] = sum_{j=0}^{n-1} X[j] e^{- 2 \pi j k / n  i} for forward FFT
  Y[k] = sum_{j=0}^{n-1} X[j] e^{+ 2 \pi j k / n  i} for backward FFT

  @param[in]   mu   the direction to be FFTed
  @param[in/out]  Float fpin[_nl][_n][_nr][2] 
  @param[in]   _nl, _n, _nr  
  @param[in] fft_forward, =0 for backward fourier transformation 
  @param[in] flag_dist_back is zero if the result of fft should not be redistribute

  */
  FFT_one_dir(const int _mu,
	      Float* fpin,
	      const int _nl, const int _n, const int _nr,
	      int fft_forward=1,
	      int flag_dist_back=1 ):
    mu(_mu), N(GJP.Nodes(mu)), n(_n),  nl(_nl), nr(_nr)   
  {
    cname="FFT_one_dir"; 
    char* fname="FFT_one_dir(...)";

    nlr    = nl*nr;
    nlr_N  = nlr/N;
    nnlr_N = n*nlr_N;

    if(nlr%N)
      ERR.NotImplemented(cname, fname, "nl*nr is not a multiple of N.\nPad DATA, or modify code (should be easy), sorry.\n\n");

    
    //Gets the grid coordinate of this node in the direction.
    P = GJP.NodeCoor(mu);
      

    // allocate temporary vector
    Float* fptmp = (Float *) smalloc(n*nlr * sizeof(Float));
    if(fptmp == 0) ERR.Pointer(cname,fname, "fptmp");
    VRB.Smalloc(cname,fname,"fptmp",fptmp, n*nlr * sizeof(Float));

    /*
     * reorder and communicate input data
     * then do the local fft,
     *  and distribute back if needed
     */
    skew(  fpin  , fptmp );
    shift( fptmp, fpin   );  //  the first arg is the output
    local_fft( fptmp, fpin, fft_forward  ); //  the first arg is the output
    
    if(! flag_dist_back ) {

      moveMem( fpin,  fptmp, n*nlr*sizeof(IFloat) );

    } else {

      unshift( fptmp, fpin); // the first arg is the output
      unskew( fptmp, fpin);

    } 

    // Free temporary vector
    VRB.Sfree(cname,fname, "fptmp",fptmp);
    sfree(fptmp);
  }
  
};


CPS_END_NAMESPACE
#endif

#endif //#ifdef USE_FFTW
