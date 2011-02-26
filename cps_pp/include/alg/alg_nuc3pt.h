//------------------------------------------------------------------
//
// alg_nuc3pt.h
//
// Header file for all alg classes relevant to Wilson-type fermion
// nuc3pt' spectrum. The type of glue or fermion is given as
// an argument of type Lattice& to the constructor. If the 
// fermion type is not F_CLASS_WILSON or F_CLASS_CLOVER or 
// F_CLASS_DWF the QPropW constructors exit with a general error.
//
//------------------------------------------------------------------


#ifndef INCLUDED_ALG_NUC3PT_H
#define INCLUDED_ALG_NUC3PT_H

#include <util/lattice.h>
#include <alg/alg_base.h>

#include <alg/qpropw.h>

#include <alg/corrfunc.h>
#include <util/momentum.h>
#include <alg/nuc3pt_arg.h>
#include <alg/nuc2pt.h>
#include <alg/nuc3pt.h>
#include <alg/nuc3pt_gammar.h>

#if TARGET == QCDOC
#include <comms/sysfunc_cps.h>
#endif


CPS_START_NAMESPACE


enum DIR {
 X =  0, 
 Y =  1,
 Z =  2,
 T =  3,
 G5 = -5
} ;


//------------------------------------------------------------------
//
// AlgNuc3pt is derived from Alg and is relevant to  
// meson three point functions with Wilson type 
// fermions (i.e. Wilson, Clover, Dwf). 
// The type of fermion is determined by the argument to the 
// constructor. If the fermion type is not F_CLASS_WILSON or 
// F_CLASS_CLOVER or F_CLASS_DWF the constructors exit with a 
// general error.
//
//------------------------------------------------------------------
class AlgNuc3pt : public Alg
{
 private:
    char* cname;

    Nuc3ptArg* Nuc3pt_arg ;
        // The argument structure for the
        // three point calculation

    QPropW* q_prop ; // Stores the quark propagator
    QPropWSeqBar* u_s_prop ; // Stores the up   quark sequential propagator
    QPropWSeqBar* d_s_prop ; // Stores the down quark sequential propagator

    FILE *fp ; // The I/O file pointer

    void GetThePropagator(int,Float) ;
    void GetTheSeqPropagator(int, Float, SourceType, int*, ProjectType) ;
    
    //void QIO_SaveProp(QPropW& prop, char*, char*, int) ;
    //void QIO_ReadProp(QPropW& prop, char*) ;

    void OpenFile();
    void CloseFile() 
      {
//	Fprintf(stdout,"====closing nuc3pt file\n");
	if(fp != NULL)
	  Fclose(fp);
      }

   
    /*!
      Computes the Scalar insertion 
      \f[
      {\cal O} = \overline{q}  q
      \f]
     */
    void calc_Scalar()
      {     
	Gamma S;
	Nuc3ptGamma Scal(S) ;
	Scal.Calc3pt(*u_s_prop,*q_prop);
	Scal.Calc3pt(*d_s_prop,*q_prop);
	OpenFile(); 
	Scal.Print(fp) ;
	CloseFile();
      }
    


#ifdef _NEDM

    //! Integration trick for NEDM
    void calc_NEDM(){ 
      Gamma Gt(T);
      Nuc3ptGammaR NEDM(Gt,Z) ;
      NEDM.Calc3pt(*u_s_prop,*q_prop);
      NEDM.Calc3pt(*d_s_prop,*q_prop);
      OpenFile(); 
      NEDM.Print(fp) ;
      CloseFile();
    }

    //! Integration trick for NEDM
    void calc_NMM(){
      Gamma Gx(X);
      Nuc3ptGammaR NMM(Gx,Y) ;
      NMM.Calc3pt(*u_s_prop,*q_prop);
      NMM.Calc3pt(*d_s_prop,*q_prop);
      OpenFile(); 
      NMM.Print(fp) ;
      CloseFile();
    }

#endif

    /*!
      Computes the vector current using the
      \f[
      {\cal O}_\mu = \overline{q} \gamma_\mu q
      \f]

      all possible mometa are inserted
     */
    void calc_Vector(const ThreeMom& q) ;

    /*! 
      Computes the Axial vector current using the
      \f[
      {\cal O}_{5\mu} = i \overline{q} \gamma_5 \gamma_\mu q
      \f]  
      
      all possible mometa are inserted
    */
    void calc_Axial(const ThreeMom& q) ;


    /*! 
      Computes the Axial vector current using the
      \f[
      {\cal O}_{5} = \overline{q} \gamma_5 q
      \f]  
      
      all possible mometa are inserted
    */
    void calc_PScalar(const ThreeMom& q) ;

    /*!
      Computes the energy momentum form factor

      \f[
      {\cal O}_{14} = \overline{q} \left[ \gamma_1 
                       \stackrel{\displaystyle \leftrightarrow}{D}_4
		       +\gamma_4 
                       \stackrel{\displaystyle \leftrightarrow}{D}_1
                       \right] q
      \f]

      with polarized projector
    */
    void calc_EnergyMomentum(const ThreeMom& mom) ;
    
    //UnPolarized
    /*!
      Computes the vector current using the
      \f[
      {\cal O}_4 = \overline{q} \gamma_4 q
      \f]
     */
    void calc_Vector()
      {     
	Gamma Gt(T);
	Nuc3ptGamma VectCurr(Gt) ;
	VectCurr.Calc3pt(*u_s_prop,*q_prop);
	VectCurr.Calc3pt(*d_s_prop,*q_prop);
	OpenFile(); 
	VectCurr.Print(fp) ;
	CloseFile();
      }

    /*!
      Computes the first unpolarized moment
      \f[
              \langle x\rangle_q^{(a)} 
      \f]

      \f[
      {\cal O}_{14} = \overline{q} \left[ \gamma_1 
                       \stackrel{\displaystyle \leftrightarrow}{D}_4
		       +\gamma_4 
                       \stackrel{\displaystyle \leftrightarrow}{D}_1
                       \right] q
      \f]

      The operator needs non-zero momentum \f$P_1\$ in order to be computed.
    */
    void calc_X_q_a()
      {
	Gamma Gx(X);
	Gamma Gt(T);
	Derivative Der_t(T);
	Derivative Der_x(T);

	Nuc3ptStru Xq_xt(Gx, Der_t);
	Xq_xt.Calc3pt(*u_s_prop, *q_prop);
	Xq_xt.Calc3pt(*d_s_prop, *q_prop);

	Nuc3ptStru Xq_tx(Gt,Der_x);
	Xq_tx.Calc3pt(*u_s_prop, *q_prop);
	Xq_tx.Calc3pt(*d_s_prop, *q_prop);
	
	Xq_xt += Xq_tx ;

	OpenFile(); 
	Fprintf(fp,"The next is: 14 + 41\n");
	Xq_xt.Print(fp) ; 
	CloseFile();
      }

    /*!
      Computes the first unpolarized moment
      \f[
              \langle x\rangle_q^{(b)} 
      \f]

      \f[
      {\cal O}_{44} = \overline{q} \left[ \gamma_4 
                       \stackrel{\displaystyle \leftrightarrow}{D}_4 
		       -\frac{1}{3} 
                        \sum_{k=1}^3 \gamma_k 
                        \stackrel{\displaystyle \leftrightarrow}{D}_k
                       \right] q
      \f]

      The operator does not need non-zero momentum to by computed.
    */
    void calc_X_q_b()
      {
	Gamma Gt(T) ;
	Derivative Der_t(T) ;
	Nuc3ptStru Xq_tt(Gt,Der_t) ;
	Xq_tt.Calc3pt(*u_s_prop,*q_prop);
	Xq_tt.Calc3pt(*d_s_prop,*q_prop);
	 
	for(int k(X);k<T;k++)
	  {
	    Gamma Gk(k) ;
	    Derivative Der_k(k) ;
	    Nuc3ptStru tmp(Complex(-1.0/3.0,0.0), Gk, Der_k) ;
	    tmp.Calc3pt(*u_s_prop,*q_prop);
	    tmp.Calc3pt(*d_s_prop,*q_prop);
	    Xq_tt+=tmp ;
	  }
	OpenFile(); 
	Fprintf(fp,"The next is: 44 - 1/3 (11 + 22 + 33)\n");
	Xq_tt.Print(fp) ; 
	CloseFile();
      }

    /*!
      Computes the second unpolarized moment
      \f[
              \langle x^2\rangle_q 
      \f]

      \f[
      {\cal O}_{411} = \overline{q} \left[ \gamma_4 
                       \stackrel{\displaystyle \leftrightarrow}{D}_1 
		       \stackrel{\displaystyle \leftrightarrow}{D}_1 
		       -\frac{1}{2} 
                        \sum_{k=2}^3 \gamma_4 
                        \stackrel{\displaystyle \leftrightarrow}{D}_k
			\stackrel{\displaystyle \leftrightarrow}{D}_k
                       \right] q
      \f]

      The operator needs non-zero momentum (\$P_1\f$) to by computed.
    */
    void calc_X2_q()
      {
	Gamma Gt(T);
	Derivative D_xx(X,X);
	Nuc3ptStru X2q_txx(Gt,D_xx) ;
	X2q_txx.Calc3pt(*u_s_prop,*q_prop);
	X2q_txx.Calc3pt(*d_s_prop,*q_prop);
	 
	for(int k(Y);k<T;k++)
	  {
	    Derivative D_kk(k,k);
	    Nuc3ptStru tmp(Complex(-0.5,0.0),Gt,D_kk) ;
	    tmp.Calc3pt(*u_s_prop,*q_prop);
	    tmp.Calc3pt(*d_s_prop,*q_prop);
	    X2q_txx+=tmp ;
	  }
	OpenFile(); 
	Fprintf(fp,"The next is: 411 - 1/2(422+433)\n");
	X2q_txx.Print(fp) ; 
	CloseFile();	
      }

     /*!
      Computes the third unpolarized moment
      \f[
              \langle x^2\rangle_q 
      \f]

      \f[
      {\cal O}_{1144} = \overline{q} \left[ \gamma_1 
                       \stackrel{\displaystyle \leftrightarrow}{D}_1 
		       \stackrel{\displaystyle \leftrightarrow}{D}_4
		       \stackrel{\displaystyle \leftrightarrow}{D}_4
		       +\gamma_2 
                       \stackrel{\displaystyle \leftrightarrow}{D}_2 
		       \stackrel{\displaystyle \leftrightarrow}{D}_3
		       \stackrel{\displaystyle \leftrightarrow}{D}_3
		       -\gamma_1 
                       \stackrel{\displaystyle \leftrightarrow}{D}_1 
		       \stackrel{\displaystyle \leftrightarrow}{D}_3
		       \stackrel{\displaystyle \leftrightarrow}{D}_3
		       -\gamma_2 
                       \stackrel{\displaystyle \leftrightarrow}{D}_2 
		       \stackrel{\displaystyle \leftrightarrow}{D}_4
		       \stackrel{\displaystyle \leftrightarrow}{D}_4
                       \right] q
      \f]

      The operator needs non-zero momentum (\$P_1\f$) to by computed.
    */
    void calc_X3_q()
      {
	Gamma Gx(X);
	Gamma Gy(Y);

	Derivative D_xtt(X,T,T);
	Nuc3ptStru X3q_xxtt(Gx,D_xtt) ;
	X3q_xxtt.Calc3pt(*u_s_prop,*q_prop);
	X3q_xxtt.Calc3pt(*d_s_prop,*q_prop);
	
	{
	  Derivative D_yzz(Y,Z,Z) ;
	  Nuc3ptStru tmp(Gy,D_yzz) ;
	  tmp.Calc3pt(*u_s_prop,*q_prop);
	  tmp.Calc3pt(*d_s_prop,*q_prop);
	  X3q_xxtt += tmp ;
	}
	
	{
	  Derivative D_xzz(X,Z,Z) ;
	  Nuc3ptStru tmp(Complex(-1.0,0.0),Gx,D_xzz) ;
	  tmp.Calc3pt(*u_s_prop,*q_prop);
	  tmp.Calc3pt(*d_s_prop,*q_prop);
	  X3q_xxtt += tmp ;
	}
	
	{
	  Derivative D_ytt(Y,T,T) ;
	  Nuc3ptStru tmp(Complex(-1.0,0.0),Gy, D_ytt) ;
	  tmp.Calc3pt(*u_s_prop,*q_prop);
	  tmp.Calc3pt(*d_s_prop,*q_prop);
	  X3q_xxtt += tmp ;
	}

	
	OpenFile(); 
	Fprintf(fp,"The next is: 1144 + 2233 - 1133 - 2244\n");
	X3q_xxtt.Print(fp) ; 
	CloseFile();
      }

    //Polarized

    /*!
      Computes the Axial vector current using the
      \f[
      {\cal O}_{53} = i \overline{q} \gamma_5 \gamma_3 q
      \f]

    */
    void calc_Axial()
      {
	Gamma G5z(G5,Z);
	Nuc3ptGamma AxialCurr(Complex(0.0,1.0),G5z) ;
	AxialCurr.Calc3pt(*u_s_prop,*q_prop);
	AxialCurr.Calc3pt(*d_s_prop,*q_prop);
	
	OpenFile();
	AxialCurr.Print(fp) ;
	CloseFile();
      }


    /*!
      Computes the first polarized moment
      \f[
               \langle x\rangle_{\Delta q}^{(a)} 
      \f]

      \f[
      {\cal O}_{13} = i \overline{q} \left[ \gamma_5 \gamma_1  
                       \stackrel{\displaystyle \leftrightarrow}{D}_3
		       + \gamma_5 \gamma_3 
                       \stackrel{\displaystyle \leftrightarrow}{D}_1
                       \right] q
      \f]

      The operator needs non-zero momentum \f$P_1\f$ in order to be computed.

      We use the i since the projector used is \f$i\gamma_5\gamma_3\f$.
    */
    void calc_X_Dq_a()
      {
	Nuc3ptStru XDq_xz(Complex(0.0,1.0),Gamma(G5,X),Derivative(Z)) ;
	XDq_xz.Calc3pt(*u_s_prop,*q_prop);
	XDq_xz.Calc3pt(*d_s_prop,*q_prop);
	 
	Nuc3ptStru XDq_zx(Complex(0.0,1.0),Gamma(G5,Z),Derivative(X)) ;
	XDq_zx.Calc3pt(*u_s_prop,*q_prop);
	XDq_zx.Calc3pt(*d_s_prop,*q_prop);

	XDq_xz+=XDq_xz ;

	OpenFile();
	Fprintf(fp,"The next is: 13 + 31\n");
	XDq_zx.Print(fp);
	CloseFile();
      }

    /*!
      Computes the first polarized moment
      \f[
               \langle x\rangle_{\Delta q}^{(b)} 
      \f]

      \f[
      {\cal O}_{34} = i \overline{q} \left[ \gamma_5 \gamma_3  
                       \stackrel{\displaystyle \leftrightarrow}{D}_4
		       + \gamma_5 \gamma_4 
                       \stackrel{\displaystyle \leftrightarrow}{D}_3
                       \right] q
      \f]

      The operator does not need non-zero momentum.

      We use the i since the projector used is \f$i\gamma_5\gamma_3\f$.
    */
    void calc_X_Dq_b()
      {
	Nuc3ptStru XDq_zt(Complex(0.0,1.0),Gamma(G5,Z),Derivative(T)) ;
	XDq_zt.Calc3pt(*u_s_prop,*q_prop);
	XDq_zt.Calc3pt(*d_s_prop,*q_prop);
	 
	Nuc3ptStru XDq_tz(Complex(0.0,1.0),Gamma(G5,T),Derivative(Z)) ;
	XDq_tz.Calc3pt(*u_s_prop,*q_prop);
	XDq_tz.Calc3pt(*d_s_prop,*q_prop);

	XDq_zt+=XDq_tz ;

	XDq_zt.setTag("p") ;
	OpenFile();
	Fprintf(fp,"The next is: 34 + 43\n");
	XDq_zt.Print(fp);
	CloseFile();

      }


    void calc_X2_Dq()
      {
      }

    /*!
      Computes twist three term (\f$ d_1\f$) of the first polarized moment

      \f[
      {\cal O}_{34} = i \overline{q} \left[ \gamma_5 \gamma_3  
                       \stackrel{\displaystyle \leftrightarrow}{D}_4
		       - \gamma_5 \gamma_4 
                       \stackrel{\displaystyle \leftrightarrow}{D}_3
                       \right] q
      \f]

      The operator does not need non-zero momentum.

      We use the i since the projector used is \f$i\gamma_5\gamma_3\f$.
    */
    void calc_d1() // Could be compined with X_Dq_a
      {
	Nuc3ptStru d1_zt(Complex(0.0,1.0),Gamma(G5,Z),Derivative(T)) ;
	d1_zt.Calc3pt(*u_s_prop,*q_prop);
	d1_zt.Calc3pt(*d_s_prop,*q_prop);
	 
	Nuc3ptStru XDq_tz(Complex(0.0,-1.0),Gamma(G5,T),Derivative(Z)) ;
	XDq_tz.Calc3pt(*u_s_prop,*q_prop);
	XDq_tz.Calc3pt(*d_s_prop,*q_prop);

	d1_zt+=XDq_tz ;
	d1_zt.setTag("m") ;

	OpenFile();
	Fprintf(fp,"The next is: 34 - 43\n");
	d1_zt.Print(fp);
	CloseFile();
      }
    void calc_d2()
      {
      }

    
    //Transversity

    /*!
      Computes the Axial vector current using the
      \f[
      {\cal O}_{534} = i \overline{q} \gamma_5 \gamma_3 \gamma_4q
      \f]

    */
    void calc_Tensor()
      {
	Nuc3ptGamma Tensor(Complex(0.0,1.0),Gamma(G5,Z,T)) ;
	Tensor.Calc3pt(*u_s_prop,*q_prop);
	Tensor.Calc3pt(*d_s_prop,*q_prop);

	OpenFile();
	Tensor.Print(fp) ;
	CloseFile();
      }
    
    /*!
      Computes the first Transversity moment
      \f[
               \langle x\rangle_{\delta q} 
      \f]

      \f[
      {\cal O}_{341} = i \overline{q} \left[ \gamma_5 \gamma_3 \gamma_4  
                       \stackrel{\displaystyle \leftrightarrow}{D}_1
		       + \gamma_5 \gamma_3 \gamma_1
                       \stackrel{\displaystyle \leftrightarrow}{D}_4
                       \right] q
      \f]

      The operator needs non-zero momentum \f$P_1\f$ in order to be computed.

      The first term is proportional to \f$S_3 M_N P_1\f$ the second
      to \f$S_3 P_1 M_N\f$. So we expect both to be non-zero.
      For that reason we sum them up befor we print them

      We use the i since the projector used is \f$i\gamma_5\gamma_3\f$.
    */
    void calc_X_dq()
      {
	Nuc3ptStru Xdq_ztx(Complex(0.0,1.0),Gamma(G5,Z,T),Derivative(X)) ;
	Xdq_ztx.Calc3pt(*u_s_prop, *q_prop);
	Xdq_ztx.Calc3pt(*d_s_prop, *q_prop);
	 
	Nuc3ptStru Xdq_zxt(Complex(0.0,1.0),Gamma(G5,Z,X),Derivative(T)) ;
	Xdq_zxt.Calc3pt(*u_s_prop, *q_prop);
	Xdq_zxt.Calc3pt(*d_s_prop, *q_prop);
	
	Xdq_ztx += Xdq_zxt ;

	OpenFile();
	Fprintf(fp,"The next is: 341 + 314\n");
	Xdq_ztx.Print(fp) ;
	CloseFile();
      }

    /*!
      Computes the conserved vector current using the
      \f[
      {\cal O}_\mu = \overline{q} \gamma_\mu q
      \f]

      all possible mometa are inserted
     */
    void calc_Cons_Vector(int i, ThreeMom*) ;

    /*! 
      Computes the conserved Axial vector and vector currents using the
      \f[
      {\cal O}_{5\mu} = i \overline{q} \gamma_5 \gamma_\mu q
      \f]  
      
      all possible mometa are inserted
    */
    void calc_Cons_Axial_Vector(int i, ThreeMom*) ;


 public:
    AlgNuc3pt(Lattice & latt, CommonArg* c_arg, Nuc3ptArg* arg);

    virtual ~AlgNuc3pt();

    void run();

    void test_run();

    void run_hbd();

    void run_2pt();
};

CPS_END_NAMESPACE

#endif
