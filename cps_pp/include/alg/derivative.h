//------------------------------------------------------------------
//  
//   derivative.h
//  
//   Header file for the  Shift class and Derivative classes
//  
//   March 2001
//  
//   Kostas Orginos
//  
//------------------------------------------------------------------*/
#ifndef INCLUDED_DERIVATIVE_H
#define INCLUDED_DERIVATIVE_H
//#include <stdio.h>
#include <util/data_types.h>
#include <util/site.h>
#include <util/qcdio.h>

CPS_START_NAMESPACE

enum Dir {FRWD=+1,
	  BKWD=-1 };

enum Quark {QUARK =+1,
	    AQUARK=-1 };


/*!------------------------------------------------------------------
  
  The Shift class is used to represent Derivative operators for quark fields
  
   \f[
   \stackrel{\displaystyle \leftrightarrow}{D}_\mu = 
   \frac{1}{2}\left( \stackrel{\displaystyle \rightarrow}{D}_\mu -
                     \stackrel{\displaystyle \leftarrow}{D}_\mu \right)
   \f]
 
   Where \f$\stackrel{\displaystyle \rightarrow}{D}\f$ acts on the quark
   and \f$\stackrel{\displaystyle \leftarrow}{D}\f$ acts on the anti-quark

   Each derivative is the difference of a forward hopping
   (\f$D^{+}_\mu\f$) term and a backward hopping term (\f$D^{-}_\mu\f$)

   i.e.

   \f[
   \stackrel{\displaystyle \rightarrow}{D}_\mu = 
   \frac{1}{2}\left( \stackrel{\displaystyle \rightarrow}{D}^+_\mu -
                     \stackrel{\displaystyle \rightarrow}{D}^-_\mu \right)
   \f]

   \f[
   \stackrel{\displaystyle \rightarrow}{D}^+_\mu = U_\mu(x) 
           \stackrel{\displaystyle \rightarrow}{\delta}_{y,x+\mu}
   \f]
   \f[
   \stackrel{\displaystyle \rightarrow}{D}^-_\mu = U^{\dagger}_\mu(x-\mu) 
           \stackrel{\displaystyle \rightarrow}{\delta}_{y,x-\mu}
   \f]
   and

   \f[
   \stackrel{\displaystyle \leftarrow}{D}_\mu = 
   \frac{1}{2}\left( \stackrel{\displaystyle \leftarrow}{D}^+_\mu -
                     \stackrel{\displaystyle \leftarrow}{D}^-_\mu \right)
   \f]
  
   \f[
   \stackrel{\displaystyle \leftarrow}{D}^+_\mu = U^{\dagger}_\mu(x) 
           \stackrel{\displaystyle \leftarrow}{\delta}_{y,x+\mu}
   \f]
   \f[
   \stackrel{\displaystyle \leftarrow}{D}^-_\mu = U_\mu(x-\mu) 
           \stackrel{\displaystyle \leftarrow}{\delta}_{y,x-\mu}
   \f]
   
//------------------------------------------------------------------*/
class Shift
{
  int  d     ; // Forward or backward
  int  q     ; // Act on quark -> or on anti-quark <-

  bool more ; // Used to iterate over all possible shifts
public:
  Shift(){}
  Shift(int d_, int q_): d(d_),q(q_){}
  
  /*! Forward or backward */
  int D() const {return d;}
  /*! Act on quark or anti-quark */
  int Q() const {return q;}
  /*! The sign of the term represented 

    Forward  quark: +1 | Forward  anti-quark: -1

    Backward quark: -1 | Backward anti-quark: +1

    if \f$d=+1\f$ for forward and \f$d=-1\f$ for backward 

    and if \f$q=+1\f$ for quark and \f$d=-1\f$ for anti-quark
    
    then
    \f[
    \mbox{Float sign()} = d*q 
    \f]
  */
  Float sign(){return (Float)(d*q);} 

  int bit(){return -d*q;} 

  void Modify(int d_, int q_){d=d_;q=q_;}

  /*! used to iteraterate over all possible shifts */
  void Next() ;

  /*! used to iteraterate over all possible shifts */
  bool NotEnd(){return more ; }
  /*! used to iteraterate over all possible shifts */
  void Start(){q=QUARK;d=FRWD;more=true;}
  /*! used to iteraterate over all possible shifts with Q fixed to QUARK*/
  void NextQuark(){d-=2;more=(d>=-1);}


  //For debugging
  //void print(){printf("shift: d=%2i q=%2i | %4g\n",d,q,sign());}

  ~Shift(){}
};

/*!
  This class represents the derivative operators for the structure functions
  Look at the documentation of the Shift class for the definitions of
  the derivatives.

  In order to implement the derivative we break it up in to terms which
  are a string of Shifts. The Shift is either a hop forward or backward
  on a quark or on an anti-quark. A forward hop is a + sign. A backward
  hop is a - sign. Action on a quark is a + sign, on an anti-quark is
  a - sign.  Each term comes with an overall factor of \f$\frac{1}{4}\f$
  (she definition in the Shift documentation).

  The Derivative class works as following: The member  int* indx 
  is a list of indices with lentgh \f$n=\mbox{Num}\f$.
  \f[
  \mbox{\rm int  *indx} = \left[\mu_1\,\mu_2\,\mu_3\,...\,\mu_n\right]
  \f]

  This list represents the operator:
  \f[
  {\cal O}_{\mu_1\,\mu_2\,\mu_3\,...\,\mu_n} = 
  \stackrel{\displaystyle \leftrightarrow}{D}_{\mu_1}\,
  \stackrel{\displaystyle \leftrightarrow}{D}_{\mu_2}\,
  \stackrel{\displaystyle \leftrightarrow}{D}_{\mu_3}\,...\,
  \stackrel{\displaystyle \leftrightarrow}{D}_{\mu_n}
  \f]

  The Derivative class has an iterator that loops over all possible terms
  (each term is a Shift), and stores the current Shift in  Shift* sh.

  Example of a Shift:
  \f[
  \mbox{*sh} = \stackrel{\displaystyle \rightarrow}{D}^+_{\mu_1}\,
  \stackrel{\displaystyle \rightarrow}{D}^-_{\mu_2}\,
  \stackrel{\displaystyle \leftarrow}{D}^-_{\mu_3}\, ...
  \stackrel{\displaystyle \rightarrow}{D}^+_{\mu_n}
                
  \f]


  void CalcEndPoints (int *aquark, int *quark, Site& s) calculates the
  position of the anti-quark and the quark corresponding to the current
  Shift.

  void CalcListOfLinks (int *lnk) returns the list of directions (lnk)
  that correspond to the shift. The directions are described in the
  form needed by Lattice:: PathOrdPlus i.e. 
  (x y z t -x -y -z -t) --->  (0 1 2 3 4 5 6 7)

  Float Fact () returns the factor the current Shift is multiplied
  by
  \f[
  \mbox{Float Fact ()} = \frac{1}{4^n} (-1)^{N_b} (-1)^{N_a}
  \f]
  Where \f$N_b\f$ the number of backward quark hops and
        \f$N_a\f$ the number of forward anti-quark hops 
*/
class Derivative
{
  int* indx ;
  int Num ;

  bool more ;
  Shift* sh ;

  char *tag ;

  /*! Computes the operator Tag, a string to be used in 3pt IO */
  void Tag() ;

 public:
  Derivative():Num(0){Tag();}

  Derivative(int mu) ;

  Derivative(int mu,int nu) ;
  
  Derivative(int mu,int nu,int rho) ;

  Derivative(int *mu,int N) ;
  
  /*! Copy constructor */
  Derivative(const Derivative& tt) ;

  
  /*! used to iteraterate over all possible terms in the derivatice */
  void Start() ;

  /*! used to iteraterate over all possible terms in the derivative */
  void Next() ;

  /*! used to iteraterate over all possible terms in the derivative 
    that acts only on quarks
   */
  void NextQuark() ;

  /*! used to iteraterate over all possible shifts */
  bool NotEnd(){return more ; }

  /*! 
    It caclulates the location of the quark and the antiquark corresponding
    to the current shift. Site& s is the current Site.
  */
  void CalcEndPoints(int *aquark, int *quark, Site& s) ;

  /*! 
    It caclulates the list of dirs the Lattice:: PathOrdPlus needs to construct
    the gauge connection needed by the derivative
  */
  void CalcListOfLinks(int *lnk) ;

  /*! 
    It caclulates the sign associated with the current shift and
    the normalization factor.
    \f[
    \mbox{Float Fact ()} = \frac{1}{4^n} (-1)^{N_b} (-1)^{N_a}
    \f]
    Where \f$N_b\f$ the number of backward quark hops and
    \f$N_a\f$ the number of forward anti-quark hops 
  */
  Float Fact() ;


  /*! Returns the index representing the DTerm.
    It is a number from 0 to \f$2^{Num}-1\f$ corresponding to the term
    from the list of DTerms. The mapping is the following:
    foreach shift the associate a bit b=d*q (-1-->0 1-->1)
    
    \li b=-1 the link is forward (there is a minus in the definition of bit)
    \li b=1  the link is backward
   
    
    Example: 
    \verbatim
    Path F ----> Indx 0
    Path B ----> Indx 1

    Path FF ---> Indx 0
    Path BF ---> Indx 1
    Path FB ---> Indx 2
    Path BB ---> Indx 3

    Path FFF ---> Indx 0
    Path FBF ---> Indx 1
    Path FFB ---> Indx 2
    Path FBB ---> Indx 3
    Path BFF ---> Indx 4
    Path BBF ---> Indx 5
    Path BFB ---> Indx 6
    Path BBB ---> Indx 7    
    \endverbatim
   */
  int DTermIndx() ;

  /*! Returns the number of derivative indices */
  int N(){return Num;}

  /*! prints the tag */
  void printTag(FILE *fp){Fprintf(fp,"%s",tag);}

  // For debugging 
  //void sh_print(){for(int i(0);i<Num; i++) sh[i].print(); printf("\n"); }

  /*! Returns the number of DTerms (derivatives only on quarks): 
    \f[
    \mbox{NDTerms()} = 2^{\mbox{Num}}
    \f]
   */
  int NDTerms(){return 1<<Num ;}
  
  ~Derivative() ;

} ;

CPS_END_NAMESPACE

#endif // !INCLUDED_DERIVATIVE_H  


