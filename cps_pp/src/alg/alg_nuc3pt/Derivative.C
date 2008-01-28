#include <alg/derivative.h>

CPS_START_NAMESPACE

/*! used to iteraterate over all possible shifts */
void Shift::Next()
{
  if(q==QUARK)
    {
      q-=2 ;
      return ;
    }
  if(d==FRWD)
    {
      q=QUARK;
      d-=2 ;
      return ;
    } 
  more=false ;
  
}

/*! Returns the operator Tag, a string to be used in 3pt IO */
void Derivative::Tag()
{
  char *cname="Derivative" ;
  char *fname="Tag" ;
  
  //check the shift and indx pointers
  if(Num>0){
    if(indx == 0) ERR.Pointer(cname, fname, "indx");
    if(sh   == 0) ERR.Pointer(cname, fname, "sh");
  }

  tag = new char[2*Num+1]; 
  if(tag==0) ERR.Pointer(cname, fname, "tag");

  tag[0]='\0' ; // Start with an empty string
  for(int i(0);i<Num;i++)
    sprintf(tag,"%sD%i",tag,indx[i]+1) ;
}

Derivative::Derivative(int mu):Num(1)
{
  indx = new int[Num] ;
  sh = new Shift[Num] ;
  indx[0]=mu;
  Tag();
}

Derivative::Derivative(int mu,int nu):Num(2)
{
  indx = new int[Num] ;
  sh = new Shift[Num] ;
  indx[0]=mu;
  indx[1]=nu;
  Tag();
}

Derivative::Derivative(int mu,int nu,int rho):Num(3)
{
  indx = new int[Num] ;
  sh = new Shift[Num] ;
  indx[0]=mu;
  indx[1]=nu;
  indx[2]=rho;
  Tag();
}

Derivative::Derivative(int *mu,int N):Num(N)
{
  indx = new int[Num] ;
  sh = new Shift[Num] ;
  for(int i=0;i<N;i++) indx[i]=mu[i];
  Tag();
}

/*! Copy constructor */
Derivative::Derivative(const Derivative& tt)
{
  Num = tt.Num ;
  indx = new int[Num] ;
  sh = new Shift[Num] ;
  for(int i(0);i<Num;i++) 
    indx[i]=tt.indx[i];
  Tag();
}

void Derivative::Start()
{
  more=true ;
  for (int i(0) ; i<Num; i++)
    sh[i].Start() ;
  if(Num==0)
    more=false ;
}

/*! used to iteraterate over all possible terms in the derivative */
void Derivative::Next()
{
  for(int i(0);i<Num;i++)
    {
      sh[i].Next();
      if(sh[i].NotEnd())
	return ;
      else
	sh[i].Start() ;
    }
  more=false ;
}

/*! used to iteraterate over all possible terms in the derivative 
  that acts only on quarks
*/
void Derivative::NextQuark()
{
  for(int i(0);i<Num;i++)
    {
      sh[i].NextQuark();
      if(sh[i].NotEnd())
	return ;
      else
	sh[i].Start() ;
    }
  more=false ;
}

/*! 
  It caclulates the location of the quark and the antiquark corresponding
  to the current shift. Site& s is the current Site.
*/
void Derivative::CalcEndPoints(int *aquark, int *quark, Site& s)
{
  // Start from the current site
  for(int mu(0);mu<4;mu++)
    quark[mu]  = aquark[mu] = s.Coor(mu) ;
  
  for(int i(0);i<Num; i++)
    if(sh[i].Q()==QUARK)
      quark[indx[i]] += sh[i].D() ;
    else 
      aquark[indx[i]] += sh[i].D() ;
}

/*! 
  It caclulates the list of dirs the Lattice:: PathOrdPlus needs to construct
  the gauge connection needed by the derivative
*/
void Derivative::CalcListOfLinks(int *lnk)
{
  for(int i(0);i<Num; i++)
    if( ((sh[i].D()==FRWD)&&(sh[i].Q()==QUARK)) || 
	((sh[i].D()==BKWD)&&(sh[i].Q()==AQUARK)) )
      lnk[i]=indx[i] ; // 0 1 2 3 : x y z t
    else
      lnk[i]=indx[i]+4 ; // 4 5 6 7 : -x -y -z -t 
}

/*! 
  It caclulates the sign associated with the current shift and
  the normalization factor.
  \f[
  \mbox{Float Fact ()} = \frac{1}{4^n} (-1)^{N_b} (-1)^{N_a}
  \f]
  Where \f$N_b\f$ the number of backward quark hops and
  \f$N_a\f$ the number of forward anti-quark hops 
*/
Float Derivative::Fact()
{
  Float s(1.0) ;
  for(int i(0);i<Num; i++)
    s *= 0.25*sh[i].sign() ;
  
  return s ;
}

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
int Derivative::DTermIndx()
{
  int r(0) ;
  for(int i(0) ; i<Num ; i++)
    {
      int t(0);
      if(sh[i].bit()>0)
	t = 1 << i ;
      r+=t ;
    }
  return r ;
}

Derivative::~Derivative()
{
  if (Num>0)
    { 
      delete [] indx;
      delete [] sh ;
    }
  delete [] tag ;
}

CPS_END_NAMESPACE
