/*! Returns the operator Tag, a string to be used in 3pt IO */
#include <alg/gamma.h>

CPS_START_NAMESPACE

void Gamma::Tag()
{

  if(_n==0){
    tag = new char[4] ;
    sprintf(tag,"G0");
  }
  else{

    //check if the _mu pointer is allocated 
    if(_mu==0) ERR.Pointer("Gamma", "Tag", "_mu");
      
    int *tt = new int[_n] ;
    //check if the tt pointer is allocated 
    if(tt==0) ERR.Pointer("Gamma", "Tag", "tt");
  
    for(int i=0;i<_n;i++)
      {
	if(_mu[i]==-5) 
	  tt[i]=5 ;
	else
	  tt[i]=_mu[i]+1 ;
      }
    
    switch (_n)
      {
      case 1: 
	tag = new char[4] ;
	sprintf(tag,"G%i",tt[0]); 
	break ;
      case 2:
	tag = new char[6] ;
	sprintf(tag,"G%iG%i",tt[0],tt[1]); 
	break ;
      default:
	tag = new char[8] ;
	sprintf(tag,"G%iG%iG%i",tt[0],tt[1],tt[2]);       
      }
    
    delete [] tt ;
  }

  //check if the tag pointer is allocated 
  if(tag==0) ERR.Pointer("Gamma", "Tag", "tag");
  
}

Gamma::Gamma(                   ): // scalar 
  _n(0){Tag() ;} ;

Gamma::Gamma(int m              ): 
  _n(1){_mu = new int[1]; _mu[0]= m; Tag();}

Gamma::Gamma(int m, int n       ): 
  _n(2){_mu = new int[2]; _mu[0]= m; _mu[1]= n; Tag();}

Gamma::Gamma(int m, int n, int r): 
  _n(3){_mu = new int[3]; _mu[0]= m; _mu[1]= n; _mu[2]= r; Tag();} 

/*! Copy constructor */
Gamma::Gamma(const Gamma& tt)
{
  _n=tt._n ;
  if(_n>0){
    _mu = new int [_n] ;
    setIndices(tt._mu) ;
  }
  Tag() ;
}

Gamma::~Gamma(){delete [] tag ; if(_n>0) delete [] _mu ;} 

CPS_END_NAMESPACE
