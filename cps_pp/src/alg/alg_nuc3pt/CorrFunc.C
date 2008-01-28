//------------------------------------------------------------------
//
// CorrFunc.C
//
// Implementation of class CorrFunc
// 
// July 2001
//
// Kostas Orginos
//
//------------------------------------------------------------------

#include <alg/corrfunc.h>

CPS_START_NAMESPACE

CorrFunc::CorrFunc()
{
  char *cname = "CorrFunc";
  char *fname = "CorrFunc()";

  VRB.Func(cname, fname);

  Nt=GJP.Tnodes()*GJP.TnodeSites();
  func=(Complex*)smalloc(Nt*sizeof(Complex));
  if(func == 0) ERR.Pointer(cname,fname, "trace");
  VRB.Smalloc(cname,fname, "func", func, Nt*sizeof(Complex));
  // Zero the correlation function
  for(int t = 0 ; t<Nt; t++) func[t] = 0.0 ;
}

//Copy constuctor 
CorrFunc::CorrFunc(const CorrFunc& rhs)
{
  char *cname = "CorrFunc";
  char *fname = "CorrFunc(const CorrFunc&)";
  VRB.Func(cname, fname);
  Nt=rhs.Nt ;
  func=(Complex*)smalloc(Nt*sizeof(Complex));
  if(func == 0) ERR.Pointer(cname,fname, "trace");
  VRB.Smalloc(cname,fname, "func", func, Nt*sizeof(Complex));
  for(int t = 0 ; t<Nt; t++) func[t] = rhs.func[t] ; ;
}

CorrFunc& CorrFunc::operator=(const CorrFunc& rhs)
{
  char *cname = "CorrFunc";
  char *fname = "operator=()";
  
  VRB.Func(cname, fname);
  if(Nt != rhs.Nt )
    {
      Nt = rhs.Nt ;
      sfree(func) ;
      func=(Complex*)smalloc(Nt*sizeof(Complex));
      if(func == 0) ERR.Pointer(cname,fname, "trace");
      VRB.Smalloc(cname,fname, "func", func, Nt*sizeof(Complex));
    }
  
  // copy correlation function
  for(int t = 0 ; t<Nt; t++) func[t] = rhs.func[t] ;
  return *this ;
}

CorrFunc& CorrFunc::operator+=(const CorrFunc& rhs)
{
  char *cname = "CorrFunc";
  char *fname = "operator+=()";
  
  VRB.Func(cname, fname);
  if(Nt != rhs.Nt )
    {
      ERR.General(cname,fname,"Correlation can not be added!\n") ;
    }
  
  // add correlation function
  for(int t = 0 ; t<Nt; t++) func[t] += rhs.func[t] ;

  return *this ;
}

CorrFunc& CorrFunc::operator-=(const CorrFunc& rhs)
{
  char *cname = "CorrFunc";
  char *fname = "operator-=()";
  
  VRB.Func(cname, fname);
  if(Nt != rhs.Nt )
    {
      ERR.General(cname,fname,"Correlation can not be subtracted!\n") ;
    }
  
  // subtract correlation function
  for(int t = 0 ; t<Nt; t++) func[t] -= rhs.func[t] ;
  
  return *this ;
}

CorrFunc& CorrFunc::operator*=(const CorrFunc& rhs)
{
  char *cname = "CorrFunc";
  char *fname = "operator*()";
  
  VRB.Func(cname, fname);
  if(Nt != rhs.Nt )
    {
      ERR.General(cname,fname,"Correlation can not be multiplied!\n") ;
    }
  
  // multiply correlation function
  for(int t = 0 ; t<Nt; t++) func[t] *= rhs.func[t] ;

  return *this ;
}

CorrFunc& CorrFunc::operator*=(const Float& r)
{
  char *cname = "CorrFunc";
  char *fname = "operator*(Float)";
  
  VRB.Func(cname, fname);
  
  // multiply correlation function
  for(int t = 0 ; t<Nt; t++) func[t] *= r ;

  return *this ;
}

CorrFunc& CorrFunc::operator*=(const Complex& c)
{
  char *cname = "CorrFunc";
  char *fname = "operator*(Float)";
  
  VRB.Func(cname, fname);
  
  // multiply correlation function
  for(int t = 0 ; t<Nt; t++) func[t] *= c ;

  return *this ;
}

CPS_END_NAMESPACE



