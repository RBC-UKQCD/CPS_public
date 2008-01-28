//------------------------------------------------------------------
//
// Nuc3pt_GAMMAR.C
//
// Implementation of the Nucl3pt_GAMMAR  class 
//
// August 2003
//
// Federico Berruto
//
//------------------------------------------------------------------

#include <alg/nuc3pt.h>
#include <alg/nuc3pt_gammar.h>

CPS_START_NAMESPACE

Nuc3ptGammaR::Nuc3ptGammaR(Gamma op, int dd):Nuc3pt(),dir(dd),G(op)
{
  cname = "Nuc3ptGammaR";    
  //char *fname = "Nuc3ptGammaR(Gamma,int)";
  //VRB.Func(cname,fname);
}

/*! 
  Computes
  \f[
  Trace(seqQ \, 
  \gamma_{\mu}\,r_{\nu} \, Quark)
  \f]
**/
void Nuc3ptGammaR::InsertOp(CorrFunc& tmp,QPropW& seqQ, QPropW& Quark)
{
  char *fname = "InsertOp()";
  
  VRB.Func(cname, fname);
  
  Site s ;
  for(s.Begin();s.End();s.nextSite())
    {
      int t(s.physT()) ;
      WilsonMatrix sq(seqQ[s.Index()]) ;
      for(int mu(0);mu<G.N();mu++)
	sq.gr(G[mu]) ; // Multiply by gamma_mu 
      // Multiply by x,y,z (dir=0,1,2)
      tmp[t] += Float(s.physCoor(dir)) * Trace(sq,Quark[s.Index()]) ;
    }
}

void Nuc3ptGammaR::PrintTheTag(FILE *fp)
{
  G.printTag(fp);
  Fprintf(fp,"_R%i",dir+1) ;
}


CPS_END_NAMESPACE
