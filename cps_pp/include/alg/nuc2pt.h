//------------------------------------------------------------------
//
// nuc2pt.h
//
// Header file for the Nucl2pt  class 
//
// April 2001
//
// Kostas Orginos
//
//------------------------------------------------------------------

#ifndef INCLUDED_NUC2PT_H
#define INCLUDED_NUC2PT_H
#include <alg/corrfunc.h>
#include <alg/diquark.h>
#include <alg/qpropw.h>
#include <string.h>

CPS_START_NAMESPACE

enum NucOp {NUC_C, NUC_G5C};

class Nuc2pt
{
  char *cname ;
  char *source ;
  char *sink ;
  char mom[12] ;
  int srcdesc[4] ; // source descriptor
  int HalfFerm ; // Flag indicating Non-Relativistic sources and sinks

  NucOp nuc ;
  SourceType SnkType ;
  Float quark_mass ;
  CorrFunc plus_parity ;
  CorrFunc minus_parity ;

  ProjectType pos_par ;
  ProjectType neg_par ;

  // Warning calcNuclceon overwrites Q2 only
  void calc_nuc_(WilsonMatrix& Q1, WilsonMatrix& Q2, WilsonMatrix& Q3, int t)
    //Q1 and Q2 form the diquark. Q3 is the single quark.
    {

      switch (nuc)
	{
	case NUC_G5C:
	  Q2.ccl(5).ccr(5);
	  Q2.diq(Q1); // compute the diquark!
	  break ;
	case NUC_C:
	  Q2.ccl(1).ccr(-1);
	  Q2.diq(Q1); // compute the diquark!
	  Q3.gl(-5).gr(-5); // This is the correct operator: hep-lat/0202022
	  break ;
	default:
	  ERR.General(cname,"calc_nuc_","Unknown Operator\n") ;
	}	          
      plus_parity[t]  += prop_nucleon(pos_par,Q3,Q2) ;
      minus_parity[t] += prop_nucleon(neg_par,Q3,Q2) ;
    }

 // Warning calcNuclceon overwrites Q2 only
  void calc_nuc_(WilsonMatrix& Q1, WilsonMatrix& Q2, WilsonMatrix& Q3, int t, int tsrc)
    //Q1 and Q2 form the diquark. Q3 is the single quark.
    {

      int Nt(plus_parity.TimeSize());

      switch (nuc)
        {
        case NUC_G5C:
          Q2.ccl(5).ccr(5);
          Q2.diq(Q1); // compute the diquark!
          break ;
        case NUC_C:
          Q2.ccl(1).ccr(-1);
          Q2.diq(Q1); // compute the diquark!
          Q3.gl(-5).gr(-5); // This is the correct operator: hep-lat/0202022
          break ;
        default:
          ERR.General(cname,"calc_nuc_","Unknown Operator\n") ;
        }          
      plus_parity[(t+Nt-tsrc)%Nt]  += prop_nucleon(pos_par,Q3,Q2) ;
      minus_parity[(t+Nt-tsrc)%Nt] += prop_nucleon(neg_par,Q3,Q2) ;
    }

  // Warning calcNuclceon overwrites Q2 only
  void calc_nuc_(WilsonMatrix& Q1, WilsonMatrix& Q2, WilsonMatrix& Q3, WilsonMatrix& Q4, WilsonMatrix& Q5, int t, int dohalf)
    //Q1 and Q2 form the diquark. Q3 is the single quark.
    {

      Diquark diq;
      WilsonVector S;
      WilsonMatrix Dq, Uq;
      int Nspin(4);
      if(dohalf) Nspin=2;
      int nspn2;

      switch (nuc)
	{
	case NUC_G5C:
	  Q2.ccl(5).ccr(5);
	  Q3.ccl(5);
	  Q4.ccr(5);
	  for(int spin=0; spin<Nspin; spin++)
	  for(int color=0; color<3; color++){
	  diq.D_diquark(Q2, Q1, Q3, Q4, spin, color);
	  diq.Project(S,PPAR);
	  Dq.load_vec(spin,color,S);
	  if(dohalf){
	    nspn2=spin+2;
	    Dq.load_vec(nspn2,color,S);
	  }	  
	  diq.U_diquark(Q2, Q1, spin, color);
	  diq.Project(S,PPAR);
	  Uq.load_vec(spin,color,S);
	  if(dohalf){
	    nspn2=spin+2;
	    Uq.load_vec(nspn2,color,S);
	  }	  
	  }
	  break ;
	case NUC_C:
	  Q2.ccl(1).ccr(-1);
	  Q2.diq(Q1); // compute the diquark!
	  Q3.gl(-5).gr(-5); // This is the correct operator: hep-lat/0202022
	  break ;
	default:
	  ERR.General(cname,"calc_nuc_","Unknown Operator\n") ;
	}	          
      plus_parity[t]  += Trace(Q5,Uq)/2;
      minus_parity[t] += Trace(Q5,Dq);
    }

public:
  //Zero Momentum nucleon 
  Nuc2pt(NucOp op, SourceType snk) ; 
  Nuc2pt(NucOp op, SourceType snk, ProjectType) ; 
 
  void calcNucleon(QPropW& quark) ;
  //for non-degenerate quarks
  void calcNucleon(QPropW& quark, QPropW& quark2) ;

  //Non-Zero Momentum nucleon from wall sources
  void calcWallMomNucleon(QPropW& quark, QPropWMomSrc& mom_quark) ;

  //Non-Zero Momentum nucleon from any source but wall
  void calcMomNucleon(QPropW& quark, ThreeMom& mom) ;


  //Zero Momentum nucleon with HalfFerm
  void calcNucleon(QPropW& quark, int dohalf) ;
  //Non-Zero Momentum nucleon from any source but wall with HalfFerm
  void calcMomNucleon(QPropW& quark, ThreeMom& mom, int dohalf) ;


  Complex prop_nucleon(int dd, const WilsonMatrix& p1, const WilsonMatrix& p2);
  Complex prop_nucleon(int dd, const WilsonMatrix& p1, const WilsonMatrix& p2, const WilsonMatrix& p3);

  void Zero()
    {
      plus_parity.Zero();
      minus_parity.Zero();
    }

  ~Nuc2pt(){}
  
  void Print(FILE *fp) ;

};

CPS_END_NAMESPACE

#endif //!INCLUDED_NUC2PT_H





