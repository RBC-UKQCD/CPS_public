#include <alg/nuc2pt.h>

CPS_START_NAMESPACE

Nuc2pt::Nuc2pt(NucOp op, SourceType snk) :  HalfFerm(0), nuc(op), SnkType(snk),
					    pos_par(PPAR),neg_par(NPAR)
{
  cname = "Nuc2pt";    
  char *fname = "Nuc2pt()";
    
  VRB.Func(cname,fname);

}

Nuc2pt::Nuc2pt(NucOp op, SourceType snk, ProjectType proj) : 
  HalfFerm(0), nuc(op), SnkType(snk), pos_par(proj)
{
  cname = "Nuc2pt";    
  char *fname = "Nuc2pt()";
    
  VRB.Func(cname,fname);

  if((pos_par!=PPAR   )&&(pos_par!=PPAR_5X)&&
     (pos_par!=PPAR_5Y)&&(pos_par!=PPAR_5Z)&&
     (pos_par!=PPAR_5 ))
    ERR.General(cname,fname,"Unknown proj.\n");

  neg_par = (ProjectType)((int) pos_par + 1) ;

}

/*!
  Calculate zero momentum nucleon two point functions
*/
void Nuc2pt::calcNucleon(QPropW& quark){
  char *cname = "Nuc2pt" ;
  char *fname = "calcNucleon(QPropW&)" ;

  strcpy(mom,"0 0 0");  

  if(quark.DoHalfFermion())
    HalfFerm=1;

  srcdesc[0]=srcdesc[1]=srcdesc[2]=srcdesc[3]=-1 ;
  srcdesc[4]=srcdesc[5]=-1 ; // TIZB
  switch (quark.SrcType())
    {
    case WALL:
      srcdesc[0] = quark.SourceTime() ;
      source="WALL" ;
      break ;
    case BOX:
      srcdesc[0] = quark.BoxSrcStart() ; 
      srcdesc[1] = quark.BoxSrcEnd() ;
      srcdesc[2] = quark.SourceTime() ;
      if( quark. BoxSrcUseXYZOffset() ){ // use QPropWArg.{x,y,z} as offset
	srcdesc[3] = quark.PointSrcX() ; 
	srcdesc[4] = quark.PointSrcX() ; 
	srcdesc[5] = quark.PointSrcX() ; 
      }
      source="BOX" ;
      break ;
    case POINT:
      srcdesc[0] = quark.PointSrcX() ; 
      srcdesc[1] = quark.PointSrcY() ;
      srcdesc[2] = quark.PointSrcZ() ;
      srcdesc[3] = quark.SourceTime() ;
      source = "POINT" ;
      break ;
    case GAUSS_GAUGE_INV:
      srcdesc[0] = quark.Gauss_N() ; 
      srcdesc[1] = int(100*quark.Gauss_W()) ;
      srcdesc[2] = quark.SourceTime() ;
      source = "GAUSS" ;
      break ;
    default:
      ERR.General(cname,fname,"Bad source type\n") ;
    }
  quark_mass = quark.Mass() ;
    
  WilsonMatrix Q1;  
  WilsonMatrix Q2;
  
    switch (SnkType)
      {
      case POINT:
      case GAUSS_GAUGE_INV:	
	{
	  if(SnkType==GAUSS_GAUGE_INV)
	    sink="GAUSS" ;
	  else
	    sink="POINT" ;
	  Site s ;
	  for(s.Begin();s.End();s.nextSite())
	    {
	      int t = s.physT() ;
	      Q1 = quark[s.Index()] ;
	      if(HalfFerm) Q1.PParProjectSink() ;
	      Q2 = Q1 ; 
	      calc_nuc_(Q1, Q2, Q1, t) ;	  
	    }
	  plus_parity.GlobalSum() ;
	  minus_parity.GlobalSum() ;
	}
	break ;
      case WALL:
	{
	  sink="WALL" ;

	  if(!quark.GFixedSnk())
	    ERR.General(cname,fname,"Wall sink needs gauge fixing\n") ;

	  for(int t = 0 ; t<plus_parity.TimeSize();t++)
	    {
	      Q1 =  quark.WallSinkProp(t) ;
	      if(HalfFerm) Q1.PParProjectSink() ;
	      Q2 = Q1 ;
	      calc_nuc_(Q1, Q2, Q1, t) ;
	    }
	}
	break ;
      default:
	ERR.General(cname,fname,"Unknown Sink Type\n") ;
      } 
}

/*!
  Calculate zero momentum nucleon two point functions
  for non-degenerate quarks.
*/
void Nuc2pt::calcNucleon(QPropW& upquark, QPropW& downquark){
  char *cname = "Nuc2pt" ;
  char *fname = "calcNucleon(QPropW&, QPropW&)" ;

  strcpy(mom,"0 0 0");

  if(upquark.DoHalfFermion())
    HalfFerm=1;

  srcdesc[0]=srcdesc[1]=srcdesc[2]=srcdesc[3]=-1 ;
  srcdesc[4]=srcdesc[5]=-1 ; // TIZB
  switch (downquark.SrcType())
    {
    case WALL:
      srcdesc[0] = downquark.SourceTime() ;
      source="WALL" ;
      break ;
    case BOX:
      srcdesc[0] = downquark.BoxSrcStart() ; 
      srcdesc[1] = downquark.BoxSrcEnd() ;
      srcdesc[2] = downquark.SourceTime() ;
      if( downquark. BoxSrcUseXYZOffset() ){ // use QPropWArg.{x,y,z} as offset
	srcdesc[3] = downquark.PointSrcX() ; 
	srcdesc[4] = downquark.PointSrcX() ; 
	srcdesc[5] = downquark.PointSrcX() ; 
      }
      source="BOX" ;
      break ;
    case POINT:
      srcdesc[0] = downquark.PointSrcX() ; 
      srcdesc[1] = downquark.PointSrcY() ;
      srcdesc[2] = downquark.PointSrcZ() ;
      srcdesc[3] = downquark.SourceTime() ;
      source = "POINT" ;
      break ;
    case GAUSS_GAUGE_INV:
      srcdesc[0] = downquark.Gauss_N() ; 
      srcdesc[1] = int(100*downquark.Gauss_W()) ;
      srcdesc[2] = downquark.SourceTime() ;
      source = "GAUSS" ;
      break ;
    default:
      ERR.General(cname,fname,"Bad source type\n") ;
    }

  quark_mass = downquark.Mass() ;
  int tsrc = downquark.SourceTime();

  if(downquark.SrcType() != upquark.SrcType())
    ERR.General(cname,fname,"Up and down source types differ\n") ;
  if(upquark.SourceTime() != downquark.SourceTime())
    ERR.General(cname,fname,"Up and down source times differ\n") ;

  WilsonMatrix Q1;  
  WilsonMatrix Q2;
  
    switch (SnkType)
      {
      case POINT:
        {
          sink="POINT" ;
          Site s ;
          for(s.Begin();s.End();s.nextSite())
            {
              int t = s.physT() ;
              Q1 = upquark[s.Index()] ;
              if(HalfFerm) Q1.PParProjectSink() ;
              Q2 = downquark[s.Index()] ;
              if(HalfFerm) Q2.PParProjectSink() ;
              calc_nuc_(Q1, Q2, Q1, t, tsrc) ;    
            } 
          plus_parity.GlobalSum() ;
          minus_parity.GlobalSum() ;
        }
        break ;
      case WALL:
        {
          sink="WALL" ;

          if(!upquark.GFixedSnk() || !downquark.GFixedSnk())
            ERR.General(cname,fname,"Wall sink needs gauge fixing\n") ;

          for(int t = 0 ; t<plus_parity.TimeSize();t++)
            {
              Q1 =  upquark.WallSinkProp(t) ;
              if(HalfFerm) Q1.PParProjectSink() ;
              Q2 =  downquark.WallSinkProp(t) ;
              if(HalfFerm) Q2.PParProjectSink() ;
              calc_nuc_(Q1, Q2, Q1, t, tsrc) ;
            }
        }
        break ;
      default:
        ERR.General(cname,fname,"Unknown Sink Type\n") ;
      } 
}

/*!
  Non-Zero Momentum nucleon  from wall sources
*/
void Nuc2pt::calcWallMomNucleon(QPropW& quark, QPropWMomSrc& mom_quark)
{
  char *fname = "Nuc2pt()";
  ThreeMom sink_mom(mom_quark.Mom()) ;
  sink_mom.conj();
  //The momentum on the sink has to be complex conjugated!!
  // p --> - p 
  
  srcdesc[0]=srcdesc[1]=srcdesc[2]=srcdesc[3]=-1 ;
  srcdesc[4]=srcdesc[5]=-1 ; // TIZB
  switch (quark.SrcType())
    {
    case WALL:
      srcdesc[0] = quark.SourceTime() ;
      source="WALL" ;
      break ;
    default:
      ERR.General(cname,fname,"Bad source type\n") ;
    }
  
  quark_mass = quark.Mass() ;
    
  WilsonMatrix Q1;  
  WilsonMatrix Q2;
  WilsonMatrix Q3;  
  
  SnkType=POINT;
  sprintf(mom,"%i %i %i", mom_quark.Mom().cmp(0), 
	  mom_quark.Mom().cmp(1), mom_quark.Mom().cmp(2) );
  Site s ;
  for(s.Begin();s.End();s.nextSite())
    {
      int t = s.physT() ;
      Q1 =  quark[s.Index()] ;  
      if(HalfFerm) Q1.PParProjectSink() ;
      Q2 = Q1 ;
      // here is the sink so I need exp(ipx) 
      // I did conjugate the momentum at the construction
      Q3 = sink_mom.Fact(s)*mom_quark[s.Index()]  ; 
      if(HalfFerm) Q3.PParProjectSink() ;
      calc_nuc_(Q1, Q2, Q3,t) ;	  
    }
  plus_parity.GlobalSum() ;
  minus_parity.GlobalSum() ;
}

/*!
  Non-Zero Momentum nucleon from box sources or point
*/
void Nuc2pt::calcMomNucleon(QPropW& quark, ThreeMom& Mom)
{
  char *fname = "Nuc2pt()";
  ThreeMom sink_mom(Mom) ;
  sink_mom.conj();
  //The momentum on the sink has to be complex conjugated!!
  // p --> - p 
  
  srcdesc[0]=srcdesc[1]=srcdesc[2]=srcdesc[3]=-1 ;
  srcdesc[4]=srcdesc[5]=-1 ; // TIZB
  switch (quark.SrcType())
    {
    case BOX:
      srcdesc[0] = quark.BoxSrcStart() ; 
      srcdesc[1] = quark.BoxSrcEnd() ;
      srcdesc[2] = quark.SourceTime() ;
      source="BOX" ;
      if( quark. BoxSrcUseXYZOffset() ){ // use QPropWArg.{x,y,z} as offset
	srcdesc[3] = quark.PointSrcX() ; 
	srcdesc[4] = quark.PointSrcX() ; 
	srcdesc[5] = quark.PointSrcX() ; 
      }
      break ;
    case POINT:
      srcdesc[0] = quark.PointSrcX() ; 
      srcdesc[1] = quark.PointSrcY() ;
      srcdesc[2] = quark.PointSrcZ() ;
      srcdesc[3] = quark.SourceTime() ;
      source="POINT" ;
      break ;
    case GAUSS_GAUGE_INV:
      srcdesc[0] = quark.Gauss_N() ; 
      srcdesc[1] = int(100*quark.Gauss_W()) ;
      srcdesc[2] = quark.SourceTime() ;
      source = "GAUSS" ;
      break ;
    default:
      ERR.General(cname,fname,"Bad source type\n") ;
    }
  
  quark_mass = quark.Mass() ;
  
  WilsonMatrix Q1;  
  WilsonMatrix Q2;
  WilsonMatrix Q3;  
  
  //SnkType=POINT;
  sprintf(mom,"%i %i %i", Mom.cmp(0), Mom.cmp(1), Mom.cmp(2) );

  Site s ;

  for(s.Begin();s.End();s.nextSite())
    {
      int t = s.physT() ;
      Q1 = quark[s.Index()] ;
      if(HalfFerm) Q1.PParProjectSink() ;
      Q2 = Q1 ; 
      // here is the sink so I need exp(ipx) 
      // I did conjugate the momentum at the construction
      Q3 = Mom.Fact(s)*Q2 ;
      if(HalfFerm) Q3.PParProjectSink() ;
      calc_nuc_(Q1, Q2, Q3, t) ;	  
    }
  plus_parity.GlobalSum() ;
  minus_parity.GlobalSum() ;
}




// TY add START
/*!
  Calculate zero momentum nucleon two point functions
*/
void Nuc2pt::calcNucleon(QPropW& quark, int dohalf){
  char *cname = "Nuc2pt" ;
  char *fname = "calcNucleon(QPropW&)" ;

  strcpy(mom,"0 0 0");  

  HalfFerm=dohalf;

  srcdesc[0]=srcdesc[1]=srcdesc[2]=srcdesc[3]=-1 ;
  srcdesc[4]=srcdesc[5]=-1 ; // TIZB
  switch (quark.SrcType())
    {
    case WALL:
      srcdesc[0] = quark.SourceTime() ;
      source="WALL" ;
      break ;
    case BOX:
      srcdesc[0] = quark.BoxSrcStart() ; 
      srcdesc[1] = quark.BoxSrcEnd() ;
      srcdesc[2] = quark.SourceTime() ;
      if( quark. BoxSrcUseXYZOffset() ){ // use QPropWArg.{x,y,z} as offset
	srcdesc[3] = quark.PointSrcX() ; 
	srcdesc[4] = quark.PointSrcX() ; 
	srcdesc[5] = quark.PointSrcX() ; 
      }
      source="BOX" ;
      break ;
    case POINT:
      srcdesc[0] = quark.PointSrcX() ; 
      srcdesc[1] = quark.PointSrcY() ;
      srcdesc[2] = quark.PointSrcZ() ;
      srcdesc[3] = quark.SourceTime() ;
      source = "POINT" ;
      break ;
    case GAUSS_GAUGE_INV:
      srcdesc[0] = quark.Gauss_N() ; 
      srcdesc[1] = int(100*quark.Gauss_W()) ;
      srcdesc[2] = quark.SourceTime() ;
      source = "GAUSS" ;
      break ;
    default:
      ERR.General(cname,fname,"Bad source type\n") ;
    }
  quark_mass = quark.Mass() ;
    
  WilsonMatrix Q1;  
  WilsonMatrix Q2;
  WilsonMatrix Q3;
  WilsonMatrix Q4;
  WilsonMatrix Q5;
  
    switch (SnkType)
      {
      case POINT:
      case GAUSS_GAUGE_INV:	
	{
	  if(SnkType==GAUSS_GAUGE_INV)
	    sink="GAUSS" ;
	  else
	    sink="POINT" ;
	  Site s ;
	  for(s.Begin();s.End();s.nextSite())
	    {
	      int t = s.physT() ;
	      Q1 = quark[s.Index()] ;
	      if(HalfFerm) Q1.PParProjectSink() ;
	      Q2 = Q1 ; 
	      Q3 = Q1 ; 
	      Q4 = Q1 ; 
	      Q5 = Q1 ; 
	      calc_nuc_(Q1, Q2, Q3, Q4, Q5, t, dohalf) ;
	    }
	  plus_parity.GlobalSum() ;
	  minus_parity.GlobalSum() ;
	}
	break ;
      case WALL:
	{
	  sink="WALL" ;

	  if(!quark.GFixedSnk())
	    ERR.General(cname,fname,"Wall sink needs gauge fixing\n") ;

	  for(int t = 0 ; t<plus_parity.TimeSize();t++)
	    {
	      Q1 =  quark.WallSinkProp(t) ;
	      if(HalfFerm) Q1.PParProjectSink() ;
	      Q2 = Q1 ;
	      calc_nuc_(Q1, Q2, Q1, t) ;
	    }
	}
	break ;
      default:
	ERR.General(cname,fname,"Unknown Sink Type\n") ;
      } 
}

/*!
  Non-Zero Momentum nucleon from box sources or point
*/
void Nuc2pt::calcMomNucleon(QPropW& quark, ThreeMom& Mom, int dohalf)
{
  char *fname = "Nuc2pt()";
  ThreeMom sink_mom(Mom) ;
  sink_mom.conj();
  //The momentum on the sink has to be complex conjugated!!
  // p --> - p 

  HalfFerm=dohalf;
  
  srcdesc[0]=srcdesc[1]=srcdesc[2]=srcdesc[3]=-1 ;
  srcdesc[4]=srcdesc[5]=-1 ; // TIZB
  switch (quark.SrcType())
    {
    case BOX:
      srcdesc[0] = quark.BoxSrcStart() ; 
      srcdesc[1] = quark.BoxSrcEnd() ;
      srcdesc[2] = quark.SourceTime() ;
      if( quark. BoxSrcUseXYZOffset() ){ // use QPropWArg.{x,y,z} as offset
	srcdesc[3] = quark.PointSrcX() ; 
	srcdesc[4] = quark.PointSrcX() ; 
	srcdesc[5] = quark.PointSrcX() ; 
      }
      source="BOX" ;
      break ;
    case POINT:
      srcdesc[0] = quark.PointSrcX() ; 
      srcdesc[1] = quark.PointSrcY() ;
      srcdesc[2] = quark.PointSrcZ() ;
      srcdesc[3] = quark.SourceTime() ;
      source="POINT" ;
      break ;
    case GAUSS_GAUGE_INV:
      srcdesc[0] = quark.Gauss_N() ; 
      srcdesc[1] = int(100*quark.Gauss_W()) ;
      srcdesc[2] = quark.SourceTime() ;
      source = "GAUSS" ;
      break ;
    default:
      ERR.General(cname,fname,"Bad source type\n") ;
    }
  
  quark_mass = quark.Mass() ;
  
  WilsonMatrix Q1;  
  WilsonMatrix Q2;
  WilsonMatrix Q3;  
  WilsonMatrix Q4;
  WilsonMatrix Q5;  
  
  //SnkType=POINT;
  sprintf(mom,"%i %i %i", Mom.cmp(0), Mom.cmp(1), Mom.cmp(2) );

  Site s ;

  for(s.Begin();s.End();s.nextSite())
    {
      int t = s.physT() ;
      Q1 = quark[s.Index()] ;
      if(HalfFerm) Q1.PParProjectSink() ;
      Q2 = Q1 ; 
      Q3 = Q1 ; 
      Q4 = Q1 ; 
      // here is the sink so I need exp(ipx) 
      // I did conjugate the momentum at the construction
      Q5 = Mom.Fact(s)*Q2 ;
      calc_nuc_(Q1, Q2, Q3, Q4, Q5, t, dohalf) ;
    }
  plus_parity.GlobalSum() ;
  minus_parity.GlobalSum() ;
}
// TY add END


void Nuc2pt::Print(FILE *fp)
{
  if(fp==NULL) return ;
  char *oper ;
  char *proj ;
  char *sink ;

  switch (nuc)
    {
    case NUC_G5C:
      oper="NUC_G5C" ;
      break ;
    case NUC_C:
      oper="NUC_C"  ;
      break ;
    default:
      oper="U" ; //just keep gcc happy
    }

  switch (pos_par)
    {
    case PPAR_5X:
      proj="5X" ;
      break ;
    case PPAR_5Y:
      proj="5Y"  ;
      break ;
    case PPAR_5Z:
      proj="5Z"  ;
      break ;
    case PPAR_5:
      proj="5"  ;
      break ;
    default:
      proj="" ; //if its PPAR don't print anything
    }

  switch (SnkType)
    {
    case POINT:
      sink="POINT" ;
      break ;
    case WALL:
      sink="WALL" ;
      break ;
    case GAUSS_GAUGE_INV:
      sink = "GAUSS" ;
      break ;
    default:
      sink="" ; // Will never happen!
    } 


  Fprintf(fp,"STARTPROP\n") ;
  Fprintf(fp,"MASSES:  %e %e %e\n",
	  (float)quark_mass,(float)quark_mass,(float)quark_mass) ;
  Fprintf(fp,"SOURCE: %s ",source) ;
  for(int i(0) ; i<6;i++) // TIZB
    if(srcdesc[i]>-1)
      Fprintf(fp,"%i ",srcdesc[i]) ;
  Fprintf(fp,"\n" ) ;
  Fprintf(fp,"SINK: %s\n",sink ) ;
  if(HalfFerm) Fprintf(fp,"NR-SPINORS\n") ;
  Fprintf(fp,"MOM: %s\n",mom ) ;
  Fprintf(fp,"OPER: %s_PP%s %s_NP%s\n",oper,proj,oper,proj) ;
  for(int t=0; t<plus_parity.TimeSize();t++)
    {
      Fprintf(fp,"%i  ",t) ;
      plus_parity.print(fp,t) ;
      minus_parity.print(fp,t) ;
      Fprintf(fp,"\n") ;
    }
  Fprintf(fp,"ENDPROP\n");
}

/*!
  Nucleon state: 
  B(i,j) = M_u(i,a;j,a') * D(k,a';k,a) + M_u(i,a;k,a') * D(k,a';j,a)
*/
Complex Nuc2pt::prop_nucleon(int dd, const WilsonMatrix& p1, 
			     const WilsonMatrix& p2)
{

  SpinMatrix tr=(Float)0.0;
  Complex res(0.0,0.0);
        
  for(int s3=0;s3<4;++s3){
    for(int s1=0;s1<4;++s1){
      for(int c1=0;c1<3;++c1){
	for(int s2=0;s2<4;++s2){
	  for(int c2=0;c2<3;++c2){
	    tr(s1,s3) 
	      +=  p1.wmat().d[s1].c[c1].d[s3].c[c2]
                * p2.wmat().d[s2].c[c2].d[s2].c[c1]
		+ p1.wmat().d[s1].c[c1].d[s2].c[c2]
		* p2.wmat().d[s2].c[c2].d[s3].c[c1];
	  }
	}
      }
    }
  }
// Projection of baryon_propagator or vector and scalar 3pt_Function
//  Tr{ (1 + gamma4) x Baryon_Propagator } for
//    dd = 0 i.e. decay direction along t
//  Tr{ (1 - gamma4) x Baryon_Propagator } for
//    dd = 1 i.e. decay direction along t
// Projection for axial and tensor 3pt_function
//  Tr{ (1 + gamma4) x (i*gamma5*gammaZ) x 3pt_Function } for
//    dd = 2 i.e. decay direction along t
//  Tr{ (1 - gamma4) x (i*gamma5*gammaZ) x 3pt_Function } for
//    dd = 3 i.e. decay direction along t
// Projection for axial and tensor 3pt_function
//  Tr{ (1 + gamma4) x (gamma5*gammaY) x 3pt_Function } for
//    dd = 4 i.e. decay direction along t
//  Tr{ (1 - gamma4) x (gamma5*gammaY) x 3pt_Function } for
//    dd = 5 i.e. decay direction along t
// Projection for axial and tensor 3pt_function
//  Tr{ (1 + gamma4) x (i*gamma5*gammaX) x 3pt_Function } for
//    dd = 6 i.e. decay direction along t
//  Tr{ (1 - gamma4) x (i*gamma5*gammaX) x 3pt_Function } for
//    dd = 7 i.e. decay direction along t
//  Tr{ (1 + gamma4) gamma5 x Baryon_Propagator } for
//    dd = 10 i.e. decay direction along t
//  Tr{ (1 - gamma4) gamma5 x Baryon_Propagator } for
//    dd = 11 i.e. decay direction along t
  switch (dd) {
  case 0:
    res =  tr(0,0)+tr(1,1)+tr(2,2)+tr(3,3)
          +tr(3,1)+tr(2,0)+tr(1,3)+tr(0,2);
    break;
  case 1:
    res =  tr(0,0)+tr(1,1)+tr(2,2)+tr(3,3)
          -tr(3,1)-tr(2,0)-tr(1,3)-tr(0,2);
    break;
  case 2:
    res = -tr(0,2)+tr(1,3)-tr(2,0)+tr(3,1)
          -tr(0,0)+tr(1,1)-tr(2,2)+tr(3,3);
    break;
  case 3:
    res = -tr(0,2)+tr(1,3)-tr(2,0)+tr(3,1)
          +tr(0,0)-tr(1,1)+tr(2,2)-tr(3,3);
    break;
  case 4:
    res =  tr(0,3)-tr(1,2)+tr(2,1)-tr(3,0)
          +tr(0,1)-tr(1,0)+tr(2,3)-tr(3,2);
    break;
  case 5:
    res =  tr(0,3)-tr(1,2)+tr(2,1)-tr(3,0)
          -tr(0,1)+tr(1,0)-tr(2,3)+tr(3,2);
    break;
  case 6:
    res = -tr(0,3)-tr(1,2)-tr(2,1)-tr(3,0)
          -tr(0,1)-tr(1,0)-tr(2,3)-tr(3,2);
    break;
  case 7:
    res = -tr(0,3)-tr(1,2)-tr(2,1)-tr(3,0)
          +tr(0,1)+tr(1,0)+tr(2,3)+tr(3,2);
    break;
  case 10:
    res =  tr(0,0)+tr(1,1)-tr(2,2)-tr(3,3)
          -tr(3,1)-tr(2,0)+tr(1,3)+tr(0,2);
    break ;
  case 11:
    res =  tr(0,0)+tr(1,1)-tr(2,2)-tr(3,3)
          +tr(3,1)+tr(2,0)-tr(1,3)-tr(0,2);
    break ;
  default: //never reached
    ERR.General(cname, "prop_nuc", "invalid number\n");
  }
  return res;
}

/*!
  Nucleon state: 
  B(i,j) = M_u(i,a;j,a') * D(k,a';k,a) + M_u(i,a;k,a') * D(k,a';j,a)
*/
Complex Nuc2pt::prop_nucleon(int dd, const WilsonMatrix& p1, 
			     const WilsonMatrix& p2, const WilsonMatrix& p3)
{

  SpinMatrix tr=(Float)0.0;
  Complex res(0.0,0.0);
        
  for(int s3=0;s3<4;++s3){
    for(int s1=0;s1<4;++s1){
      for(int c1=0;c1<3;++c1){
	for(int s2=0;s2<4;++s2){
	  for(int c2=0;c2<3;++c2){
	    tr(s1,s3) 
	      +=  p1.wmat().d[s1].c[c1].d[s3].c[c2]
                * p2.wmat().d[s2].c[c2].d[s2].c[c1]
		+ p1.wmat().d[s1].c[c1].d[s2].c[c2]
		* p3.wmat().d[s2].c[c2].d[s3].c[c1];
	  }
	}
      }
    }
  }
// Projection of baryon_propagator or vector and scalar 3pt_Function
//  Tr{ (1 + gamma4) x Baryon_Propagator } for
//    dd = 0 i.e. decay direction along t
//  Tr{ (1 - gamma4) x Baryon_Propagator } for
//    dd = 1 i.e. decay direction along t
// Projection for axial and tensor 3pt_function
//  Tr{ (1 + gamma4) x (i*gamma5*gammaZ) x 3pt_Function } for
//    dd = 2 i.e. decay direction along t
//  Tr{ (1 - gamma4) x (i*gamma5*gammaZ) x 3pt_Function } for
//    dd = 3 i.e. decay direction along t
// Projection for axial and tensor 3pt_function
//  Tr{ (1 + gamma4) x (gamma5*gammaY) x 3pt_Function } for
//    dd = 4 i.e. decay direction along t
//  Tr{ (1 - gamma4) x (gamma5*gammaY) x 3pt_Function } for
//    dd = 5 i.e. decay direction along t
// Projection for axial and tensor 3pt_function
//  Tr{ (1 + gamma4) x (i*gamma5*gammaX) x 3pt_Function } for
//    dd = 6 i.e. decay direction along t
//  Tr{ (1 - gamma4) x (i*gamma5*gammaX) x 3pt_Function } for
//    dd = 7 i.e. decay direction along t
//  Tr{ (1 + gamma4) gamma5 x Baryon_Propagator } for
//    dd = 10 i.e. decay direction along t
//  Tr{ (1 - gamma4) gamma5 x Baryon_Propagator } for
//    dd = 11 i.e. decay direction along t
  switch (dd) {
  case 0:
    res =  tr(0,0)+tr(1,1)+tr(2,2)+tr(3,3)
          +tr(3,1)+tr(2,0)+tr(1,3)+tr(0,2);
    break;
  case 1:
    res =  tr(0,0)+tr(1,1)+tr(2,2)+tr(3,3)
          -tr(3,1)-tr(2,0)-tr(1,3)-tr(0,2);
    break;
  case 2:
    res = -tr(0,2)+tr(1,3)-tr(2,0)+tr(3,1)
          -tr(0,0)+tr(1,1)-tr(2,2)+tr(3,3);
    break;
  case 3:
    res = -tr(0,2)+tr(1,3)-tr(2,0)+tr(3,1)
          +tr(0,0)-tr(1,1)+tr(2,2)-tr(3,3);
    break;
  case 4:
    res =  tr(0,3)-tr(1,2)+tr(2,1)-tr(3,0)
          +tr(0,1)-tr(1,0)+tr(2,3)-tr(3,2);
    break;
  case 5:
    res =  tr(0,3)-tr(1,2)+tr(2,1)-tr(3,0)
          -tr(0,1)+tr(1,0)-tr(2,3)+tr(3,2);
    break;
  case 6:
    res = -tr(0,3)-tr(1,2)-tr(2,1)-tr(3,0)
          -tr(0,1)-tr(1,0)-tr(2,3)-tr(3,2);
    break;
  case 7:
    res = -tr(0,3)-tr(1,2)-tr(2,1)-tr(3,0)
          +tr(0,1)+tr(1,0)+tr(2,3)+tr(3,2);
    break;
  case 10:
    res =  tr(0,0)+tr(1,1)-tr(2,2)-tr(3,3)
          -tr(3,1)-tr(2,0)+tr(1,3)+tr(0,2);
    break ;
  case 11:
    res =  tr(0,0)+tr(1,1)-tr(2,2)-tr(3,3)
          +tr(3,1)+tr(2,0)-tr(1,3)-tr(0,2);
    break ;
  default: //never reached
    ERR.General(cname, "prop_nuc", "invalid number\n");
  }
  return res;
}

CPS_END_NAMESPACE
