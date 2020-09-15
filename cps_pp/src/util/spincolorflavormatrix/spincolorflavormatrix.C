#include <stdio.h>
#include <alg/lanc_arg.h>
#include <util/spincolorflavormatrix.h>

CPS_START_NAMESPACE

//We can use the G-parity prop-conj relation to flip the source momentum a posteriori providing the source has the form  e^{-ipx} \eta    where \eta obeys  \sigma_2 C\gamma^5 \eta^* = \eta\sigma_2 C\gamma^5
//This applies to pretty much all standard source types. However *this code cannot check this is true, it just assumes it!!*
void SpinColorFlavorMatrix::flipSourceMomentum(){
  //Note   A.ccl(-1) = C A
  this->cconj();
  this->ccl(-1).gl(-5).pl(sigma2);
  this->ccr(-1).gr(-5).pr(sigma2);
}

// void SpinColorFlavorMatrix::generate(QPropWcontainer &from_f0, QPropWcontainer &from_f1, Lattice &lattice, const int &site, const PropSplane &splane){
//   //Generate from a pair of propagators with different source flavor but otherwise identical    
//   wmat[0][0] = getSite(site,0,from_f0,lattice,splane);
//   wmat[1][0] = getSite(site,1,from_f0,lattice,splane);
//   wmat[0][1] = getSite(site,0,from_f1,lattice,splane);
//   wmat[1][1] = getSite(site,1,from_f1,lattice,splane);
// }

void SpinColorFlavorMatrix::generate(QPropWcontainer &from_f0, QPropWcontainer &from_f1, Lattice &lattice, const int site, const PropSplane splane){
  generate(from_f0.getProp(lattice), from_f1.getProp(lattice), site, splane);
}

void SpinColorFlavorMatrix::generate(QPropW &from_f0, QPropW &from_f1, const int site, const PropSplane splane){
  //Generate from a pair of propagators with different source flavor but otherwise identical    
  wmat[0][0] = getSite(site,0,from_f0,splane);
  wmat[1][0] = getSite(site,1,from_f0,splane);
  wmat[0][1] = getSite(site,0,from_f1,splane);
  wmat[1][1] = getSite(site,1,from_f1,splane);
}



void SpinColorFlavorMatrix::generate_from_cconj_pair(QPropWcontainer &from, QPropWcontainer &from_conj, Lattice &lattice, const int site, const PropSplane splane){
  int flav = from.flavor();
  if(from_conj.flavor()!=flav){
    ERR.General(cname,"generate_from_cconj_pair","Requires partner 'from_conj' to have same flavor");
  }
  generate_from_cconj_pair(from.getProp(lattice), from_conj.getProp(lattice), flav, site, splane);
}

void SpinColorFlavorMatrix::generate_from_cconj_pair(QPropW &from, QPropW &from_conj, const int from_flav, const int site, const PropSplane splane){
  //Use the propagator conjugate relation to compute the 2x2 flavor matrix propagator using a  propagator from a single flavor and one with the complex conjugate of the source used for the first    

  if(from_flav == 0){
    wmat[0][0] = getSite(site,0,from,splane);
    wmat[1][0] = getSite(site,1,from,splane);

    wmat[0][1] = getSite(site,1,from_conj,splane);
    wmat[1][1] = getSite(site,0,from_conj,splane);
      
    wmat[0][1].cconj();
    wmat[1][1].cconj();
      
    wmat[0][1].ccl(-1).gl(-5).ccr(1).gr(-5); //ccl(-1) mults by C from left, ccr(1) mults by C=-C^dag from right
    wmat[1][1].ccl(-1).gl(-5).ccr(-1).gr(-5); //has opposite sign to the above 
  }else{
    wmat[0][1] = getSite(site,0,from,splane);
    wmat[1][1] = getSite(site,1,from,splane);

    wmat[1][0] = getSite(site,0,from_conj,splane);
    wmat[0][0] = getSite(site,1,from_conj,splane);

    wmat[1][0].cconj();
    wmat[0][0].cconj();

    wmat[1][0].ccl(-1).gl(-5).ccr(1).gr(-5); 
    wmat[0][0].ccl(-1).gl(-5).ccr(-1).gr(-5);
  }

}






// void SpinColorFlavorMatrix::generate_from_cconj_pair(QPropWcontainer &from, QPropWcontainer &from_conj, Lattice &lattice, const int &site, const PropSplane &splane){
//   //Use the propagator conjugate relation to compute the 2x2 flavor matrix propagator using a  propagator from a single flavor and one with the complex conjugate of the source used for the first    
//   int flav = from.flavor();
//   if(from_conj.flavor()!=flav){
//     ERR.General(cname,"generate_from_cconj_pair","Requires partner 'from_conj' to have same flavor");
//   }

//   if(flav == 0){
//     wmat[0][0] = getSite(site,0,from,lattice,splane);
//     wmat[1][0] = getSite(site,1,from,lattice,splane);

//     wmat[0][1] = getSite(site,1,from_conj,lattice,splane);
//     wmat[1][1] = getSite(site,0,from_conj,lattice,splane);
      
//     wmat[0][1].cconj();
//     wmat[1][1].cconj();
      
//     wmat[0][1].ccl(-1).gl(-5).ccr(1).gr(-5); //ccl(-1) mults by C from left, ccr(1) mults by C=-C^dag from right
//     wmat[1][1].ccl(-1).gl(-5).ccr(-1).gr(-5); //has opposite sign to the above 
//   }else{
//     wmat[0][1] = getSite(site,0,from,lattice,splane);
//     wmat[1][1] = getSite(site,1,from,lattice,splane);

//     wmat[1][0] = getSite(site,0,from_conj,lattice,splane);
//     wmat[0][0] = getSite(site,1,from_conj,lattice,splane);

//     wmat[1][0].cconj();
//     wmat[0][0].cconj();

//     wmat[1][0].ccl(-1).gl(-5).ccr(1).gr(-5); 
//     wmat[0][0].ccl(-1).gl(-5).ccr(-1).gr(-5);
//   }

// }
  
void SpinColorFlavorMatrix::generate_from_real_source(QPropWcontainer &from, Lattice &lattice, const int &site, const PropSplane &splane){
  //Use the prop conj relation for a real source to compute the full 2x2 propagator

  int flav = from.flavor();
  if(flav == 0){
    wmat[0][0] = getSite(site,0,from,lattice,splane);
    wmat[1][0] = getSite(site,1,from,lattice,splane);
    wmat[0][1] = wmat[1][0];
    wmat[1][1] = wmat[0][0];

    wmat[0][1].cconj();
    wmat[1][1].cconj();

    wmat[0][1].ccl(-1).gl(-5).ccr(1).gr(-5); //ccl(-1) mults by C from left, ccr(1) mults by C=-C^dag from right
    wmat[1][1].ccl(-1).gl(-5).ccr(-1).gr(-5); //has opposite sign to the above 
  }else{ //flavour 1 source
    wmat[0][1] = getSite(site,0,from,lattice,splane);
    wmat[1][1] = getSite(site,1,from,lattice,splane);
    wmat[1][0] = wmat[0][1];
    wmat[0][0] = wmat[1][1];

    wmat[1][0].cconj();
    wmat[0][0].cconj();

    wmat[1][0].ccl(-1).gl(-5).ccr(1).gr(-5); 
    wmat[0][0].ccl(-1).gl(-5).ccr(-1).gr(-5);
  }
}

void SpinColorFlavorMatrix::generate(QPropWcontainer &from, Lattice &lattice, const int &site, const PropSplane &splane){
  const char* fname = "generate(PropagatorContainer &from, const int &site)";
    
  if(!GJP.Gparity()){
    ERR.General(cname,fname,"Require 2f G-parity BCs to be active");
  }
  //Test to see if we have a propagator with the other flavor
  GparityOtherFlavPropAttrArg* otherfarg; //if a propagator with the same properties but the other flavor exists then use both to generate the matrix
  if(from.getAttr(otherfarg)){
    QPropWcontainer &otherfprop = PropManager::getProp(otherfarg->tag).convert<QPropWcontainer>();
    if(otherfprop.flavor() == from.flavor()) ERR.General(cname,fname,"Found a propagator %s with supposedly the other flavor to propagator %s, but in fact the flavors are identical!",otherfarg->tag,from.tag());
    QPropWcontainer *f0prop; QPropWcontainer *f1prop;
    if(from.flavor() == 0){ f0prop = &from; f1prop = &otherfprop; }
    else { f1prop = &from; f0prop = &otherfprop; }
      
    return generate(*f0prop,*f1prop,lattice,site,splane);
  }

  //Test to see if we have a propagator that has the complex conjugate of the source
  GparityComplexConjSourcePartnerPropAttrArg *cconjsrcarg;
  if(from.getAttr(cconjsrcarg)){
    QPropWcontainer &cconjsrcprop = PropManager::getProp(cconjsrcarg->tag).convert<QPropWcontainer>();
    if(cconjsrcprop.flavor() != from.flavor()) ERR.General(cname,fname,"Found a propagator %s with supposedly the complex conjugate source to propagator %s, but in fact they have different flavors!",cconjsrcarg->tag,from.tag());
    return generate_from_cconj_pair(from, cconjsrcprop, lattice, site, splane);
  }

  //Test for real sources: cos source (point or wall) or zero momentum point source
  PointSourceAttrArg *pt;
  MomentumAttrArg *mom;
  MomCosAttrArg *cos;
  WallSourceAttrArg *wall;
  if( ( from.getAttr(mom) && from.getAttr(cos) ) || ( from.getAttr(pt) && !from.getAttr(mom) ) || ( from.getAttr(wall) && !from.getAttr(mom) ) ){
    return generate_from_real_source(from,lattice,site,splane);
  }

  ERR.General(cname,fname,"Cannot generate prop elements with source %s",from.tag());
}




Complex Trace(const SpinColorFlavorMatrix& a, const SpinColorFlavorMatrix& b){
  Complex trace(0.0,0.0);
  for(int f2=0;f2<2;f2++){
    for(int f1=0;f1<2;f1++){
      for(int s2=0;s2<4;++s2){
	for(int c2=0;c2<3;++c2){
	  for(int s1=0;s1<4;++s1){
	    for(int c1=0;c1<3;++c1){
	      trace += a(s1,c1,f1,s2,c2,f2)*b(s2,c2,f2,s1,c1,f1);
	    }
	  }
	}
      }
    }
  }
  return trace;
}

Matrix SpinFlavorTrace(const SpinColorFlavorMatrix& a, const SpinColorFlavorMatrix& b){
  Matrix trace(0.0);

  for(int c3=0;c3<3;c3++){
    for(int c1=0;c1<3;++c1){
      for(int f1=0;f1<2;f1++){
	for(int s1=0;s1<4;++s1){
	  for(int f2=0;f2<2;f2++){
	    for(int s2=0;s2<4;++s2){
	      for(int c2=0;c2<3;++c2){
		trace(c1,c3) += a(s1,c1,f1,s2,c2,f2) * b(s2,c2,f2,s1,c3,f1);
	      }
	    }
	  }
	}
      }
    }
  }
  return trace;
}

SpinFlavorMatrix ColorTrace(const SpinColorFlavorMatrix& a, const SpinColorFlavorMatrix& b){
  SpinFlavorMatrix trace(0.0);

  for(int s1=0;s1<4;++s1){
    for(int c1=0;c1<3;++c1){
      for(int f1=0;f1<2;f1++){
	for(int f3=0;f3<2;f3++){
	  for(int s3=0;s3<4;++s3){
	    for(int s2=0;s2<4;++s2){
	      for(int c2=0;c2<3;++c2){
		for(int f2=0;f2<2;f2++){
		  trace(s1,f1,s3,f3) += a(s1,c1,f1,s2,c2,f2) * b(s2,c2,f2,s3,c1,f3);
		}
	      }
	    }
	  }
	}
      }
    }
  }
  return trace;
}

FlavorSpinMatrix ColorTraceFS(const SpinColorFlavorMatrix& a, const SpinColorFlavorMatrix& b){
  FlavorSpinMatrix trace(0.0);

  for(int s1=0;s1<4;++s1){
    for(int c1=0;c1<3;++c1){
      for(int f1=0;f1<2;f1++){
	for(int f3=0;f3<2;f3++){
	  for(int s3=0;s3<4;++s3){
	    for(int s2=0;s2<4;++s2){
	      for(int c2=0;c2<3;++c2){
		for(int f2=0;f2<2;f2++){
		  trace(f1,f3)(s1,s3) += a(s1,c1,f1,s2,c2,f2) * b(s2,c2,f2,s3,c1,f3);
		}
	      }
	    }
	  }
	}
      }
    }
  }
  return trace;
}





void FlavorSpinColorMatrix::generate(QPropWcontainer &from, Lattice &lattice, const int &site, const PropSplane &splane){
  const static char* cname = "FlavorSpinColorMatrix";
  const static char* fname = "generate(PropagatorContainer &from, const int &site)";
    
  if(!GJP.Gparity()){
    ERR.General(cname,fname,"Require 2f G-parity BCs to be active");
  }
  PointSourceAttrArg *pt;
  MomentumAttrArg *mom;
  MomCosAttrArg *cos;
  WallSourceAttrArg *wall;
  if( ( from.getAttr(mom) && from.getAttr(cos) ) || ( from.getAttr(pt) && !from.getAttr(mom) ) || ( from.getAttr(wall) && !from.getAttr(mom) ) ){
    //cos source (point or wall) or zero momentum point source
  }else{
    GparityOtherFlavPropAttrArg* otherfarg; //if a propagator with the same properties but the other flavor exists then use both to generate the matrix
    if(from.getAttr(otherfarg)){
      QPropWcontainer &otherfprop = PropManager::getProp(otherfarg->tag).convert<QPropWcontainer>();
      if(otherfprop.flavor() == from.flavor()) ERR.General(cname,fname,"Found a propagator %s with supposedly the other flavor to propagator %s, but in fact the flavors are identical!",otherfarg->tag,from.tag());
      QPropWcontainer *f0prop; QPropWcontainer *f1prop;
      if(from.flavor() == 0){ f0prop = &from; f1prop = &otherfprop; }
      else { f1prop = &from; f0prop = &otherfprop; }
	
      return generate(*f0prop,*f1prop,lattice,site,splane);
    }else ERR.General(cname,fname,"Cannot generate prop elements without a real source, e.g. a Cos (Point or Wall) or Zero-Momentum Point or Wall source, or else two props with the same attributes but different flavors and a GparityOtherFlavPropAttrArg given");
  }
    
  int flav = from.flavor();
  WilsonMatrix wmat[2][2];

  if(flav == 0){
    wmat[0][0] = getSite(site,0,from,lattice,splane);
    wmat[1][0] = getSite(site,1,from,lattice,splane);
    wmat[0][1] = wmat[1][0];
    wmat[1][1] = wmat[0][0];

    wmat[0][1].cconj();
    wmat[1][1].cconj();

    wmat[0][1].ccl(-1).gl(-5).ccr(1).gr(-5); //ccl(-1) mults by C from left, ccr(1) mults by C=-C^dag from right
    wmat[1][1].ccl(-1).gl(-5).ccr(-1).gr(-5); //has opposite sign to the above 
  }else{ //flavour 1 source
    wmat[0][1] = getSite(site,0,from,lattice,splane);
    wmat[1][1] = getSite(site,1,from,lattice,splane);
    wmat[1][0] = wmat[0][1];
    wmat[0][0] = wmat[1][1];

    wmat[1][0].cconj();
    wmat[0][0].cconj();

    wmat[1][0].ccl(-1).gl(-5).ccr(1).gr(-5); 
    wmat[0][0].ccl(-1).gl(-5).ccr(-1).gr(-5);
  }
    
  for(int f1=0;f1<2;++f1) for(int s1=0;s1<4;++s1) for(int c1=0;c1<3;++c1)
						    for(int f2=0;f2<2;++f2) for(int s2=0;s2<4;++s2) for(int c2=0;c2<3;++c2)
												      p.f[f1].d[s1].c[c1].f[f2].d[s2].c[c2] = wmat[f1][f2](s1,c1,s2,c2);
}

//multiply on left by a flavor matrix
FlavorSpinColorMatrix & FlavorSpinColorMatrix::pl(const FlavorMatrixType &type){
  const static char* cname = "FlavorSpinColorMatrix";

  if(type == F0 || type == F1){
    int f1 = (type == F0 ? 1 : 0);
    for(int s1=0;s1<4;++s1) for(int c1=0;c1<3;++c1)
			      for(int f2=0;f2<2;++f2) for(int s2=0;s2<4;++s2) for(int c2=0;c2<3;++c2)
										p.f[f1].d[s1].c[c1].f[f2].d[s2].c[c2] = 0.0;
    return *this;

  }else if(type == Fud){
    for(int s1=0;s1<4;++s1) for(int c1=0;c1<3;++c1)
			      for(int f2=0;f2<2;++f2) for(int s2=0;s2<4;++s2) for(int c2=0;c2<3;++c2){
				    Complex swp = p.f[0].d[s1].c[c1].f[f2].d[s2].c[c2];
				    p.f[0].d[s1].c[c1].f[f2].d[s2].c[c2] = p.f[1].d[s1].c[c1].f[f2].d[s2].c[c2];
				    p.f[1].d[s1].c[c1].f[f2].d[s2].c[c2] = swp;
				  }
    return *this;

  }else if(type == sigma0){
    return *this;
  }else if(type == sigma1){
    return pl(Fud);
  }else if(type == sigma2){
    for(int s1=0;s1<4;++s1) for(int c1=0;c1<3;++c1)
			      for(int f2=0;f2<2;++f2) for(int s2=0;s2<4;++s2) for(int c2=0;c2<3;++c2){
				    Complex swp = Complex(0,1)*p.f[0].d[s1].c[c1].f[f2].d[s2].c[c2];
				    p.f[0].d[s1].c[c1].f[f2].d[s2].c[c2] = Complex(0,-1)*p.f[1].d[s1].c[c1].f[f2].d[s2].c[c2];
				    p.f[1].d[s1].c[c1].f[f2].d[s2].c[c2] = swp;
				  }
    return *this;
  }else if(type == sigma3){
    for(int s1=0;s1<4;++s1) for(int c1=0;c1<3;++c1)
    for(int f2=0;f2<2;++f2) for(int s2=0;s2<4;++s2) for(int c2=0;c2<3;++c2)
	p.f[1].d[s1].c[c1].f[f2].d[s2].c[c2] *= -1;
    return *this;
  }
  ERR.General(cname,"pl(const FlavorMatrixType &type)","Unknown FlavorMatrixType");
}
//multiply on right by a flavor matrix
FlavorSpinColorMatrix & FlavorSpinColorMatrix::pr(const FlavorMatrixType &type){
  const static char* cname = "FlavorSpinColorMatrix";

  if(type == F0 || type == F1){
    int f2 = (type == F0 ? 1 : 0);
    for(int f1=0;f1<2;++f1) for(int s1=0;s1<4;++s1) for(int c1=0;c1<3;++c1)
						      for(int s2=0;s2<4;++s2) for(int c2=0;c2<3;++c2)
										p.f[f1].d[s1].c[c1].f[f2].d[s2].c[c2] = 0.0;
    return *this;

  }else if(type == Fud){
    for(int f1=0;f1<2;++f1) for(int s1=0;s1<4;++s1) for(int c1=0;c1<3;++c1)
						      for(int s2=0;s2<4;++s2) for(int c2=0;c2<3;++c2){
							  Complex swp = p.f[f1].d[s1].c[c1].f[0].d[s2].c[c2];
							  p.f[f1].d[s1].c[c1].f[0].d[s2].c[c2] = p.f[f1].d[s1].c[c1].f[1].d[s2].c[c2];
							  p.f[f1].d[s1].c[c1].f[1].d[s2].c[c2] = swp;
							}
    return *this;
  }else if(type == sigma0){
    return *this;
  }else if(type == sigma1){
    return pr(Fud);
  }else if(type == sigma2){
    for(int f1=0;f1<2;++f1) for(int s1=0;s1<4;++s1) for(int c1=0;c1<3;++c1)
    for(int s2=0;s2<4;++s2) for(int c2=0;c2<3;++c2){
      Complex swp = Complex(0,-1)*p.f[f1].d[s1].c[c1].f[0].d[s2].c[c2];
      p.f[f1].d[s1].c[c1].f[0].d[s2].c[c2] = Complex(0,1)*p.f[f1].d[s1].c[c1].f[1].d[s2].c[c2];
      p.f[f1].d[s1].c[c1].f[1].d[s2].c[c2] = swp;
    }
    return *this;
  }else if(type == sigma3){
    for(int f1=0;f1<2;++f1) for(int s1=0;s1<4;++s1) for(int c1=0;c1<3;++c1)
    for(int s2=0;s2<4;++s2) for(int c2=0;c2<3;++c2)
       p.f[f1].d[s1].c[c1].f[1].d[s2].c[c2] *= -1;
    return *this;
  }
  ERR.General(cname,"pr(const FlavorMatrixType &type)","Unknown FlavorMatrixType");
}







#define TIMESPLUSONE(a,b) { b=a; }
#define TIMESMINUSONE(a,b) { b=-a; }
#define TIMESPLUSI(a,b) { b=Complex(-a.imag(),a.real()); }
#define TIMESMINUSI(a,b) { b=Complex(a.imag(),-a.real()); }

//CK: The below are the same as for WilsonMatrix only with the additional loop over flavour components

FlavorSpinColorMatrix& FlavorSpinColorMatrix::gl(int dir){
  int i; /*color*/
  int c2,s2;    /* column indices, color and spin */
  spin_color_flavor_matrix src=p;
  int f1,f2;

  switch(dir){
    case 0:
        for(f1=0;f1<2;++f1) for(f2=0;f2<2;++f2)
        for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
            TIMESPLUSI(  src.f[f1].d[3].c[i].f[f2].d[s2].c[c2],
                p.f[f1].d[0].c[i].f[f2].d[s2].c[c2] );
            TIMESPLUSI(  src.f[f1].d[2].c[i].f[f2].d[s2].c[c2],
                p.f[f1].d[1].c[i].f[f2].d[s2].c[c2] );
            TIMESMINUSI( src.f[f1].d[1].c[i].f[f2].d[s2].c[c2],
                p.f[f1].d[2].c[i].f[f2].d[s2].c[c2] );
            TIMESMINUSI( src.f[f1].d[0].c[i].f[f2].d[s2].c[c2],
                p.f[f1].d[3].c[i].f[f2].d[s2].c[c2] );
        }
        break;
    case 1:
        for(f1=0;f1<2;++f1) for(f2=0;f2<2;++f2)
        for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
            TIMESMINUSONE( src.f[f1].d[3].c[i].f[f2].d[s2].c[c2],
                p.f[f1].d[0].c[i].f[f2].d[s2].c[c2] );
            TIMESPLUSONE(  src.f[f1].d[2].c[i].f[f2].d[s2].c[c2],
                p.f[f1].d[1].c[i].f[f2].d[s2].c[c2] );
            TIMESPLUSONE(  src.f[f1].d[1].c[i].f[f2].d[s2].c[c2],
                p.f[f1].d[2].c[i].f[f2].d[s2].c[c2] );
            TIMESMINUSONE( src.f[f1].d[0].c[i].f[f2].d[s2].c[c2],
                p.f[f1].d[3].c[i].f[f2].d[s2].c[c2] );
        }
        break;
    case 2:
        for(f1=0;f1<2;++f1) for(f2=0;f2<2;++f2)
        for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
            TIMESPLUSI(  src.f[f1].d[2].c[i].f[f2].d[s2].c[c2],
                p.f[f1].d[0].c[i].f[f2].d[s2].c[c2] );
            TIMESMINUSI( src.f[f1].d[3].c[i].f[f2].d[s2].c[c2],
                p.f[f1].d[1].c[i].f[f2].d[s2].c[c2] );
            TIMESMINUSI( src.f[f1].d[0].c[i].f[f2].d[s2].c[c2],
                p.f[f1].d[2].c[i].f[f2].d[s2].c[c2] );
            TIMESPLUSI(  src.f[f1].d[1].c[i].f[f2].d[s2].c[c2],
                p.f[f1].d[3].c[i].f[f2].d[s2].c[c2] );
        }
	break;
    case 3:
        for(f1=0;f1<2;++f1) for(f2=0;f2<2;++f2)
        for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
            TIMESPLUSONE( src.f[f1].d[2].c[i].f[f2].d[s2].c[c2],
                p.f[f1].d[0].c[i].f[f2].d[s2].c[c2] );
            TIMESPLUSONE( src.f[f1].d[3].c[i].f[f2].d[s2].c[c2],
                p.f[f1].d[1].c[i].f[f2].d[s2].c[c2] );
            TIMESPLUSONE( src.f[f1].d[0].c[i].f[f2].d[s2].c[c2],
                p.f[f1].d[2].c[i].f[f2].d[s2].c[c2] );
            TIMESPLUSONE( src.f[f1].d[1].c[i].f[f2].d[s2].c[c2],
                p.f[f1].d[3].c[i].f[f2].d[s2].c[c2] );
        }
        break;
    case -5:
        for(f1=0;f1<2;++f1) for(f2=0;f2<2;++f2)
        for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
            TIMESPLUSONE(  src.f[f1].d[0].c[i].f[f2].d[s2].c[c2],
                p.f[f1].d[0].c[i].f[f2].d[s2].c[c2] );
            TIMESPLUSONE(  src.f[f1].d[1].c[i].f[f2].d[s2].c[c2],
                p.f[f1].d[1].c[i].f[f2].d[s2].c[c2] );
            TIMESMINUSONE( src.f[f1].d[2].c[i].f[f2].d[s2].c[c2],
                p.f[f1].d[2].c[i].f[f2].d[s2].c[c2] );
            TIMESMINUSONE( src.f[f1].d[3].c[i].f[f2].d[s2].c[c2],
                p.f[f1].d[3].c[i].f[f2].d[s2].c[c2] );
        }
        break;
    default:
	break;
  }
  return *this;
}


FlavorSpinColorMatrix& FlavorSpinColorMatrix::gr(int dir){

  int i; /*color*/
  int c1, s1;    /* row indices, color and spin */
  int f1,f2;
  spin_color_flavor_matrix src=p;

  switch(dir){
    case 0:
        for(f1=0;f1<2;++f1) for(f2=0;f2<2;++f2)
        for(i=0;i<3;i++)for(s1=0;s1<4;s1++)for(c1=0;c1<3;c1++){
            TIMESMINUSI( src.f[f1].d[s1].c[c1].f[f2].d[3].c[i],
                p.f[f1].d[s1].c[c1].f[f2].d[0].c[i] );
            TIMESMINUSI( src.f[f1].d[s1].c[c1].f[f2].d[2].c[i],
                p.f[f1].d[s1].c[c1].f[f2].d[1].c[i] );
            TIMESPLUSI(  src.f[f1].d[s1].c[c1].f[f2].d[1].c[i],
                p.f[f1].d[s1].c[c1].f[f2].d[2].c[i] );
            TIMESPLUSI(  src.f[f1].d[s1].c[c1].f[f2].d[0].c[i],
                p.f[f1].d[s1].c[c1].f[f2].d[3].c[i] );
        }
        break;
    case 1:
        for(f1=0;f1<2;++f1) for(f2=0;f2<2;++f2)
        for(i=0;i<3;i++)for(s1=0;s1<4;s1++)for(c1=0;c1<3;c1++){
            TIMESMINUSONE( src.f[f1].d[s1].c[c1].f[f2].d[3].c[i],
                p.f[f1].d[s1].c[c1].f[f2].d[0].c[i] );
            TIMESPLUSONE(  src.f[f1].d[s1].c[c1].f[f2].d[2].c[i],
                p.f[f1].d[s1].c[c1].f[f2].d[1].c[i] );
            TIMESPLUSONE(  src.f[f1].d[s1].c[c1].f[f2].d[1].c[i],
                p.f[f1].d[s1].c[c1].f[f2].d[2].c[i] );
            TIMESMINUSONE( src.f[f1].d[s1].c[c1].f[f2].d[0].c[i],
                p.f[f1].d[s1].c[c1].f[f2].d[3].c[i] );
        }
        break;
    case 2:
        for(f1=0;f1<2;++f1) for(f2=0;f2<2;++f2)
        for(i=0;i<3;i++)for(s1=0;s1<4;s1++)for(c1=0;c1<3;c1++){
            TIMESMINUSI( src.f[f1].d[s1].c[c1].f[f2].d[2].c[i],
                p.f[f1].d[s1].c[c1].f[f2].d[0].c[i] );
            TIMESPLUSI(  src.f[f1].d[s1].c[c1].f[f2].d[3].c[i],
                p.f[f1].d[s1].c[c1].f[f2].d[1].c[i] );
            TIMESPLUSI(  src.f[f1].d[s1].c[c1].f[f2].d[0].c[i],
                p.f[f1].d[s1].c[c1].f[f2].d[2].c[i] );
            TIMESMINUSI( src.f[f1].d[s1].c[c1].f[f2].d[1].c[i],
                p.f[f1].d[s1].c[c1].f[f2].d[3].c[i] );
        }
        break;
    case 3:
        for(f1=0;f1<2;++f1) for(f2=0;f2<2;++f2)
        for(i=0;i<3;i++)for(s1=0;s1<4;s1++)for(c1=0;c1<3;c1++){
            TIMESPLUSONE( src.f[f1].d[s1].c[c1].f[f2].d[2].c[i],
                p.f[f1].d[s1].c[c1].f[f2].d[0].c[i] );
            TIMESPLUSONE( src.f[f1].d[s1].c[c1].f[f2].d[3].c[i],
                p.f[f1].d[s1].c[c1].f[f2].d[1].c[i] );
            TIMESPLUSONE( src.f[f1].d[s1].c[c1].f[f2].d[0].c[i],
                p.f[f1].d[s1].c[c1].f[f2].d[2].c[i] );
            TIMESPLUSONE( src.f[f1].d[s1].c[c1].f[f2].d[1].c[i],
                p.f[f1].d[s1].c[c1].f[f2].d[3].c[i] );
        }
        break;
    case -5:
        for(f1=0;f1<2;++f1) for(f2=0;f2<2;++f2)
        for(i=0;i<3;i++)for(s1=0;s1<4;s1++)for(c1=0;c1<3;c1++){
            TIMESPLUSONE(  src.f[f1].d[s1].c[c1].f[f2].d[0].c[i],
                p.f[f1].d[s1].c[c1].f[f2].d[0].c[i] );
            TIMESPLUSONE(  src.f[f1].d[s1].c[c1].f[f2].d[1].c[i],
                p.f[f1].d[s1].c[c1].f[f2].d[1].c[i] );
            TIMESMINUSONE( src.f[f1].d[s1].c[c1].f[f2].d[2].c[i],
                p.f[f1].d[s1].c[c1].f[f2].d[2].c[i] );
            TIMESMINUSONE( src.f[f1].d[s1].c[c1].f[f2].d[3].c[i],
                p.f[f1].d[s1].c[c1].f[f2].d[3].c[i] );
        }
        break;
    default:
        //VRB.Result(cname,fname,"BAD CALL TO gl()\n");
        break;
  }
  return *this;
}

//Note, this has the same incorrect ordering as for WilsonMatrix:   ccl(1)  =  C^{-1}M   ccl(-1) = CM
FlavorSpinColorMatrix& FlavorSpinColorMatrix::ccl(int dir){
  int i; /*color*/
  int c2,s2;    /* column indices, color and spin */
  int f1,f2;
  spin_color_flavor_matrix src=p;

  switch(dir){
    case 1:
        for(f1=0;f1<2;++f1) for(f2=0;f2<2;++f2)
        for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
           TIMESPLUSONE(  src.f[f1].d[3].c[i].f[f2].d[s2].c[c2],
                p.f[f1].d[2].c[i].f[f2].d[s2].c[c2] );
	   TIMESMINUSONE( src.f[f1].d[2].c[i].f[f2].d[s2].c[c2],
                p.f[f1].d[3].c[i].f[f2].d[s2].c[c2] );
	   TIMESMINUSONE( src.f[f1].d[1].c[i].f[f2].d[s2].c[c2],
                p.f[f1].d[0].c[i].f[f2].d[s2].c[c2] );
	   TIMESPLUSONE(  src.f[f1].d[0].c[i].f[f2].d[s2].c[c2],
                p.f[f1].d[1].c[i].f[f2].d[s2].c[c2] );        	
        }
        break;
    case -1:
        for(f1=0;f1<2;++f1) for(f2=0;f2<2;++f2)
        for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
           TIMESMINUSONE(  src.f[f1].d[3].c[i].f[f2].d[s2].c[c2],
                p.f[f1].d[2].c[i].f[f2].d[s2].c[c2] );
	   TIMESPLUSONE( src.f[f1].d[2].c[i].f[f2].d[s2].c[c2],
                p.f[f1].d[3].c[i].f[f2].d[s2].c[c2] );
	   TIMESPLUSONE( src.f[f1].d[1].c[i].f[f2].d[s2].c[c2],
                p.f[f1].d[0].c[i].f[f2].d[s2].c[c2] );
	   TIMESMINUSONE(  src.f[f1].d[0].c[i].f[f2].d[s2].c[c2],
                p.f[f1].d[1].c[i].f[f2].d[s2].c[c2] );        	
        }
        break;
    case 5:
        for(f1=0;f1<2;++f1) for(f2=0;f2<2;++f2)
        for(i=0;i<3;i++)for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
            TIMESMINUSONE( src.f[f1].d[3].c[i].f[f2].d[s2].c[c2],
                p.f[f1].d[2].c[i].f[f2].d[s2].c[c2] );
	    TIMESPLUSONE(  src.f[f1].d[2].c[i].f[f2].d[s2].c[c2],
                p.f[f1].d[3].c[i].f[f2].d[s2].c[c2] );
            TIMESMINUSONE( src.f[f1].d[1].c[i].f[f2].d[s2].c[c2],
                p.f[f1].d[0].c[i].f[f2].d[s2].c[c2] );
            TIMESPLUSONE(  src.f[f1].d[0].c[i].f[f2].d[s2].c[c2],
                p.f[f1].d[1].c[i].f[f2].d[s2].c[c2] );
        }
        break;
    default:
	//VRB.Result(cname,fname,"BAD CALL TO ccl()\n");
	break;
  }
  return *this;
}

FlavorSpinColorMatrix& FlavorSpinColorMatrix::ccr(int dir){
  int i; /*color*/
  int c1,s1;    /* column indices, color and spin */
  spin_color_flavor_matrix src=p;
  int f1,f2;

  switch(dir){
    case 1:
        for(f1=0;f1<2;++f1) for(f2=0;f2<2;++f2)
        for(i=0;i<3;i++)for(s1=0;s1<4;s1++)for(c1=0;c1<3;c1++){
            TIMESPLUSONE(  src.f[f1].d[s1].c[c1].f[f2].d[3].c[i],
                p.f[f1].d[s1].c[c1].f[f2].d[2].c[i] );
            TIMESMINUSONE( src.f[f1].d[s1].c[c1].f[f2].d[2].c[i],
                p.f[f1].d[s1].c[c1].f[f2].d[3].c[i] );
            TIMESMINUSONE( src.f[f1].d[s1].c[c1].f[f2].d[1].c[i],
                p.f[f1].d[s1].c[c1].f[f2].d[0].c[i] );
            TIMESPLUSONE(  src.f[f1].d[s1].c[c1].f[f2].d[0].c[i],
                p.f[f1].d[s1].c[c1].f[f2].d[1].c[i] );
        }
        break;
    case -1:
        for(f1=0;f1<2;++f1) for(f2=0;f2<2;++f2)
        for(i=0;i<3;i++)for(s1=0;s1<4;s1++)for(c1=0;c1<3;c1++){
            TIMESMINUSONE(  src.f[f1].d[s1].c[c1].f[f2].d[3].c[i],
                p.f[f1].d[s1].c[c1].f[f2].d[2].c[i] );
            TIMESPLUSONE( src.f[f1].d[s1].c[c1].f[f2].d[2].c[i],
                p.f[f1].d[s1].c[c1].f[f2].d[3].c[i] );
            TIMESPLUSONE( src.f[f1].d[s1].c[c1].f[f2].d[1].c[i],
                p.f[f1].d[s1].c[c1].f[f2].d[0].c[i] );
            TIMESMINUSONE(  src.f[f1].d[s1].c[c1].f[f2].d[0].c[i],
                p.f[f1].d[s1].c[c1].f[f2].d[1].c[i] );
        }
        break;
    case 5:
        for(f1=0;f1<2;++f1) for(f2=0;f2<2;++f2)
        for(i=0;i<3;i++)for(s1=0;s1<4;s1++)for(c1=0;c1<3;c1++){
            TIMESMINUSONE( src.f[f1].d[s1].c[c1].f[f2].d[3].c[i],
                p.f[f1].d[s1].c[c1].f[f2].d[2].c[i] );
            TIMESPLUSONE(  src.f[f1].d[s1].c[c1].f[f2].d[2].c[i],
                p.f[f1].d[s1].c[c1].f[f2].d[3].c[i] );
            TIMESMINUSONE( src.f[f1].d[s1].c[c1].f[f2].d[1].c[i],
                p.f[f1].d[s1].c[c1].f[f2].d[0].c[i] );
            TIMESPLUSONE(  src.f[f1].d[s1].c[c1].f[f2].d[0].c[i],
                p.f[f1].d[s1].c[c1].f[f2].d[1].c[i] );
        }
        break;
    default:
	//VRB.Result(cname,fname,"BAD CALL TO ccr()\n");
	break;
  }
  return *this;
}


//! left multiply by 1/2(1-gamma_5)
FlavorSpinColorMatrix& FlavorSpinColorMatrix::glPL(){
  for(int f1=0;f1<2;++f1) for(int f2=0;f2<2;++f2)
  for(int srow =0;srow<2;srow++)
    for(int crow=0;crow<3;crow++)
      for(int scol=0;scol<4;scol++)
	for(int ccol=0;ccol<3;ccol++)
	  p.f[f1].d[srow].c[crow].f[f2].d[scol].c[ccol] = 0.0;
  return *this;
}


//! left multiply by 1/2(1+gamma_5)
FlavorSpinColorMatrix& FlavorSpinColorMatrix::glPR(){
  for(int f1=0;f1<2;++f1) for(int f2=0;f2<2;++f2)
  for(int srow =2;srow<4;srow++)
    for(int crow=0;crow<3;crow++)
      for(int scol=0;scol<4;scol++)
	for(int ccol=0;ccol<3;ccol++)
	  p.f[f1].d[srow].c[crow].f[f2].d[scol].c[ccol] = 0.0;
  return *this;
}

//! right multiply by 1/2(1-gamma_5)
FlavorSpinColorMatrix& FlavorSpinColorMatrix::grPL(){
  for(int f1=0;f1<2;++f1) for(int f2=0;f2<2;++f2)
  for(int srow =0;srow<4;srow++)
    for(int crow=0;crow<3;crow++)
      for(int scol=0;scol<2;scol++)
	for(int ccol=0;ccol<3;ccol++)
	  p.f[f1].d[srow].c[crow].f[f2].d[scol].c[ccol] = 0.0;
  return *this;
}

//! right multiply by 1/2(1+gamma_5)
FlavorSpinColorMatrix& FlavorSpinColorMatrix::grPR(){
  for(int f1=0;f1<2;++f1) for(int f2=0;f2<2;++f2)
  for(int srow =0;srow<4;srow++)
    for(int crow=0;crow<3;crow++)
      for(int scol=2;scol<4;scol++)
	for(int ccol=0;ccol<3;ccol++)
	  p.f[f1].d[srow].c[crow].f[f2].d[scol].c[ccol] = 0.0;
  return *this;
}












CPS_END_NAMESPACE
