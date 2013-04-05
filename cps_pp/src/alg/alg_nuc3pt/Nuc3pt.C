//------------------------------------------------------------------
//
// Nuc3pt.C
//
// Implementation of the Nucl3pt  class 
//
// April 2002
//
// Kostas Orginos
//
//------------------------------------------------------------------

#include <alg/nuc3pt.h>

CPS_START_NAMESPACE

//#define DEBUG


Nuc3pt::Nuc3pt() : factor(1.0),src_t(-1),snk_t(-1),mom()  
{
  cname = "Nuc3pt";    
  //char *fname = "Nuc3pt()";
  //VRB.Func(cname,fname);
  
  srcdesc[0]=srcdesc[1]=srcdesc[2]=-1 ;
}

Nuc3pt::Nuc3pt(ThreeMom m) : factor(1.0),src_t(-1),snk_t(-1),mom(m)  
{
  cname = "Nuc3pt";    
  //char *fname = "Nuc3pt()";
  //VRB.Func(cname,fname);
  
  srcdesc[0]=srcdesc[1]=srcdesc[2]=-1 ;  
}
  
Nuc3pt::Nuc3pt(Complex cc) : factor(cc),src_t(-1),snk_t(-1),mom()  
{
  cname = "Nuc3pt";    
  //char *fname = "Nuc3pt()";
  //VRB.Func(cname,fname);
  
  srcdesc[0]=srcdesc[1]=srcdesc[2]=-1 ;
}

Nuc3pt::Nuc3pt(ThreeMom m, Complex cc) : factor(cc),src_t(-1),snk_t(-1),mom(m)
{
  cname = "Nuc3pt";    
  //char *fname = "Nuc3pt()";    
  //VRB.Func(cname,fname);
  
  srcdesc[0]=srcdesc[1]=srcdesc[2]=-1 ;
}



/*! IO stuff */
void Nuc3pt::Print(FILE *fp) 
{
  char *fname = "Print(FILE *fp)";

  VRB.Func(cname,fname);

  if(fp==NULL) return ;

  Fprintf(fp,"START_NUC3PT\n") ;
  Fprintf(fp,"MASSES:  %e %e %e\n",
	  (float)quark_mass,(float)quark_mass,(float)quark_mass) ;
  Fprintf(fp,"SOURCE: %s ",source) ;
  for(int i(0) ; i<3;i++)
    if(srcdesc[i]>-1)
      Fprintf(fp,"%i ",srcdesc[i]) ;
  Fprintf(fp,"%i\n",src_t ) ;
  Fprintf(fp,"SINK: %s %i\n",sink,snk_t) ;
  Fprintf(fp,"SNK_MOM: %i %i %i\n",
	  snk_mom.cmp(0),snk_mom.cmp(1),snk_mom.cmp(2)) ;
  //  Fprintf(fp,"OPER: %s\n",oper) ;
  Fprintf(fp,"OPER: "); PrintTheTag(fp) ; Fprintf(fp,"\n");
  Fprintf(fp,"OP_MOM: %i %i %i\n",mom.cmp(0),mom.cmp(1),mom.cmp(2)) ;
  Fprintf(fp,"FACT: %e %e\n",factor.real(),factor.imag()) ;
  Fprintf(fp,"PROJ: %s\n",project) ;
  Fprintf(fp,"QUARKS:    up      down\n") ;
  for(int t=0; t<u_cf.TimeSize();t++)
    {
      Fprintf(fp,"%i  ",t) ;
      u_cf.print(fp,t) ;
      d_cf.print(fp,t) ;
      Fprintf(fp,"\n") ;
    }
  Fprintf(fp,"END_NUC3PT\n");
}

/*!
  Calculates the 3pt function given the sequential propagator and the 
  quark propagator.
 */
void Nuc3pt::Calc3pt(QPropWSeqBar& seqQ, QPropW& Quark)
{
   char *fname = "Calc3pt(QPropWSeqBar&, QPropW&)";

  VRB.Func(cname,fname);

  if(snk_t<0) // Sink time slice not set
    snk_t = seqQ.SourceTime() ;
  if(src_t<0) // Source time slice not set
    src_t = Quark.SourceTime() ;

  if((src_t != Quark.SourceTime()) || (snk_t != seqQ.SourceTime()))
    ERR.General(cname,fname,"OOPS! source/sink times don't match\n") ;

  //I only do degenerate quark masses
  quark_mass = Quark.Mass() ;
  if (quark_mass != seqQ.Mass())
    ERR.General(cname,fname,"OOPS! up and down masses differ!\n") ;
  
  // Check what kind of projection we are doing
  switch (seqQ.Projection())
    {
    case PPAR_5X:
      project="PPAR_5X" ;
      break ;
    case PPAR_5Y:
      project="PPAR_5Y" ;
      break ;
    case PPAR_5Z:
      project="PPAR_5Z" ;
      break ;
    case PPAR:
      project="PPAR" ;
      break ;
    case PPAR_XY:
      project="PPAR_XY" ;
      break ;
    default:
      ERR.General(cname,fname,"OOPS! Unknown projection\n") ;
    }

  char *msg = "OOPS! source changed!\n" ;
  //Setup the source Tag
  switch (Quark.SrcType())
    {
    case POINT:
      if(srcdesc[0]<0) srcdesc[0] = Quark.PointSrcX() ; 
      if(srcdesc[1]<0) srcdesc[1] = Quark.PointSrcY() ;
      if(srcdesc[2]<0) srcdesc[2] = Quark.PointSrcZ() ;
      if((srcdesc[0]!=Quark.PointSrcX())||
	 (srcdesc[1]!=Quark.PointSrcY())||
	 (srcdesc[2]!=Quark.PointSrcY()))
	ERR.General(cname,fname,msg) ;
      source="POINT" ;
      break ;
    case WALL:
      source="WALL" ;
      break ;
    case BOX:
      if(srcdesc[0]<0) srcdesc[0] = Quark.BoxSrcStart() ; 
      if(srcdesc[1]<0) srcdesc[1] = Quark.BoxSrcEnd() ;
      if((srcdesc[0]!=Quark.BoxSrcStart())||(srcdesc[1]!=Quark.BoxSrcEnd()))
	ERR.General(cname,fname,msg) ;
      source="BOX" ;
      break ;
    case GAUSS_GAUGE_INV:
      if(srcdesc[0]<0) srcdesc[0] = Quark.Gauss_N() ; 
      if(srcdesc[1]<0) srcdesc[1] = int(100*Quark.Gauss_W()) ;
      if((srcdesc[0]!=Quark.Gauss_N())||(srcdesc[1]!=int(100*Quark.Gauss_W())))
	ERR.General(cname,fname,msg) ;
      source="GAUSS" ;
      break ;
    default:
      ERR.General(cname,fname,"OOPS! Unknown source.\n") ;
    }
  
  //Setup the sink Tag
  if((seqQ.SrcType()==PROT_U_SEQ)||(seqQ.SrcType()==PROT_D_SEQ))
    {
      if(seqQ.SeqSmearSource()==GAUSS_GAUGE_INV)
	sink="GAUSS" ;
      else
        sink="POINT" ;
    }

  // Copy the sink momentum (just for IO reasons)
  snk_mom = seqQ.Mom() ; 

  CorrFunc tmp ; // tmp is zeroed at construction

  InsertOp(tmp,seqQ,Quark) ;

  //Global sum
  tmp.GlobalSum();
  
  //Multiply by the factor...
  tmp *=  factor ;

  switch (seqQ.SrcType())
    {
    case PROT_U_SEQ:
      u_cf = tmp ;
      break ;
    case PROT_D_SEQ:
      d_cf = tmp ;
      break ;
    default:
      ERR.General(cname,fname,"Unknown quark...\n") ;
    }
}


/*! IO stuff */
void Nuc3pt::Print(FILE *fp, int Nmom, ThreeMom* momo) 
{
  char *fname = "Print(FILE *fp)";

  VRB.Func(cname,fname);

  if(fp==NULL) return ;

  for(int ip=0;ip<Nmom;ip++) {
  Fprintf(fp,"START_NUC3PT\n") ;
  Fprintf(fp,"MASSES:  %e %e %e\n",
	  (float)quark_mass,(float)quark_mass,(float)quark_mass) ;
  Fprintf(fp,"SOURCE: %s ",source) ;
  for(int i(0) ; i<3;i++)
    if(srcdesc[i]>-1)
      Fprintf(fp,"%i ",srcdesc[i]) ;
  Fprintf(fp,"%i\n",src_t ) ;
  Fprintf(fp,"SINK: %s %i\n",sink,snk_t) ;
  Fprintf(fp,"SNK_MOM: %i %i %i\n",
	  snk_mom.cmp(0),snk_mom.cmp(1),snk_mom.cmp(2)) ;
  //  Fprintf(fp,"OPER: %s\n",oper) ;
  Fprintf(fp,"OPER: "); PrintTheTag(fp) ; Fprintf(fp,"\n");
  Fprintf(fp,"OP_MOM: %i %i %i\n",momo[ip].cmp(0),momo[ip].cmp(1),momo[ip].cmp(2)) ;
  Fprintf(fp,"FACT: %e %e\n",factor.real(),factor.imag()) ;
  Fprintf(fp,"PROJ: %s\n",project) ;
  Fprintf(fp,"QUARKS:    up      down\n") ;
  for(int t=0; t<u_cf.TimeSize();t++)
    {
      Fprintf(fp,"%i  ",t) ;
      u_cfm[ip]->print(fp,t) ;
      d_cfm[ip]->print(fp,t) ;
      Fprintf(fp,"\n") ;
    }
  Fprintf(fp,"END_NUC3PT\n");

  delete u_cfm[ip];
  delete d_cfm[ip];
  }
}

/*!
  Calculates the 3pt function given the sequential propagator and the 
  quark propagator.
 */
void Nuc3pt::Calc3pt(QPropWSeqBar& seqQ, QPropW& Quark, int Nmom, ThreeMom* momo)
{
   char *fname = "Calc3pt(QPropWSeqBar&, QPropW&)";

  VRB.Func(cname,fname);

  if(snk_t<0) // Sink time slice not set
    snk_t = seqQ.SourceTime() ;
  if(src_t<0) // Source time slice not set
    src_t = Quark.SourceTime() ;

  if((src_t != Quark.SourceTime()) || (snk_t != seqQ.SourceTime()))
    ERR.General(cname,fname,"OOPS! source/sink times don't match\n") ;

  //I only do degenerate quark masses
  quark_mass = Quark.Mass() ;
  if (quark_mass != seqQ.Mass())
    ERR.General(cname,fname,"OOPS! up and down masses differ!\n") ;
  
  // Check what kind of projection we are doing
  switch (seqQ.Projection())
    {
    case PPAR_5X:
      project="PPAR_5X" ;
      break ;
    case PPAR_5Y:
      project="PPAR_5Y" ;
      break ;
    case PPAR_5Z:
      project="PPAR_5Z" ;
      break ;
    case PPAR:
      project="PPAR" ;
      break ;
    case PPAR_XY:
      project="PPAR_XY" ;
      break ;
    default:
      ERR.General(cname,fname,"OOPS! Unknown projection\n") ;
    }

  char *msg = "OOPS! source changed!\n" ;
  //Setup the source Tag
  switch (Quark.SrcType())
    {
    case POINT:
      if(srcdesc[0]<0) srcdesc[0] = Quark.PointSrcX() ; 
      if(srcdesc[1]<0) srcdesc[1] = Quark.PointSrcY() ;
      if(srcdesc[2]<0) srcdesc[2] = Quark.PointSrcZ() ;
      if((srcdesc[0]!=Quark.PointSrcX())||
	 (srcdesc[1]!=Quark.PointSrcY())||
	 (srcdesc[2]!=Quark.PointSrcY()))
	ERR.General(cname,fname,msg) ;
      source="POINT" ;
      break ;
    case WALL:
      source="WALL" ;
      break ;
    case BOX:
      if(srcdesc[0]<0) srcdesc[0] = Quark.BoxSrcStart() ; 
      if(srcdesc[1]<0) srcdesc[1] = Quark.BoxSrcEnd() ;
      if((srcdesc[0]!=Quark.BoxSrcStart())||(srcdesc[1]!=Quark.BoxSrcEnd()))
	ERR.General(cname,fname,msg) ;
      source="BOX" ;
      break ;
    case GAUSS_GAUGE_INV:
      if(srcdesc[0]<0) srcdesc[0] = Quark.Gauss_N() ; 
      if(srcdesc[1]<0) srcdesc[1] = int(100*Quark.Gauss_W()) ;
      if((srcdesc[0]!=Quark.Gauss_N())||(srcdesc[1]!=int(100*Quark.Gauss_W())))
	ERR.General(cname,fname,msg) ;
      source="GAUSS" ;
      break ;
    default:
      ERR.General(cname,fname,"OOPS! Unknown source.\n") ;
    }
  
  //Setup the sink Tag
  if((seqQ.SrcType()==PROT_U_SEQ)||(seqQ.SrcType()==PROT_D_SEQ))
    {
      if(seqQ.SeqSmearSource()==GAUSS_GAUGE_INV)
	sink="GAUSS" ;
      else
        sink="POINT" ;
    }

  // Copy the sink momentum (just for IO reasons)
  snk_mom = seqQ.Mom() ; 

//  CorrFunc tmp[Nmom] ; // tmp is zeroed at construction
  CorrFunc *tmp = new CorrFunc[Nmom];
  if(tmp == 0) ERR.Pointer(cname, fname, "tmp");

  InsertOp(tmp,seqQ,Quark,Nmom,momo) ;

  //Global sum
  for(int ip=0;ip<Nmom;ip++) tmp[ip].GlobalSum();
  
  //Multiply by the factor...
  for(int ip=0;ip<Nmom;ip++) tmp[ip] *=  factor ;

  for(int ip=0;ip<Nmom;ip++) {
  switch (seqQ.SrcType())
    {
    case PROT_U_SEQ:
      *u_cfm[ip] = tmp[ip] ;
      break ;
    case PROT_D_SEQ:
      *d_cfm[ip] = tmp[ip] ;
      break ;
    default:
      ERR.General(cname,fname,"Unknown quark...\n") ;
    }
  }
  delete[] tmp;
}


Nuc3ptGamma::Nuc3ptGamma(Gamma op):Nuc3pt(),G(op)
{
  cname = "Nuc3ptGamma";    
  //char *fname = "Nuc3ptGamma(Gamma)";
  //VRB.Func(cname,fname);
}

Nuc3ptGamma::Nuc3ptGamma(ThreeMom m, Gamma op):Nuc3pt(m),G(op)
{
  cname = "Nuc3ptGamma";    
  //char *fname = "Nuc3ptGamma(ThreeMom,Gamma)";
  //VRB.Func(cname,fname);
}

Nuc3ptGamma::Nuc3ptGamma(Complex cc, Gamma op):Nuc3pt(cc),G(op)
{
  cname = "Nuc3ptGamma";    
  //char *fname = "Nuc3ptGamma(Complex,Gamma)";
  //VRB.Func(cname,fname);
}

Nuc3ptGamma::Nuc3ptGamma(ThreeMom m, Complex cc, Gamma op):Nuc3pt(m,cc),G(op)
{
  cname = "Nuc3ptGamma";    
  //char *fname = "Nuc3ptGamma(ThreeMom,Complex,Gamma)";
  //VRB.Func(cname,fname);
}
  

/*! 
  Computes
  \f[
   Trace(seqQ \, 
     \gamma_{\mu_1} \, \gamma_{\mu_2} \, \gamma_{\mu_3}\, 
     e^{-i mom . x}\, Quark)
  \f]
  where \f$\mu\f$ is the member int gamma and ThreeMom mom is
  the momentum injected at the operator. 
**/
void Nuc3ptGamma::InsertOp(CorrFunc& tmp,QPropW& seqQ, QPropW& Quark)
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
      Complex cc(mom.Fact(s)) ; // The momentum factor exp(-ipx)
      cc *= Trace(sq,Quark[s.Index()]) ;
      tmp[t] += cc ;
    }
}

void Nuc3ptGamma::InsertOp(CorrFunc* tmp,QPropW& seqQ, QPropW& Quark, int Nmom, ThreeMom* momo)
{
  char *fname = "InsertOp()";
}


Nuc3ptStru::Nuc3ptStru(Gamma gg, Derivative dd):Nuc3pt(),G(gg),D(dd)
{
  cname = "Nuc3ptStru";    
  //char *fname = "Nuc3ptStru(G,D)";
  //VRB.Func(cname,fname);
  tag[0]='\0' ;
}

Nuc3ptStru::Nuc3ptStru(Complex cc, Gamma gg, Derivative dd):Nuc3pt(cc),
							    G(gg),D(dd)
{
  cname = "Nuc3ptStru";    
  //char *fname = "Nuc3ptStru(C,G,D)";
  //VRB.Func(cname,fname); 
  tag[0]='\0' ;
}

Nuc3ptStru::Nuc3ptStru(ThreeMom m, Gamma gg, Derivative dd):Nuc3pt(m),
							    G(gg),D(dd)
{
  cname = "Nuc3ptStru";    
  //char *fname = "Nuc3ptStru(M,G,D)";
  //VRB.Func(cname,fname);
  tag[0]='\0' ;
}

Nuc3ptStru::Nuc3ptStru(ThreeMom m, Complex cc, Gamma gg, Derivative dd):
  Nuc3pt(m,cc),G(gg),D(dd)
{
  cname = "Nuc3ptStru";    
  //char *fname = "Nuc3ptStru(M,C,G,D)";
  //VRB.Func(cname,fname);
  tag[0]='\0' ;
}

/*! 
  Computes
  \f[
   \mbox{Trace}(\mbox{seqQ} \, 
     \gamma_{\mu_1} \, \gamma_{\mu_2} \, \gamma_{\mu_3}\, 
     \stackrel{\displaystyle \leftrightarrow}{D}_\nu 
     \stackrel{\displaystyle \leftrightarrow}{D}_\rho ...
     e^{-i mom . x}\, \mbox{Quark})
  \f]
  where \f$\mu\f$ is the member int gamma and ThreeMom mom is
  the momentum injected at the operator.

  The string of derivatives gets broken up in strings of forward
  and backward hops. For each term the quark and anti-quark are sitting
  on different sites and they are connected with the gauge connection U.
  The gauge connection U is calculated using the PathOrdPlus routine
  given the location of the anti-quark (starting point) and a list
  of hops. The Derivative class is used not only  to do the break-up
  of the string of derivatives in strings of hops, but also to 
  compute the locations of the quark and the anti-quark and 
  provide the list of links needed to compute the gauge connection U
  for each term.

  In order to speed up the code first we compute for all sites
  all derivative terms with the derivatives applied on the quark i.e
  the anti-quark is not displaced. 
  These terms are not summed up but they are stored individually in the
  DTerms Terms class. 
  Then I we go again through the all sites and construct the derivatives.
  If for a given term the anti-quark is displaced we just have to communicate
  a single Complex number from that site. The DTerms class handles all that.
  The routine Derivative :: DTermIndx() is the mapping of the derivative
  terms to a particular storage order. The DTerm class does not have
  to know the storage order. Only the Derivative class knows about it.
 
 NOTE: The quark is at the sink of the Quark propagator
        while the anti-quark is at to source of the sequential propagator.
	
**/
void Nuc3ptStru::InsertOp(CorrFunc& tmp,QPropW& seqQ, QPropW& Quark)
{
  char *fname = "InsertOp()";

  VRB.Func(cname, fname);
  
  if(seqQ.GFixedSnk() || Quark.GFixedSnk())
    ERR.General(cname,fname,"Gauge Fixed sinks!\n") ;

  Lattice& lat=seqQ.AlgLattice() ; // The lattice

#ifdef DEBUG  
printf("%s:: NDTerms= %i\n",cname,D.NDTerms());
#endif

  DTerms Terms(D.NDTerms()); // The derivative terms
  Matrix U ;
  WilsonMatrix aQ ; // the anti-quark
  WilsonMatrix  Q ; // quark
  WilsonMatrix aQU ; // anti-quark times Gauge path
  WilsonMatrix buff ; // Buffer needed for communication 
  int aq_vec[4]; // Anti-quark location
  int q_vec[4] ; // quark location
  int *dirs = new int[D.N()];
  if(dirs == 0) ERR.Pointer(cname, fname, "dirs");

  Site s ;
  //first compute all derivatives terms on quarks
  for(s.Begin();s.End();s.nextSite()) // Loop over sites
    {
      aQ=seqQ[s.Index()] ;
      // Apply the gamma matrices on the sequential  propagator
      // We multiply the source indices of the sequential propagator
      // which are the anti-quark
      for(int mu(0);mu<G.N();mu++)
	aQ.gr(G[mu]) ; // Multiply by gamma_mu 

      // Loop over all derivative terms on the Quark
      for(D.Start();D.NotEnd();D.NextQuark()) 
	{
	  // Calculate the endpoints of the shift
	  D.CalcEndPoints(aq_vec,q_vec,s);
	  // The sink indices of the quark propagator are the quark
	  Q=Quark.GetMatrix(q_vec,buff) ;
	  
	 
	  D.CalcListOfLinks(dirs);

#ifdef DEBUG
printf("\n%s:: dirs ",cname);
for(int i(0);i<D.N();i++) printf("%i ",dirs[i]);
printf(" aq_vec: ");
for(int i(0);i<4;i++) printf("%i ",aq_vec[i]);
printf(" q_vec: ") ; 
for(int i(0);i<4;i++) printf("%i ",q_vec[i]);
#endif

	  // Start from the anti-quark and move along the path described
	  // by the list of dirs multiplying all the encountered links
	  // at the end of the day you should be at the position of the quark
          if(D.N()>1)  
	    {
	      U.ZeroMatrix();
	      lat.PathOrdProdPlus(U,aq_vec,dirs,D.N()) ;
	      //lat.PathOrdProd(U,aq_vec,dirs,D.N()) ;
	    }
	  else //if D.N()==1 PathOrdPlus does not work
	    {
	      int abs_dir = dirs[0]&3; //abs_dir is {0,1,2,3}
	      if(abs_dir == dirs[0]) // dir is forward
		U = *(lat.GetBufferedLink(aq_vec, abs_dir));
	      else // dir is backward. q_vec is now the site where the link is
		U.Dagger(*(lat.GetBufferedLink(q_vec, abs_dir)));
	    }

#ifdef DEBUG
Complex cboo(U.Tr()) ; 
printf("  TrU= (%g %g)\n",cboo.real(),cboo.imag());       
#endif

	  //Now that the gauge connection is constructed multiply it on
	  //the quark i.e the source indices of the anti-quark propagator
	  aQU.UMultSource(U,aQ) ;

	  Terms(D.DTermIndx(),s.Index()) = Trace(aQU,Q) ;

#ifdef DEBUG
printf("   (anti-Q, Q) aQU = (%6.3f %6.3f), (%6.3f %6.3f) (%6.3f %6.3f)\n",
       aQ.Trace().real(),aQ.Trace().imag(),
       Q.Trace().real(),Q.Trace().imag(),
       aQU.Trace().real(),aQU.Trace().imag());
printf("%s:: DTermIndx(): %i site: %i Value: ( %6.3f %6.3f ) \n",
       cname,D.DTermIndx(),s.Index(),
       Terms(D.DTermIndx(),s.Index()).real(), 
       Terms(D.DTermIndx(),s.Index()).imag());
#endif


	}
    }

  Complex tt ;
  // Loop over sites and construct the operator
  for(s.Begin();s.End();s.nextSite()) 
    {
      int t(s.physT()) ;
      for(D.Start();D.NotEnd();D.Next()) // Loop over all derivative terms
	{
	  //Calculate the endpoints of the shift
	  //aq_vec : anti quark location
	  // q_vec :      quark location
	  D.CalcEndPoints(aq_vec,q_vec,s);
	  Complex cc(Terms.GetTerm(D.DTermIndx(),aq_vec,tt)) ;

	  //Finaly the momentum factors
	  cc *= mom.Fact(s) ; // The momentum factor exp(-ipx)	
	  cc *= D.Fact()    ; /* The sign and the normalization factor 
				 for the current derivative term         */

	  tmp[t] += cc ;

#ifdef DEBUG
printf("%s:: DTermIndx(): %i site: %i Fact: %9.6f Value: ( %6.3f %6.3f ) \n",
       cname, D.DTermIndx(),s.Index(),D.Fact(),cc.real(),cc.imag());
#endif
	
	}
    }
  delete [] dirs ;
}

void Nuc3ptStru::InsertOp(CorrFunc* tmp,QPropW& seqQ, QPropW& Quark, int Nmom, ThreeMom* momo)
{
  char *fname = "InsertOp()";
}


/*!
  It is used to add up symmetry variations of a given opperator.
  I have no good way of updating the factor and tag for now
  so the factor and tag of the left had side are kept.
  The programmer at the alg_nuc3pt level has to be carefull of
  what gets added up and what the printed tag and factor mean.
 */
Nuc3pt& Nuc3pt::operator+=(Nuc3pt& rhs)
{
  char *fname = "+=(Nuc3pt&)";
  char *error = "operators do not match\n" ;

  VRB.Func(cname,fname);

  if((quark_mass != rhs.quark_mass)||(src_t!=rhs.src_t)||(snk_t!=rhs.snk_t))
    ERR.General(cname,fname,error);

  if((strcmp(source, rhs.source ))||
     (strcmp(sink  , rhs.sink   ))||
     (strcmp(project,rhs.project)))
    ERR.General(cname,fname,error);
  
  for(int k(0);k<3;k++)
    if((mom.cmp(k)!=rhs.mom.cmp(k))||(snk_mom.cmp(k)!=rhs.snk_mom.cmp(k)))
      ERR.General(cname,fname,error);

  u_cf+=rhs.u_cf ;
  d_cf+=rhs.d_cf ;

  return *this ;
}

/*!
  It is used to add up symmetry variations of a given opperator.
  I have no good way of updating the factor and tag for now
  so the factor and tag of the left had side are kept.
  The programmer at the alg_nuc3pt level has to be carefull of
  what gets subtructed and what the printed tag and factor mean.
 */
Nuc3pt& Nuc3pt::operator-=(Nuc3pt& rhs)
{
  char *fname = "+=(Nuc3pt&)";
  char *error = "operators do not match\n" ;

  VRB.Func(cname,fname);

  if((quark_mass != rhs.quark_mass)||(src_t!=rhs.src_t)||(snk_t!=rhs.snk_t))
    ERR.General(cname,fname,error);

  if((strcmp(source, rhs.source ))||
     (strcmp(sink  , rhs.sink   ))||
     (strcmp(project,rhs.project)))
    ERR.General(cname,fname,error);
  
  for(int k(0);k<3;k++)
    if((mom.cmp(k)!=rhs.mom.cmp(k))||(snk_mom.cmp(k)!=rhs.snk_mom.cmp(k)))
      ERR.General(cname,fname,error);

  u_cf-=rhs.u_cf ;
  d_cf-=rhs.d_cf ;

  return *this ;
}


void Nuc3ptStru::PrintTheTag(FILE *fp)
{
  G.printTag(fp);
  Fprintf(fp,"_") ;
  D.printTag(fp) ;
  if(strlen(tag)>0) Fprintf(fp,"_%s",tag) ;
} 



Nuc3ptCons::Nuc3ptCons(Gamma op):Nuc3pt(),G(op)
{
  cname = "Nuc3ptCons";    
  //char *fname = "Nuc3ptCons(Gamma)";
  //VRB.Func(cname,fname);
}

Nuc3ptCons::Nuc3ptCons(ThreeMom m, Gamma op):Nuc3pt(m),G(op)
{
  cname = "Nuc3ptCons";    
  //char *fname = "Nuc3ptCons(ThreeMom,Gamma)";
  //VRB.Func(cname,fname);
}

Nuc3ptCons::Nuc3ptCons(Complex cc, Gamma op):Nuc3pt(cc),G(op)
{
  cname = "Nuc3ptCons";    
  //char *fname = "Nuc3ptCons(Complex,Gamma)";
  //VRB.Func(cname,fname);
}

Nuc3ptCons::Nuc3ptCons(ThreeMom m, Complex cc, Gamma op):Nuc3pt(m,cc),G(op)
{
  cname = "Nuc3ptCons";    
  //char *fname = "Nuc3ptCons(ThreeMom,Complex,Gamma)";
  //VRB.Func(cname,fname);
}
  
Nuc3ptCons::Nuc3ptCons(int Nmom, Gamma op):Nuc3pt(),G(op)
{
  cname = "Nuc3ptCons";    
  //char *fname = "Nuc3ptCons(ThreeMom,Gamma)";
  //VRB.Func(cname,fname);

  for(int ip=0;ip<Nmom;ip++) {
    u_cfm[ip] = new CorrFunc();
    d_cfm[ip] = new CorrFunc();
  }
}

Nuc3ptCons::Nuc3ptCons(int Nmom, Complex cc, Gamma op):Nuc3pt(cc),G(op)
{
  cname = "Nuc3ptCons";    
  //char *fname = "Nuc3ptCons(ThreeMom,Complex,Gamma)";
  //VRB.Func(cname,fname);

  for(int ip=0;ip<Nmom;ip++) {
    u_cfm[ip] = new CorrFunc();
    d_cfm[ip] = new CorrFunc();
  }
}
  

//int siteOffset(const int lcl[], const int lcl_sites[]);
  //-----------------------------------------------------------------
  // TY Add Start
int Nuc3ptCons::siteOffset(const int lcl[], const int lcl_sites[]) const
{
  int l = 4 - 1;
  int offset = lcl[l];
  while (l-- > 0) {
      offset *= lcl_sites[l];
      offset += lcl[l];
  }
  return offset;
}
  // TY Add End
  //-----------------------------------------------------------------

/*! 
  Computes
  \f[
   Trace(seqQ \, 
     \gamma_{\mu_1} \, \gamma_{\mu_2} \, \gamma_{\mu_3}\, 
     e^{-i mom . x}\, Quark)
  \f]
  where \f$\mu\f$ is the member int gamma and ThreeMom mom is
  the momentum injected at the operator. 
**/
void Nuc3ptCons::InsertOp(CorrFunc& tmp,QPropW& seqQ, QPropW& Quark)
{
  char *fname = "InsertOp()";

  VRB.Func(cname, fname);

  WilsonMatrix v0, v1, v0_next, v1_next;
  WilsonMatrix temp_v0D, temp_v0D_next;
  WilsonMatrix temp_v0_1pD, temp_v0_1mD_next;
  WilsonMatrix temp_v1_u;
  WilsonMatrix temp_v1_next_u;


  int shift_t=GJP.TnodeCoor()*GJP.TnodeSites();
  int vol= (GJP.VolNodeSites()/GJP.TnodeSites());

  int dir = G[G.N()-1];


  const int LORENTZs(4); 
  int lcl_sites[LORENTZs]; 
  lcl_sites[0] = GJP.XnodeSites();
  lcl_sites[1] = GJP.YnodeSites();
  lcl_sites[2] = GJP.ZnodeSites();
  lcl_sites[3] = GJP.TnodeSites();

  int lcl_node[LORENTZs];
  lcl_node[0] = GJP.XnodeCoor();
  lcl_node[1] = GJP.YnodeCoor();
  lcl_node[2] = GJP.ZnodeCoor();
  lcl_node[3] = GJP.TnodeCoor();

  int lcl2glb_offset[LORENTZs];
  for (int i = 0; i < LORENTZs; ++i) 
    lcl2glb_offset[i] = lcl_sites[i] * lcl_node[i];

 
  Lattice& lat = seqQ.AlgLattice();
  Matrix *gauge_field = lat.GaugeField();  

  int ls_glb = GJP.SnodeSites() * GJP.Snodes();
  int ls_ini = 0;
  int ls_end = ls_glb;
  if(GJP.Snodes()==2) {
    ls_ini = GJP.SnodeCoor() * GJP.SnodeSites();
    ls_end = ( GJP.SnodeCoor() + 1 ) * GJP.SnodeSites();
  }

  char dummy[50];
  Site site;

  for(int s=ls_ini;s<ls_end;s++){
    //for(int s=0;s<ls_glb;s++){ 
    Quark.RestoreQPropLs(dummy, s);
    if(GJP.Snodes()==2) {
    seqQ.RestoreQPropLs(dummy, ls_end-1-s+ls_ini);
    } else {
    seqQ.RestoreQPropLs(dummy, ls_end-1-s);
    }

    // This operations are included in QPropWSeqProtU(D)Src.
    // When we restore 5d prop the operations are nessesary.
    //Multiply by gamma5 and take the dagger to make it in to quark.
    for (site.Begin();site.End();site.nextSite()) {
      seqQ[site.Index()].gl(-5);
      seqQ[site.Index()].hconj();
    }


    int lcl[LORENTZs];
    int lcl_next[LORENTZs]; // Next site along propagation direction


    VRB.Clock(cname,fname,"RO\n");
    for(lcl[3]=0; lcl[3]<=GJP.TnodeSites()-1; lcl[3]++) {
    int t=lcl[3]+shift_t;

    for(lcl[0]=0; lcl[0]<=GJP.XnodeSites()-1; lcl[0]++) 
    for(lcl[1]=0; lcl[1]<=GJP.YnodeSites()-1; lcl[1]++) 
    for(lcl[2]=0; lcl[2]<=GJP.ZnodeSites()-1; lcl[2]++) {

      Float coeff = 1.0;
      if( G.N() == 2 ) {
      coeff = -1.0;
      if(GJP.Snodes()==2) {
	if(GJP.SnodeCoor()==1) coeff = -coeff;
      } else {
	if( s > ls_glb/2 - 1 ) coeff = -coeff;
      }
      }

      int lcl_offset = siteOffset(lcl,lcl_sites);

      // coordinates and offset for lcl_next
      for (int i  = 0; i < LORENTZs; i++ ) 
	lcl_next[i] = ( (i == dir) ? (lcl[i]+1)%lcl_sites[i]
			  : lcl[i] );

      int lcl_next_offset = siteOffset(lcl_next,lcl_sites);

      // U_mu(x) where mu = dir
      Matrix * link = gauge_field + siteOffset(lcl,lcl_sites) * 4 + dir ;



      // S_F(x, s)
      v0=Quark[lcl_offset];
      // S_F^+(x, ls_glb-1-s)
      v1=seqQ[lcl_offset];


      // v0_next = S_F(x+dir, s)
      // v1_next = S_F^+(x+dir, ls_glb-1-s)
      if ((lcl[dir]+1) == lcl_sites[dir]) {
	getPlusData( (IFloat *)&v0_next,
		     (IFloat *)(&Quark[lcl_next_offset]), sizeof(WilsonMatrix)/sizeof(Float),
		     dir) ;
	getPlusData( (IFloat *)&v1_next,
		     (IFloat *)(&seqQ[lcl_next_offset]), sizeof(WilsonMatrix)/sizeof(Float), 
		     dir) ;
	
	// fix boundary condition
	switch( dir ) {
	case 0:
	  if (GJP.XnodeBc()==BND_CND_APRD) coeff = -coeff ;
	  break;
	case 1:
	  if (GJP.YnodeBc()==BND_CND_APRD) coeff = -coeff ;
	  break;
	case 2:
	  if (GJP.ZnodeBc()==BND_CND_APRD) coeff = -coeff ;
	  break;
	case 3:
	  if (GJP.TnodeBc()==BND_CND_APRD) coeff = -coeff ;
	  break;
	} // end switch

      } else {
	v0_next = Quark[lcl_next_offset] ;
	v1_next = seqQ[lcl_next_offset] ;
      }

      // Gamma^{dir} S_F(x, s)
      temp_v0D = v0; temp_v0D.gl(dir);
      // Gamma^{dir} S_F(x+dir, s)
      temp_v0D_next = v0_next; temp_v0D_next.gl(dir);

      // ( 1 + Gamma^{dir} ) S_F(x, s)
      temp_v0_1pD = v0 + temp_v0D;
      // ( 1 - Gamma^{dir} ) S_F(x+dir, s)
      temp_v0_1mD_next = v0_next - temp_v0D_next;

      // S_F^+(x, ls_glb-1-s) U
      temp_v1_u.UMultSource(*link, v1);

      // S_F^+(x+dir, ls_glb-1-s) U^+
      temp_v1_next_u.UdagMultSource(*link, v1_next);

      Complex cc(0);
      //-Tr( ( 1 - Gamma^{dir} ) S_F(x+dir, s) S_F^+(x, ls_glb-1-s) U )
      cc-=coeff*( Trace( temp_v0_1mD_next, temp_v1_u ) );
      //+Tr( ( 1 + Gamma^{dir} ) S_F(x, s) S_F^+(x+dir, ls_glb-1-s) U^+ )
      cc+=coeff*( Trace( temp_v0_1pD, temp_v1_next_u ) );
      Site st(lcl[0],lcl[1],lcl[2],lcl[3]);
      cc *=  mom.Fact(st); // The momentum factor exp(-ipx)
      tmp[t] += cc ;
    } // xyz loop

    } // t loop

  } // ls loop

  for(int tt=0; tt<=GJP.TnodeSites()-1; tt++) {
    int t=tt+shift_t;
    if(GJP.Snodes()==2) {
      Float rsum, isum;
      rsum = tmp[t].real();
      glb_sum_dir(&rsum, 4);
      isum = tmp[t].imag();
      glb_sum_dir(&isum, 4);
      tmp[t] = Complex( rsum, isum );
    }
    tmp[t]=tmp[t]/2;
  }
}


/*! 
  Computes
  \f[
   Trace(seqQ \, 
     \gamma_{\mu_1} \, \gamma_{\mu_2} \, \gamma_{\mu_3}\, 
     e^{-i mom . x}\, Quark)
  \f]
  where \f$\mu\f$ is the member int gamma and ThreeMom mom is
  the momentum injected at the operator. 
**/
void Nuc3ptCons::InsertOp(CorrFunc* tmp,QPropW& seqQ, QPropW& Quark, int Nmom, ThreeMom* momo)
{
  char *fname = "InsertOp()";

  VRB.Func(cname, fname);

  WilsonMatrix v0, v1, v0_next, v1_next;
  WilsonMatrix temp_v0D, temp_v0D_next;
  WilsonMatrix temp_v0_1pD, temp_v0_1mD_next;
  WilsonMatrix temp_v1_u;
  WilsonMatrix temp_v1_next_u;


  int shift_t=GJP.TnodeCoor()*GJP.TnodeSites();
  int vol= (GJP.VolNodeSites()/GJP.TnodeSites());

  int dir = G[G.N()-1];


  const int LORENTZs(4); 
  int lcl_sites[LORENTZs]; 
  lcl_sites[0] = GJP.XnodeSites();
  lcl_sites[1] = GJP.YnodeSites();
  lcl_sites[2] = GJP.ZnodeSites();
  lcl_sites[3] = GJP.TnodeSites();

  int lcl_node[LORENTZs];
  lcl_node[0] = GJP.XnodeCoor();
  lcl_node[1] = GJP.YnodeCoor();
  lcl_node[2] = GJP.ZnodeCoor();
  lcl_node[3] = GJP.TnodeCoor();

  int lcl2glb_offset[LORENTZs];
  for (int i = 0; i < LORENTZs; ++i) 
    lcl2glb_offset[i] = lcl_sites[i] * lcl_node[i];

 
  Lattice& lat = seqQ.AlgLattice();
  Matrix *gauge_field = lat.GaugeField();  

  int ls_glb = GJP.SnodeSites() * GJP.Snodes();
  int ls_ini = 0;
  int ls_end = ls_glb;
  if(GJP.Snodes()==2) {
    ls_ini = GJP.SnodeCoor() * GJP.SnodeSites();
    ls_end = ( GJP.SnodeCoor() + 1 ) * GJP.SnodeSites();
  }

  char dummy[50];
  Site site;

  for(int s=ls_ini;s<ls_end;s++){ 
    //for(int s=0;s<ls_glb;s++){ 
    Quark.RestoreQPropLs(dummy, s);
    if(GJP.Snodes()==2) {
    seqQ.RestoreQPropLs(dummy, ls_end-1-s+ls_ini);
    } else {
    seqQ.RestoreQPropLs(dummy, ls_end-1-s);
    }
    //seqQ.RestoreQPropLs(dummy, ls_glb-1-s);

    // This operations are included in QPropWSeqProtU(D)Src.
    // When we restore 5d prop the operations are nessesary.
    //Multiply by gamma5 and take the dagger to make it in to quark.
    for (site.Begin();site.End();site.nextSite()) {
      seqQ[site.Index()].gl(-5);
      seqQ[site.Index()].hconj();
    }


    int lcl[LORENTZs];
    int lcl_next[LORENTZs]; // Next site along propagation direction


    VRB.Clock(cname,fname,"RO\n");
    for(lcl[3]=0; lcl[3]<=GJP.TnodeSites()-1; lcl[3]++) {
    int t=lcl[3]+shift_t;

    for(lcl[0]=0; lcl[0]<=GJP.XnodeSites()-1; lcl[0]++) 
    for(lcl[1]=0; lcl[1]<=GJP.YnodeSites()-1; lcl[1]++) 
    for(lcl[2]=0; lcl[2]<=GJP.ZnodeSites()-1; lcl[2]++) {

      Float coeff = 1.0;
      if( G.N() == 2 ) {
      coeff = -1.0;
      if(GJP.Snodes()==2) {
	if(GJP.SnodeCoor()==1) coeff = -coeff;
      } else {
	if( s > ls_glb/2 - 1 ) coeff = -coeff;
      }
      }

      int lcl_offset = siteOffset(lcl,lcl_sites);

      // coordinates and offset for lcl_next
      for (int i  = 0; i < LORENTZs; i++ ) 
	lcl_next[i] = ( (i == dir) ? (lcl[i]+1)%lcl_sites[i]
			  : lcl[i] );

      int lcl_next_offset = siteOffset(lcl_next,lcl_sites);

      // U_mu(x) where mu = dir
      Matrix * link = gauge_field + siteOffset(lcl,lcl_sites) * 4 + dir ;



      // S_F(x, s)
      v0=Quark[lcl_offset];
      // S_F^+(x, ls_glb-1-s)
      v1=seqQ[lcl_offset];


      // v0_next = S_F(x+dir, s)
      // v1_next = S_F^+(x+dir, ls_glb-1-s)
      if ((lcl[dir]+1) == lcl_sites[dir]) {
	getPlusData( (IFloat *)&v0_next,
		     (IFloat *)(&Quark[lcl_next_offset]), sizeof(WilsonMatrix)/sizeof(Float),
		     dir) ;
	getPlusData( (IFloat *)&v1_next,
		     (IFloat *)(&seqQ[lcl_next_offset]), sizeof(WilsonMatrix)/sizeof(Float), 
		     dir) ;
	
	// fix boundary condition
	switch( dir ) {
	case 0:
	  if (GJP.XnodeBc()==BND_CND_APRD) coeff = -coeff ;
	  break;
	case 1:
	  if (GJP.YnodeBc()==BND_CND_APRD) coeff = -coeff ;
	  break;
	case 2:
	  if (GJP.ZnodeBc()==BND_CND_APRD) coeff = -coeff ;
	  break;
	case 3:
	  if (GJP.TnodeBc()==BND_CND_APRD) coeff = -coeff ;
	  break;
	} // end switch

      } else {
	v0_next = Quark[lcl_next_offset] ;
	v1_next = seqQ[lcl_next_offset] ;
      }

      // Gamma^{dir} S_F(x, s)
      temp_v0D = v0; temp_v0D.gl(dir);
      // Gamma^{dir} S_F(x+dir, s)
      temp_v0D_next = v0_next; temp_v0D_next.gl(dir);

      // ( 1 + Gamma^{dir} ) S_F(x, s)
      temp_v0_1pD = v0 + temp_v0D;
      // ( 1 - Gamma^{dir} ) S_F(x+dir, s)
      temp_v0_1mD_next = v0_next - temp_v0D_next;

      // S_F^+(x, ls_glb-1-s) U
      temp_v1_u.UMultSource(*link, v1);

      // S_F^+(x+dir, ls_glb-1-s) U^+
      temp_v1_next_u.UdagMultSource(*link, v1_next);

      Complex cc(0);
      //-Tr( ( 1 - Gamma^{dir} ) S_F(x+dir, s) S_F^+(x, ls_glb-1-s) U )
      cc-=coeff*( Trace( temp_v0_1mD_next, temp_v1_u ) );
      //+Tr( ( 1 + Gamma^{dir} ) S_F(x, s) S_F^+(x+dir, ls_glb-1-s) U^+ )
      cc+=coeff*( Trace( temp_v0_1pD, temp_v1_next_u ) );
      Site st(lcl[0],lcl[1],lcl[2],lcl[3]);
      for(int ip=0;ip<Nmom;ip++) {
      Complex cd = cc * momo[ip].Fact(st); // The momentum factor exp(-ipx)
      tmp[ip][t] += cd ;
      } // ip loop
    } // xyz loop

    } // t loop

  } // ls loop

  for(int ip=0;ip<Nmom;ip++) {
    for(int tt=0; tt<=GJP.TnodeSites()-1; tt++) {
      int t=tt+shift_t;
      if(GJP.Snodes()==2) {
	Float rsum, isum;
	rsum = tmp[ip][t].real();
	glb_sum_dir(&rsum, 4);
	isum = tmp[ip][t].imag();
	glb_sum_dir(&isum, 4);
	tmp[ip][t] = Complex( rsum, isum );
      }
      tmp[ip][t]=tmp[ip][t]/2;
    }
  }

}


CPS_END_NAMESPACE
