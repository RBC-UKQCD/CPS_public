// To: Tom
// From: Taku
//
// Here I commented out "#include <u1.h>" and
// defines `twist_links' with `cTimesVec' imported from U1Lattice class.
// So there is some redundancy between this code and "U1.C".
//
// Likewise I changed to use `twist_links' instead of `multiply_u1_link'.
// I won't check if I was doing this correctly, especially about
// `twist_links'. Be careful.
//



//------------------------------------------------------------------
//
// alg_nuc3pt.C
//
// AlgNuc3pt is derived from Alg and is relevant to  
// nuc3pt correlation functions with Wilson-type 
// fermions (i.e. Wilson, Clover, Dwf). 
// The type of fermion is determined by the argument to the 
// constructor but it should not be a non Wilson type lattice.
//
//------------------------------------------------------------------

#include <stdlib.h>	// exit()
#include <stdio.h>
#include <alg/alg_nuc3pt.h>
#include <alg/nuc3pt.h>
#include <alg/meson.h>
//#include <util/u1.h>

CPS_START_NAMESPACE



//-----------------------------------------------------------------
//  Twisted B.C stuff (similar redundant stuff in U1.C, U1.h)
//---------------------------------------------------------------

// A = b*C
void cTimesVec(IFloat *a, 
	       IFloat re, 
	       IFloat im, 
	       const IFloat *c,
	       int len)
{
  for(int i = 0; i < len; i += 2, c += 2) 
    {
      *a++ = re * *c     - im * *(c+1);   // real part
      *a++ = re * *(c+1) + im * *c;       // imag part
    }
}

// multiply exp( - i Q ) in `mu' direction. 
void twist_links(Lattice &lat, const Float Q, const int mu)
{
  // multiply su3 links by u1 links
  // assumes both are in canonical order!
  // link =  exp(-i A_mu Q) where Q is the 
  // quark electric charge or twist angle
  int x[4];
  Complex cc(cos(Q), -sin(Q)); 
  
  Matrix temp;
  for(x[3]=0; x[3]<GJP.TnodeSites();x[3]++){
    for(x[2]=0; x[2]<GJP.ZnodeSites();x[2]++){
      for(x[1]=0; x[1]<GJP.YnodeSites();x[1]++){
	for(x[0]=0; x[0]<GJP.XnodeSites();x[0]++){
	  int offset_x = lat.GsiteOffset(x);
	  // the link matrix at current site, in direction mu:
	  Matrix *su3_link = lat.GaugeField() + offset_x + mu;
	  //cTimesVec((float*)&temp, *(float*)phase, *((float*)phase+1), 
	  cTimesVec((IFloat*)&temp, (IFloat)cc.real(), (IFloat)cc.imag(), 
		    (IFloat*)su3_link, 18);
	  moveMem((IFloat*)su3_link, (IFloat*)&temp, 18);
	}
      }
    }
  }
}




//------------------------------------------------------------------
// Constructor 
//------------------------------------------------------------------
AlgNuc3pt::AlgNuc3pt(Lattice& latt, CommonArg *c_arg, Nuc3ptArg *arg): 
  Alg(latt, c_arg),fp(NULL)
{
  cname = "AlgNuc3pt";
  char *fname = "AlgNuc3pt(L&,CommonArg*,Nuc3ptArg*)";
  VRB.Func(cname,fname);

  arg->check_args() ;

  // Initialize the argument pointer
  //----------------------------------------------------------------
  if(arg == 0) ERR.Pointer(cname,fname, "arg");
  Nuc3pt_arg = arg;
  
  if((Nuc3pt_arg->src_type != BOX)&&(Nuc3pt_arg->src_type != POINT)&&(Nuc3pt_arg->src_type != GAUSS_GAUGE_INV))
    ERR.General(cname,fname, "Unsupported source type!\n");

  // Print out input parameters --- This is irrelevant. Has to go.
  //----------------------------------------------------------------
  VRB.Input(cname,fname,
	    "stop_rsd = %g\n",float (Nuc3pt_arg->cg.stop_rsd));
  VRB.Input(cname,fname,
	    "max_num_iter = %d\n", Nuc3pt_arg->cg.max_num_iter);
  
  //-------------------Added by M.F.Lin start----------------------------// 
  // array of quark propagators for the multi-sink calculations
  // array of size 1 reduces to the ordinary one sink calculation
  
  if ( Nuc3pt_arg-> calc_seqQ == MULT_SEQ 
	|| Nuc3pt_arg -> calc_seqQ == WRITE_MULT_SEQ 
	|| Nuc3pt_arg->calc_seqQ == READ_MULT_SEQ) 
     num_qprop = Nuc3pt_arg->num_src;
  else
     num_qprop = 1;
  
  q_prop = (QPropW**) smalloc(cname,fname,"q_prop",num_qprop*sizeof(QPropW *)) ;
  for ( int n = 0; n < num_qprop; n++ ) q_prop[n] = NULL;
  //-------------------Added by M.F.Lin end------------------------------//
 
  u_s_prop=NULL ;
  d_s_prop=NULL ;
  //???
}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
AlgNuc3pt::~AlgNuc3pt() {
  char *fname = "~AlgNuc3pt()";
  VRB.Func(cname,fname);
  //delete prop if it is still alive
  for ( int n = 0; n < num_qprop; n++ ) 
  {
      if(q_prop[n]!=NULL)
         delete q_prop[n] ;
  }
  sfree(q_prop);
  
  if(u_s_prop!=NULL)
    delete u_s_prop ;
  if(d_s_prop!=NULL)
    delete d_s_prop ;
  //???
}


//------------------------------------------------------------------
//
//------------------------------------------------------------------
void AlgNuc3pt::run()
{
  char *fname = "run()";
  VRB.Func(cname,fname);

  if( Nuc3pt_arg->DoHalfFermion == 2 && Nuc3pt_arg->DoConserved == 1 ) {
    ERR.General(cname,fname,"DoConserved should be 2 when DoHalfFermion = 2!!!\n");
  }

  // Always do point sink...
  Nuc2pt Nuc_c5(NUC_G5C,POINT) ;
  // do a few mesons for kicks
  Meson pion;
  pion.setGamma(-5);
  Meson vector[3];
  vector[0].setGamma(0);
  vector[1].setGamma(1);
  vector[2].setGamma(2);
  char *MesonSrcTag ;
  switch (Nuc3pt_arg->src_type){
  case  POINT:
    MesonSrcTag = "POINT" ;
    break; 
  case  BOX:
    MesonSrcTag = "GFBOX" ;
    break; 
  case  GAUSS_GAUGE_INV:
    MesonSrcTag = "GSmear" ;
    break; 
  default: // It should  never get here! 
    MesonSrcTag = "" ;// make GCC happy!
    ERR.General(cname,fname,"OOPS!!!\n");
    break ;
  }
  pion.setSrc(MesonSrcTag);
  vector[0].setSrc(MesonSrcTag);
  vector[1].setSrc(MesonSrcTag);
  vector[2].setSrc(MesonSrcTag);

  int ZeroMom[3] ;
  ZeroMom[0]=ZeroMom[1]=ZeroMom[2]=0 ;
  
  ProjectType ptype[4] ;
  ptype[0] = PPAR_5X ;
  ptype[1] = PPAR_5Y ;
  ptype[2] = PPAR_5Z ;
  ptype[3] = PPAR_5  ;


  OpenFile();
  Fprintf(fp,"List of momenta:\n");
  int Nmom(56) ;
  ThreeMom sink_mom[56] ;
  int count(0);
  for(int p1=-2;p1<3;p1++)
    for(int p2=-2;p2<3;p2++)
      for(int p3=-2;p3<3;p3++)
	if((p1*p1+p2*p2+p3*p3)<=Nuc3pt_arg->MaxMom2)
	  if((p1*p1+p2*p2+p3*p3)!=0)// eliminate the p=0
	    {
	      Fprintf(fp,"\t\t\t%i: %i %i %i\n",count,p1,p2,p3);
	      sink_mom[count] = ThreeMom(p1,p2,p3);
	      count++ ;
	    }
  CloseFile();
  if(count>56)
    ERR.General(cname,fname,"No of momenta > 56\n");

  //set the number of momenta
  Nmom=count ; 
  // non-zero momentum give bad signal with BOX sources...
  // do not ever do it!!
  if( Nuc3pt_arg->src_type == BOX && (Nuc3pt_arg->BoxEnd - Nuc3pt_arg->BoxStart) != GJP.Xnodes()*GJP.XnodeSites()-1)  Nmom=0 ; 

  Float qmass(Nuc3pt_arg->Mass(0)) ;
  // Loop over source times
  OpenFile();
  Fprintf(fp,"Doing quark Mass: %g\n", qmass); 
  Fprintf(fp,"Doing %d sources inc= %d\n", 
	  Nuc3pt_arg->num_src, Nuc3pt_arg->source_inc); 
  CloseFile();


  // start: ts=t_source, then increment by "source_inc"
  int ts=Nuc3pt_arg->t_source;
  int mt[5];
  
  if ( Nuc3pt_arg->num_mult > 1) ts=Nuc3pt_arg->mt[0]; // use the t_source, source_inc counter if not doing MultGauss. -MFL

  for(int i_source=0; i_source < Nuc3pt_arg->num_src; i_source++)
    {
      int t_sink = (ts + Nuc3pt_arg->t_sink)%(GJP.Tnodes()*GJP.TnodeSites());
      OpenFile();
      Fprintf(fp,"Doing source/sink time slices: %d %d\n", ts, t_sink);
      CloseFile();

      // First calculate the needed two point functions
      Nuc_c5.Zero() ;
      pion.Zero() ;
      vector[0].Zero() ;
      vector[1].Zero() ;
      vector[2].Zero() ;
      char out_prop[200];

      // If we don't do the coherent sink, we only allocate 1 propagator
      // hence we reuse the memory for the propagator for different source locations.
      // Otherwise we need to have all the propagators present at the same time
      // -- MFL
      int n;
      if (Nuc3pt_arg->num_src == num_qprop) n = i_source;
      else n = 0;

      GetThePropagator(n, ts, qmass);

	Nuc_c5.calcNucleon(*q_prop[n]) ;
	pion.setMass(qmass,qmass) ;
	pion.calcMeson(*q_prop[n],*q_prop[n]) ;
	vector[0].setMass(qmass,qmass) ;
	vector[0].calcMeson(*q_prop[n],*q_prop[n]) ;
	vector[1].setMass(qmass,qmass) ;
	vector[1].calcMeson(*q_prop[n],*q_prop[n]) ;
	vector[2].setMass(qmass,qmass) ;
	vector[2].calcMeson(*q_prop[n],*q_prop[n]) ;

	// Print out 2pt function results
	//----------------------------------------------------------------
	OpenFile();
	Nuc_c5.Print(fp) ;
	pion.Print(fp) ;
	vector[0].Print(fp) ;
	vector[1].Print(fp) ;
	vector[2].Print(fp) ;
	CloseFile();
	
	//Do the projections needed for disconnected Ga
	for(int p(0);p<4;p++)
	  {
	    Nuc2pt Nuc_c5_p(NUC_G5C,POINT,ptype[p]) ;
//	    Nuc_c5_p.Zero();
	    
	    if(Nuc3pt_arg->DoGa1Proj||(p>2))
	      {
		Nuc_c5_p.calcNucleon(*q_prop[n]) ;
		
		OpenFile();
		Nuc_c5_p.Print(fp) ;
		CloseFile();
	      }
	    
	  }
	
	//do some non-zero momenta
	{
	  Nuc2pt Nuc_c5_p(NUC_G5C,POINT,PPAR_5) ;
	  for(int m(0);m<Nmom;m++)
	    {
	      
	      Nuc_c5.Zero() ;
	      Nuc_c5.calcMomNucleon(*q_prop[n],sink_mom[m]) ;
	      // Print out 2pt function results
	      //----------------------------------------------------------------
	      OpenFile();
	      Nuc_c5.Print(fp) ;
	      CloseFile();
	      
	      // Calculate the smeared smeared momentum two point function
	      // to extract alpha for  NEDM reweighting	    
	      Nuc_c5_p.Zero() ;
	      Nuc_c5_p.calcMomNucleon(*q_prop[n],sink_mom[m]) ;
	      // Print out 2pt function result
	      //----------------------------------------------------------------
	      OpenFile();
	      Nuc_c5_p.Print(fp) ;
	      CloseFile();


	    }// End momentum
	}

	//Do the smeared sink stuff
	if((Nuc3pt_arg->src_type==GAUSS_GAUGE_INV)&&(Nuc3pt_arg->DoSS2ptF)){
	  
	  Nuc2pt Nuc_c5_ss(NUC_G5C,GAUSS_GAUGE_INV) ;
	  
	  QPropWGaussSrc smr(*q_prop[n]) ;
	  QPropWGaussArg gauss_arg = q_prop[n]->GaussArg();
	  smr.GaussSmearSinkProp(gauss_arg);
	  
	  // First calculate the needed two point functions
	  Nuc_c5_ss.Zero() ;
	  Nuc_c5_ss.calcNucleon(smr) ;
	  
	  pion.Zero() ;
	  pion.calcMeson(smr,smr) ;
	  
	  
	  // Print out 2pt function results
	  //----------------------------------------------------------------
	  OpenFile();
	  Nuc_c5_ss.Print(fp) ;
	  pion.Print(fp) ;
	  CloseFile();
	  
	  //Do the projections needed for disconnected Ga
	  for(int p(0);p<4;p++)
	    {
	      Nuc2pt Nuc_c5_ss_p(NUC_G5C,GAUSS_GAUGE_INV,ptype[p]) ;
	      if(Nuc3pt_arg->DoGa1Proj||(p>2)){
		Nuc_c5_ss_p.calcNucleon(smr) ;
		OpenFile();
		Nuc_c5_ss_p.Print(fp) ;
		CloseFile();
	      }
	    }
	  
	  //do some non-zero momenta
	  {
	    Nuc2pt Nuc_c5_ss_p(NUC_G5C,GAUSS_GAUGE_INV,PPAR_5) ;
	    for(int m(0);m<Nmom;m++)
	      {
		Nuc_c5_ss.Zero() ;
		Nuc_c5_ss.calcMomNucleon(smr,sink_mom[m]) ;
		// Print out 2pt function results
		//----------------------------------------------------------------
		OpenFile();
		Nuc_c5_ss.Print(fp) ;
		CloseFile();
		
		// Calculate the smeared smeared momentum two point function
		// to extract alpha for  NEDM reweighting
		Nuc_c5_ss_p.Zero() ;
		Nuc_c5_ss_p.calcMomNucleon(smr,sink_mom[m]) ;
		// Print out 2pt function results
		//----------------------------------------------------------------
		OpenFile();
		Nuc_c5_ss_p.Print(fp) ;
		CloseFile();
	      }// End momentum
	    
	    
	    if(Nuc3pt_arg->DoHalfFermion == 2) {
	      OpenFile();
	      Fprintf(fp,"DoHalfFermion=%d non-rel nucleon 2-pt\n",Nuc3pt_arg->DoHalfFermion);
	      CloseFile();
	      
	      q_prop[n]->NonRelProp( 1 );
	      
	      {
		Nuc2pt Nuc_c5_p_nr(NUC_G5C,POINT) ;
		Nuc_c5_p_nr.Zero() ;
		Nuc_c5_p_nr.calcNucleon(*q_prop[n]) ;
		OpenFile();
		Nuc_c5_p_nr.Print(fp) ;
		CloseFile();
		
		for(int m(0);m<Nmom;m++) {
		  Nuc_c5_p_nr.Zero() ;
		  Nuc_c5_p_nr.calcMomNucleon(*q_prop[n],sink_mom[m]) ;
		  OpenFile();
		  Nuc_c5_p_nr.Print(fp);
		  CloseFile();
		}
	      }
	      
	      smr.NonRelProp( 0 );
	      
	      int dohalf = Nuc3pt_arg->DoHalfFermion;
	      // First calculate the needed two point functions
	      {
		Nuc2pt Nuc_c5_g_nr(NUC_G5C,GAUSS_GAUGE_INV) ;
		Nuc_c5_g_nr.Zero() ;
		Nuc_c5_g_nr.calcNucleon(smr) ;
		
		// Print out 2pt function results
		//----------------------------------------------------------------
		OpenFile();
		Nuc_c5_g_nr.Print(fp) ;
		CloseFile();
		
		//do some non-zero momenta
		for(int m(0);m<Nmom;m++)
		  {
		    Nuc_c5_g_nr.Zero() ;
		    Nuc_c5_g_nr.calcMomNucleon(smr,sink_mom[m]) ;
		    // Print out 2pt function results
		    //----------------------------------------------------------------
		    OpenFile();
		    Nuc_c5_g_nr.Print(fp) ;
		    CloseFile();
		  }// End momentum
	      }
	    } // DoHalfF = 2
 //         } // End loop over qprop
	  
	  // Finally, smear the sink at time t_sink for 3-pt functions (below) 
	  for(int nt=0; nt<Nuc3pt_arg->num_mult; nt++){
	    // Multi Gauss t_sink set
	    mt[nt]=Nuc3pt_arg->mt[nt]; //save
	    Nuc3pt_arg->mt[nt]+=Nuc3pt_arg->t_sink; // locations of the proton sinks
            Nuc3pt_arg->mt[nt]=Nuc3pt_arg->mt[nt]%(GJP.Tnodes()*GJP.TnodeSites());
	    q_prop[n]->GaussSmearSinkProp(Nuc3pt_arg->mt[nt],q_prop[n]->GaussArg());
          }
	  
	} //end smeared sink
      }


      //Now do the 3pt functions

      // If doing coherent sinks, don't calculate the 3pt functions until all the 
      // forward propagators have been calculated. --MFL
      int do_seq = 0;
//      int t_sink;

      if(Nuc3pt_arg->calc_seqQ != MULT_SEQ 
	&& Nuc3pt_arg->calc_seqQ != WRITE_MULT_SEQ 
	&& Nuc3pt_arg->calc_seqQ != READ_MULT_SEQ) {
	do_seq = 1;
        t_sink = ts + Nuc3pt_arg->t_sink;
      }
      // once all the forward propagators have been calculated,
      // do the sequential propagaotrs. --MFL
      else if (i_source == num_qprop-1) {
	do_seq = 1;
	t_sink = Nuc3pt_arg->t_source + Nuc3pt_arg->t_sink;
      }


      if(Nuc3pt_arg->DoUnPolarized && do_seq)
	{
	  // for conserved vector currents
	  Gamma Gt(T);
	  Nuc3ptCons VectCurr(Gt);
	  Nuc3ptCons *VectCurrp[4];
	  if(Nuc3pt_arg->DoConserved) {
	  for (int i(X);i<4;i++){
	    DIR d = DIR(i) ;
	    Gamma G(d);
	    VectCurrp[i] = new Nuc3ptCons(Nmom,G) ;
	  }
	  }

	  OpenFile();
	  Fprintf(fp,"UnPolarized Zero mom.\n");
	  CloseFile();

  	  //first do the zero momentum un-polarized stuff

	  //up-quark
	  GetTheSeqPropagator(t_sink,qmass,
			      PROT_U_SEQ,ZeroMom,PPAR);

	  //down-quark
	  GetTheSeqPropagator(t_sink,qmass,
			      PROT_D_SEQ,ZeroMom,PPAR);

	  calc_Scalar();         //needed for the sigma term
	  calc_Vector();	 //Vector current
	  calc_X_q_b();          //<x>_q (b) does not need non-zero momentum
	  for(int i(0) ; i<Nmom ; i++) calc_Vector(sink_mom[i]);

	  //conserved current
	  if(Nuc3pt_arg->DoConserved) {
	    if(GJP.Snodes()==2) u_s_prop->SwapQPropLs();

	    for ( int nt = 0; nt < num_qprop; nt++ ) {
	      VectCurr.Calc3pt(*u_s_prop,*q_prop[nt]);
	      for (int i(X);i<4;i++)
		VectCurrp[i]->Calc3pt(*u_s_prop,*q_prop[nt],Nmom,sink_mom);
	    }
	    char* dummy;
	    u_s_prop->RestoreOrgProp(dummy,1);
	    u_s_prop->DeleteQPropLs();

	    if(GJP.Snodes()==2) d_s_prop->SwapQPropLs();
	    
	    for ( int nt = 0; nt < num_qprop; nt++ ) {
	      VectCurr.Calc3pt(*d_s_prop,*q_prop[nt]);
	      for (int i(X);i<4;i++)
	      VectCurrp[i]->Calc3pt(*d_s_prop,*q_prop[nt],Nmom,sink_mom);
	    }
	    
	    d_s_prop->RestoreOrgProp(dummy,1);
	    d_s_prop->DeleteQPropLs();

	    OpenFile();
	    Fprintf(fp,"UnPolarized Zero mom Conserved Vector\n");
	    VectCurr.Print(fp) ;
	    for (int i(X);i<4;i++)
	      VectCurrp[i]->Print(fp,Nmom,sink_mom) ;
	    CloseFile();
	  }

	  if(Nuc3pt_arg->DoConserved) {
	  for (int i(X);i<4;i++)
	    delete VectCurrp[i];
	  }
	}

      if(Nuc3pt_arg->DoUnPolarizedMom && do_seq)
	{
	  OpenFile();
	  Fprintf(fp,"UnPolarized mom.\n");
	  CloseFile();
	  int UnitMom[3];
	  UnitMom[1]=1 ; // check the conventions
	  UnitMom[0]=UnitMom[2]=0 ;
	  
	  //up-quark
	  GetTheSeqPropagator(t_sink,qmass,
			      PROT_U_SEQ,UnitMom,PPAR);
	  //down-quark
	  GetTheSeqPropagator(t_sink,qmass,
			      PROT_D_SEQ,UnitMom,PPAR);

	  calc_Scalar();         //needed for the sigma term
	  calc_Vector();	 //Vector current
	  calc_X_q_b();          //<x>_q (b) does not need non-zero momentum
	  for(int i(0) ; i<Nmom ; i++) calc_Vector(sink_mom[i]);

	  calc_X_q_a(); //<x>_q (a) needs non-zero momentum
	  calc_X2_q(); //<x^2>_q    needs non-zero momentum
	  calc_X3_q(); //<x^3>_q    needs non-zero momentum

	}

      // Polarized 
      if(Nuc3pt_arg->DoPolarized && do_seq)
	{
	  // for conserved axial and vector currents
	  Gamma G5z(G5,Z);
	  Nuc3ptCons AxialCurr(Complex(0.0,1.0),G5z);
	  Nuc3ptCons *AxialCurrp[4];
	  Nuc3ptCons *VectCurrp[4];
	  if(Nuc3pt_arg->DoConserved) {
	  for (int i(X);i<4;i++){
	    DIR d = DIR(i);
	    Gamma G5d(G5,d);
	    AxialCurrp[i] = new Nuc3ptCons(Nmom,Complex(0.0,1.0),G5d);

	    d = DIR(i) ;
	    Gamma G(d);
	    VectCurrp[i] = new Nuc3ptCons(Nmom,G) ;
	  }
	  }

	  OpenFile();
	  Fprintf(fp,"Polarized Zero mom.\n");
	  CloseFile();
	  //
	  //up-quark
	  GetTheSeqPropagator(t_sink,qmass,
			      PROT_U_SEQ,ZeroMom,PPAR_5Z);

	  //down-quark
	  GetTheSeqPropagator(t_sink,qmass,
			      PROT_D_SEQ,ZeroMom,PPAR_5Z);
	  
	  calc_Axial(); //Axial-Vector current i Gamma_5 Gamma_z
	  calc_Tensor() ;//Tensor current i Gamma_5 Gamma_z Gamma_t
	  calc_X_Dq_b() ; 
	  calc_d1() ;
	  
	  for(int i(0) ; i<Nmom ; i++){
	    calc_Vector(sink_mom[i]);//neutron EDM and mangetic moment
	    calc_Axial(sink_mom[i]);//axial form factors
	    calc_EnergyMomentum(sink_mom[i]) ; // proton spin
	    calc_PScalar(sink_mom[i]);//pseudoscalar form factors
	  }

	  if(Nuc3pt_arg->DoConserved) {

	    if(GJP.Snodes()==2) {
	      u_s_prop->SwapQPropLs();
	      d_s_prop->SwapQPropLs();
	    }

	    for ( int nt = 0; nt < num_qprop; nt++ ) {

	      AxialCurr.Calc3pt(*u_s_prop,*q_prop[nt]);
	      for(int i(X);i<4;i++){
		AxialCurrp[i]->Calc3pt(*u_s_prop,*q_prop[nt],Nmom,sink_mom);
		VectCurrp[i]->Calc3pt(*u_s_prop,*q_prop[nt],Nmom,sink_mom);
	      }
	      
	      AxialCurr.Calc3pt(*d_s_prop,*q_prop[nt]);
	      for(int i(X);i<4;i++){
		AxialCurrp[i]->Calc3pt(*d_s_prop,*q_prop[nt],Nmom,sink_mom);
		VectCurrp[i]->Calc3pt(*d_s_prop,*q_prop[nt],Nmom,sink_mom);
	      }

	      OpenFile();
	      Fprintf(fp,"Polarized Zero mom Conserved Axial and Vector\n");
	      AxialCurr.Print(fp) ;
	      for (int i(X);i<4;i++) {
		VectCurrp[i]->Print(fp,Nmom,sink_mom) ;
	      }
	      for (int i(X);i<4;i++) {
		AxialCurrp[i]->Print(fp,Nmom,sink_mom) ;
	      }
	      CloseFile();
	    }
	    
	    char* dummy;
	    u_s_prop->RestoreOrgProp(dummy,1);
	    u_s_prop->DeleteQPropLs();
	    d_s_prop->RestoreOrgProp(dummy,1);
	    d_s_prop->DeleteQPropLs();
	  }
	  	  
	  if(Nuc3pt_arg->DoConserved) {
	    for(int i(X);i<4;i++) {
	      delete AxialCurrp[i];
	      delete VectCurrp[i];
	    }
	  }
	}

      
      if(Nuc3pt_arg->DoPolarizedMom && do_seq)
	{
	  OpenFile();
	  Fprintf(fp,"Polarized with mom.\n");
	  CloseFile();
	  int UnitMom[3];
	  UnitMom[1]=1 ; // check the conventions
	  UnitMom[0]=UnitMom[2]=0 ;
	  
	  //up-quark
	  GetTheSeqPropagator(t_sink,qmass,
			      PROT_U_SEQ,UnitMom,PPAR_5Z);
	  //down-quark
	  GetTheSeqPropagator(t_sink,qmass,
			      PROT_D_SEQ,UnitMom,PPAR_5Z);
	  
	  for(int i(0) ; i<Nmom ; i++){
	    calc_Vector(sink_mom[i]);//neutron EDM and magnetic moment
	    calc_Axial(sink_mom[i]);//axial form factors
	    calc_EnergyMomentum(sink_mom[i]) ; // proton spin
	  }
	  calc_X_Dq_a() ;
	  calc_X2_Dq() ;
	  calc_X_dq() ;
	  calc_d2() ;
	}

      if(Nuc3pt_arg->DoConserved) q_prop[n]->DeleteQPropLs();

      ts+=Nuc3pt_arg->source_inc;
      for(int nt=0; nt<Nuc3pt_arg->num_mult; nt++){
	Nuc3pt_arg->mt[nt]=mt[nt];
	Nuc3pt_arg->mt[nt]+=Nuc3pt_arg->source_inc;
      }

/*
      if(i_source==1){ // obsolete. commented out  --MFL
	ts+=Nuc3pt_arg->mt[4];
	for(int nt=0; nt<Nuc3pt_arg->num_mult; nt++){
	  Nuc3pt_arg->mt[nt]+=Nuc3pt_arg->mt[4];
	}
      }
*/
    } // end loop over source timeslices
} // end run


/*!
  Computes the Quark propagator needed for both the 2pt and the 3pt functions

  \li int time is the source time slice
  \li Float mass is the quark mass

 */
void AlgNuc3pt::GetThePropagator(int n, int time, Float mass){

  char *fname = "GetThePropagator(i,f)";
  VRB.Func(cname,fname);

  QPropWArg qp_arg ;
  qp_arg.cg = (Nuc3pt_arg->cg);
  qp_arg.cg.mass = mass ;
  qp_arg.t = time ;

  qp_arg.StartSrcSpin = 0 ;
  qp_arg.EndSrcSpin = 4 ;
  qp_arg.StartSrcColor = 0 ;
  qp_arg.EndSrcColor = 3 ;

//  qp_arg.StartSrcSpin  = (Nuc3pt_arg->StartSrcSpin);
//  qp_arg.EndSrcSpin    = (Nuc3pt_arg->EndSrcSpin);
//  qp_arg.StartSrcColor = (Nuc3pt_arg->StartSrcColor);
//  qp_arg.EndSrcColor   = (Nuc3pt_arg->EndSrcColor);

  if (qp_arg.StartSrcSpin >= qp_arg.EndSrcSpin 
	|| qp_arg.StartSrcColor >= qp_arg.EndSrcColor)   
  ERR.General(cname,fname,"spin and color indices out of range, [%d, %d], [%d, %d]\n", 
	qp_arg.StartSrcSpin, qp_arg.EndSrcSpin, qp_arg.StartSrcColor, qp_arg.EndSrcColor);

  //No gaugefixing at the sink. Derivative operators will fail if gauge fixed
  qp_arg.gauge_fix_snk=0; 

  //Do half fermions ?
  qp_arg.do_half_fermion = Nuc3pt_arg->DoHalfFermion ;
  if(Nuc3pt_arg->DoHalfFermion == 2) qp_arg.do_half_fermion = 0;
   // for both
  if(Nuc3pt_arg->DoHalfFermion == 1) qp_arg.EndSrcSpin = 2;
 
  //Save 5d prop in memory
  if(Nuc3pt_arg->DoConserved) qp_arg.save_ls_prop = 2;

  // delete q_prop if it is allocated
  if(q_prop[n] != NULL)
      delete q_prop[n] ;
  
  // save prop
  qp_arg.save_prop=0;
  if(Nuc3pt_arg->calc_QProp==WRITE_QPROP) {
    qp_arg.save_prop=1;
    char chartmp[200];
    sprintf(chartmp,"Nuc3pt");
    qp_arg.ensemble_id = chartmp;
    qp_arg.ensemble_label = Nuc3pt_arg->ensemble_label;
    qp_arg.seqNum=Nuc3pt_arg->ensemble_id;
  }

  int mu;

  switch (Nuc3pt_arg->src_type){
  case BOX:
    {
      QPropWBoxArg box_arg;
      qp_arg.gauge_fix_src=1; // Gauge fixing is needed for Box source
      box_arg.box_start = Nuc3pt_arg->BoxStart;
      box_arg.box_end   = Nuc3pt_arg->BoxEnd;
      if(Nuc3pt_arg->DoPerPlusAper){
	  BndCndType save_bc = GJP.Tbc() ; // Save current BC
	  GJP.Tbc(BND_CND_PRD);
	  q_prop[n]             =new QPropWBoxSrc(AlgLattice(),&qp_arg, &box_arg,common_arg);
	  GJP.Tbc(BND_CND_APRD);
	  QPropWBoxSrc *aprop=new QPropWBoxSrc(AlgLattice(),&qp_arg, &box_arg,common_arg);
	  q_prop[n]->Average(*aprop) ;
	  delete aprop ;
	  GJP.Tbc(save_bc); // Restore original BC
      }else{
	// impose twisted b.c. if requested
	if(Nuc3pt_arg->theta != 0){
	  
	  for(mu=0;mu<3;mu++){
	    twist_links(AlgLattice(), 
			Nuc3pt_arg->theta, mu);
	  }
	}
	  q_prop[n] = new QPropWBoxSrc(AlgLattice(),&qp_arg, &box_arg, common_arg);
	  q_prop[n]->Run();
      }
    }
    break ;
  case POINT:
    {
      qp_arg.gauge_fix_src=0; // No Gauge fixing is needed for Point source
      qp_arg.x = Nuc3pt_arg->x[0];
      qp_arg.y = Nuc3pt_arg->x[1];
      qp_arg.z = Nuc3pt_arg->x[2];
      if(Nuc3pt_arg->DoPerPlusAper){
	  BndCndType save_bc = GJP.Tbc() ; // Save current BC
	  GJP.Tbc(BND_CND_PRD);
	  q_prop[n]            =new QPropWPointSrc(AlgLattice(),&qp_arg,common_arg);
	  GJP.Tbc(BND_CND_APRD);
	  QPropWPointSrc *apr=new QPropWPointSrc(AlgLattice(),&qp_arg,common_arg);
	  q_prop[n]->Average(*apr) ;
	  delete apr ;
	  GJP.Tbc(save_bc); // Restore original BC
      }else{
	// impose twisted b.c. if requested
	if(Nuc3pt_arg->theta != 0){
	  for(mu=0;mu<3;mu++){
	    twist_links(AlgLattice(), 
			     Nuc3pt_arg->theta, mu);
	  }
	}
	  q_prop[n] = new QPropWPointSrc(AlgLattice(),&qp_arg,common_arg);
	  q_prop[n]->Run();
      }
    }
    break ;
  case GAUSS_GAUGE_INV:
    {
      QPropWGaussArg gauss_arg;
      qp_arg.gauge_fix_src=0; // No Gauge fixing is needed for Point source
      qp_arg.x = Nuc3pt_arg->x[0];
      qp_arg.y = Nuc3pt_arg->x[1];
      qp_arg.z = Nuc3pt_arg->x[2];
      gauss_arg.gauss_N  =   Nuc3pt_arg->gauss_N ;
      gauss_arg.gauss_W  =   Nuc3pt_arg->gauss_W ;
      // Multi Gauss
      gauss_arg.nt = Nuc3pt_arg->num_mult ;
      for(int nt=0; nt<gauss_arg.nt; nt++) gauss_arg.mt[nt] = Nuc3pt_arg->mt[nt] ;
      //Ape Smearing
      gauss_arg.gauss_link_smear_type = Nuc3pt_arg->gauss_link_smear_type;
      gauss_arg.gauss_link_smear_coeff = Nuc3pt_arg->gauss_link_smear_coeff;
      gauss_arg.gauss_link_smear_N = Nuc3pt_arg->gauss_link_smear_N;

      if(Nuc3pt_arg->DoPerPlusAper){
	  BndCndType save_bc = GJP.Tbc() ; // Save current BC
	  GJP.Tbc(BND_CND_PRD);
	  q_prop[n]              =new QPropWMultGaussSrc(AlgLattice(),&qp_arg,&gauss_arg,common_arg);
	  GJP.Tbc(BND_CND_APRD);
	
	  QPropWMultGaussSrc *apr=new QPropWMultGaussSrc(AlgLattice(),&qp_arg,&gauss_arg, common_arg);
	  q_prop[n]->Average(*apr) ;
	  delete apr ;
	  GJP.Tbc(save_bc); // Restore original BC
      }else{
	// impose twisted b.c. if requested
	if(Nuc3pt_arg->theta != 0){
	  //printf("Made it here\n");
	  for(mu=0;mu<3;mu++){
	    twist_links(AlgLattice(), 
			Nuc3pt_arg->theta, mu);
	  }
	}

	char out_prop[200];
	
	qp_arg.file = out_prop;

	 sprintf(out_prop,"%s%i",Nuc3pt_arg->prop_file,time);
	  
	 Fprintf(stdout, "prop outfile = %s\n", qp_arg.file);
	
	  if(Nuc3pt_arg->calc_QProp != READ_QPROP){
	    q_prop[n] = new QPropWMultGaussSrc(AlgLattice(),&qp_arg,&gauss_arg,common_arg);
	  } else {
	    q_prop[n] = new QPropWGaussSrc(AlgLattice(),&qp_arg,&gauss_arg,common_arg,qp_arg.file);
	    q_prop[n]->Allocate(0);
//	    q_prop[n]->ReLoad(qp_arg.file); //TODO: need different filenames for different propagators (MFL)
	    q_prop[n]->RestoreQProp(out_prop,0);
	  }
      }
    }
    break ;
/*
  case SUM_MOM:
    {
      qp_arg.gauge_fix_src=1;
      // make sources with these momenta, and sum
      int **mom;
      int nsrc = 3;
      mom = (int**) smalloc(nsrc * sizeof(int*) );
      for(int i=0;i<nsrc;i++) mom[i] = (int*) smalloc(3*sizeof(int));
            
      mom[0][0]=1; mom[0][1]=0; mom[0][2]=0;
      mom[1][0]=0; mom[1][1]=1; mom[1][2]=0;
      mom[2][0]=0; mom[2][1]=0; mom[2][2]=1;
   
      // impose twisted b.c. if requested
      if(Nuc3pt_arg->theta != 0){
	//printf("Made it here\n");
	for(mu=0;mu<3;mu++){
	  twist_links(AlgLattice(), Nuc3pt_arg->theta, mu);
	}
      }
      q_prop = new QPropWSumMomSrc(AlgLattice(),&qp_arg,nsrc,mom,common_arg);
    }
    break ;
 */
 default: // It should  never get here! 
    ERR.General(cname,fname,"Invalid source type\n");
    break ;
  }
}

/*!
  Computes the sequential propagator using the quark propagator 
  stored in the q_prop.

  \li int time is the sink time slice
  \li Float mass is the quark mass
  \li SourceType type specifies up quark or down quark sequential source 
      propagators.
  \li int *mom specifies the momentum at the sink 
  \li ProjectType P specifies the type of projection to be used.
*/
void  AlgNuc3pt::GetTheSeqPropagator(int time, Float mass, SourceType type, 
				     int *mom, ProjectType P)
{
  char *fname = "GetTheSeqPropagator(i,f,...)";
  VRB.Func(cname,fname);

  QPropWArg qp_arg ;
  qp_arg.cg = (Nuc3pt_arg->cg);
  qp_arg.cg.mass = mass ;
  qp_arg.t = time ;

  qp_arg.StartSrcSpin  = 0;
  qp_arg.EndSrcSpin    = 4;
  qp_arg.StartSrcColor = 0;
  qp_arg.EndSrcColor   = 3;

//  qp_arg.StartSrcSpin  = (Nuc3pt_arg->StartSrcSpin);
//  qp_arg.EndSrcSpin    = (Nuc3pt_arg->EndSrcSpin);
//  qp_arg.StartSrcColor = (Nuc3pt_arg->StartSrcColor);
//  qp_arg.EndSrcColor   = (Nuc3pt_arg->EndSrcColor);

  // for multi-sink calculations. -MFL
  int *t_sink=(int *) smalloc(cname,fname,"t_sink",num_qprop*sizeof(int));
  
  for ( int t=0; t<num_qprop; t++ ) {
	t_sink[t] = time + Nuc3pt_arg->source_inc * t;
	Fprintf(stdout,"=====sink at %d\n",t_sink[t]);
  }
 
  //No gaugefixing at the sink. Derivative operators will fail if gauge fixed
  qp_arg.gauge_fix_snk=0; 
  //The Sequential source code figures out if the source need GF
  //It only need if the sink of q_prop is gauge fixed 
  //(which is not in our case)
 
  //Do half fermions ?
  qp_arg.do_half_fermion = Nuc3pt_arg->DoHalfFermion ;
 
  //Multi Gauss
  QPropWGaussArg gauss_arg;
  gauss_arg.nt = Nuc3pt_arg->num_mult;
  for(int nt=0; nt<gauss_arg.nt; nt++) gauss_arg.mt[nt] = Nuc3pt_arg->mt[nt] ;
  //Ape Smearing
  gauss_arg.gauss_link_smear_type = Nuc3pt_arg->gauss_link_smear_type;
  gauss_arg.gauss_link_smear_coeff = Nuc3pt_arg->gauss_link_smear_coeff;
  gauss_arg.gauss_link_smear_N = Nuc3pt_arg->gauss_link_smear_N;

  qp_arg.save_prop=0;
  if(Nuc3pt_arg->calc_seqQ==WRITE_SEQ || Nuc3pt_arg->calc_seqQ==WRITE_MULT_SEQ) {
    qp_arg.save_prop=1;
    char chartmp[200];
    sprintf(chartmp,"Nuc3pt");
    qp_arg.ensemble_id = chartmp;
    qp_arg.ensemble_label = Nuc3pt_arg->ensemble_label;
    qp_arg.seqNum=Nuc3pt_arg->ensemble_id;
  }

  char out_prop[200];
  char prjct[10];
  switch (P)
    {
    case PPAR:
      { sprintf( prjct, "PPAR" );
      }
      break;
    case PPAR_5Z:
      { sprintf( prjct, "PPAR5Z" );
      }
      break;
    default: // It should  never get here! 
      ERR.General(cname,fname,"OOOPS!!!\n");
    }

  if(q_prop==NULL)
    ERR.General(cname,fname,"OOPS! q_prop not allocated!!\n");
  else
    {
      for(int nt=0; nt<num_qprop; nt++) 
	if(q_prop[nt]==NULL)
	ERR.General(cname,fname,"OOPS! q_prop[%d] not allocated!!\n",nt);
    }

  switch (type)
    {
    case PROT_U_SEQ: // U-quark sequential source 
      {
	// delete u_s_prop if it is allocated
	if(u_s_prop!=NULL)
	  delete u_s_prop ;
	if(Nuc3pt_arg->DoPerPlusAper){
	    BndCndType save_bc = GJP.Tbc() ; // Save current BC
	    GJP.Tbc(BND_CND_PRD);
	    u_s_prop = new QPropWMultSeqProtUSrc(AlgLattice(),num_qprop, q_prop, mom, P, 
					     &qp_arg,&gauss_arg,common_arg,t_sink);
	    GJP.Tbc(BND_CND_APRD);
	    QPropWMultSeqProtUSrc *aprop = 
	      new QPropWMultSeqProtUSrc(AlgLattice(),num_qprop, q_prop,mom,P,&qp_arg,&gauss_arg,common_arg,t_sink);
	    u_s_prop->Average(*aprop) ;
	    delete aprop ;
	    GJP.Tbc(save_bc); // Restore original BC

	}
	else{
	  //Save 5d prop in memory
	  if(Nuc3pt_arg->DoConserved) qp_arg.save_ls_prop = 2;
          qp_arg.file=out_prop;

	  //Assuming we place the sources and sinks evenly on the lattice.
	  //There is little reason why we don't want to do that. --MFL
	  sprintf( out_prop, "prop_seq_u_m%g_tsrc%d_tsnk%d_GS_w%g_n%d_%s_%s_%d",mass,q_prop[0]->SourceTime(),time,Nuc3pt_arg->gauss_W,Nuc3pt_arg->gauss_N,prjct,Nuc3pt_arg->ensemble_label,Nuc3pt_arg->ensemble_id );

	 //commented out for the moment. --MFL 
	 //if(Nuc3pt_arg->num_mult==2) sprintf( out_prop, "prop_seq_u_m%g_tsrc%d_%d_tsnk%d_%d_GS_w%g_n%d_%s_%s_%d",mass,q_prop->SourceTime(),Nuc3pt_arg->mt[1]-Nuc3pt_arg->t_sink,Nuc3pt_arg->mt[0],Nuc3pt_arg->mt[1],Nuc3pt_arg->gauss_W,Nuc3pt_arg->gauss_N,prjct,Nuc3pt_arg->ensemble_label,Nuc3pt_arg->ensemble_id );
	  if(Nuc3pt_arg->calc_seqQ != READ_SEQ){
	    u_s_prop = new QPropWMultSeqProtUSrc(AlgLattice(), num_qprop, q_prop, mom, P, 
					     &qp_arg,&gauss_arg, common_arg,t_sink);
	  } else {
	    u_s_prop = new QPropWMultSeqProtUSrc(AlgLattice(), 1, q_prop, mom, P, &qp_arg, &gauss_arg, common_arg, out_prop);
	    u_s_prop->Allocate(0);
//	    u_s_prop->RestoreQProp(out_prop,0);
	    u_s_prop->ReLoad(out_prop);	
	  }
	}
      }
      break;
    case PROT_D_SEQ: // D-quark sequential source 
      {
	// delete d_s_prop if it is allocated
	if(d_s_prop!=NULL)
	  delete d_s_prop ;
	if(Nuc3pt_arg->DoPerPlusAper){
	  BndCndType save_bc = GJP.Tbc() ; // Save current BC
	  GJP.Tbc(BND_CND_PRD);
	  d_s_prop = new QPropWMultSeqProtDSrc(AlgLattice(),num_qprop, q_prop, mom, P, 
					   &qp_arg,&gauss_arg, common_arg,t_sink);
	  GJP.Tbc(BND_CND_APRD);
	  QPropWMultSeqProtDSrc *aprop = 
	    new QPropWMultSeqProtDSrc(AlgLattice(),num_qprop, q_prop, mom,P,&qp_arg,&gauss_arg, common_arg,t_sink);
	  d_s_prop->Average(*aprop) ;
	  delete aprop ;
	  GJP.Tbc(save_bc); // Restore original BC
	}
	else{
	  //Save 5d prop in file
          if(Nuc3pt_arg->DoConserved == 1) {
	    qp_arg.save_ls_prop = 1;
	    if (qp_arg.file == NULL) 
	      qp_arg.file = (char*)smalloc(100*sizeof(char));
	    sprintf(qp_arg.file,"prop_d");
	  }
	  if(Nuc3pt_arg->DoConserved == 2) qp_arg.save_ls_prop = 2;
          qp_arg.file=out_prop;

	  sprintf( out_prop, "prop_seq_d_m%g_tsrc%d_tsnk%d_GS_w%g_n%d_%s_%s_%d",mass,q_prop[0]->SourceTime(),time,Nuc3pt_arg->gauss_W,Nuc3pt_arg->gauss_N,prjct,Nuc3pt_arg->ensemble_label,Nuc3pt_arg->ensemble_id );
	  //	  if(Nuc3pt_arg->num_mult==2) sprintf( out_prop, "prop_seq_d_m%g_tsrc%d_%d_tsnk%d_%d_GS_w%g_n%d_%s_%s_%d",mass,q_prop->SourceTime(),Nuc3pt_arg->mt[1]-Nuc3pt_arg->t_sink,Nuc3pt_arg->mt[0],Nuc3pt_arg->mt[1],Nuc3pt_arg->gauss_W,Nuc3pt_arg->gauss_N,prjct,Nuc3pt_arg->ensemble_label,Nuc3pt_arg->ensemble_id );

	  if(Nuc3pt_arg->calc_seqQ != READ_SEQ){
	    d_s_prop = new QPropWMultSeqProtDSrc(AlgLattice(),num_qprop, q_prop, mom, P, 
					   &qp_arg, &gauss_arg,common_arg,t_sink);
	  } else {
	    d_s_prop = new QPropWMultSeqProtDSrc(AlgLattice(), 1, q_prop, mom, P, &qp_arg, &gauss_arg, common_arg, out_prop);
	    d_s_prop->Allocate(0);
	    d_s_prop->ReLoad(out_prop);
	  }
	}
      }
      break ;
    default: // It should  never get here! 
      ERR.General(cname,fname,"OOOPS!!!\n");
    }
 sfree(t_sink);
}


void AlgNuc3pt::OpenFile()
{
  char *fname="OpenFile()" ;
  if(common_arg->results != 0)
    {
      if((fp = Fopen((char *)common_arg->results, "a")) == NULL )
	ERR.FileA(cname,fname, (char *)common_arg->results);
    }
}


/*!
  Computes the vector current using the
  \f[
  {\cal O}_\mu = \overline{q} \gamma_\mu q
  \f]
  
  all possible mometa are inserted
*/
void AlgNuc3pt::calc_Vector(const ThreeMom& q)
{     
  for (int i(X);i<4;i++){
    for ( int n = 0; n < num_qprop; n++ ) {
      DIR d = DIR(i) ;
      Gamma G(d);
      Nuc3ptGamma VectCurr(q,G) ;
      VectCurr.Calc3pt(*u_s_prop,*q_prop[n]);
      VectCurr.Calc3pt(*d_s_prop,*q_prop[n]);
      OpenFile(); 
      VectCurr.Print(fp) ;
      CloseFile();
    }
  }
}

/*!
  Computes the Axial vector current using the
  \f[
  {\cal O}_{5\mu} = i \overline{q} \gamma_5 \gamma_\mu q
  \f]  

  all possible mometa are inserted
*/
void AlgNuc3pt::calc_Axial(const ThreeMom& q)
{ 
  for ( int n = 0; n < num_qprop; n++ ) {
    for (int i(X);i<4;i++){
      DIR d = DIR(i) ;
      Gamma G5d(G5,d);
      Nuc3ptGamma AxialCurr(q,Complex(0.0,1.0),G5d) ;
      AxialCurr.Calc3pt(*u_s_prop,*q_prop[n]);
      AxialCurr.Calc3pt(*d_s_prop,*q_prop[n]);
      
      OpenFile();
    AxialCurr.Print(fp) ;
    CloseFile();
    }
  }
}

/*!
  Computes the Axial vector current using the
  \f[
  {\cal O}_{5} = \overline{q} \gamma_5 q
  \f]  

  all possible mometa are inserted
*/
void AlgNuc3pt::calc_PScalar(const ThreeMom& q)
{ 
  for ( int n = 0; n < num_qprop; n++ ) {
    Gamma G(G5);
    Nuc3ptGamma PScalarCurr(q,G) ;
    PScalarCurr.Calc3pt(*u_s_prop,*q_prop[n]);
    PScalarCurr.Calc3pt(*d_s_prop,*q_prop[n]);
    
    OpenFile();
    PScalarCurr.Print(fp) ;
    CloseFile();
  }
}


/*!
  Computes the first unpolarized moment
  
  \f[
  {\cal O}_{14} = \overline{q} \left[ \gamma_1 
  \stackrel{\displaystyle \leftrightarrow}{D}_4
  +\gamma_4 
  \stackrel{\displaystyle \leftrightarrow}{D}_1
  \right] q
  \f]
  
  with polarized projector
*/
void AlgNuc3pt::calc_EnergyMomentum(const ThreeMom& mom)
{
  OpenFile(); 
  Fprintf(fp,"The next is: Energy momentum k4 + 4k\n");
  CloseFile();
  for ( int n = 0; n < num_qprop; n++ ) {
    for (int i(X);i<4;i++){
      DIR d = DIR(i) ;
      Gamma Gx(d);
      Gamma Gt(T);
      Derivative Der_t(T);
      Derivative Der_x(d);
      
      Nuc3ptStru Xq_xt(mom,Gx, Der_t);
      Xq_xt.Calc3pt(*u_s_prop, *q_prop[n]);
      Xq_xt.Calc3pt(*d_s_prop, *q_prop[n]);
      
      Nuc3ptStru Xq_tx(mom,Gt,Der_x);
      Xq_tx.Calc3pt(*u_s_prop, *q_prop[n]);
      Xq_tx.Calc3pt(*d_s_prop, *q_prop[n]);
      
      Xq_xt += Xq_tx ;
      OpenFile(); 
      Xq_xt.Print(fp) ; 
      CloseFile();
    }
  }
}

/*!
  Computes the conserved vector current using the
  \f[
  {\cal O}_\mu = \overline{q} \gamma_\mu q
  \f]
  
  all possible mometa are inserted
*/
void AlgNuc3pt::calc_Cons_Vector(int Nmom, ThreeMom* mom)
{     
  Gamma Gt(T);
  Nuc3ptCons VectCurr(Gt);
  const int MaxNmom=50;
  if(Nmom>MaxNmom)
    ERR.General(cname,"calc_Cons_Vector","Nmom(%d)>MaxNmom(%d)",Nmom,MaxNmom);
  Nuc3ptCons *VectCurrp[MaxNmom][4];
  for(int ip(0);ip<Nmom;ip++) for (int i(X);i<4;i++){
      DIR d = DIR(i) ;
      Gamma G(d);
      VectCurrp[ip][i] = new Nuc3ptCons(mom[ip],G) ;
    }
  
  for ( int n = 0; n < num_qprop; n++ ) {
    QPropW* quark = new QPropWGaussSrc(*q_prop[n]);
    
    VectCurr.Calc3pt(*u_s_prop,*quark);
    for(int ip(0);ip<Nmom;ip++) 
      for (int i(X);i<4;i++)
	VectCurrp[ip][i]->Calc3pt(*u_s_prop,*quark);
    
    u_s_prop->DeleteQPropLs();
    
    if(Nuc3pt_arg->DoConserved == 1) {
      char dummy[30];
      d_s_prop->RestoreQPropLs_ftom(dummy);
    }
    
    VectCurr.Calc3pt(*d_s_prop,*quark);
    for(int ip(0);ip<Nmom;ip++) for (int i(X);i<4;i++)
				  VectCurrp[ip][i]->Calc3pt(*d_s_prop,*quark);
    
    d_s_prop->DeleteQPropLs();
    
    OpenFile(); 
    VectCurr.Print(fp) ;
    for(int ip(0);ip<Nmom;ip++) for (int i(X);i<4;i++)
				  VectCurrp[ip][i]->Print(fp) ;
    CloseFile();
    
    delete quark;
  }

  for(int ip(0);ip<Nmom;ip++) 
    for (int i(X);i<4;i++)
      delete VectCurrp[ip][i];

}

/*!
  Computes the conserved Axial vector current using the
  \f[
  {\cal O}_{5\mu} = i \overline{q} \gamma_5 \gamma_\mu q
  \f]  

  all possible mometa are inserted
*/
void AlgNuc3pt::calc_Cons_Axial_Vector(int Nmom, ThreeMom* mom)
{ 
  Gamma G5z(G5,Z);
  Nuc3ptCons AxialCurr(Complex(0.0,1.0),G5z);
  const int MaxNmom=50;
  if(Nmom>MaxNmom)
    ERR.General(cname,"calc_Cons_Vector","Nmom(%d)>MaxNmom(%d)",Nmom,MaxNmom);
  Nuc3ptCons *AxialCurrp[MaxNmom][4];
  for(int ip(0);ip<Nmom;ip++) for (int i(X);i<4;i++){
    DIR d = DIR(i);
    Gamma G5d(G5,d);
    AxialCurrp[ip][i] = new Nuc3ptCons(mom[ip],Complex(0.0,1.0),G5d);
  }

  Nuc3ptCons *VectCurrp[MaxNmom][4];
  for(int ip(0);ip<Nmom;ip++) for (int i(X);i<4;i++){
    DIR d = DIR(i) ;
    Gamma G(d);
    VectCurrp[ip][i] = new Nuc3ptCons(mom[ip],G) ;
  }

  for ( int n = 0; n < num_qprop; n++ ) {
    QPropW* quark = new QPropWGaussSrc(*q_prop[n]);
    
    AxialCurr.Calc3pt(*u_s_prop,*quark);
    for(int ip(0);ip<Nmom;ip++) for (int i(X);i<4;i++){
	AxialCurrp[ip][i]->Calc3pt(*u_s_prop,*quark);
	VectCurrp[ip][i]->Calc3pt(*u_s_prop,*quark);
      }
    
    u_s_prop->DeleteQPropLs();
    
    if(Nuc3pt_arg->DoConserved == 1) {
      char dummy[30];
      d_s_prop->RestoreQPropLs_ftom(dummy);
    }
    
    AxialCurr.Calc3pt(*d_s_prop,*quark);
    for(int ip(0);ip<Nmom;ip++) for (int i(X);i<4;i++) {
	AxialCurrp[ip][i]->Calc3pt(*d_s_prop,*quark);
	VectCurrp[ip][i]->Calc3pt(*d_s_prop,*quark);
      }

    d_s_prop->DeleteQPropLs();
    
    OpenFile();
    AxialCurr.Print(fp) ;
    for(int ip(0);ip<Nmom;ip++) {
      for (int i(X);i<4;i++) {
	VectCurrp[ip][i]->Print(fp) ;
      }
      for (int i(X);i<4;i++) {
	AxialCurrp[ip][i]->Print(fp) ;
      }
    }
    CloseFile();

    delete quark;
  }

  for(int ip(0);ip<Nmom;ip++) for (int i(X);i<4;i++) {
      delete AxialCurrp[ip][i];
      delete VectCurrp[ip][i];
    }

}

CPS_END_NAMESPACE
