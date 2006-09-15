//------------------------------------------------------------------
//
// QPropW.C
//
// Kostas Orginos  (February 2002)
//
// The class functions for QPropW.
//
// For now this is specific to three colors. The constructor
// will exit if the number of colors is not equal to three.
//
// In AlgThreept the QPropWRandWallSrc has to be replaced with 
// QPropWRandSlabSrc
//
//
//------------------------------------------------------------------

#include <stdlib.h>     // exit()
#include <stdio.h>
#include <string.h>
#include <alg/common_arg.h>
#include <comms/glb.h>
#include <comms/scu.h>
// #include <util/data_io.h>
#include <comms/sysfunc.h>

#include <fcntl.h>      // read and write control flags,
#include <unistd.h>     // close(). These are needed for io parts to
                        // compile on PCs
#include <alg/qpropw.h>
#include <util/qcdio.h>

CPS_START_NAMESPACE

// Allocate space for propagators
void QPropW::Allocate(int mid) {

  char *fname = "Allocate(int)";
  VRB.Func(cname, fname);

  if (!mid) {
    if (prop == NULL) { // Allocate only if needed
	prop = (WilsonMatrix*)smalloc(GJP.VolNodeSites()*sizeof(WilsonMatrix));
	if (prop == 0) ERR.Pointer(cname, fname, "prop");
	VRB.Smalloc(cname, fname, "prop", prop,
				GJP.VolNodeSites() * sizeof(WilsonMatrix));
	}
  } else {
    if (midprop == NULL) { // Allocate only if needed
	  midprop = (WilsonMatrix*)smalloc(GJP.VolNodeSites()*sizeof(WilsonMatrix));
	  if (midprop == 0) ERR.Pointer(cname, fname, "midprop");
	  VRB.Smalloc(cname, fname, "midprop", midprop,
				  GJP.VolNodeSites() * sizeof(WilsonMatrix));
	}
  }
}
// Free space for propagators
void QPropW::Delete(int mid) {

  char *fname = "Delete(int)";
  VRB.Func(cname, fname);

  if (!mid) {
    if (prop != NULL) {
	  VRB.Sfree(cname, fname, "prop", prop);
	  sfree(prop);
	  prop = NULL;
	}
  } else {
	if (midprop != NULL) {
	  VRB.Sfree(cname, fname, "midprop", midprop);
	  sfree(midprop);
	  midprop = NULL;
	}
  }
}

// Constructor without and with QPropWArg
QPropW::QPropW(Lattice& lat, CommonArg* c_arg): Alg(lat, c_arg) {

  char *fname = "QPropW(L&, ComArg*)";
  cname = "QPropW";
  VRB.Func(cname, fname);
  
  prop = NULL;
  midprop = NULL;
}
QPropW::QPropW(Lattice& lat, QPropWArg* arg, CommonArg* c_arg): 
  Alg(lat, c_arg) {

  char *fname = "QPropW(L&, QPropWArg*, ComArg*)";
  cname = "QPropW";
  VRB.Func(cname, fname);

  qp_arg = *arg;
  prop = NULL;
  midprop = NULL;
}
// copy constructor
QPropW::QPropW(const QPropW& rhs):Alg(rhs),prop(NULL),midprop(NULL) {

  char *fname = "QPropW(const QPropW&)";
  cname = "QPropW";
  VRB.Func(cname, fname);
   
  Allocate(PROP);
  for (int i=0;i<GJP.VolNodeSites();i++)
    prop[i] = rhs.prop[i];
  
  if (rhs.StoreMidprop()) {
	Allocate(MIDPROP);
	for (int i=0;i<GJP.VolNodeSites();i++)
	  midprop[i] = rhs.midprop[i];
  }

  qp_arg = rhs.qp_arg;
}

// asignment operator
QPropW& QPropW::operator=(const QPropW& rhs) { 

  char *fname = "operator=(const QPropW& rhs)";
  VRB.Func(cname, fname);

  if (this != &rhs) {

    Allocate(PROP);
    
    for (int i=0;i<GJP.VolNodeSites();i++)
      prop[i] = rhs.prop[i];

    if (rhs.StoreMidprop()) {
	  Allocate(MIDPROP);
	  for (int i=0;i<GJP.VolNodeSites();i++)
	  midprop[i]=rhs.midprop[i];
	}
  
    qp_arg = rhs.qp_arg;
  }
  
  return *this;
}

// averaging constructor
QPropW::QPropW(QPropW& prop1, QPropW& prop2):Alg(prop1)
{
   cname = "QPropW"; 
   char *fname = "QPropW(QPropW&,QPropW&)";
   VRB.Func(cname, fname);

   prop = NULL;
   midprop = NULL;

   Allocate(PROP);
   for (int i=0; i<GJP.VolNodeSites(); i++)
     prop[i] = ((Float)0.5)*(prop1.prop[i]+prop2.prop[i]);

   qp_arg = prop1.qp_arg;
}

void QPropW::Run() {

   char *fname = "Run()";
   VRB.Func(cname, fname);

   // Set the node size of the full (non-checkerboarded) fermion field
   //----------------------------------------------------------------
   //int f_size = GJP.VolNodeSites() * Lat.FsiteSize()/GJP.SnodeSites();
   int iter;
   Float true_res;

   int Nspins = 4; // Number of spin components to be done
   
   // Flag set if sequential propagator 
   int seq_src = ((SrcType()==PROT_U_SEQ)||
				  (SrcType()==PROT_D_SEQ)||
				  (SrcType()==MESSEQ)      );

   if (DoHalfFermion()) Nspins = 2;

   // does prop exist? Assume it does not.
   int do_cg = 1;

   /****************************************************************
     The code below is temporarily isolated for purposes of merging
     with CPS main branch,  12/09/04, Oleg Loktik
   -------------------- Quarantine starts --------------------------

   if (pfs_file_exists(Arg.file)){
     // read data into prop
     RestoreQProp(Arg.file,PROP); // Only restores stuff into prop 

     // multiply source by 1/2(1+gamma_t). If the propagator  
     // on disk is half fermion it does nothing. Otherwise it gets
     // converted to the half fermion propagator. 
     if (Arg.DoHalfFermion) 
       for (int s(0);s<GJP.VolNodeSites();s++)
	 prop[s].PParProjectSink() ;
      
     do_cg = 0;
   }
   -------------------- End of quarantine -------------------------*/ 

   if (do_cg) {

     Allocate(PROP); 
     if (StoreMidprop()) Allocate(MIDPROP);
     
     FermionVectorTp src;
     FermionVectorTp sol;
     FermionVectorTp midsol;
     
     for (int spn=0; spn < Nspins; spn++)
       for (int col=0; col < GJP.Colors(); col++) {
		 
		 // initial guess (Zero)
		 sol.ZeroSource();
		 SetSource(src,spn,col);
		 
		 if ((DoHalfFermion())&&(!seq_src)) // Rotate to chiral basis
		   src.DiracToChiral();
		 
		 // Get the prop
		 printf("Before CG in QpropW.Run() \n");
		 CG(src, sol, midsol, iter, true_res);
		 //gauge fix solution
		 FixSol(sol);
		 if (StoreMidprop()) FixSol(midsol);
		 
		 // Collect solutions in propagator.
		 LoadRow(spn,col,sol,midsol);
		 
		 if (DoHalfFermion()) {// copy spin 0 to spin 1 and spin 2 to spin 3
		   int spn2 = spn + 2;
		   if (seq_src) {
			 LoadRow(spn2,col,sol,midsol);
		   } else { // Regular propagator zero the extra components
			 src.ZeroSource();
			 LoadRow(spn2,col,src,src);
		   }
		 }
		 
		 if (common_arg->results != 0) {
		   FILE *fp;
		   if ((fp = Fopen((char *)common_arg->results, "a")) == NULL) {
			 ERR.FileA(cname,fname, (char *)common_arg->results);
		   }
		   Fprintf(fp, "Cg iters = %d true residual = %e\n",
				   iter, (float)true_res);
		   Fclose(fp);
		 }
		 
       } // End spin-color loop
	 
     // Rotate the source indices to Chiral basis if needed
     if ((DoHalfFermion())&&(!seq_src)) {	
	   for (int s=0;s<GJP.VolNodeSites();s++)
		 prop[s].SinkChiralToDirac(); // multiply by V^\dagger
	   
	   if (StoreMidprop())
		 for (int s=0;s<GJP.VolNodeSites();s++)
		   midprop[s].SinkChiralToDirac(); // multiply by V^\dagger
	 }
   }
   // save prop
   if (do_cg && qp_arg.save_prop) {
     SaveQProp(qp_arg.file,PROP); 
   }
}


// Do conjugate gradient
void QPropW::CG(FermionVectorTp& source, FermionVectorTp& sol, 
		FermionVectorTp& midsol, int& iter, Float& true_res) {

  char *fname = "CG(source&, sol&, midsol&, int&, Float&)";
  VRB.Func(cname, fname);

  /*
  CgArg cg_arg;
  cg_arg.mass = 0.03;
  cg_arg.max_num_iter = 1000;
  cg_arg.stop_rsd = 1.0e-8;
  */

  Lattice& Lat = AlgLattice();

  // Set the node size of the full (non-checkerboarded) fermion field
  //----------------------------------------------------------------
  int ls = GJP.SnodeSites();
  int ls_glb = GJP.SnodeSites()*GJP.Snodes();
  int f_size = GJP.VolNodeSites() * Lat.FsiteSize()/GJP.SnodeSites();
  int f_size_5d = f_size * ls;

  // Do inversion
  //----------------------------------------------------------------
  if (Lat.Fclass() == F_CLASS_DWF) {
    Vector *src_4d    = (Vector*)source.data();
    Vector *sol_4d    = (Vector*)sol.data();
    Vector *midsol_4d = (Vector*)midsol.data();
    Vector *src_5d    = (Vector*)smalloc(f_size_5d * sizeof(IFloat));
    Vector *sol_5d    = (Vector*)smalloc(f_size_5d * sizeof(IFloat));
    if (src_5d == 0)
      ERR.Pointer(cname,fname, "src_5d");
    VRB.Smalloc(cname,fname, "src_5d", src_5d, f_size_5d * sizeof(IFloat));
    if (sol_5d == 0)
      ERR.Pointer(cname,fname, "sol_5d");
    VRB.Smalloc(cname,fname, "sol_5d", sol_5d, f_size_5d * sizeof(IFloat));

    //Get the lattice form the Alg base class
    Lattice& Lat = this->AlgLattice() ;

    Lat.Ffour2five(src_5d, src_4d, 0, ls_glb-1);
    Lat.Ffour2five(sol_5d, sol_4d, ls_glb-1, 0);

	iter = Lat.FmatInv(sol_5d, src_5d, &(qp_arg.cg), &true_res,
					   CNV_FRM_YES, PRESERVE_NO);
        
   /****************************************************************
     The code below is temporarily isolated for purposes of merging
     with CPS main branch,  12/09/04, Oleg Loktik 
   -------------------- Quarantine starts --------------------------

     iter = Lat.FmatInv4dSrc(sol_5d, src_5d, 0, ls_glb-1, &(Arg.cg), &true_res, 
    	                     CNV_FRM_YES, PRESERVE_NO);
   -------------------- End of quarantine -------------------------*/ 

    // prop on walls
    Lat.Ffive2four(sol_4d, sol_5d, ls_glb-1, 0);
    // midpoint prop
    if (StoreMidprop())
      Lat.Ffive2four(midsol_4d, sol_5d, ls_glb/2, ls_glb/2-1);

    VRB.Sfree(cname,fname, "sol_5d", sol_5d);
    sfree(sol_5d);
    VRB.Sfree(cname,fname, "src_5d", src_5d);
    sfree(src_5d);

  } else {
    iter = Lat.FmatInv((Vector*)sol.data(),(Vector*)source.data(),
					   &(qp_arg.cg), &true_res, CNV_FRM_YES, PRESERVE_NO);
  }

}

/*!
  gauge fix solution - only works for coulomb
  guage in time, or Landau gauge
*/
void QPropW::FixSol(FermionVectorTp& sol) {

   char *fname = "FixSol()";
   VRB.Func(cname, fname);

   // replace with check for FIX_GAUGE_NONE ??
   if (GFixedSnk()) {
	 Lattice& latt(AlgLattice());
	 switch ( latt.FixGaugeKind() ) {
	 case ( FIX_GAUGE_NONE ):
	   ERR.General(cname,fname,"Gauge fixing matrices not calculated\n");
	   break;
	 case ( FIX_GAUGE_COULOMB_T ):
	   // GaugeFixSink seems to be broken for anything
	   // except timeslice fixing ( dir isn't even used 
	   // in the function )
	   sol.GaugeFixSink(latt, 3);
	   break;
	   
	 case ( FIX_GAUGE_LANDAU ):
	   sol.LandauGaugeFixSink(latt);
	   break;
	   
	 default:
	   // should never be reached
	   ERR.General(cname,fname,"Unimplemented gauge fixing\n");
	   break;
	 }
   }
}

/*!
  Collect the solutions in the propagator
 */
void QPropW::LoadRow(int spin, int color, FermionVectorTp& sol, 
					 FermionVectorTp& midsol) {
  int i;
  for (int s=0; s<GJP.VolNodeSites(); s++) {
    i = s*SPINOR_SIZE; // FermionVector index
    prop[s].load_row(spin, color, (wilson_vector &)sol[i]);
  }
  if (StoreMidprop()) { // Collect solutions in midpoint propagator.
    for (int s=0; s<GJP.VolNodeSites(); s++) {
      i = s*SPINOR_SIZE; // lattice site
      midprop[s].load_row(spin, color, (wilson_vector &)midsol[i]);
    }
  }
} 

//sets the filename for IO
//void QPropW::SetFileName(char *nm)
//{
//  char *fname = "SetFileName(char*)";
//
//  if (Arg.file!=NULL) {
//    VRB.Sfree(cname, fname, "Arg.file", qp_arg.file);
//    sfree(Arg.file);
//  }
//
//  qp_arg.file = (char*) smalloc(100) ; // string is 100 characters long
//  if (Arg.file == 0) ERR.Pointer(cname, fname, "Arg.file");
//  VRB.Smalloc(cname, fname, "Arg.file", qp_arg.file, 100 );
//
//  strcpy(Arg.file,nm);
//}

// Needed by alg_threept. Could be replaced by the disk system.
void QPropW::ShiftPropForward(int n) {

   char *fname = "ShiftPropForward(int)";

   Float* recv_buf;
   Float* send_buf;
   int len =  12 * 12 * 2;
   // size of transfer in words
   recv_buf = (Float*)smalloc(len*sizeof(Float));
   if (recv_buf == 0) ERR.Pointer(cname,fname, "recv_buf");
   VRB.Smalloc(cname, fname, "recv_buf", recv_buf, 
			   12 * 12 * sizeof(Float));
 
   for (int j=0; j<n; j++) {
     // shift 1 node in t-dir.  prop -> prop
     for (int i=0; i<GJP.VolNodeSites(); i++) {
       send_buf = (Float*)&prop[i];
       getMinusData((IFloat*)recv_buf, (IFloat*)send_buf, len, 3);
       moveMem((IFloat *)&prop[i], (IFloat*)recv_buf, len*sizeof(IFloat) );
     }
   }

   VRB.Sfree(cname, fname, "recv_buf", recv_buf);
   sfree(recv_buf);

}
// Needed by alg_threept. Could be replaced by the disk system.
void QPropW::ShiftPropBackward(int n) {

   char *fname = "ShiftPropBack()";

   Float* recv_buf;
   Float* send_buf;
   int len = 12 * 12 * 2;
       // size of transfer in words
   recv_buf = (Float*)smalloc(len*sizeof(Float));
   if (recv_buf == 0) ERR.Pointer(cname,fname, "recv_buf");
   VRB.Smalloc(cname, fname, "recv_buf", recv_buf, 
	       12 * 12 * sizeof(Float));
 
   for (int j=0; j<n; j++) {
     // shift 1 node in t-dir.  prop -> prop
     for (int i=0; i<GJP.VolNodeSites(); i++) {
       send_buf = (Float*)&prop[i];
       getPlusData((IFloat*)recv_buf, (IFloat*)send_buf, len, 3);
       moveMem((IFloat*)&prop[i], (IFloat*)recv_buf, len*sizeof(IFloat));
     }
   }

   VRB.Sfree(cname, fname, "recv_buf", recv_buf);
   sfree(recv_buf);

}

/*!
  Compute the average of the current propagator and the propagaror Q
  \f[
  prop[i] =\frac{1}{2}(prop[i] + Q.prop[i])
  \f]
  It does this for both prop and midprop

  NOTE: This way of doing things saves an extra propagator in storage
  the averaging constructor should be removed since it is now obsolete.
 */
void QPropW::Average(QPropW& Q) {

  if ((Q.prop != NULL)&&(prop !=NULL))
    for (int i=0; i< GJP.VolNodeSites(); i++)
      prop[i] = ((Float)0.5)*(prop[i]+Q.prop[i]);

  if ((Q.midprop != NULL)&&(midprop !=NULL))
    for (int i=0; i< GJP.VolNodeSites(); i++)
      midprop[i] = ((Float)0.5)*(midprop[i]+Q.midprop[i]);
}

/*!
 Purpose:
   get a WilsonMatrix at specified coordinates (vec).
    can deal with vec being off node.

 Arguments:
\li   vec:     coordinates [x,y,z,t] of the WilsonMatrix we want to fetch. 
               These coordinates are relative to the [0,0,0,0]
               site of the node. They  could be out-of-range, i.e., 
               located off-node.
\li   tmp:     Buffer for the Matrix if communication is needed.
\li   return:  a reference to the Matrix. If off-node, it points to tmp.

   WARNING: It only works on prop and not on midprop
**/
WilsonMatrix& QPropW::GetMatrix(const int *vec, WilsonMatrix& tmp) const {

  // offset out-of-range coordinates site[] into on_node_site[]
  // in order to locate the Matrix
  //------------------------------------------------------------------------
  int on_node_site[4],site[4];
  int on_node = 1;
  WilsonMatrix *on_node_wmat;
  {
    for (int i=0; i<4; i++) {
      site[i] = on_node_site[i] = vec[i] ;
      while (on_node_site[i] < 0) {
        on_node_site[i] += GJP.NodeSites(i) ;
      }
      on_node_site[i] %= GJP.NodeSites(i) ;
      if (on_node_site[i] != site[i]) {
        on_node = 0;
      }
    }
    on_node_wmat = prop + (on_node_site[0] + GJP.XnodeSites()*(
						   on_node_site[1] + GJP.YnodeSites()*(
                           on_node_site[2] + GJP.ZnodeSites()*
                           on_node_site[3])));
  }

#ifndef PARALLEL
//VRB.FuncEnd(cname, fname) ;
  return *on_node_wmat;
#endif

  // send to the destination node if the site is off-node
  //------------------------------------------------------------------------
  if (on_node) {
	//VRB.FuncEnd(cname, fname);
    return *on_node_wmat;
  } else {
    WilsonMatrix send = *on_node_wmat;
    WilsonMatrix &recv = tmp;
    for (int i=0; i<4; i++) {
      while (site[i] != on_node_site[i]) {
        if (site[i] < 0) {
	  // the WilsonMatrix has 288 number of floats 
          getMinusData((IFloat*)&recv, (IFloat*)&send, sizeof(WilsonMatrix),i);
          on_node_site[i] -= GJP.NodeSites(i);
        } else {
	  // the WilsonMatrix has 288 number of floats 
          getPlusData ((IFloat*)&recv, (IFloat*)&send, sizeof(WilsonMatrix),i);
          on_node_site[i] += GJP.NodeSites(i);
        }
        send = recv;
      }
    }
//  VRB.FuncEnd(cname, fname) ;
    return recv ;
  }
}



void QPropW::SetSource(FermionVectorTp& src, int spin, int color) {

   char *fname = "SetSource()";
   VRB.Func(cname, fname);

   // Do nothing
   // This function should be overridden by specific QPropW types
   
}

Complex& QPropW::rand_src(int i) const { 
  //Do nothing ....
  ERR.General("QPropW","rand_src","No random source\n");
  // This is just to keep the compiler happy
  return *((Complex*)prop); 
}

WilsonMatrix QPropW::WallSinkProp(int t_sink) {

    WilsonMatrix wmat = (Float)0.0;

    for (int i=0; i<GJP.VolNodeSites(); i++) {
	  int t = i/(GJP.VolNodeSites()/GJP.TnodeSites());
	  t += GJP.TnodeCoor()*GJP.TnodeSites();
	  if (t != t_sink) continue;
	  wmat += prop[i];	  
    }
	
#ifdef PARALLEL
    slice_sum((Float*)&wmat, 288, 99);
#endif
	
    return wmat;
}


QPropW::~QPropW() {

  char *fname = "~QPropW()";
  VRB.Func(cname, fname);

  Delete(PROP);
  Delete(MIDPROP);
  //  if (Arg.file!=NULL) {
  //  VRB.Sfree(cname, fname, "Arg.file", qp_arg.file);
  //  sfree(Arg.file);
  // }
}

void QPropW::SaveQProp(char* name, int mid) {

  char *fname = "SaveQProp()";
  
  VRB.Func(cname, fname);

   /****************************************************************
     The code below is temporarily isolated for purposes of merging
     with CPS main branch,  12/09/04, Oleg Loktik 
   -------------------- Quarantine starts --------------------------

  if (! mid) {
    unsigned int* data;
    data = (unsigned int*)prop;
    save_data(name, data, GJP.VolNodeSites() * sizeof(WilsonMatrix));
    
    if (common_arg->results != 0) {
      FILE *fp;
      if ( (fp = Fopen((char *)common_arg->results, "a")) == NULL ) {
	ERR.FileA(cname,fname, (char *)common_arg->results);
      }
      Fprintf(fp, "Saved prop in file %s\n", name);
      Fclose(fp);
    }
  } else {
    unsigned int* data;
    data = (unsigned int*)midprop;
    save_data(name, data, GJP.VolNodeSites() * sizeof(WilsonMatrix));
    if (common_arg->results != 0) {
      FILE *fp;
      if ( (fp = Fopen((char *)common_arg->results, "a")) == NULL ) {
	ERR.FileA(cname,fname, (char *)common_arg->results);
      }
      Fprintf(fp, "Saved mid-point prop in file %s\n", name);
      Fclose(fp);
    }
  }
   -------------------- Quarantine ends ----------------------------*/
}

// Restore prop
void QPropW::RestoreQProp(char* name, int mid) {

  char *fname = "RestoreQProp()";
  VRB.Func(cname, fname);

   /****************************************************************
     The code below is temporarily isolated for purposes of merging
     with CPS main branch,  12/09/04, Oleg Loktik 
   -------------------- Quarantine starts --------------------------

  if (! mid) {

    Allocate(PROP) ;
    
    unsigned int* data;
    data = (unsigned int*)prop;
    read_data(name, data, GJP.VolNodeSites() * sizeof(WilsonMatrix),0);

    if (common_arg->results != 0) {
      FILE *fp;
      if ( (fp = Fopen((char *)common_arg->results, "a")) == NULL ) {
	ERR.FileA(cname,fname, (char *)common_arg->results);
      }
      Fprintf(fp, "Read prop from file %s\n", name);
      Fclose(fp);
    }
  } else {
    
    Allocate(MIDPROP) ;
    unsigned int* data;
    data = (unsigned int*)midprop;
    read_data(name, data, GJP.VolNodeSites() * sizeof(WilsonMatrix),0);

    if (common_arg->results != 0) {
      FILE *fp;
      if ( (fp = Fopen((char *)common_arg->results, "a")) == NULL ) {
	ERR.FileA(cname,fname, (char *)common_arg->results);
      }
      Fprintf(fp, "Read mid-point prop from file %s\n", name);
      Fclose(fp);
    }
  }
  -------------------- Quarantine ends ----------------------------*/

}


//HueyWen and Peter
QPropWGFLfuncSrc::QPropWGFLfuncSrc(Lattice &lat, 
				   CgArg *arg,
				   CommonArg *common_arg,
				   Float (*fn) (int,int,int,int)
				   ):QPropW(lat, common_arg)
{
   func = fn;

   char *fname = "QPropWGFLfuncSrc(L&, CgArg*, blah)";
   char *cname = "QPropWGFLfuncSrc";
   VRB.Func(cname, fname);

   // Set the node size of the full (non-checkerboarded) fermion field
   //----------------------------------------------------------------
   int f_size = GJP.VolNodeSites() * lat.FsiteSize()/GJP.SnodeSites();
   int iter;
   Float true_res;

   // allocate space for the quark propagator
   //----------------------------------------------------------------
   prop = (WilsonMatrix*) smalloc(GJP.VolNodeSites() * sizeof(WilsonMatrix));
   if (prop == 0) ERR.Pointer(cname, fname, "prop");
   VRB.Smalloc(cname, fname, "prop", prop, 12 * f_size * sizeof(Float));
   FermionVectorTp src;
   FermionVectorTp sol;

   for (int spin = 0; spin < 4; spin++)
     for (int color = 0; color < GJP.Colors(); color++){

        // initial guess
        sol.SetVolSource(color, spin);
        // set the source

        src.SetGFLfuncSource(lat, color, spin, fn);

        // Get the prop
        CG(lat, arg, src, sol, iter, true_res);
	
	// HueyWen 
	//sol.LandauGaugeFixSink(lat, 3);
	sol.LandauGaugeFixSink(lat);

        // Collect solutions in propagator.
        int j;
        for(int i=0; i<f_size; i+=SPINOR_SIZE){
           j=i/SPINOR_SIZE; // lattice site
           prop[j].load_vec(spin, color, (wilson_vector &)sol[i]);
        }

    	if(common_arg->results != 0){
      	  FILE *fp;
      	  if( (fp = Fopen((char *)common_arg->results, "a")) == NULL ) {
       	     ERR.FileA(cname,fname, (char *)common_arg->results);
      	  }
      	  Fprintf(fp, "Cg iters = %d true residual = %e\n",
			iter, true_res);
      	  Fclose(fp);
    	}


   } // End spin-color loop
}

QPropWGFLfuncSrc::~QPropWGFLfuncSrc()
{
  char *fname = "~QPropWGFLfuncSrc()";
  char *cname = "QPropWGFLfuncSrc";
  VRB.Func(cname, fname);
  VRB.Sfree(cname, fname, "prop", prop);
  sfree(prop);
}


//------------------------------------------------------------------
// Quark Propagator (Wilson type) with Wall Source
//------------------------------------------------------------------
QPropWWallSrc::QPropWWallSrc(Lattice& lat, CommonArg* c_arg):QPropW(lat, c_arg) {
 
  char *fname = "QPropWWallSrc(L&, ComArg*)";
  cname = "QPropWWallSrc";
  VRB.Func(cname, fname);
}
QPropWWallSrc::QPropWWallSrc(Lattice& lat,  QPropWArg* arg, CommonArg* c_arg):
  QPropW(lat, arg, c_arg) { 

  char *fname = "QPropWWallSrc(L&, QPropWArg*, ComArg*)";
  cname = "QPropWWallSrc";
  VRB.Func(cname, fname);

  // get the propagator
  Run();
}
QPropWWallSrc::QPropWWallSrc(QPropWWallSrc& prop1, QPropWWallSrc& prop2): 
  QPropW(prop1,prop2) {
  
  char *fname = "QPropWWallSrc(prop&, prop&)";
  cname = "QPropWWallSrc";
  VRB.Func(cname, fname);
}

//set wall source
void QPropWWallSrc::SetSource(FermionVectorTp& src, int spin, int color) {

  char *fname = "SetSource()";
  VRB.Func(cname, fname);
  
  src.ZeroSource() ;
  src.SetWallSource(color, spin, qp_arg.t);
  if (GFixedSrc())
    src.GFWallSource(AlgLattice(), spin, 3, qp_arg.t);
}

// === for exponetial smeared source ===================================

// exponential smeared source
QPropWExpSrc::QPropWExpSrc(Lattice& lat, CommonArg* c_arg):QPropW(lat, c_arg)
{
  char *fname = "QPropWExpSrc(L&, ComArg*)";
  cname = "QPropWExpSrc";

  VRB.Func(cname, fname);
}

QPropWExpSrc::QPropWExpSrc(Lattice& lat,  QPropWArg* arg, QPropWExpArg *e_arg,
  CommonArg* c_arg): QPropW(lat, arg, c_arg), exp_arg (*e_arg)
{
  char *fname = "QPropWExpSrc(L&, QPropWArg*, ComArg*)";
  cname = "QPropWExpSrc";

  VRB.Func(cname, fname);

  // get the propagator
  Run();
}


QPropWExpSrc::QPropWExpSrc(QPropWExpSrc* prop1, QPropWExpSrc* prop2) :
  QPropW(*prop1, *prop2)
{
  //char *fname = "QPropWExpSrc(prop*, prop*)";
  //cname = "QPropWExpSrc";
  //VRB.Func(cname, fname);

}

QPropWExpSrc::QPropWExpSrc(QPropWExpSrc& prop1, QPropWExpSrc& prop2) : 
QPropW(prop1,prop2)
{
  //char *fname = "QPropWExpSrc(prop&, prop&)";
  //cname = "QPropWExpSrc";
  //VRB.Func(cname, fname);

}

QPropWExpSrc::QPropWExpSrc(QPropW& prop1) : 
QPropW(prop1)
{
  //CJ: Have to set exp_arg: how?
  //char *fname = "QPropW(prop&)";
  //cname = "QPropWExpSrc";
  //VRB.Func(cname, fname);
  exp_arg.exp_A=1.2;
  exp_arg.exp_B=0.1;
  exp_arg.exp_C=8;
}

// Set exponential smeared source
void QPropWExpSrc::SetSource(FermionVectorTp& src, int spin, int color)
{
  char *fname = "SetSource()";
  VRB.Func(cname, fname);

  src.ZeroSource() ;
  src.SetExpSource(color, spin, qp_arg.x, qp_arg.y, qp_arg.z, qp_arg.t,
                   exp_arg.exp_A, exp_arg.exp_B, exp_arg.exp_C);
  if(qp_arg.gauge_fix_src)
    src.GFWallSource(AlgLattice(), spin, 3, qp_arg.t);
}

//------------------------------------------------------------------
// Quark Propagator (Wilson type) with Momentum Source
//------------------------------------------------------------------
QPropWMomSrc::QPropWMomSrc(Lattice& lat, CommonArg* c_arg):
  QPropWWallSrc(lat, c_arg) {
  
  char *fname = "QPropWMomSrc(L&, ComArg*)";
  cname = "QPropWMomSrc";
  VRB.Func(cname, fname);
}
QPropWMomSrc::QPropWMomSrc(Lattice& lat,  QPropWArg* arg, 
			   int* p, CommonArg* c_arg): 
  QPropWWallSrc(lat, c_arg), mom(p) {
 
  char *fname = "QPropWMomSrc(L&, QPropWArg*, ComArg*)";
  cname = "QPropWMomSrc";
  VRB.Func(cname, fname);

  Run();
}
// copy constructor
QPropWMomSrc::QPropWMomSrc(const QPropWMomSrc& rhs) : 
  QPropWWallSrc(rhs), mom(rhs.mom) {

  char *fname = "QPropW(const QPropW&)";
  cname = "QPropW";
  VRB.Func(cname, fname);
}

void QPropWMomSrc::SetSource(FermionVectorTp& src, int spin, int color) {
  char *fname = "SetSource()";
  VRB.Func(cname, fname);
  
  src.ZeroSource();
  src.SetMomSource(color, spin, qp_arg.t, mom);
  if (GFixedSrc())
    src.GFWallSource(AlgLattice(), spin, 3, qp_arg.t);
}

//------------------------------------------------------------------
// Quark Propagator (Wilson type) with Volume Source
//------------------------------------------------------------------
QPropWVolSrc::QPropWVolSrc(Lattice& lat, CommonArg* c_arg) : QPropW(lat, c_arg) {

  char *fname = "QPropWVolSrc(L&, ComArg*)";
  cname = "QPropWVolSrc";
  VRB.Func(cname, fname);
}
QPropWVolSrc::QPropWVolSrc(Lattice& lat,  QPropWArg* arg, CommonArg* c_arg) : 
  QPropW(lat, arg, c_arg) {
 
  char *fname = "QPropWVolSrc(L&, QPropWArg*, ComArg*)";
  cname = "QPropWVolSrc";
  VRB.Func(cname, fname);

  Run();
}

void QPropWVolSrc::SetSource(FermionVectorTp& src, int spin, int color) {

  char *fname = "SetSource()";
  VRB.Func(cname, fname);
  
  src.SetVolSource(color, spin);
  if (GFixedSrc())
    for (int t=0; t<GJP.Tnodes()*GJP.TnodeSites(); t++)
      src.GFWallSource(AlgLattice(), spin, 3, t);
}

//------------------------------------------------------------------
// Quark Propagator (Wilson type) with Point Source
//------------------------------------------------------------------
QPropWPointSrc::QPropWPointSrc(Lattice& lat, CommonArg* c_arg) : 
  QPropW(lat, c_arg) {
 
  char *fname = "QPropWPointSrc(L&, ComArg*)";
  cname = "QPropWPointSrc";
  VRB.Func(cname, fname);
}
QPropWPointSrc::QPropWPointSrc(Lattice& lat,  QPropWArg* arg,
							   CommonArg* c_arg) : 
  QPropW(lat, arg, c_arg) {
 
  char *fname = "QPropWPointSrc(L&, ComArg*)";
  cname = "QPropWPointSrc";
  VRB.Func(cname, fname);

  Run();
}

void QPropWPointSrc::SetSource(FermionVectorTp& src, int spin, int color) {

  char *fname = "SetSource()";
  VRB.Func(cname, fname);

  src.SetPointSource(color,spin,qp_arg.x,qp_arg.y,qp_arg.z,qp_arg.t);
}

//------------------------------------------------------------------
// Quark Propagator (Wilson type) with Box Source
//------------------------------------------------------------------
QPropWBoxSrc::QPropWBoxSrc(Lattice& lat, CommonArg* c_arg) :
  QPropW(lat, c_arg) { 

  char *fname = "QPropWBoxSrc(L&, ComArg*)";
  cname = "QPropWBoxSrc";
  VRB.Func(cname, fname);
}
QPropWBoxSrc::QPropWBoxSrc(Lattice& lat,  QPropWArg* arg, CommonArg* c_arg) : 
  QPropW(lat, arg, c_arg) {
 
  char *fname = "QPropWBoxSrc(L&, QPropWArg*, ComArg*)";
  cname = "QPropWBoxSrc";
  VRB.Func(cname, fname);

  Run();
}

void QPropWBoxSrc::SetSource(FermionVectorTp& src, int spin, int color) {

  char *fname = "SetSource()";
  VRB.Func(cname, fname);
  
  src.SetBoxSource(color, spin, box_arg.box_start, box_arg.box_end, qp_arg.t );
  if (GFixedSrc()) 
    src.GFWallSource(AlgLattice(), spin, 3, qp_arg.t);
}

//------------------------------------------------------------------
// Quark Propagator (Wilson type) with Random Source
//------------------------------------------------------------------
QPropWRand::QPropWRand(Lattice& lat, CommonArg* c_arg) : QPropW(lat, c_arg) {

  char *fname = "QPropWRand(L&, ComArg*)";
  cname = "QPropWRand";
  VRB.Func(cname, fname);

  rsrc = NULL;
}

// Allocate and Delete space for the random source
void QPropWRand::AllocateRsrc() {

  char *fname = "AllocateRsrc()";
  VRB.Func(cname, fname);

  if (rsrc == NULL) {
	int rsrc_size = 2 * GJP.VolNodeSites();
	rsrc = (Float*)smalloc(rsrc_size * sizeof(Float));
	if (rsrc == 0) ERR.Pointer(cname, fname, "rsrc");
	VRB.Smalloc(cname, fname, "rsrc", rsrc, rsrc_size * sizeof(Float));
  }
}
void QPropWRand::DeleteRsrc() {

  char *fname = "DeleteRsrc()";
  VRB.Func(cname, fname);

  if (rsrc != NULL) {
	VRB.Sfree(cname, fname, "rsrc", rsrc);
	sfree(rsrc);
	rsrc = NULL;
  }
}

QPropWRand::QPropWRand(Lattice& lat,  QPropWArg* arg, QPropWRandArg *r_arg, 
CommonArg* c_arg) : QPropW(lat, arg, c_arg),rand_arg(*r_arg) 
{
  char *fname = "QPropWRand(L&, QPropWArg*, ComArg*)";
  cname = "QPropWRand";
  VRB.Func(cname, fname);

  rsrc = NULL;
  AllocateRsrc();

  int rsrc_size = 2*GJP.VolNodeSites();

  if (rand_arg.rng == GAUSS) {
	GaussianRandomGenerator rng(0.5);
	rng.Reset(rand_arg.seed);
	for (int i=0; i<rsrc_size; i++)
	  rsrc[i] = rng.Rand();
  }
  if (rand_arg.rng == UONE) {
	UniformRandomGenerator rng(0,6.2831853);
	rng.Reset(rand_arg.seed);
	for (int i=0; i<rsrc_size/2; i++) {
	  Float theta(rng.Rand());
	  rsrc[2*i  ] = cos(theta); // real part
	  rsrc[2*i+1] = sin(theta); // imaginary part
	}
  }
  if (rand_arg.rng == ZTWO) {
	UniformRandomGenerator rng(0.0,1.0);
	rng.Reset(rand_arg.seed);
	for (int i=0; i<rsrc_size/2; i++) {
	  if (rng.Rand()>.5) {
	    rsrc[2*i] =  1.0;
	  } else {
	    rsrc[2*i] = -1.0;
	  }
	  rsrc[2*i+1] = 0.0; // source is purely real
	}
  }
  
}

QPropWRand::QPropWRand(const QPropWRand& rhs) : QPropW(rhs),rsrc(NULL) {

  char *fname = "QPropW(const QPropW&)";
  cname = "QPropW";
  VRB.Func(cname, fname);

  AllocateRsrc();
  for (int i=0; i<2*GJP.VolNodeSites(); i++)
    rsrc[i] = rhs.rsrc[i];
}

Complex& QPropWRand::rand_src(int i) const   
{ return ((Complex*)rsrc)[i]; }

QPropWRand& QPropWRand::operator=(const QPropWRand& rhs) {

  char *fname = "operator=(const QPropWRand& rhs)";
  VRB.Func(cname, fname);

  if (this != &rhs) {
	QPropW::operator=(rhs) ; // This copies the QPropW stuff...
      
	AllocateRsrc();    
	for (int i=0; i<2*GJP.VolNodeSites(); i++)
	  rsrc[i] = rhs.rsrc[i];
  }
  
  return *this;
}

QPropWRand::~QPropWRand() {
  char *fname = "~QPropWRand()";
  VRB.Func(cname, fname);

  DeleteRsrc();
}

void QPropWRand::ShiftPropForward(int n) {

  char *fname = "ShiftPropForward()";
  VRB.Func(cname, fname);

  QPropW::ShiftPropForward(n) ;

  Float* recv_buf;
  Float* send_buf;
  int len = 12 * 12 * 2;
  int len2 = 4;
  // size of transfers in words
  recv_buf = (Float*)smalloc(len*sizeof(Float));
  if (recv_buf == 0) ERR.Pointer(cname,fname, "recv_buf");
  VRB.Smalloc(cname, fname, "recv_buf", recv_buf,
			  len * sizeof(Float));
  
  for (int j=0; j<n; j++) {
    // shift 1 node in t-dir.  prop -> prop
    for (int i=0; i < 2 * GJP.VolNodeSites(); i+=len2) {
      send_buf = &rsrc[i];
      getMinusData((IFloat*)recv_buf, (IFloat*)send_buf, len2, 3);
      moveMem((IFloat*)&rsrc[i], (IFloat*)recv_buf, len2*sizeof(IFloat));
    }
    
  }
  
  VRB.Sfree(cname, fname, "recv_buf", recv_buf);
  sfree(recv_buf);
}
void QPropWRand::ShiftPropBackward(int n) {

  char *fname = "ShiftPropBackward()";
  VRB.Func(cname, fname);

  QPropW::ShiftPropBackward(n) ;

  Float* recv_buf;
  Float* send_buf;
  int len = 12 * 12 * 2;
  int len2 = 4;
  // size of transfer in words
  recv_buf = (Float*)smalloc(len * sizeof(Float));
  if (recv_buf == 0) ERR.Pointer(cname,fname, "recv_buf");
  VRB.Smalloc(cname, fname, "recv_buf", recv_buf,
	      len * sizeof(Float));
  
  for (int j=0; j<n; j++) {
    // shift 1 node in t-dir.  prop -> prop
    for (int i=0; i < 2 * GJP.VolNodeSites(); i+=len2) {
      send_buf = &rsrc[i];
      getPlusData((IFloat*)recv_buf, (IFloat*)send_buf, len2, 3);
      moveMem((IFloat*)&rsrc[i], (IFloat*)recv_buf, len2*sizeof(IFloat));
    }
    
  }
  
  VRB.Sfree(cname, fname, "recv_buf", recv_buf);
  sfree(recv_buf);
}

// Restore prop
void QPropWRand::RestoreQProp(char* name, int mid) {

   char *fname = "RestoreQProp()";
   VRB.Func(cname, fname);
   /****************************************************************
     The code below is temporarily isolated for purposes of merging
     with CPS main branch,  12/09/04, Oleg Loktik 
   -------------------- Quarantine starts --------------------------

   
   QPropW::RestoreQProp(name,mid) ;

   AllocateRsrc();
   
   unsigned int* data;
   data = (unsigned int*)rsrc;
   read_data(name, data, 
	     GJP.VolNodeSites() * sizeof(Complex),
	     GJP.VolNodeSites() * sizeof(WilsonMatrix));
   printf("Read rsrc from file %s\n", name);
  -------------------- Quarantine ends ---------------------------*/
}
// Save prop
void QPropWRand::SaveQProp(char* name, int mid) {
  char *fname = "SaveQProp()";
  VRB.Func(cname, fname);
   /****************************************************************
     The code below is temporarily isolated for purposes of merging
     with CPS main branch,  12/09/04, Oleg Loktik 
   -------------------- Quarantine starts --------------------------
  

  QPropW::SaveQProp(name,mid) ;
  
  unsigned int* data;
  data = (unsigned int*)rsrc;
  append_data(name, data, GJP.VolNodeSites() * sizeof(Complex));
  printf("Appended rsrc to file %s\n", name);
  DeleteRsrc();
  -------------------- Quarantine ends ---------------------------*/
}

//------------------------------------------------------------------
// Quark Propagator (Wilson type) with Random Wall Source
//------------------------------------------------------------------
QPropWRandWallSrc::QPropWRandWallSrc(Lattice& lat, CommonArg* c_arg) : 
  QPropWRand(lat, c_arg) {
  
  char *fname = "QPropWRandWallSrc(L&, ComArg*)";
  cname = "QPropWRandWallSrc";
  VRB.Func(cname, fname);
}
QPropWRandWallSrc::QPropWRandWallSrc(Lattice& lat,  QPropWArg* arg, 
  QPropWRandArg *r_arg, CommonArg* c_arg) 
  : QPropWRand(lat, arg, r_arg, c_arg) {
 
  char *fname = "QPropWRandWallSrc(L&, ComArg*)";
  cname = "QPropWRandWallSrc";
  VRB.Func(cname, fname);

  Run();
}

void QPropWRandWallSrc::SetSource(FermionVectorTp& src, int spin, int color) {

  char *fname = "SetSource()"; 
  VRB.Func(cname, fname);

  if (rsrc==NULL) {//need random numbers may implemented later
    ERR.General(cname,fname,"No randrom numbers found!\n") ;
  }

  src.ZeroSource();
  src.SetWallSource(color, spin, qp_arg.t, rsrc);
  if (GFixedSrc())
    src.GFWallSource(AlgLattice(), spin, 3, qp_arg.t);
}


//------------------------------------------------------------------
// Quark Propagator (Wilson type) with Random Volume Source
//------------------------------------------------------------------
QPropWRandVolSrc::QPropWRandVolSrc(Lattice& lat, CommonArg* c_arg) : 
  QPropWRand(lat, c_arg) {

  char *fname = "QPropWRandVolSrc(L&, ComArg*)";
  cname = "QPropWRandVolSrc";
  VRB.Func(cname, fname);
}
QPropWRandVolSrc::QPropWRandVolSrc(Lattice& lat,  QPropWArg* arg, 
  QPropWRandArg *r_arg, CommonArg* c_arg) 
  : QPropWRand(lat, arg, r_arg, c_arg) {
 
  char *fname = "QPropWRandVolSrc(L&, ComArg*)";
  cname = "QPropWRandVolSrc";
  VRB.Func(cname, fname);

  Run();
}

//set the random source
void QPropWRandVolSrc::SetSource(FermionVectorTp& src, int spin, int color) {

  char *fname = "SetSource()"; 
  VRB.Func(cname, fname);

  if (rsrc==NULL) {//need random numbers may implemented later
    ERR.General(cname,fname,"No randrom numbers found!\n") ;
  }
  
  src.SetVolSource(color, spin, rsrc);
  if (GFixedSrc()) 
    for (int t=0;t<GJP.Tnodes()*GJP.TnodeSites(); t++)
      src.GFWallSource(AlgLattice(), spin, 3, t);
}


//------------------------------------------------------------------
// Quark Propagator (Wilson type) with Random Slab Source
//------------------------------------------------------------------
QPropWRandSlabSrc::QPropWRandSlabSrc(Lattice& lat, CommonArg* c_arg) : 
  QPropWRand(lat, c_arg) {

  char *fname = "QPropWRandSlabSrc(L&, ComArg*)";
  cname = "QPropWRandSlabSrc";
  VRB.Func(cname, fname);
}

QPropWRandSlabSrc::QPropWRandSlabSrc(Lattice& lat,  QPropWArg* arg, 
  QPropWSlabArg *s_arg, CommonArg* c_arg) 
  : QPropWRand(lat, arg, &(s_arg->rand_arg), c_arg),slab_width(s_arg->slab_width) {
 
  char *fname = "QPropWRandSlabSrc(L&, ComArg*)";
  cname = "QPropWRandSlabSrc";
  VRB.Func(cname, fname);

  Run();
}

QPropWRandSlabSrc::QPropWRandSlabSrc(Lattice &lat, QPropWArg* arg, 
				     Float* src, CommonArg* c_arg): 
  QPropWRand(lat, c_arg) {

  char *fname = "QPropWRandSlabSrc(L&, QPropWArg*, Float*, CommonArg*)";
  cname = "QPropWRandSlabSrc";
  VRB.Func(cname, fname);
   
  for (int i=0; i<2*GJP.VolNodeSites(); i++)
    rsrc[i] = src[i];

  Run();
}

void QPropWRandSlabSrc::SetSource(FermionVectorTp& src, int spin, int color) {

  char *fname = "SetSource()"; 
  VRB.Func(cname, fname);

  if (rsrc==NULL) {//need random numbers may implemented later
    ERR.General(cname,fname,"No randrom numbers found!\n") ;
  }
  src.ZeroSource();
  
  for (int t=qp_arg.t; t < qp_arg.t+slab_width; t++) {
    src.SetWallSource(color, spin, t, rsrc);
    if (GFixedSrc()) 
      src.GFWallSource(AlgLattice(), spin, 3, t);
  }
}

//------------------------------------------------------------------
// Quark Propagator (Wilson type) with Sequential Sources
//------------------------------------------------------------------
QPropWSeq::QPropWSeq(Lattice& lat, QPropW& q, int *p, 
		     QPropWArg* q_arg, CommonArg* c_arg) : 
  QPropW(lat, q_arg, c_arg), quark(q), mom(p) {

  char *fname = "QPropWSeq(L&, ComArg*)";
  cname = "QPropWSeq";
  VRB.Func(cname, fname);
  
  // Stores the quark mass of the source propagator
  // Needed by Yasumichi's not degenerate mass runs
  quark_mass = quark.Mass() ; 
  
  //if the QPropW used for constructing the sequential source
  //propagator is done using HalfFerion set the DoHalfFermion 
  //in case the user forgot to do so.
  if (q.DoHalfFermion()) qp_arg.do_half_fermion = 1;
}

QPropWSeqMesSrc::QPropWSeqMesSrc(Lattice& lat, QPropW& quark,  int *p, 
				 int g, QPropWArg* q_arg, CommonArg* c_arg) : 
  QPropWSeq(lat, quark, p, q_arg, c_arg),gamma(g) {

  char *fname = "QPropWSeq(L&,...)";
  cname = "QPropWSeq";
  VRB.Func(cname, fname);

  Run();
}

void QPropWSeqMesSrc::SetSource(FermionVectorTp& src, int spin, int color) {

  char *fname = "SetSource()";
  cname = "QPropWSeqMesSrc";
  VRB.Func(cname, fname);

  if (color < 0 || color >= GJP.Colors())
    ERR.General(cname, fname,
				"Color index out of range: color = %d\n", color);
  
  if (spin < 0 || spin > 3)
    ERR.General(cname, fname,
				"Spin index out of range: spin = %d\n", spin);

  src.ZeroSource();
  Site s;
  WilsonMatrix tmp;
  for (s.Begin();s.End();s.nextSite())
    if (qp_arg.t == s.physT()) {
      tmp = quark[s.Index()];
      tmp.gl(gamma); // Multiply by the meson operator
      // multiply by the mommentum factor exp(ipx)
      Complex tt(conj(mom.Fact(s)));
      tmp*=tt;
      src.CopyWilsonMatSink(s.Index(),spin,color,tmp);
    }

  // Gauge fix the source. If QPropW sink is gauge fixed
  // this has to be done!
  if (quark.GFixedSnk()) {
	for (int ss=0;ss<4;ss++)
	  src.GFWallSource(AlgLattice(), ss, 3, qp_arg.t);
  }
}

QPropWSeqBar::QPropWSeqBar(Lattice& lat, QPropW& quark,  int *p, 
						   ProjectType pp, QPropWArg* q_arg, 
						   CommonArg* c_arg) : 
  QPropWSeq(lat, quark, p, q_arg, c_arg), proj(pp) {

  char *fname = "QPropWSeqBar(L&,...)";
  cname = "QPropWSeq";
  VRB.Func(cname, fname);
}

QPropWSeqProtDSrc::QPropWSeqProtDSrc(Lattice& lat, QPropW& quark,  int *p, 
									 ProjectType pp, QPropWArg* q_arg,
									 CommonArg* c_arg):
  QPropWSeqBar(lat, quark, p, pp, q_arg, c_arg) {

  char *fname = "QPropWSeqProtDSrc(L&,...)";
  cname = "QPropWSeq";
  VRB.Func(cname, fname);

  Run();

  //Multiply by gamma5 and take the dagger to make it in to quark.
  Site s;
  for (s.Begin();s.End();s.nextSite()) {
	QPropW::operator[](s.Index()).gl(-5);
	QPropW::operator[](s.Index()).hconj();
  }
}

void QPropWSeqProtDSrc::SetSource(FermionVectorTp& src, int spin, int color) {

  char *fname = "SetSource()";  
  VRB.Func(cname, fname);

  if (color < 0 || color >= GJP.Colors())
    ERR.General(cname, fname,
				"Color index out of range: color = %d\n", color);
  
  if (spin < 0 || spin > 3)
    ERR.General(cname, fname,
				"Spin index out of range: spin = %d\n", spin);
  
  src.ZeroSource();
  Site s;
  Diquark diq;
  WilsonVector S;

  WilsonMatrix q;
  WilsonMatrix OqO;
  WilsonMatrix qO ;
  WilsonMatrix Oq ;
  for (s.Begin();s.End();s.nextSite())
    if (qp_arg.t == s.physT()) {
	  int i = s.Index();
	  q = quark[i];
	  // If DoHalfFermion is on we have non-relativistic sources
	  // Multiply the sink by 1/2(1+gamma_t) when we do Half Spinors
	  // By doing this we implement the non-relativistic sink
	  if (DoHalfFermion()) q.PParProjectSink(); 
	  Oq = qO = q;
	  // multiply C*gamma_5 left  C is the charge conjugation
	  Oq.ccl(5);
	  // multiply C*gamma_5 right C is the charge conjugation
	  qO.ccr(5);
	  OqO = Oq;
	  // multiply C*gamma_5 right C is the charge conjugation
	  // Shoichi's code misses a minus sign here (or in ccl)
	  OqO.ccr(5);
	  // spin is denoted as delta in notes
	  // color is denoted as d in notes
	  diq.D_diquark(OqO, q, Oq, qO, spin, color);
	  diq.Project(S,proj);
	  
	  //multiply by the momentum factor exp(ipx)
	  Complex tt(conj(mom.Fact(s)));
	  S*=tt;

	  S.conj();

	  //if non-relativistic sink is needed multiply  by 1/2(1+gamma_t)
	  if (DoHalfFermion()) S.PParProject(); 
	  //Note the order first project then multiply by gamma5 
	  //multily by gamma5
	  S.gamma(-5);	
	  
	  src.CopyWilsonVec(i,S);
	}
  
  // Gauge fix the source. If quark sink is gauge fixed
  // this has to be done!
  if (quark.GFixedSnk()) {
	for (int ss=0;ss<4;ss++)
	  src.GFWallSource(AlgLattice(), ss, 3, qp_arg.t);
  }
}

QPropWSeqProtUSrc::QPropWSeqProtUSrc(Lattice& lat, QPropW& quark,  int *p, 
									 ProjectType pp, QPropWArg* q_arg,
									 CommonArg* c_arg):
  QPropWSeqBar(lat, quark, p, pp, q_arg, c_arg) {

  char *fname = "QPropWSeqProtUSrc(L&,...)";
  cname = "QPropWSeq";
  VRB.Func(cname, fname);

  Run();

  //Multiply by gamma5 and take the dagger to make it in to quark.
  Site s;
  for (s.Begin();s.End();s.nextSite()) {
	QPropW::operator[](s.Index()).gl(-5);
	QPropW::operator[](s.Index()).hconj(); 
  }
}

void QPropWSeqProtUSrc::SetSource(FermionVectorTp& src, int spin, int color) {

  char *fname = "SetSource()";  
  VRB.Func(cname, fname);

  if (color < 0 || color >= GJP.Colors())
    ERR.General(cname, fname,
				"Color index out of range: color = %d\n", color);
  
  if (spin < 0 || spin > 3)
    ERR.General(cname, fname,
				"Spin index out of range: spin = %d\n", spin);
  
  src.ZeroSource();
  Site s;
  Diquark diq;
  WilsonVector S;
  WilsonMatrix q;
  WilsonMatrix OqO;

  for (s.Begin();s.End();s.nextSite())
    if (qp_arg.t == s.physT()) {
	  int i(s.Index());
	  q = quark[i];
	  // If DoHalfFermion is on we have non-relativistic sources
	  // Multiply the sink by 1/2(1+gamma_t) when we do Half Spinors
	  // By doing this we implement the non-relativistic sink
	  if (DoHalfFermion()) q.PParProjectSink();
	  OqO = q;
	  
	  // multiply C*gamma_5 left  C is the charge conjugation
	  OqO.ccl(5);
	  // multiply C*gamma_5 right C is the charge conjugation
	  OqO.ccr(5);
	  // spin is denoted as delta in notes
	  // color is denoted as d in notes
	  
	  diq.U_diquark(OqO, q, spin, color);
	  diq.Project(S,proj);
	  
	  //multiply by the momentum factor exp(ipx)
	  Complex tt(conj(mom.Fact(s)));
	  S*=tt;
	  
	  S.conj();
	  
	  //if non-relativistic sink is needed multiply  by 1/2(1+gamma_t)
	  if (DoHalfFermion()) S.PParProject(); 
	  //Note the order first project then multiply by gamma5 
	  //multily by gamma5
	  S.gamma(-5);
	  
	  src.CopyWilsonVec(i,S);
	}// Loop over sites
  
  // Gauge fix the source. If quark sink is gauge fixed
  // this has to be done!
  if (quark.GFixedSnk()) {
	for (int ss=0;ss<4;ss++)
	  src.GFWallSource(AlgLattice(), ss, 3, qp_arg.t);
  }
}

CPS_END_NAMESPACE

