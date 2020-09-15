// -*- mode:c++; c-basic-offset:2 -*-
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

#include <config.h>
#include <stdlib.h>		// exit()
#include <stdio.h>
#include <string.h>
#include <cassert>
#include <alg/common_arg.h>
#include <comms/glb.h>
#include <util/vector.h>
#include <comms/scu.h>
#include <comms/sysfunc_cps.h>

#include <fcntl.h>		// read and write control flags,
#include <unistd.h>		// close(). These are needed for io parts to
			// compile on PCs
#include <cassert>		// for assert()

#include <alg/qpropw.h>
#include <util/qcdio.h>

#include <util/qioarg.h>
#include <util/sproj_tr.h>
#include <util/time_cps.h>
#include <util/omp_wrapper.h>

//YA
#include <alg/alg_plaq.h>
#include <alg/alg_smear.h>
#include <alg/no_arg.h>

#ifdef USE_QMP
#include <qmp.h>
#endif
#include <alg/mobius_arg.h>
//CK
#include <util/fpconv.h>
#include <util/checksum.h>

#ifdef USE_BFM

//CK: these are redefined by BFM (to the same values)
#undef ND
#undef SPINOR_SIZE
#undef HALF_SPINOR_SIZE
#undef GAUGE_SIZE

#include <util/lattice/fbfm.h>
#endif

#define VOLFMT QIO_VOLFMT

CPS_START_NAMESPACE
// Allocate space for propagators
void QPropW::Allocate (int mid)
{
  const char *fname = "Allocate(int)";
  VRB.Func (cname, fname);

  int sz =
    (GJP.Gparity ()? 2 : 1) * GJP.VolNodeSites () * sizeof (WilsonMatrix);

  switch (mid) {
  case 0:
    if (prop == NULL) {		// Allocate only if needed
      prop = (WilsonMatrix *) smalloc (cname, fname, "prop", sz);
    }
    break;
  case 1:
    if (midprop == NULL) {	// Allocate only if needed
      midprop = (WilsonMatrix *) smalloc (cname, fname, "midprop", sz);
      VRB.Result (cname, fname, "midprop=%p\n", midprop);
    }
    break;
  case 2:
    if (prop5d == NULL) {	// Allocate only if needed
      int ls = GJP.SnodeSites ();
      prop5d =
	(WilsonMatrix *) smalloc (cname, fname, "prop5d",
				  ls * GJP.VolNodeSites () *
				  sizeof (WilsonMatrix));
    }
    break;
  default:
    ERR.General (cname, fname, "Bad prop index in Allocate()\n");
    break;
  }

}

// Free space for propagators
void QPropW::Delete (int mid)
{

  const char *fname = "Delete(int)";
  VRB.Func (cname, fname);
  switch (mid) {

  case 0:
    if (prop != NULL) {
      VRB.Sfree (cname, fname, "prop", prop);
      sfree (prop);
      prop = NULL;
    }
    break;
  case 1:
    if (midprop != NULL) {
      VRB.Sfree (cname, fname, "midprop", midprop);
      sfree (midprop);
      midprop = NULL;
    }
    break;
  case 2:
    if (prop5d != NULL) {
      VRB.Sfree (cname, fname, "prop5d", prop5d);
      sfree (prop5d);
      prop5d = NULL;
    }
    break;
  default:
    ERR.General (cname, fname, "Bad prop index in Allocate()\n");
    break;
  }

}

// Constructor without and with QPropWArg
QPropW::QPropW (Lattice & lat, CommonArg * c_arg)
:Alg (lat, c_arg)
{

  char *fname = "QPropW(L&, ComArg*)";
  cname = "QPropW";
  VRB.Func (cname, fname);

  prop = NULL;
  midprop = NULL;
  prop5d = NULL;
  // EES
  propls = NULL;

  // YA
  lat_back = NULL;
  link_status_smeared = false;
  sink_type = POINT;
}

QPropW::QPropW (Lattice & lat, QPropWArg * arg, CommonArg * c_arg)
:Alg (lat, c_arg)
{
  char *fname = "QPropW(L&, QPropWArg*, ComArg*)";
  cname = "QPropW";
  VRB.Func (cname, fname);

  qp_arg = *arg;
  prop = NULL;
  midprop = NULL;
  propls = NULL;
  prop5d = NULL;

  VRB.Result(cname,fname,"mass=%g\n",arg->cg.mass);
  lat.SetMassArg(arg->cg.mass);


  // YA
  lat_back = NULL;
  link_status_smeared = false;
  sink_type = POINT;

  //-----------------------------------------------------------------
  // TY Add Start
  if (qp_arg.save_ls_prop == 1) {
    char sname[100];
    for (int nls (0); nls < GJP.SnodeSites (); nls++) {
#if TARGET == QCDOC
      sprintf (sname, "/pfs/%s.m%0.3f.l%d.id%d.dat",
	       qp_arg.file, qp_arg.cg.mass, nls, UniqueID ());
#else
      sprintf (sname, "pfs/%s.m%0.3f.l%d.id%d.dat",
	       qp_arg.file, qp_arg.cg.mass, nls, UniqueID ());
#endif
      FILE *fp;
      if ((fp = fopen (sname, "w")) != NULL) {
	VRB.Flow (cname, fname, "Remove file %s\n", sname);
      } else {
	ERR.FileA (cname, fname, sname);
      }
      fclose (fp);
    }
  }
  if (qp_arg.save_ls_prop == 2) {
    VRB.Debug (cname, fname, "Before allocate porpls\n");
    if (propls == NULL) {
      propls =
	(WilsonMatrix *) smalloc (cname, fname, "propls",
				  GJP.VolNodeSites () * GJP.SnodeSites () *
				  sizeof (WilsonMatrix));
      VRB.Debug (cname, fname, "Allocate porpls\n");
    }
  }
  // TY Add End
  //-----------------------------------------------------------------
}


// read in a prop from a file
QPropWRead::QPropWRead (Lattice & lat, QPropWArg * arg, CommonArg * c_arg,
			SourceType type)
:QPropW (lat, arg, c_arg)
{
  char *fname = "QPropWRead(L&, QPropWArg*, ComArg*, SourceType)";
  cname = "QPropWRead";
  VRB.Func (cname, fname);

  src_type = type;
  qp_arg = *arg;

  Allocate (0);			//space for 4d prop
  RestoreQProp (qp_arg.file, 0);
}



// copy constructor
QPropW::QPropW (const QPropW & rhs):Alg (rhs), midprop (NULL), prop (NULL),
prop5d (NULL)
{

  char *fname = "QPropW(const QPropW&)";
  cname = "QPropW";
  VRB.Func (cname, fname);

  Allocate (PROP);

  int sz = GJP.VolNodeSites ();
  if (GJP.Gparity ())
    sz *= 2;

  for (int i = 0; i < sz; i++)
    prop[i] = rhs.prop[i];

  if (rhs.StoreMidprop ()) {
    Allocate (MIDPROP);
    for (int i = 0; i < sz; i++)
      midprop[i] = rhs.midprop[i];
  }
  // YA
  lat_back = NULL;
  link_status_smeared = false;
  sink_type = rhs.sink_type;

  //-----------------------------------------------------------------
  // TY Add Start
  propls = rhs.propls;
  // TY Add End
  //-----------------------------------------------------------------
  qp_arg = rhs.qp_arg;
}

// asignment operator
QPropW & QPropW::operator= (const QPropW & rhs)
{

  char *fname = "operator=(const QPropW& rhs)";
  VRB.Func (cname, fname);

  if (this != &rhs) {

    Allocate (PROP);

    int sz = GJP.VolNodeSites ();
    if (GJP.Gparity ())
      sz *= 2;

    for (int i = 0; i < sz; i++)
      prop[i] = rhs.prop[i];

    if (rhs.StoreMidprop ()) {
      Allocate (MIDPROP);
      for (int i = 0; i < sz; i++)
	midprop[i] = rhs.midprop[i];
    }

    qp_arg = rhs.qp_arg;

    // YA
    lat_back = NULL;
    link_status_smeared = false;
    sink_type = rhs.sink_type;

  }

  return *this;
}

// averaging constructor
QPropW::QPropW (QPropW & prop1, QPropW & prop2):Alg (prop1)
{
  cname = "QPropW";
  char *fname = "QPropW(QPropW&,QPropW&)";
  VRB.Func (cname, fname);

  prop = NULL;
  midprop = NULL;
  prop5d = NULL;

  // YA
  lat_back = NULL;
  link_status_smeared = false;
  sink_type = prop1.sink_type;

  Allocate (PROP);

  int sz = GJP.VolNodeSites ();
  if (GJP.Gparity ())
    sz *= 2;

  for (int i = 0; i < sz; i++)
    prop[i] = ((Float) 0.5) * (prop1.prop[i] + prop2.prop[i]);

  qp_arg = prop1.qp_arg;
}


void QPropW::MultiplyProp (WilsonMatrix * to, WilsonMatrix * from,
			   const Float & renFac)
{
  int sz = GJP.VolNodeSites ();
  if (GJP.Gparity ())
    sz *= 2;

  for (int i = 0; i < sz; i++)
    to[i] = renFac * from[i];
}


void QPropW::GetGaugeFixInfo (char *gfixInfo)
{
  switch (AlgLattice ().FixGaugeKind ()) {

  case FIX_GAUGE_NONE:
    sprintf (gfixInfo, "no GF");
    break;

  case FIX_GAUGE_LANDAU:
    sprintf (gfixInfo, "Landau GF, StpCnd=%0.0E",
	     AlgLattice ().FixGaugeStopCond ());
    break;

  case FIX_GAUGE_COULOMB_T:
    sprintf (gfixInfo, "Coulomb(T) GF, StpCnd=%0.0E",
	     AlgLattice ().FixGaugeStopCond ());
    break;

  default:
    sprintf (gfixInfo, "UNKNOWN GF");

  }
}

void QPropW::GetFermionInfo (char *fermionInfo)
{
  switch (AlgLattice ().Fclass ()) {

  case F_CLASS_DWF:
    sprintf (fermionInfo, "DWF, Ls=%i, M5=%0.2f", GJP.Sites (4),
	     GJP.DwfHeight ());
    break;

  case F_CLASS_NONE:
    sprintf (fermionInfo, "NO FERMION TYPE");
    break;

  case F_CLASS_STAG:
    sprintf (fermionInfo, "staggered fermion");
    break;

  case F_CLASS_WILSON:
    sprintf (fermionInfo, "Wilson fermion");
    break;

  case F_CLASS_CLOVER:
    sprintf (fermionInfo, "Clover fermion");
    break;

  case F_CLASS_ASQTAD:
    sprintf (fermionInfo, "aSqTad fermion");
    break;

  case F_CLASS_P4:
    sprintf (fermionInfo, "P4 fermion");
    break;

  default:
    sprintf (fermionInfo, "UNKNOWN FERMION TYPE");

  }
}

// EES merged with ReRun()
void QPropW::Run (const int do_rerun, const Float precision)
{
  char *fname = "Run()";
  VRB.Func (cname, fname);
  //CJ: make it skip running CG when EigCG is intended
//  if (qp_arg.cg.Inverter == EIGCG)
//    return;

  Float dtime0 = dclock ();

  // Set the node size of the full (non-checkerboarded) fermion field
  //----------------------------------------------------------------
  //size_t f_size = GJP.VolNodeSites() * Lat.FsiteSize()/GJP.SnodeSites();
  int iter;
  Float true_res;

  int Nspins = 4;		// Number of spin components to be done
  // Flag set if sequential propagator 
  int seq_src = ((SrcType () == PROT_U_SEQ) ||
		 (SrcType () == PROT_D_SEQ) || (SrcType () == MESSEQ));
  if (DoHalfFermion ())
    Nspins = 2;

  // does prop exist? Assume it does not.
  int do_cg = 1;

  int StartSpin = 0;
  int EndSpin = 4;
  int StartColor = 0;
  int EndColor = 3;
  Lattice & lat = AlgLattice ();

  if (qp_arg.save_prop == 2) {
    StartSpin = qp_arg.StartSrcSpin;
    EndSpin = qp_arg.EndSrcSpin;
    StartColor = qp_arg.StartSrcColor;
    EndColor = qp_arg.EndSrcColor;
  }
  //-----------------------------------------------------------------
  // TY Add Start
  // we need to store the source
  Float *save_source = NULL;
  // TY Add End
  //-----------------------------------------------------------------
  WilsonMatrix *read_prop = NULL;
  WilsonMatrix *save_prop = NULL;
  int glb_walls = GJP.TnodeSites () * GJP.Tnodes ();

if (do_cg) {

    Allocate (PROP);
    // zero the prop (in case we are doing QED)
    for (int i = 0; i < GJP.VolNodeSites (); i++) {
      prop[i] = 0.0;
    }
    VRB.Result (cname, fname, "Fclass()=%d\n", lat.Fclass ());
    if (lat.F5D () || lat.Fclass () == F_CLASS_BFM) {
      Allocate (PROP5D);
      if (StoreMidprop ())
	Allocate (MIDPROP);

      // zero the 5d prop (in case we are doing QED, or 1 color)
      for (int s = 0; s < GJP.SnodeSites (); s++) {
	int vol = GJP.VolNodeSites ();
	for (int i = 0; i < vol; i++) {
	  int site5d = (i + vol * s);
	  prop5d[site5d] = 0.0;
	}
      }
      // zero the midprop (in case we are doing QED)
      if (StoreMidprop ()) {
	for (int i = 0; i < GJP.VolNodeSites (); i++) {
	  midprop[i] = 0.0;
	}
      }
    }

    FermionVectorTp src;
    FermionVectorTp sol;
    //if(AlgLattice().Fclass() == F_CLASS_DWF || AlgLattice().Fclass() == F_CLASS_MOBIUS )
    FermionVectorTp midsol;

    //-----------------------------------------------------------------
    // TY Add Start
    // For conserved axial current
    if (AlgLattice ().F5D ()) {
      conserved = (Float *) smalloc (cname, fname, "d_conserved_p",
				     glb_walls * sizeof (Float));
      for (int i = 0; i < glb_walls; i++)
	conserved[i] = 0.0;
    } else {
      conserved = NULL;
    }

    spnclr_cnt = 0;

    // we need to store the source
    if (qp_arg.save_prop || do_rerun) {
      int sz = GJP.VolNodeSites () * 288 * sizeof (Float);
      if (GJP.Gparity ())
	sz *= 2;
      save_source = (Float *) smalloc (cname, fname, "save_source", sz);
    }
    // TY Add End
    //-----------------------------------------------------------------

    // in case we do a rerun, we also need to store a propagator
    if (do_rerun) {
      int sz = GJP.VolNodeSites () * sizeof (WilsonMatrix);
      if (GJP.Gparity ())
	sz *= 2;
      read_prop = (WilsonMatrix *) smalloc (cname, fname, "read_prop", sz);

#ifdef USE_QIO
      qio_readPropagator readPropQio (qp_arg.file, QIO_FULL_SOURCE, read_prop,
				      save_source, GJP.argc (), GJP.argv (),
				      VOLFMT);
#endif //USE_QIO

      if (AlgLattice ().Fclass () == F_CLASS_DWF){
	MultiplyProp (read_prop, read_prop, 5.0 - GJP.DwfHeight ());

    } else if (AlgLattice ().Fclass () == F_CLASS_MOBIUS) {

      Float renFac = GJP.Mobius_b () * (4 - GJP.DwfHeight ()) + GJP.DwfA5Inv ();

      for (int ii (0); ii < GJP.VolNodeSites (); ++ii)
	*(read_prop + ii) *= renFac;

    }
  }
  //-----------------------------------------------------------------
  // M. Lightman
  // For m_res
  if (AlgLattice ().F5D ()) {
    j5q_pion =
      (Float *) smalloc (cname, fname, "d_j5q_pion_p",
			 glb_walls * sizeof (Float));
//       j5q_pion = (Float *) smalloc(fsize);

    Float *flt_p = (Float *) j5q_pion;
    for (int i = 0; i < glb_walls; i++)
      *flt_p++ = 0.0;
  } else {
    j5q_pion = NULL;
  }
  // End M. Lightman
  //-----------------------------------------------------------------

  for (int spn = StartSpin; spn < EndSpin; spn++)
    for (int col = StartColor; col < EndColor; col++) {
      VRB.Result (cname, fname, "Starting inversion for spin %d color %d\n",
		  spn, col);

      Float dt_src = -dclock ();
      // initial guess (Zero)
      sol.ZeroSource ();

      if (!do_rerun) {
	SetSource (src, spn, col);

	// store the source
	if (qp_arg.save_prop) {
	  int nmat = GJP.VolNodeSites ();
	  if (GJP.Gparity ())
	    nmat *= 2;
	  for (int index (0); index < nmat; ++index)
	    for (int mm (0); mm < 4; ++mm)
	      for (int cc (0); cc < GJP.Colors (); ++cc) {
		// now same ordering as propagator [volume][spin][color][solution_spin][solution_color][ReIm]
		*(save_source + 288 * index + 72 * mm + 24 * cc + 6 * spn +
		  2 * col) = src[24 * index + 6 * mm + 2 * cc];
		*(save_source + 288 * index + 72 * mm + 24 * cc + 6 * spn +
		  2 * col + 1) = src[24 * index + 6 * mm + 2 * cc + 1];
	      }
	}
      } else {			// rerun
	int nmat = GJP.VolNodeSites ();
	if (GJP.Gparity ())
	  nmat *= 2;
	for (int index (0); index < nmat; ++index) {
	  WilsonMatrix *tmp_mat = (WilsonMatrix *) save_source + index;
	  src.CopyWilsonMatSink (index, spn, col, *tmp_mat);
	}
      }

      if ((DoHalfFermion ()) && (!seq_src))	// Rotate to chiral basis
	src.DiracToChiral ();
      dt_src += dclock ();
      VRB.Result (cname, fname,
		  "Time taken to fix source,etc: %17.10e seconds.\n", dt_src);

      dt_src = -dclock ();
      // Get the prop
      VRB.Debug (cname, fname, "Before CG in QpropW.Run() \n");
      //CG(src, sol, midsol, iter, true_res);
      CG (spn, col, src, sol, midsol, iter, true_res);

      //gauge fix solution
      dt_src += dclock ();
      VRB.Result (cname, fname, "Time taken to CG: %17.10e seconds.\n", dt_src);
      dt_src = -dclock ();
      FixSol (sol);
      if (StoreMidprop ())
	FixSol (midsol);

      // Collect solutions in propagator.
      LoadRow (spn, col, sol, midsol);

      if (DoHalfFermion ()) {	// copy spin 0 to spin 1 and spin 2 to spin 3
	int spn2 = spn + 2;
	if (seq_src) {
	  LoadRow (spn2, col, sol, midsol);
	} else {		// Regular propagator zero the extra components
	  src.ZeroSource ();
	  LoadRow (spn2, col, src, src);
	}
      }
      dt_src += dclock ();
      VRB.Result (cname, fname,
		  "Time taken to fix sink,etc: %17.10e seconds.\n", dt_src);

      if (common_arg->results != 0) {
	FILE *fp;
	if ((fp = Fopen ((char *) common_arg->results, "a")) == NULL) {
	  ERR.FileA (cname, fname, (char *) common_arg->results);
	}
	Fprintf (fp, "Cg iters = %d true residual = %e\n",
		 iter, (Float) true_res);
	Fclose (fp);
      }
    }				// End spin-color loop

  // Rotate the source indices to Chiral basis if needed
  if ((DoHalfFermion ()) && (!seq_src)) {
    for (int s = 0; s < GJP.VolNodeSites (); s++)
      prop[s].SinkChiralToDirac ();	// multiply by V^\dagger

    if (StoreMidprop ())
      for (int s = 0; s < GJP.VolNodeSites (); s++)
	midprop[s].SinkChiralToDirac ();	// multiply by V^\dagger
  }
}

  //Print out time taken to invert
Float dtime1 = dclock ();
VRB.Result (cname, fname, "Time taken to invert: %17.10e seconds.\n",
	    dtime1 - dtime0);

  //-----------------------------------------------------------------
  // TY Add Start
  // Print out conserved axial results
if (AlgLattice ().Fclass () == F_CLASS_DWF
    || AlgLattice ().Fclass () == F_CLASS_MOBIUS) {
  int time_size = GJP.TnodeSites () * GJP.Tnodes ();
  for (int t (0); t < time_size; t++)
    slice_sum ((Float *) & conserved[t], 1, 99);
  if (common_arg->results != 0) {
    FILE *fp;
    if ((fp = Fopen ((char *) common_arg->results, "a")) == NULL) {
      ERR.FileA (cname, fname, (char *) common_arg->results);
    }
    Fprintf (fp, "Conserved Axial w_spect\n");
    for (int t = 0; t < time_size; t++) {
      Fprintf (fp, "%d = %.16e\n", t, conserved[t]);
    }
    Fclose (fp);
  }
  sfree (cname, fname, "conserved", conserved);
}
  // TY Add End
  //-----------------------------------------------------------------

  //-----------------------------------------------------------------
  // M. Lightman
  // Print out J5q Pion contraction
if (AlgLattice ().Fclass () == F_CLASS_DWF
    || AlgLattice ().Fclass () == F_CLASS_MOBIUS) {
  int time_size = GJP.TnodeSites () * GJP.Tnodes ();
  for (int t (0); t < time_size; t++)
    slice_sum ((Float *) & j5q_pion[t], 1, 99);
  if (common_arg->results != 0) {
    FILE *fp1;
    if ((fp1 = Fopen ((char *) common_arg->results, "a")) == NULL) {
      ERR.FileA (cname, fname, (char *) common_arg->results);
    }
    Fprintf (fp1, "J5q Pion Contraction\n");
    for (int t = 0; t < time_size; t++) {
      Fprintf (fp1, "%d = %.16e\n", t, j5q_pion[t]);
    }
    Fclose (fp1);
  }
  sfree (cname, fname, "j5q_pion", j5q_pion);
}
  // End M. Lightman
  //-----------------------------------------------------------------

  //Print out time taken for midpoint calculations.
Float dtime2 = dclock ();
VRB.Result (cname, fname,
	    "Time spent calculating midpoint and conserved axial current: %17.10e seconds.\n",
	    dtime2 - dtime1);

  // save prop
// from master branch
#if 1
   if (do_cg && qp_arg.save_prop) {

     char propType[256], sourceType[256], propOutfile[256];
     char gfixInfo[256];

     switch ( AlgLattice().FixGaugeKind() ){

     case FIX_GAUGE_NONE:
       sprintf(gfixInfo,"no GF");
       break;

     case FIX_GAUGE_LANDAU:
       sprintf(gfixInfo,"Landau GF, StpCnd=%0.0E", AlgLattice().FixGaugeStopCond());
       break;

     case FIX_GAUGE_COULOMB_T:
       sprintf(gfixInfo,"Coulomb(T) GF, StpCnd=%0.0E", AlgLattice().FixGaugeStopCond());
       break;

     default:
       sprintf(gfixInfo,"UNKNOWN GF");

     }


     char fermionInfo[256];
     
     switch (AlgLattice().Fclass() ){

     case F_CLASS_DWF:
       sprintf(fermionInfo,"DWF, Ls=%i, M5=%0.2f",GJP.Sites(4),GJP.DwfHeight());
       break;

     case F_CLASS_MOBIUS:
       sprintf(fermionInfo,"MOBIUS, Ls=%i, M5=%0.2f, B=%0.2f, C=%0.2f",
	       GJP.Sites(4),GJP.DwfHeight(),GJP.Mobius_b(),GJP.Mobius_c());
       break;
	 
     case F_CLASS_NONE:
       sprintf(fermionInfo,"NO FERMION TYPE");
       break;

     case F_CLASS_STAG:
       sprintf(fermionInfo,"staggered fermion");
       break;

     case F_CLASS_WILSON:
       sprintf(fermionInfo,"Wilson fermion");
       break;

     case F_CLASS_CLOVER:
       sprintf(fermionInfo,"Clover fermion");
       break;
 	
     case F_CLASS_ASQTAD:
       sprintf(fermionInfo,"aSqTad fermion");
       break;

     case F_CLASS_P4: 
       sprintf(fermionInfo,"P4 fermion");
       break;

     case F_CLASS_NAIVE: 
       sprintf(fermionInfo,"naive fermion");
       break;

     default:
       sprintf(fermionInfo,"UNKNOWN FERMION TYPE");
      
     }

     sprintf(propType,"4D propagator, mass=%0.4f, StpCond=%0.0E,\nBC=%s%s%s%s,\n%s,\n%s", 
	     qp_arg.cg.mass, qp_arg.cg.stop_rsd,
	     ((GJP.Xbc()==BND_CND_PRD) ? "P" : "A"),
	     ((GJP.Ybc()==BND_CND_PRD) ? "P" : "A"),
	     ((GJP.Zbc()==BND_CND_PRD) ? "P" : "A"),
	     ((GJP.Tbc()==BND_CND_PRD) ? "P" : "A"),   
	     gfixInfo, fermionInfo
	     );
    
     // be a bit more sophisticated
     sprintf(sourceType,"%s-source at t=%i",SourceType_map[SrcType()].name ,SourceTime());
   
     if(!do_rerun)
       sprintf(propOutfile,"%s",qp_arg.file);
     else
       sprintf(propOutfile,"%s.rewrite",qp_arg.file);

     //in case of DWF, renormalize first
     if(AlgLattice().Fclass() == F_CLASS_DWF){
       
       Float renFac = 1./(5. - GJP.DwfHeight());
       
       save_prop = (WilsonMatrix*)smalloc(cname, fname, "save_prop",GJP.VolNodeSites()*sizeof(WilsonMatrix));
       
       for(int ii(0); ii <  GJP.VolNodeSites(); ++ii)
	 *(save_prop + ii) = renFac * prop[ii];

       
     }else if(AlgLattice().Fclass() == F_CLASS_MOBIUS){
       
       Float renFac = 1./((GJP.Mobius_b()*( 4 - GJP.DwfHeight() ) + GJP.DwfA5Inv()));
       
       save_prop = (WilsonMatrix*)smalloc(GJP.VolNodeSites()*sizeof(WilsonMatrix));
       if (save_prop == 0) ERR.Pointer(cname, fname, "pr3op");
       VRB.Smalloc(cname, fname, "save_prop", save_prop,
		   GJP.VolNodeSites() * sizeof(WilsonMatrix));
       
       for(int ii(0); ii <  GJP.VolNodeSites(); ++ii)
	 *(save_prop + ii) = renFac * prop[ii];       
     }else save_prop = &prop[0];
     
#ifdef USE_QIO
     Float qio_time = -dclock();
     
     // always writes the full 4D source
     //qio_writePropagator writePropQio(propOutfile, QIO_FULL_SOURCE, save_prop, save_source,
     //			      qp_arg.ensemble_id, qp_arg.ensemble_label, qp_arg.seqNum, propType, sourceType,
     //			      GJP.argc(), GJP.argv(), VOLFMT);
     
     // write a t-slice/slices or hypercube in some cases for the source
     qio_writePropagator writePropQio;

     writePropQio.setHeader(qp_arg.ensemble_id, qp_arg.ensemble_label, qp_arg.seqNum, propType, sourceType);
     
     // just writing one time-slice?
     if(  (!do_rerun) &&
       ( (SrcType() == POINT) || (SrcType() == VOLUME) || (SrcType() == BOX) || (SrcType() == WALL) ) ){
       VRB.Flow(cname,fname," source-type: %s only write t-slice %i to file\n",SourceType_map[SrcType()].name ,SourceTime());
       writePropQio.setSourceTslice(SourceTime());
     }

     switch(qp_arg.save_prop){
     case 1:
	writePropQio.write_12pairs(propOutfile, QIO_FULL_SOURCE, save_prop, save_source, VOLFMT);
        break;
     case 2:
        for (int spn=StartSpin; spn < EndSpin; spn++){
          for (int col=StartColor; col < EndColor; col++) {
	    char file[256];
	    sprintf(file,"%ss%dc%d",propOutfile,spn,col);
	    writePropQio.write_pair(file, QIO_FULL_SOURCE, save_prop, save_source, spn, col, VOLFMT);
          }
        }
	break;
     default: ERR.General(cname,fname,"invalid save_prop in qp arge\n");
     }
     qio_time +=dclock();
     print_time("QPropW::Run","qio_writePropagator",qio_time);


#endif // USE_QIO
     
     if(AlgLattice().Fclass() == F_CLASS_DWF || AlgLattice().Fclass() == F_CLASS_MOBIUS )
       sfree(save_prop);
     //Print out time taken to save
     Float dtime3 = dclock();
     VRB.Result(cname, fname,
                "Time taken to save: %17.10e seconds.\n",
                dtime3 - dtime2);
   }
#else
  if (do_cg && qp_arg.save_prop) {
    char propOutfile[256];
    if(!do_rerun)
      sprintf(propOutfile,qp_arg.file);
    else
      sprintf(propOutfile,"%s.rewrite",qp_arg.file);
    SaveQProp(&propOutfile[0], (WilsonMatrix*)save_source, !do_rerun);
    delta_t = Timer::relative_time();
    VRB.Result(cname,fname,"Time taken to save: %d hours %d minutes %f seconds.\n",delta_t.hours,delta_t.mins,delta_t.secs);     
  }

#endif


if (qp_arg.save_prop || do_rerun)
sfree (save_source);

  //Compare prop and read_prop if doing rerun
if (do_rerun) {
  CompareProps (prop, read_prop, precision);
  sfree (read_prop);
}

}


void QPropW::ReLoad (char *infile)
{

char *fname = "ReLoad( char *)";
int float_size = sizeof (WilsonMatrix) * GJP.VolNodeSites () / sizeof (Float);
if (GJP.Gparity ())
  float_size *= 2;

if (!prop)
  Allocate (0);			// Allocate only if needed


// if the propagator file of source with spin/color divided exists or not
int sc_part_file_exist = CheckSCfile (infile);
  //TIZB, we only read the propagator file with divided spin/color source if they exists,
  // if not, we will read the propagator with all spin/color source.
if (qp_arg.save_prop == 2 && sc_part_file_exist) {
  //        for(int spin=qpropw_arg.StartSrcSpin;spin<qpropw_arg.EndSrcSpin;spin++){
  //          for(int color=qpropw_arg.StartSrcColor;color<qpropw_arg.EndSrcColor;color++){
  for (int spin = 0; spin < 4; spin++) {
    for (int color = 0; color < 3; color++) {
      char file[256];
      sprintf (file, "%ss%dc%d", infile, spin, color);
      ReLoadSC (file, spin, color);
    }
  }
} else {
  Float *dummy_source = (Float *) smalloc (cname, fname, "dummy_source",
					   float_size * sizeof (Float));
  if (AlgLattice ().Fclass () == F_CLASS_DWF
      || AlgLattice ().Fclass () == F_CLASS_MOBIUS)
#ifdef USE_QIO
    qio_readPropagator readPropQio (infile, &prop[0], dummy_source, float_size,
				    float_size, VOLFMT);
#endif
  sfree (dummy_source);
}
if (AlgLattice ().Fclass () == F_CLASS_DWF)
  MultiplyProp (prop, prop, 5. - GJP.DwfHeight ());
}




//Compare two propagators and write out differences
void QPropW::CompareProps (WilsonMatrix * prop_A, WilsonMatrix * prop_B,
			 const Float & precision)
{
  //now compare prop and read_prop
const char *fname = "CompareProps(WilsonMatrix*, WilsonMatrix*, const Float&)";

Float errCnt (0.);
Float sumerr (0.);

int nstacked = 1;
if (GJP.Gparity ())
  nstacked = 2;
int stkoff = GJP.VolNodeSites ();

for (int stk = 0; stk < nstacked; stk++) {
  for (int index (0); index < GJP.VolNodeSites (); ++index) {

    WilsonMatrix mat_read, mat_calc;

    mat_calc = prop_A[index + stk * stkoff];
    mat_read = *(prop_B + index + stk * stkoff);

    for (int s_src (0); s_src < 4; ++s_src)
      for (int c_src (0); c_src < 3; ++c_src)
	for (int s_snk (0); s_snk < 4; ++s_snk)
	  for (int c_snk (0); c_snk < 3; ++c_snk) {

	    Complex tmp_calc = mat_calc (s_snk, c_snk, s_src, c_src);
	    Complex tmp_read = mat_read (s_snk, c_snk, s_src, c_src);

	    Float diff;
	    diff =
	      (fabs (tmp_calc.real () - tmp_read.real ()) +
	       fabs (tmp_calc.imag () -
		     tmp_read.imag ())) / sqrt ((tmp_calc.real () *
						 tmp_calc.real () +
						 tmp_calc.imag () *
						 tmp_calc.imag ()));

	    if (diff > precision) {

	      errCnt += 1.0;
	      sumerr += diff;
	      if (GJP.Gparity ()) {
		VRB.Result (cname, fname,
			    "mismatch propagator: stacked flavor idx %d index %i snk %i %i src %i %i\n %f: (%f,%f) <-> (%f,%f)\n",
			    stk, index, s_snk, c_snk, s_src, c_src, diff,
			    tmp_calc.real (), tmp_calc.imag (),
			    tmp_read.real (), tmp_read.imag ());
	      } else {
		VRB.Result (cname, fname,
			    "mismatch propagator: index %i snk %i %i src %i %i\n %f: (%f,%f) <-> (%f,%f)\n",
			    index, s_snk, c_snk, s_src, c_src, diff,
			    tmp_calc.real (), tmp_calc.imag (),
			    tmp_read.real (), tmp_read.imag ());
	      }
	    }

	  }
  }
}

glb_sum_five (&errCnt);
glb_sum_five (&sumerr);
Float averr = sumerr / errCnt;

if (fabs (errCnt) > 0.) {
  VRB.Result (cname, fname, " ReRun prop. with TOTAL NUMBER OF ERRORS: %f\n",
	      errCnt);
  VRB.Result (cname, fname, " Average error: %e\n", averr);
  VRB.Result (cname, fname, " The precision is set at: %e\n", precision);
} else {
  VRB.Result (cname, fname, " ReRun prop. successfully!\n");
  VRB.Result (cname, fname, " The precision is set at: %e\n", precision);
}

}

// reload spin-color source(s)
void QPropW::ReLoadSC (char *infile, int spin, int color)
{

char *fname = "ReLoadSC( char *)";

int float_size = sizeof (WilsonMatrix) * GJP.VolNodeSites () / sizeof (Float);


  // assumes prop is already allocated
if (prop == 0)
  ERR.General (cname, fname, "prop should be allocated already");

Float *dummy_source = (Float *) smalloc (cname, fname, "dummy_source",
					 float_size * sizeof (Float));
#ifdef USE_QIO
qio_readPropagator readPropQio (infile, spin, color, &prop[0], dummy_source,
				float_size, float_size, VOLFMT);
#endif

  // must renormalize after return in calling function
  //if(AlgLattice().Fclass() == F_CLASS_DWF){

  //Float renFac = 5.-GJP.DwfHeight();

    //for(int ii(0); ii < GJP.VolNodeSites(); ++ii)
    //prop[ii] *= renFac;

  //}

sfree (cname, fname, "dummy_source", dummy_source);

}

void QPropW::CG (Lattice & lat, CgArg * arg, FermionVectorTp & source,
	       FermionVectorTp & sol, int &iter, Float & true_res)
{
ERR.NotImplemented (cname, "CG(Lattice,CgArg,FermionVector...)");
}

// Do conjugate gradient
//void QPropW::CG(FermionVectorTp& source, FermionVectorTp& sol, 
//              FermionVectorTp& midsol, int& iter, Float& true_res) {
void QPropW::CG (int spn, int col,
	       FermionVectorTp & source, FermionVectorTp & sol,
	       FermionVectorTp & midsol, int &iter, Float & true_res)
{

char *fname = "CG(source&, sol&, midsol&, int&, Float&)";
VRB.Func (cname, fname);

Lattice & Lat = AlgLattice ();

  // Set the node size of the full (non-checkerboarded) fermion field
  //----------------------------------------------------------------
int ls = GJP.SnodeSites ();
int ls_glb = GJP.SnodeSites () * GJP.Snodes ();
size_t f_size = GJP.VolNodeSites () * Lat.FsiteSize () / GJP.SnodeSites ();
if (GJP.Gparity ())
  f_size *= 2;

size_t f_size_5d = f_size * ls;

VRB.Result(cname, fname, "f_size_5d = %d\n", f_size_5d);

  // Do inversion

  //----------------------------------------------------------------
   bool do_5d = false;
if (Lat.F5D () ) {
    do_5d = true;
}
#ifdef USE_BFM
  else if (Lat.Fclass () == F_CLASS_BFM ) {
    if (Fbfm::CurrentSolver() != WilsonFermion &&
	Fbfm::CurrentSolver() != WilsonTM)
      do_5d = true;
  }
#endif
 VRB.Result(cname,fname,"do_5d=%d\m",do_5d);

  if (do_5d) {
    Vector *src_4d = (Vector *) source.data ();
    Vector *sol_4d = (Vector *) sol.data ();
    Vector *midsol_4d = (Vector *) midsol.data ();
    Vector *src_5d =
      (Vector *) smalloc (cname, fname, "src_5d", f_size_5d * sizeof (IFloat));
    Vector *sol_5d =
      (Vector *) smalloc (cname, fname, "sol_5d", f_size_5d * sizeof (IFloat));

    //TIZB 2012-01-29
    // zero clear the sol_5d, important to avoid sys error for AMA 
    sol_5d->VecZero (f_size_5d);

    // //DEBUG - checksum the 4d source
    // {
    //   FPConv fp;
    //   enum FP_FORMAT format = FP_IEEE64LITTLE;
    //   uint32_t csum(0);

    //   Float *field_4D = (Float *)src_4d;
    //   int vol_4d = GJP.VolNodeSites();
    //   for(int x=0; x<vol_4d; x++){
    //  uint32_t csum_contrib = fp.checksum((char *)(field_4D),24,format);
    //  csum += csum_contrib;
    //  field_4D+=24;
    //   }

    //   QioControl qc;
    //   csum = qc.globalSumUint(csum);

    //   if(UniqueID()==0) printf("4D source checksum %u\n",csum);
    // }
    // //DEBUG


    //printf("CG converting 4D source to 5D\n"); //DEBUG
    Lat.Ffour2five (src_5d, src_4d, 0, ls_glb - 1);
    Lat.Ffour2five (sol_5d, sol_4d, ls_glb - 1, 0);
    VRB.Result(cname,fname,"src_5d sol_5d src_4d sol_4d=%e %e %e %e\n",
	src_5d->NormSqGlbSum(f_size_5d),
	sol_5d->NormSqGlbSum(f_size_5d),
	src_4d->NormSqGlbSum(f_size),
	sol_4d->NormSqGlbSum(f_size)
	);

{
      Float *temp_p = (Float*) src_5d;
      VRB.Result (cname, fname, "src_5d = %g %g %g %g %g %g norm=%g\n",
        *temp_p, *(temp_p+1), *(temp_p+2),
        *(temp_p+3), *(temp_p+4), *(temp_p+6),
        src_5d->NormSqGlbSum(f_size_5d));
}

    // do the MADWF or not
    if (!qp_arg.mob_arg_s) {
      VRB.Result (cname, fname, "No MADWF, Fclass()=%d\n", Lat.Fclass ());
      //   if(GJP.Gparity()) nwilson*=2;

      //   FPConv fp;
      //   enum FP_FORMAT format = FP_IEEE64LITTLE;
      //   uint32_t csum(0);

      //   Float *field_5D = (Float *)src_5d;

      //   for(int x=0; x<nwilson; x++){
      //  uint32_t csum_contrib = fp.checksum((char *)(field_5D),24,format);
      //  csum += csum_contrib;
      //  field_5D+=24;
      //   }

      //   QioControl qc;
      //   csum = qc.globalSumUint(csum);

      //   if(UniqueID()==0) printf("5D source checksum %u\n",csum);
      // }


      iter = Lat.FmatInv (sol_5d, src_5d, &(qp_arg.cg), &true_res,
			  CNV_FRM_YES, PRESERVE_NO);

    } else {
      MobiusArg *mob_arg_l = (MobiusArg *) (qp_arg.mob_arg_l);
      MobiusArg *mob_arg_s = (MobiusArg *) (qp_arg.mob_arg_s);
      //printf("TIZB entering MADWF\n");
      iter = Lat.FmatInv (sol_5d, src_5d, mob_arg_l, mob_arg_s, &true_res,
			  CNV_FRM_YES, PRESERVE_NO);
    }

      Float *temp_p = (Float*) sol_5d;
      VRB.Result (cname, fname, "sol_5d = %g %g %g %g %g %g norm=%g\n",
	*temp_p, *(temp_p+1), *(temp_p+2),
	*(temp_p+3), *(temp_p+4), *(temp_p+6),
	sol_5d->NormSqGlbSum(f_size_5d));

#define STORE5DPROP
#ifdef STORE5DPROP
    int vol = GJP.VolNodeSites ();
    if (Lat.F5D ())		//checking prop5d != NULL
      for (int s = 0; s < ls; s++) {
	for (int site = 0; site < vol; site++) {
	  int site5d = (site + vol * s);
	  for (int s1 = 0; s1 < 4; ++s1) {
	    for (int c1 = 0; c1 < 3; ++c1) {
	      int i = c1 + 3 * (s1 + 4 * site5d);
	      Rcomplex cc = *((Rcomplex *) sol_5d + i);
	      prop5d[site5d].load_elem (s1, c1, spn, col, cc);
	      //printf("Storing prop[s=%d,i=%d] spin %d color %d     %e %e\n",
	      //   s,site,s1, cl, cc.real(), cc.imag());
	    }
	  }
	}
      }
#endif

    //-----------------------------------------------------------------
    // TY Add Start
    if (qp_arg.save_ls_prop)
      for (int nls (0); nls < GJP.SnodeSites (); nls++)
	SaveQPropLs (sol_5d, qp_arg.file, nls);
    spnclr_cnt++;

    if (AlgLattice ().F5D () && !GJP.Gparity ()) {
      MeasConAxialOld (sol_5d);
    }
    // TY Add End
    //-----------------------------------------------------------------

    //-----------------------------------------------------------------
    // M. Lightman
    if (AlgLattice ().F5D ()) {
      MeasJ5qPion (sol_5d);
    }
    // End M. Lightman
    //-----------------------------------------------------------------

    // prop on walls
    Lat.Ffive2four (sol_4d, sol_5d, ls_glb - 1, 0);
    // midpoint prop
    if (StoreMidprop ())
      Lat.Ffive2four (midsol_4d, sol_5d, ls_glb / 2 - 1, ls_glb / 2);

    sfree (cname, fname, "sol_5d", sol_5d);
    sfree (cname, fname, "src_5d", src_5d);
  } else {
    //
    //! FIXME !  Don't we need zero clear for solution for AMA here ?
    //  may be none will use sloppy CG except DWF-type
VRB.Result(cname,fname,"source=%p sol=%p\n",source.data(),sol.data());
    iter = Lat.FmatInv ((Vector *) sol.data (), (Vector *) source.data (),
			&(qp_arg.cg), &true_res, CNV_FRM_YES, PRESERVE_NO);
  }

}

/*!
  gauge fix solution - only works for coulomb
  guage in time, or Landau gauge
*/
void QPropW::FixSol (FermionVectorTp & sol)
{
  char *fname = "FixSol()";
  VRB.Func (cname, fname);

  // replace with check for FIX_GAUGE_NONE ??
  if (GFixedSnk ()) {
    Lattice & latt (AlgLattice ());
    switch (latt.FixGaugeKind ()) {
    case (FIX_GAUGE_NONE):
      ERR.General (cname, fname, "Gauge fixing matrices not calculated\n");
      break;
    case (FIX_GAUGE_COULOMB_T):
      // GaugeFixSink seems to be broken for anything
      // except timeslice fixing ( dir isn't even used 
      // in the function )
      sol.GaugeFixSink (latt, 3);
      break;

    case (FIX_GAUGE_LANDAU):
      sol.LandauGaugeFixSink (latt);
      break;

    default:
      // should never be reached
      ERR.General (cname, fname, "Unimplemented gauge fixing\n");
      break;
    }
  }
}

/*!
  gauge unfix solution,needed for QPropWRestore() - only works for coulomb
  guage in time, or Landau gauge
*/
void QPropW::UnfixSol (FermionVectorTp & sol)
{

  char *fname = "FixSol()";
  VRB.Func (cname, fname);

  // replace with check for FIX_GAUGE_NONE ??
  if (GFixedSnk ()) {
    Lattice & latt (AlgLattice ());
    switch (latt.FixGaugeKind ()) {
    case (FIX_GAUGE_NONE):
      ERR.General (cname, fname, "Gauge fixing matrices not calculated\n");
      break;
    case (FIX_GAUGE_COULOMB_T):
      // GaugeFixSink seems to be broken for anything
      // except timeslice fixing ( dir isn't even used 
      // in the function )
      sol.GaugeFixSink (latt, 3, 1);
      break;

//       case ( FIX_GAUGE_LANDAU ):
//         sol.LandauGaugeFixSink(latt);
//         break;

    default:
      // should never be reached
      ERR.General (cname, fname, "Unimplemented gauge fixing\n");
      break;
    }
  }
}

/*!
  Collect the solutions in the propagator
 */
void QPropW::LoadRow (int spin, int color,
		      FermionVectorTp & sol, FermionVectorTp & midsol)
{
  VRB.Result (cname, "LoadRow()", "%d %d %p %p %p %p\n", spin, color, &sol,
	      &midsol, prop, midprop);
  int prop_sz = GJP.Gparity ()? 2 * GJP.VolNodeSites () : GJP.VolNodeSites ();

  for (int s = 0; s < prop_sz; s++) {
    int i = s * SPINOR_SIZE;
    prop[s].load_row (spin, color, (wilson_vector &) sol[i]);
  }

  // Collect solutions in midpoint propagator.
  if (StoreMidprop ()) {
#pragma omp parallel for
    for (int s = 0; s < prop_sz; s++) {
      int i = s * SPINOR_SIZE;
      midprop[s].load_row (spin, color, (wilson_vector &) midsol[i]);
    }
  }
  VRB.FuncEnd (cname, "LoadRow()");
}

/*!
  Reverse of LoadRow, needed for RestoreQProp
 */
void QPropW::SaveRow (int spin, int color, FermionVectorTp & sol,
		      FermionVectorTp & midsol)
{
  int i;
  int prop_sz = GJP.Gparity ()? 2 * GJP.VolNodeSites () : GJP.VolNodeSites ();

  for (int s = 0; s < prop_sz; s++) {
    i = s * SPINOR_SIZE;	// FermionVector index
    prop[s].save_row (spin, color, (wilson_vector &) sol[i]);
  }

  // Collect solutions in midpoint propagator.
  for (int s = 0; s < prop_sz; s++) {
    if (StoreMidprop ()) {
#pragma omp parallel for
      for (int s = 0; s < prop_sz; s++) {
	i = s * SPINOR_SIZE;	// lattice site
	midprop[s].save_row (spin, color, (wilson_vector &) midsol[i]);
      }
    }
  }
}

// Needed by alg_threept. Could be replaced by the disk system.
  void QPropW::ShiftPropForward (int n)
  {

    char *fname = "ShiftPropForward(int)";

    Float *recv_buf;
    Float *send_buf;
    int len = 12 * 12 * 2;
    // size of transfer in words
    recv_buf = (Float *) smalloc (len * sizeof (Float));
    if (recv_buf == 0)
      ERR.Pointer (cname, fname, "recv_buf");
    VRB.Smalloc (cname, fname, "recv_buf", recv_buf, 12 * 12 * sizeof (Float));

    for (int j = 0; j < n; j++) {
      // shift 1 node in t-dir.  prop -> prop
      int sz = GJP.VolNodeSites ();
      if (GJP.Gparity ())
	sz *= 2;

      for (int i = 0; i < sz; i++) {
	send_buf = (Float *) & prop[i];
	getMinusData ((IFloat *) recv_buf, (IFloat *) send_buf, len, 3);
	moveMem ((IFloat *) & prop[i], (IFloat *) recv_buf,
		 len * sizeof (IFloat));
      }
    }

    VRB.Sfree (cname, fname, "recv_buf", recv_buf);
    sfree (recv_buf);

  }
// Needed by alg_threept. Could be replaced by the disk system.
  void QPropW::ShiftPropBackward (int n)
  {

    char *fname = "ShiftPropBack()";

    Float *recv_buf;
    Float *send_buf;
    int len = 12 * 12 * 2;
    // size of transfer in words
    recv_buf = (Float *) smalloc (len * sizeof (Float));
    if (recv_buf == 0)
      ERR.Pointer (cname, fname, "recv_buf");
    VRB.Smalloc (cname, fname, "recv_buf", recv_buf, 12 * 12 * sizeof (Float));

    for (int j = 0; j < n; j++) {
      // shift 1 node in t-dir.  prop -> prop
      int sz = GJP.VolNodeSites ();
      if (GJP.Gparity ())
	sz *= 2;

      for (int i = 0; i < sz; i++) {
	send_buf = (Float *) & prop[i];
	getPlusData ((IFloat *) recv_buf, (IFloat *) send_buf, len, 3);
	moveMem ((IFloat *) & prop[i], (IFloat *) recv_buf,
		 len * sizeof (IFloat));
      }
    }

    VRB.Sfree (cname, fname, "recv_buf", recv_buf);
    sfree (recv_buf);

  }

/*!
  Compute the average of the current propagator and the propagator Q
  \f[
  prop[i] =\frac{1}{2}(prop[i] + Q.prop[i])
  \f]
  It does this for both prop and midprop

  NOTE: This way of doing things saves an extra propagator in storage
  the averaging constructor should be removed since it is now obsolete.
 */
  void QPropW::Average (QPropW & Q)
  {
    int sz = GJP.Gparity ()? 2 * GJP.VolNodeSites () : GJP.VolNodeSites ();

    if ((Q.prop != NULL) && (prop != NULL))
      for (int i = 0; i < sz; i++)
	prop[i] = ((Float) 0.5) * (prop[i] + Q.prop[i]);

    if ((Q.midprop != NULL) && (midprop != NULL))
      for (int i = 0; i < sz; i++)
	midprop[i] = ((Float) 0.5) * (midprop[i] + Q.midprop[i]);
  }

/*!
  Compute a general linear combination of the current propagator 
  and the propagator Q
  \f[
  prop[i] =a*prop[i] + b*Q.prop[i]
  \f]
  It does this for both prop and midprop

  This is needed when taking both P+A and P-A.
 */
  void QPropW::LinComb (QPropW & Q, Float a, Float b)
  {
    int sz = GJP.VolNodeSites ();
    if (GJP.Gparity ())
      sz *= 2;

    if ((Q.prop != NULL) && (prop != NULL))
      for (int i = 0; i < sz; i++)
	prop[i] = a * prop[i] + b * Q.prop[i];

    if ((Q.midprop != NULL) && (midprop != NULL))
      for (int i = 0; i < sz; i++)
	midprop[i] = a * midprop[i] + b * Q.midprop[i];
  }

/*! 
  Compute the norm^2 of the propagator.
*/
  IFloat QPropW::norm (const bool & global_sum) const
  {
    IFloat out (0.0);

    int sz = GJP.VolNodeSites ();
    if (GJP.Gparity ())
      sz *= 2;

    for (int i = 0; i < sz; i++)
      out += prop[i].norm ();

    if (global_sum)
      glb_sum (&out);
    return out;
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
  WilsonMatrix & QPropW::GetMatrix (const int *vec, WilsonMatrix & tmp) const
  {

    // offset out-of-range coordinates site[] into on_node_site[]
    // in order to locate the Matrix
    //------------------------------------------------------------------------
    int on_node_site[4], site[4];
    int on_node = 1;
    WilsonMatrix *on_node_wmat;
    {
      for (int i = 0; i < 4; i++)
      {
	site[i] = on_node_site[i] = vec[i];
	while (on_node_site[i] < 0)
	{
	  on_node_site[i] += GJP.NodeSites (i);
	}
	on_node_site[i] %= GJP.NodeSites (i);
	if (on_node_site[i] != site[i]) {
	  on_node = 0;
	}
      }
      on_node_wmat =
	prop + (on_node_site[0] +
		GJP.XnodeSites () * (on_node_site[1] +
				     GJP.YnodeSites () * (on_node_site[2] +
							  GJP.ZnodeSites () *
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
      WilsonMatrix & recv = tmp;
      for (int i = 0; i < 4; i++) {
	while (site[i] != on_node_site[i]) {
	  if (site[i] < 0) {
	    // the WilsonMatrix has 288 number of floats 
	    //getMinusData((IFloat*)&recv, (IFloat*)&send, sizeof(WilsonMatrix),i);
	    getMinusData ((IFloat *) & recv, (IFloat *) & send, 288, i);
	    on_node_site[i] -= GJP.NodeSites (i);
	  } else {
	    // the WilsonMatrix has 288 number of floats 
	    //getPlusData ((IFloat*)&recv, (IFloat*)&send, sizeof(WilsonMatrix),i);
	    getPlusData ((IFloat *) & recv, (IFloat *) & send, 288, i);
	    on_node_site[i] += GJP.NodeSites (i);
	  }
	  send = recv;
	}
      }
//  VRB.FuncEnd(cname, fname) ;
      return recv;
    }
  }



  void QPropW::SetSource (FermionVectorTp & src, int spin, int color)
  {
    const char *fname = "SetSource()";
    VRB.Func (cname, fname);
  }

  Complex & QPropW::rand_src (int i) const
  {
    //Do nothing ....
    ERR.General ("QPropW", "rand_src", "No random source\n");
    // This is just to keep the compiler happy
    return *((Complex *) prop);
  }

  WilsonMatrix QPropW::WallSinkProp (int t_sink)
  {

    WilsonMatrix wmat = (Float) 0.0;

    for (int i = 0; i < GJP.VolNodeSites (); i++) {
      int t = i / (GJP.VolNodeSites () / GJP.TnodeSites ());
      t += GJP.TnodeCoor () * GJP.TnodeSites ();
      if (t != t_sink)
	continue;
      wmat += prop[i];
    }

#ifdef PARALLEL
    slice_sum ((Float *) & wmat, 288, 99);
#endif

    return wmat;
  }

  WilsonMatrix QPropW::MomSinkProp (int t_sink, int *p)
  {

    ThreeMom mom (p);
    WilsonMatrix wmat = (Float) 0.0;

    Site s;
    for (s.Begin (); s.End (); s.nextSite ()) {
      if (s.physT () != t_sink)
	continue;
      wmat += prop[s.Index ()] * mom.Fact (s);
    }


#ifdef PARALLEL
    slice_sum ((Float *) & wmat, 288, 99);
#endif

    return wmat;
  }

  WilsonMatrix QPropW::TwistMomSinkProp (int t_sink, int *p)
  {

    ThreeMomTwist mom (p);
    WilsonMatrix wmat = (Float) 0.0;

    Site s;
    for (s.Begin (); s.End (); s.nextSite ()) {
      if (s.physT () != t_sink)
	continue;
      wmat += prop[s.Index ()] * mom.Fact (s);
    }


#ifdef PARALLEL
    slice_sum ((Float *) & wmat, 288, 99);
#endif

    return wmat;
  }

  WilsonMatrix QPropW::CosSinkProp (int t_sink, int *p)
  {

    ThreeMom mom (p);
    WilsonMatrix wmat = (Float) 0.0;

    Site s;
    for (s.Begin (); s.End (); s.nextSite ()) {
      if (s.physT () != t_sink)
	continue;
      wmat += prop[s.Index ()] * mom.FactCos (s);
    }


#ifdef PARALLEL
    slice_sum ((Float *) & wmat, 288, 99);
#endif

    return wmat;
  }

  WilsonMatrix QPropW::TwistCosSinkProp (int t_sink, int *p)
  {

    ThreeMomTwist mom (p);
    WilsonMatrix wmat = (Float) 0.0;

    Site s;
    for (s.Begin (); s.End (); s.nextSite ()) {
      if (s.physT () != t_sink)
	continue;
      wmat += prop[s.Index ()] * mom.FactCos (s);
    }


#ifdef PARALLEL
    slice_sum ((Float *) & wmat, 288, 99);
#endif

    return wmat;
  }


  QPropW::~QPropW () {

    char *fname = "~QPropW()";
    VRB.Func (cname, fname);

    Delete (PROP5D);
    Delete (PROP);
    Delete (MIDPROP);
    propls = NULL;
    //  if (Arg.file!=NULL) {
    //  VRB.Sfree(cname, fname, "Arg.file", qp_arg.file);
    //  sfree(Arg.file);
    // }

    // YA 
    //UndoLinkSmear(); // set original link back if replaced with smeared link
    if (lat_back)
      delete[]lat_back;
  }



//Save propagator stored in 'QPropW::prop'. Also save source stored at 'source'. 'allow_save_single_src_tslice' allows us to write only a single
//time slice if the source type lives only on that timeslice.
  void QPropW::SaveQProp (const char *propOutfile, WilsonMatrix * source,
			  const bool & allow_save_single_src_tslice)
  {
    const char *fname =
      "SaveQProp(WilsonMatrix *source, const bool &allow_save_single_src_tslice)";

    char propType[256], sourceType[256];

    char gfixInfo[256];
    GetGaugeFixInfo (gfixInfo);

    char fermionInfo[256];
    GetFermionInfo (fermionInfo);

    sprintf (propType,
	     "4D propagator, mass=%0.4f, StpCond=%0.0E,\nBC=%s%s%s%s,\n%s,\n%s",
	     qp_arg.cg.mass, qp_arg.cg.stop_rsd,
	     BndCndType_map[GJP.Xbc ()].name, BndCndType_map[GJP.Ybc ()].name,
	     BndCndType_map[GJP.Zbc ()].name, BndCndType_map[GJP.Tbc ()].name,
	     gfixInfo, fermionInfo);

    //sprintf(sourceType, "fullSource");
    // be a bit more sophisticated
    sprintf (sourceType, "%s-source at t=%i", SourceType_map[SrcType ()].name,
	     SourceTime ());

    WilsonMatrix *save_prop;

    //in case of DWF, renormalize first
    if (AlgLattice ().Fclass () == F_CLASS_DWF) {
      int sz = GJP.VolNodeSites () * sizeof (WilsonMatrix);
      if (GJP.Gparity ())
	sz *= 2;

      save_prop = (WilsonMatrix *) smalloc (cname, fname, "save_prop", sz);
      MultiplyProp (save_prop, prop, 1.0 / (5.0 - GJP.DwfHeight ()));
    } else {
      save_prop = &prop[0];
    }

#ifdef USE_QIO
    Float qio_time = -dclock ();

    // always writes the full 4D source
    //qio_writePropagator writePropQio(propOutfile, QIO_FULL_SOURCE, save_prop, source,
    //                          qp_arg.ensemble_id, qp_arg.ensemble_label, qp_arg.seqNum, propType, sourceType,
    //                          GJP.argc(), GJP.argv(), VOLFMT);

    // write a t-slice/slices or hypercube in some cases for the source
    qio_writePropagator writePropQio;

    writePropQio.setHeader (qp_arg.ensemble_id, qp_arg.ensemble_label,
			    qp_arg.seqNum, propType, sourceType);

    // just writing one time-slice?
    if (allow_save_single_src_tslice &&
	((SrcType () == POINT) || (SrcType () == VOLUME) || (SrcType () == BOX)
	 || (SrcType () == WALL))) {
      VRB.Flow (cname, fname,
		" source-type: %s only write t-slice %i to file\n",
		SourceType_map[SrcType ()].name, SourceTime ());
      writePropQio.setSourceTslice (SourceTime ());
    }

    if (qp_arg.save_prop == 1) {
      char *file = const_cast < char *>(propOutfile);
      writePropQio.write_12pairs (file, QIO_FULL_SOURCE, save_prop, source,
				  VOLFMT);
    } else if (qp_arg.save_prop == 2) {
      if (GJP.Gparity ())
	ERR.General (cname, fname,
		     "Saving partial G-parity propagator not yet coded\n");

      int StartSpin = qp_arg.StartSrcSpin;
      int EndSpin = qp_arg.EndSrcSpin;
      int StartColor = qp_arg.StartSrcColor;
      int EndColor = qp_arg.EndSrcColor;

      for (int spn = StartSpin; spn < EndSpin; spn++) {
	for (int col = StartColor; col < EndColor; col++) {
	  char file[256];
	  sprintf (file, "%ss%dc%d", propOutfile, spn, col);
	  writePropQio.write_pair (file, QIO_FULL_SOURCE, save_prop, source,
				   spn, col, VOLFMT);
	}
      }
    } else
      ERR.General (cname, fname, "invalid save_prop in qp arge\n");

    qio_time += dclock ();
    print_time ("QPropW::Run", "qio_writePropagator", qio_time);

#endif // USE_QIO
    if (AlgLattice ().Fclass () == F_CLASS_DWF)
      sfree (save_prop);
  }




  void QPropW::SaveQProp (char *name, int mid)
  {

    char *fname = "SaveQProp()";

    ERR.NotImplemented (cname, fname);

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
//#endif

//Greg: Added RestoreQProp, RestoreQPropLs, RestoreQPropLs_ftom
//from v5_0_18 to solve linker errors.
// Restore prop
  void QPropW::RestoreQProp (char *name, int mid)
  {

    char *fname = "RestoreQProp()";
    VRB.Func (cname, fname);

    if (prop == NULL)
      Allocate (PROP);

    // we need to store the source
    Float *read_source =
      (Float *) smalloc (GJP.VolNodeSites () * 288 * sizeof (Float));
    if (read_source == 0)
      ERR.Pointer (cname, fname, "read_source");
    VRB.Smalloc (cname, fname, "read_source", read_source,
		 GJP.VolNodeSites () * 288 * sizeof (Float));

#ifdef USE_QIO

    //char tmp_filename[256];
    // strcpy(tmp_filename, qp_arg.file);

    qio_readPropagator readPropQio (qp_arg.file, QIO_FULL_SOURCE, &prop[0],
				    read_source, GJP.argc (), GJP.argv (),
				    VOLFMT);
#endif // USE_QIO
    sfree (read_source);

    // Flag set if sequential propagator 
    int seq_src = ((SrcType () == PROT_U_SEQ) ||
		   (SrcType () == PROT_D_SEQ) || (SrcType () == MESSEQ));

    if (seq_src) {
      Site s;
      for (s.Begin (); s.End (); s.nextSite ()) {
	QPropW::operator[](s.Index ()).gl (-5);
	QPropW::operator[](s.Index ()).hconj ();
      }
    }
  }

  //-----------------------------------------------------------------
  // TY Add Start
// Save 5d prop at each ls
  void QPropW::SaveQPropLs (Vector * sol_5d, char *name, int ls)
  {

    char *fname = "SaveQPropLs()";

    VRB.Func (cname, fname);

    VRB.Flow (cname, fname, "Saving propagator to pfs...\n");

    int fv_size = GJP.Colors () * 4 * 2 * sizeof (Float) * GJP.VolNodeSites ();
    char sname[100];


    if (qp_arg.save_ls_prop == 1) {
      int skip_buf = fv_size * ls / sizeof (Vector);

#if TARGET == QCDOC
      sprintf (sname, "/pfs/%s.m%0.3f.l%d.id%d.dat",
	       qp_arg.file, qp_arg.cg.mass, ls, UniqueID ());
#else
      sprintf (sname, "pfs/%s.m%0.3f.l%d.id%d.dat",
	       qp_arg.file, qp_arg.cg.mass, ls, UniqueID ());
#endif

      FILE *fp;
      if ((fp = fopen (sname, "a")) != NULL) {
	fwrite (sol_5d + skip_buf, 1, fv_size, fp);
      } else {
	ERR.FileA (cname, fname, sname);
      }
      fclose (fp);

    }

    if (qp_arg.save_ls_prop == 2) {
      // Flag set if sequential propagator 
      int seq_src = ((SrcType () == PROT_U_SEQ) ||
		     (SrcType () == PROT_D_SEQ) || (SrcType () == MESSEQ));
      int spn = spnclr_cnt / GJP.Colors ();
      int clr = spnclr_cnt - spn * GJP.Colors ();
      int shft_buf = GJP.VolNodeSites () * ls;
      int skip_buf = shft_buf * SPINOR_SIZE * sizeof (Float) / sizeof (Vector);
      int i;
      for (int s = 0; s < GJP.VolNodeSites (); s++) {
	i = s * SPINOR_SIZE * sizeof (Float) / sizeof (Vector);	// FermionVector index
	propls[s + shft_buf].load_row (spn, clr,
				       (wilson_vector &) sol_5d[i + skip_buf]);
      }

      if (DoHalfFermion ()) {
	int spn2 = spn + 2;
	if (!seq_src) {
	  FermionVectorTp src;
	  src.ZeroSource ();
	  for (int s = 0; s < GJP.VolNodeSites (); s++) {
	    i = s * SPINOR_SIZE * sizeof (Float) / sizeof (Vector);	// FermionVector index
	    propls[s + shft_buf].load_row (spn2, clr, (wilson_vector &) src[i]);
	  }
	} else {
	  for (int s = 0; s < GJP.VolNodeSites (); s++) {
	    i = s * SPINOR_SIZE * sizeof (Float) / sizeof (Vector);	// FermionVector index
	    propls[s + shft_buf].load_row (spn2, clr,
					   (wilson_vector &) sol_5d[i +
								    skip_buf]);
	  }
	}

      }
      //printf("End propagator to pfs... spn=%d clr=%d\n",spn,clr);
    }
#ifdef PARALLEL
    QioControl sync;
    if (sync.synchronize (1) == 1)
      ERR.General (cname, fname, "Synchronize Error\n");
#endif
  }

// Restore 5d prop at ls
  void QPropW::RestoreQPropLs (char *name, int ls)
  {

    char *fname = "RestoreQPropLs()";
    VRB.Func (cname, fname);

    // Flag set if sequential propagator 
    int seq_src = ((SrcType () == PROT_U_SEQ) ||
		   (SrcType () == PROT_D_SEQ) || (SrcType () == MESSEQ));

    if (qp_arg.save_ls_prop == 1) {
      FermionVectorTp sol;
      FermionVectorTp midsol;

      int fv_size =
	GJP.Colors () * 4 * 2 * sizeof (Float) * GJP.VolNodeSites ();
      char sname[100];

#if TARGET == QCDOC
      sprintf (sname, "/pfs/%s.m%0.3f.l%d.id%d.dat",
	       qp_arg.file, qp_arg.cg.mass, ls, UniqueID ());
#else
      sprintf (sname, "pfs/%s.m%0.3f.l%d.id%d.dat",
	       qp_arg.file, qp_arg.cg.mass, ls, UniqueID ());
#endif

      int Nspin = 4;
      if (DoHalfFermion ())
	Nspin = 2;

      FILE *fp;
      if ((fp = fopen (sname, "r")) != NULL) {
	for (int spn = 0; spn < Nspin; spn++)
	  for (int col = 0; col < GJP.Colors (); col++) {
	    fread ((Vector *) sol.data (), 1, fv_size, fp);
	    LoadRow (spn, col, sol, midsol);

	    if ((DoHalfFermion ()) && (!seq_src)) {
	      int spn2 = spn + 2;
	      sol.ZeroSource ();
	      LoadRow (spn2, col, sol, midsol);
	    } else {
	      int spn2 = spn + 2;
	      LoadRow (spn2, col, sol, midsol);
	    }

	  }
      } else {
	ERR.FileA (cname, fname, name);
      }

      VRB.Flow (cname, fname, "Read prop from file %s\n", sname);
      fclose (fp);
    }

    if (qp_arg.save_ls_prop == 2) {
      size_t f_size = sizeof (WilsonMatrix) * GJP.VolNodeSites () / sizeof (Float);

      int s_local = ls % GJP.SnodeSites ();
      int s_node = ls / GJP.SnodeSites ();
      if (GJP.Snodes () > 1)
	for (int s = 0; s < GJP.VolNodeSites (); s++)
	  prop[s] = 0.0;

      if (s_node == GJP.SnodeCoor ()) {
	int shft_buf = GJP.VolNodeSites () * s_local;
	for (int s = 0; s < GJP.VolNodeSites (); s++) {
	  prop[s] = propls[s + shft_buf];
	  //printf("%d %e %e\n",s,*((Float*)&prop[s]),*((Float*)&propls[s+shft_buf]));
	}
	VRB.Debug (cname, fname, "End read propagator from memory\n");
      }
      if (GJP.Snodes () > 1 && GJP.Snodes () != 2) {
	VRB.Flow (cname, fname, "d gsum start restore\n", f_size);
	Float sum;
	Float *field_4D = (Float *) prop;
	for (int i = 0; i < f_size; i++) {
	  sum = field_4D[i];
	  glb_sum_dir (&sum, 4);
	  field_4D[i] = sum;
	}
	VRB.Flow (cname, fname, "d gsum end restore\n", f_size);
      }
    }
    // Rotate the source indices to Chiral basis if needed
    if ((DoHalfFermion ()) && (!seq_src) && (DoHalfFermion () != 2)) {
      for (int s = 0; s < GJP.VolNodeSites (); s++)
	prop[s].SinkChiralToDirac ();	// multiply by V^\dagger
    }
#ifdef PARALLEL
    QioControl sync;
    if (sync.synchronize (1) == 1)
      ERR.General (cname, fname, "Synchronize error\n");
#endif
  }

// Restore 5d prop from file to memory
  void QPropW::RestoreQPropLs_ftom (char *name)
  {

    char *fname = "RestoreQPropLs()";
    VRB.Func (cname, fname);

    // Flag set if sequential propagator 
    int seq_src = ((SrcType () == PROT_U_SEQ) ||
		   (SrcType () == PROT_D_SEQ) || (SrcType () == MESSEQ));

    if (qp_arg.save_ls_prop == 1) {
      FermionVectorTp sol;
      FermionVectorTp midsol;

      int fv_size =
	GJP.Colors () * 4 * 2 * sizeof (Float) * GJP.VolNodeSites ();
      char sname[100];

      // Allocate 5d memory
      if (propls == NULL) {
	propls =
	  (WilsonMatrix *) smalloc (GJP.VolNodeSites () * GJP.SnodeSites () *
				    sizeof (WilsonMatrix));
	if (propls == 0)
	  ERR.Pointer (cname, fname, "propls");
	VRB.Smalloc (cname, fname, "propls", propls,
		     GJP.VolNodeSites () * GJP.SnodeSites () *
		     sizeof (WilsonMatrix));
	VRB.Debug (cname, fname, "Allocate porpls\n");
      }

      FILE *fp;

      int Nspin = 4;
      if (DoHalfFermion ())
	Nspin = 2;

      for (int ls (0); ls < GJP.SnodeSites (); ls++) {
	int shft_buf = GJP.VolNodeSites () * ls;

#if TARGET == QCDOC
	sprintf (sname, "/pfs/%s.m%0.3f.l%d.id%d.dat",
		 qp_arg.file, qp_arg.cg.mass, ls, UniqueID ());
#else
	sprintf (sname, "pfs/%s.m%0.3f.l%d.id%d.dat",
		 qp_arg.file, qp_arg.cg.mass, ls, UniqueID ());
#endif
	if (fp = fopen (sname, "r")) {
	  for (int spn = 0; spn < Nspin; spn++)
	    for (int col = 0; col < GJP.Colors (); col++) {
	      fread ((Vector *) sol.data (), 1, fv_size, fp);
	      int i;
	      for (int s = 0; s < GJP.VolNodeSites (); s++) {
		i = s * SPINOR_SIZE;	// FermionVector index
		propls[s + shft_buf].load_row (spn, col,
					       (wilson_vector &) sol[i]);
	      }

	      if ((DoHalfFermion ()) && (!seq_src)) {
		int spn2 = spn + 2;
		sol.ZeroSource ();
		for (int s = 0; s < GJP.VolNodeSites (); s++) {
		  i = s * SPINOR_SIZE;	// FermionVector index
		  propls[s + shft_buf].load_row (spn2, col,
						 (wilson_vector &) sol[i]);
		}
	      } else {
		int spn2 = spn + 2;
		for (int s = 0; s < GJP.VolNodeSites (); s++) {
		  i = s * SPINOR_SIZE;	// FermionVector index
		  propls[s + shft_buf].load_row (spn2, col,
						 (wilson_vector &) sol[i]);
		}
	      }

	    }
	} else {
	  ERR.FileA (cname, fname, name);
	}

	VRB.Debug (cname, fname, "Read 5d prop from file to memory\n");
	fclose (fp);
      }				// ls loop

      // 5d prop in memory
      qp_arg.save_ls_prop = 2;
    }
#ifdef PARALLEL
    QioControl sync;
    if (sync.synchronize (1) == 1)
      ERR.General (cname, fname, "Synchronize error\n");
#endif
  }

// Swap 5d prop at ls=0 to GJP.SnodeSites()-1
  void QPropW::SwapQPropLs ()
  {

    char *fname = "RestoreQPropLs()";
    VRB.Func (cname, fname);

    if (qp_arg.save_ls_prop == 1) {
      ERR.General (cname, fname,
		   "qp_arg.save_ls_prop == 1 does not implement.\n");
    }

    if (GJP.Snodes () != 2) {
      ERR.General (cname, fname, "GJP.Snodes() should be 2.\n");
    }

    if (qp_arg.save_ls_prop == 2) {
      size_t f_size = sizeof (WilsonMatrix) * GJP.VolNodeSites () / sizeof (Float);
//    int ls;

      for (int ls = 0; ls < GJP.SnodeSites (); ls++) {
	int s_local = ls % GJP.SnodeSites ();
	int s_node = ls / GJP.SnodeSites ();
//      int l_local = ( ls + GJP.SnodeSites() ) % GJP.SnodeSites();
	int l_node = (ls + GJP.SnodeSites ()) / GJP.SnodeSites ();
	// for(int s=0; s<GJP.VolNodeSites(); s++) prop[s] = 0.0;

	VRB.Debug (cname, fname, "d gsum start swap\n", f_size);
	int shft = f_size * s_local;
	Float sum;
	Float *field_5D = (Float *) propls;
//      Float* field_4D = (Float *) prop;
	Float tmp = 0.;
	for (int i = 0; i < f_size; i++) {
	  // ls=0 to ls=GJP.SnodeCoor()
	  if (s_node == GJP.SnodeCoor ())
	    sum = field_5D[i + shft];
	  if (s_node != GJP.SnodeCoor ())
	    sum = 0;
	  glb_sum_dir (&sum, 4);
	  if (s_node != GJP.SnodeCoor ())
	    tmp = sum;

	  // ls=GJP.SnodeCoor() to ls=0
	  if (l_node == GJP.SnodeCoor ())
	    sum = field_5D[i + shft];
	  if (l_node != GJP.SnodeCoor ())
	    sum = 0;
	  glb_sum_dir (&sum, 4);

	  // propls[ls=0]=sum, propls[ls=GJP.SnodeCoor()]=tmp
	  if (s_node == GJP.SnodeCoor ())
	    field_5D[i + shft] = sum;
	  if (l_node == GJP.SnodeCoor ())
	    field_5D[i + shft] = tmp;
	}
	VRB.Debug (cname, fname, "d gsum end swap\n", f_size);
      }
    }

#ifdef PARALLEL
    QioControl sync;
    if (sync.synchronize (1) == 1)
      ERR.General (cname, fname, "Synchronize error\n");
#endif
  }


  void QPropW::DeleteQPropLs ()
  {

    char *fname = "DeleterQPropLs()";
    VRB.Func (cname, fname);

    if (propls != NULL) {
      VRB.Sfree (cname, fname, "propls", propls);
      sfree (propls);
      propls = NULL;
    }
  }


// Restore 4d prop from file
  void QPropW::RestoreOrgProp (char *name, int ls)
  {

    int ls_glb = GJP.SnodeSites () * GJP.Snodes ();

    char *fname = "RestoreOrgProp()";
    VRB.Func (cname, fname);

    // Flag set if sequential propagator 
    int seq_src = ((SrcType () == PROT_U_SEQ) ||
		   (SrcType () == PROT_D_SEQ) || (SrcType () == MESSEQ));

    FermionVectorTp sol;
    FermionVectorTp midsol;

    int org_ls = ls;

    if (qp_arg.save_ls_prop == 1) {
      int fv_size =
	GJP.Colors () * 4 * 2 * sizeof (Float) * GJP.VolNodeSites ();
      char sname[100], lname[100];

      ls = GJP.SnodeSites () - 1;
#if TARGET == QCDOC
      sprintf (sname, "/pfs/%s.m%0.3f.l%d.id%d.dat",
	       qp_arg.file, qp_arg.cg.mass, ls, UniqueID ());
#else
      sprintf (sname, "pfs/%s.m%0.3f.l%d.id%d.dat",
	       qp_arg.file, qp_arg.cg.mass, ls, UniqueID ());
#endif

      ls = 0;
#if TARGET == QCDOC
      sprintf (lname, "/pfs/%s.m%0.3f.l%d.id%d.dat",
	       qp_arg.file, qp_arg.cg.mass, ls, UniqueID ());
#else
      sprintf (lname, "pfs/%s.m%0.3f.l%d.id%d.dat",
	       qp_arg.file, qp_arg.cg.mass, ls, UniqueID ());
#endif

      int vol_4d = GJP.VolNodeSites ();

      int Nspin = 4;
      if (DoHalfFermion ())
	Nspin = 2;

      FILE *fp, *gp = NULL;
      if ((fp = fopen (sname, "r")) != NULL)
	if ((gp = fopen (lname, "r")) != NULL) {
	  for (int spn = 0; spn < Nspin; spn++)
	    for (int col = 0; col < GJP.Colors (); col++) {
	      fread ((Vector *) sol.data (), 1, fv_size, fp);
	      fread ((Vector *) midsol.data (), 1, fv_size, gp);

	      int x;
	      int i;
	      Float *field_4D;
	      Float *field_5D;
	      field_4D = (Float *) sol.data ();
	      field_5D = (Float *) midsol.data ();
	      field_4D = field_4D + 12;
	      field_5D = field_5D + 12;
	      for (x = 0; x < vol_4d; x++) {
		for (i = 0; i < 12; i++) {
		  field_4D[i] = field_5D[i];
		}
		field_4D = field_4D + 24;
		field_5D = field_5D + 24;
	      }

	      LoadRow (spn, col, sol, midsol);

	      if ((DoHalfFermion ()) && (!seq_src)) {
		int spn2 = spn + 2;
		sol.ZeroSource ();
		LoadRow (spn2, col, sol, midsol);
	      } else {
		int spn2 = spn + 2;
		LoadRow (spn2, col, sol, midsol);
	      }

	    }
	} else {
	  ERR.FileA (cname, fname, name);
	}

      fclose (fp);
      fclose (gp);
    }

    if (qp_arg.save_ls_prop == 2) {
      if (GJP.Snodes () > 1 && org_ls == 0) {
	size_t f_size =
	  sizeof (WilsonMatrix) * GJP.VolNodeSites () / sizeof (Float);

	int s_u_local = (ls_glb - 1) % GJP.SnodeSites ();
	int s_l_local = 0;
	int s_u_node = (ls_glb - 1) / GJP.SnodeSites ();
	int s_l_node = 0;
	for (int s = 0; s < GJP.VolNodeSites (); s++)
	  prop[s] = 0.0;

	if (s_u_node == GJP.SnodeCoor ()) {
	  int shft_buf_l = GJP.VolNodeSites () * s_u_local;
	  for (int s = 0; s < GJP.VolNodeSites (); s++) {
	    WilsonMatrix temp1, temp2;
	    temp1 = propls[s + shft_buf_l];
	    temp2 = 0.0;
	    for (int sr_s = 0; sr_s < 2; sr_s++)
	      for (int sr_c = 0; sr_c < 3; sr_c++) {
		wilson_vector temp3;
		temp3 = temp1.sol (sr_s, sr_c);
		temp2.load_vec (sr_s, sr_c, temp3);
	      }
	    prop[s] = temp2;
	  }
	}

	if (s_l_node == GJP.SnodeCoor ()) {
	  int shft_buf_l = GJP.VolNodeSites () * s_l_local;
	  for (int s = 0; s < GJP.VolNodeSites (); s++) {
	    WilsonMatrix temp1, temp2;
	    temp1 = propls[s + shft_buf_l];
	    temp2 = 0.0;
	    for (int sr_s = 2; sr_s < 4; sr_s++)
	      for (int sr_c = 0; sr_c < 3; sr_c++) {
		wilson_vector temp3;
		temp3 = temp1.sol (sr_s, sr_c);
		temp2.load_vec (sr_s, sr_c, temp3);
	      }
	    prop[s] = temp2;
	  }
	}

	VRB.Debug (cname, fname, "d gsum start org 0\n", f_size);
	Float sum;
	Float *field_4D = (Float *) prop;
	for (int i = 0; i < f_size; i++) {
	  sum = field_4D[i];
	  glb_sum_dir (&sum, 4);
	  field_4D[i] = sum;
	}

      } else if (GJP.Snodes () == 2 && org_ls != 0) {
	size_t f_size =
	  sizeof (WilsonMatrix) * GJP.VolNodeSites () / sizeof (Float);

	int s_u_local = GJP.SnodeSites () - 1;
	int s_l_local = 0;
	int s_u_node = 0;
	int s_l_node = 1;
	for (int s = 0; s < GJP.VolNodeSites (); s++)
	  prop[s] = 0.0;

	if (s_u_node == GJP.SnodeCoor ()) {
	  int shft_buf_l = GJP.VolNodeSites () * s_u_local;
	  for (int s = 0; s < GJP.VolNodeSites (); s++) {
	    WilsonMatrix temp1, temp2;
	    temp1 = propls[s + shft_buf_l];
	    temp2 = 0.0;
	    for (int sr_s = 0; sr_s < 2; sr_s++)
	      for (int sr_c = 0; sr_c < 3; sr_c++) {
		wilson_vector temp3;
		temp3 = temp1.sol (sr_s, sr_c);
		temp2.load_vec (sr_s, sr_c, temp3);
	      }
	    prop[s] = temp2;
	  }
	}

	if (s_l_node == GJP.SnodeCoor ()) {
	  int shft_buf_l = GJP.VolNodeSites () * s_l_local;
	  for (int s = 0; s < GJP.VolNodeSites (); s++) {
	    WilsonMatrix temp1, temp2;
	    temp1 = propls[s + shft_buf_l];
	    temp2 = 0.0;
	    for (int sr_s = 2; sr_s < 4; sr_s++)
	      for (int sr_c = 0; sr_c < 3; sr_c++) {
		wilson_vector temp3;
		temp3 = temp1.sol (sr_s, sr_c);
		temp2.load_vec (sr_s, sr_c, temp3);
	      }
	    prop[s] = temp2;
	  }
	}

	VRB.Debug (cname, fname, "d gsum start org 1\n", f_size);
	Float sum;
	Float *field_4D = (Float *) prop;
	for (int i = 0; i < f_size; i++) {
	  sum = field_4D[i];
	  glb_sum_dir (&sum, 4);
	  field_4D[i] = sum;
	}

      } else {
	ls = 0;
	int shft_buf_z = GJP.VolNodeSites () * ls;
	ls = GJP.SnodeSites () - 1;
	int shft_buf_l = GJP.VolNodeSites () * ls;
	for (int s = 0; s < GJP.VolNodeSites (); s++) {
	  WilsonMatrix temp1, temp2;
	  temp1 = propls[s + shft_buf_z];
	  temp2 = propls[s + shft_buf_l];
	  for (int sr_s = 2; sr_s < 4; sr_s++)
	    for (int sr_c = 0; sr_c < 3; sr_c++) {
	      wilson_vector temp3;
	      temp3 = temp1.sol (sr_s, sr_c);
	      temp2.load_vec (sr_s, sr_c, temp3);
	    }
	  prop[s] = temp2;
	}
      }
    }
    // Rotate the source indices to Chiral basis if needed
    if ((DoHalfFermion ()) && (!seq_src) && (DoHalfFermion () != 2)) {
      for (int s = 0; s < GJP.VolNodeSites (); s++)
	prop[s].SinkChiralToDirac ();	// multiply by V^\dagger
    }

    if (seq_src) {
      Site s;
      for (s.Begin (); s.End (); s.nextSite ()) {
	QPropW::operator[](s.Index ()).gl (-5);
	QPropW::operator[](s.Index ()).hconj ();
      }
    }
#ifdef PARALLEL
    QioControl sync;
    if (sync.synchronize (1) == 1)
      ERR.General (cname, fname, "Synchronize error\n");
#endif
  }

// Dirac to Chiral propagator
  void QPropW::NonRelProp (int lss)
  {
    if (!DoHalfFermion ()) {
      for (int s = 0; s < GJP.VolNodeSites (); s++) {
	//prop[s].PParProjectSink();
	prop[s].PParProjectSource ();
      }

      if (qp_arg.save_ls_prop == 2 && lss == 1) {
	for (int ls = 0; ls < GJP.SnodeSites (); ls++) {
	  int shft_buf = GJP.VolNodeSites () * ls;
	  for (int s = 0; s < GJP.VolNodeSites (); s++) {
	    propls[s + shft_buf].PParProjectSource ();
	  }
	}
      }
      qp_arg.do_half_fermion = 2;
    }
  }


  void (*sproj_tr[8]) (IFloat * f,
		       IFloat * v,
		       IFloat * w, int num_blk, int v_stride, int w_stride);
  int siteOffset (const int lcl[], const int lcl_sites[]);


// Measure conserved axial correlator
  void QPropW::MeasConAxialOld (Vector * sol_5d)
  {

    char *fname = "MeasConAxialOld()";
    VRB.Func (cname, fname);

    int prop_dir = 3;

//  int fv_size = GJP.Colors() * 4 * 2 * sizeof(Float) * GJP.VolNodeSites();
    int ls_glb = GJP.SnodeSites () * GJP.Snodes ();

    const int LORENTZs (4);
    int lcl_sites[LORENTZs];
    lcl_sites[0] = GJP.XnodeSites ();
    lcl_sites[1] = GJP.YnodeSites ();
    lcl_sites[2] = GJP.ZnodeSites ();
    lcl_sites[3] = GJP.TnodeSites ();

    int lclMin[LORENTZs], lclMax[LORENTZs];
    lclMin[0] = lclMin[1] = lclMin[2] = lclMin[3] = 0;
    for (int i (0); i < LORENTZs; i++)
      lclMax[i] = lcl_sites[i] - 1;

    int lcl_node[LORENTZs];
    lcl_node[0] = GJP.XnodeCoor ();
    lcl_node[1] = GJP.YnodeCoor ();
    lcl_node[2] = GJP.ZnodeCoor ();
    lcl_node[3] = GJP.TnodeCoor ();

    int lcl2glb_offset[LORENTZs];
    for (int i = 0; i < LORENTZs; ++i)
      lcl2glb_offset[i] = lcl_sites[i] * lcl_node[i];


    int SPINORs = 2 * GJP.Colors () * LORENTZs;

//  int glb_walls = lcl_sites[prop_dir];
//  int fsize = glb_walls * sizeof(Float);


    // allocate space for two 4d field 
    //-----------------------------------------------------------------------

    int d_size_4d = GJP.VolNodeSites () * SPINORs;

    Float *d_data_p1 = (IFloat *) smalloc (d_size_4d * sizeof (IFloat));
    if (d_data_p1 == 0)
      ERR.Pointer (cname, fname, "d_data_p1");
    VRB.Smalloc (cname, fname, "d_data_p1", d_data_p1,
		 d_size_4d * sizeof (IFloat));

    Float *d_data_p2 = (IFloat *) smalloc (d_size_4d * sizeof (IFloat));
    if (d_data_p2 == 0)
      ERR.Pointer (cname, fname, "d_data_p2");
    VRB.Smalloc (cname, fname, "d_data_p2", d_data_p2,
		 d_size_4d * sizeof (IFloat));

    //-----------------------------------------------------------------------
    // allocate space for two spinor fields for SCU transfer
    //-----------------------------------------------------------------------

    Float *tmp_p1 = (Float *) smalloc (SPINORs * sizeof (Float));
    if (tmp_p1 == 0)
      ERR.Pointer (cname, fname, "tmp_p1");
    VRB.Smalloc (cname, fname, "tmp_p1", tmp_p1, SPINORs * sizeof (Float));

    Float *tmp_p2 = (Float *) smalloc (SPINORs * sizeof (Float));
    if (tmp_p2 == 0)
      ERR.Pointer (cname, fname, "tmp_p2");
    VRB.Smalloc (cname, fname, "tmp_p2", tmp_p2, SPINORs * sizeof (Float));


    //-----------------------------------------------------------------------
    // allocate space for two spinor fields for gamma_5 multiplication 
    //-----------------------------------------------------------------------

    Float *v1_g5 = (Float *) smalloc (SPINORs * sizeof (Float));
    if (v1_g5 == 0)
      ERR.Pointer (cname, fname, "v1_g5");
    VRB.Smalloc (cname, fname, "v1_g5", v1_g5, SPINORs * sizeof (Float));

    Float *v1_next_g5 = (Float *) smalloc (cname, fname,
					   "v1_next_g5",
					   SPINORs * sizeof (Float));
    if (v1_next_g5 == 0)
      ERR.Pointer (cname, fname, "v1_next_g5");
    //      VRB.Smalloc(d_class_name, ctor_str, "v1_next_g5", v1_next_g5,
    //            SPINORs * sizeof(Float)) ;


    //-----------------------------------------------------------------------
    // Fill in the array of sproj_tr functions
    //-----------------------------------------------------------------------
    sproj_tr[SPROJ_XM] = sprojTrXm;
    sproj_tr[SPROJ_YM] = sprojTrYm;
    sproj_tr[SPROJ_ZM] = sprojTrZm;
    sproj_tr[SPROJ_TM] = sprojTrTm;
    sproj_tr[SPROJ_XP] = sprojTrXp;
    sproj_tr[SPROJ_YP] = sprojTrYp;
    sproj_tr[SPROJ_ZP] = sprojTrZp;
    sproj_tr[SPROJ_TP] = sprojTrTp;


    Lattice & lat = AlgLattice ();
    Matrix *gauge_field = lat.GaugeField ();


    //  To avoid duplicate work, stop at the middle of the s direction
    for (int s = 0; s < ls_glb / 2; s++) {
      lat.Ffive2four ((Vector *) d_data_p1, sol_5d, s, s);
      lat.Ffive2four ((Vector *) d_data_p2, sol_5d, ls_glb - 1 - s,
		      ls_glb - 1 - s);

      int lcl_walls = lcl_sites[prop_dir];

      for (int lclW = 0; lclW < lcl_walls; ++lclW) {
	int lcl[LORENTZs];
	int lcl_next[LORENTZs];	// Next site along propagation direction

	// Define hyperplane    
	lclMin[prop_dir] = lclMax[prop_dir] = lclW;

	for (lcl[0] = lclMin[0]; lcl[0] <= lclMax[0]; lcl[0]++)
	  for (lcl[1] = lclMin[1]; lcl[1] <= lclMax[1]; lcl[1]++)
	    for (lcl[2] = lclMin[2]; lcl[2] <= lclMax[2]; lcl[2]++)
	      for (lcl[3] = lclMin[3]; lcl[3] <= lclMax[3]; lcl[3]++) {

		int lcl_offset = siteOffset (lcl, lcl_sites) * SPINORs;

		// coordinates and offset for lcl_next
		for (int i = 0; i < LORENTZs; i++)
		  lcl_next[i] = ((i == prop_dir) ? (lcl[i] + 1) % lcl_sites[i]
				 : lcl[i]);

		int lcl_next_offset =
		  siteOffset (lcl_next, lcl_sites) * SPINORs;

		// U_mu(x) where mu = prop_dir
		Matrix *link =
		  gauge_field + siteOffset (lcl, lcl_sites) * 4 + prop_dir;

		{
		  Float *v1, *v2;
		  Float *v1_next, *v2_next;
		  Float coeff = 1.0;

		  // S_F(x, s)
		  v1 = (Float *) d_data_p1 + lcl_offset;
		  // S_F(x, ls_glb-1-s)
		  v2 = (Float *) d_data_p2 + lcl_offset;

		  // v1_next = S_F(x+prop_dir, s)
		  // v2_next = S_F(x+prop_dir, ls_glb-1-s)
		  if ((lcl[prop_dir] + 1) == lcl_sites[prop_dir]) {
		    getPlusData ((IFloat *) tmp_p1,
				 (IFloat *) d_data_p1 + lcl_next_offset,
				 SPINORs, prop_dir);
		    getPlusData ((IFloat *) tmp_p2,
				 (IFloat *) d_data_p2 + lcl_next_offset,
				 SPINORs, prop_dir);
		    v1_next = tmp_p1;
		    v2_next = tmp_p2;

		    // fix boundary condition
		    switch (prop_dir) {
		    case 0:
		      if (GJP.XnodeBc () == BND_CND_APRD)
			coeff = -coeff;
		      break;
		    case 1:
		      if (GJP.YnodeBc () == BND_CND_APRD)
			coeff = -coeff;
		      break;
		    case 2:
		      if (GJP.ZnodeBc () == BND_CND_APRD)
			coeff = -coeff;
		      break;
		    case 3:
		      if (GJP.TnodeBc () == BND_CND_APRD)
			coeff = -coeff;
		      break;
		    }		// end switch

		  } else {
		    v1_next = (Float *) d_data_p1 + lcl_next_offset;
		    v2_next = (Float *) d_data_p2 + lcl_next_offset;
		  }

		  {
		    // Gamma^5 S_F(x, s)
		    lat.Gamma5 ((Vector *) v1_g5, (Vector *) v1, 1);

		    // Gamma^5 S_F(x+prop_dir, s)
		    lat.Gamma5 ((Vector *) v1_next_g5, (Vector *) v1_next, 1);

		    Matrix tmp1, tmp2, f;
		    Float result = 0.0;

		    // tmp1 = Tr_spin ( (1 + gamma_{prop_dir} ) 
		    //       gamma_5 S_F(x+prop_dir, s) S_F^dagger(x, ls_glb-1-s))
		    sproj_tr[prop_dir + 4] ((IFloat *) & tmp1,
					    (IFloat *) v1_next_g5,
					    (IFloat *) v2, 1, 0, 0);

		    f.DotMEqual (*link, tmp1);

		    result = f.ReTr ();

		    // tmp2 = Tr_spin ( (1 - gamma_{prop_dir} ) 
		    //       gamma_5 S_F(x, s) S_F^dagger(x+prop_dir, ls_glb-1-s))
		    sproj_tr[prop_dir] ((IFloat *) & tmp2,
					(IFloat *) v1_g5,
					(IFloat *) v2_next, 1, 0, 0);

		    tmp1.Dagger (*link);
		    f.DotMEqual (tmp1, tmp2);
		    result -= f.ReTr ();

		    result *= coeff;


		    *(conserved + lclW + lcl2glb_offset[prop_dir]) += result;
		  }

		}
	      }			// for(lcl[_]..)
      }				// for (int lclW...)

      //printf("ls %d %e\n",s,conserved[0]);

    }				// for (int s ... )

    sfree (v1_next_g5);
    sfree (v1_g5);
    sfree (tmp_p2);
    sfree (tmp_p1);
    sfree (d_data_p2);
    sfree (d_data_p1);
  }
  // TY Add End
  //-----------------------------------------------------------------

//-----------------------------------------------------------------
//M. Lightman
  void QPropW::MeasJ5qPion (Vector * sol_5d)
  {

    char *fname = "MeasJ5()";
    VRB.Func (cname, fname);
    VRB.Result (cname, fname, "sol_5d=%p\n", sol_5d);

    int prop_dir = 3;
    int ls_glb = GJP.SnodeSites () * GJP.Snodes ();

    const int LORENTZs (4);
    int SPINORs = 2 * GJP.Colors () * LORENTZs;

    int lcl_sites[LORENTZs];
    lcl_sites[0] = GJP.XnodeSites ();
    lcl_sites[1] = GJP.YnodeSites ();
    lcl_sites[2] = GJP.ZnodeSites ();
    lcl_sites[3] = GJP.TnodeSites ();

    int lclMin[LORENTZs], lclMax[LORENTZs];
//  lclMin[0] = lclMin[1] = lclMin[2] = lclMin[3] = 0;
    for (int i (0); i < LORENTZs; i++) {
      lclMin[i] = 0;
      lclMax[i] = lcl_sites[i] - 1;
    }

    int lcl_node[LORENTZs];
    lcl_node[0] = GJP.XnodeCoor ();
    lcl_node[1] = GJP.YnodeCoor ();
    lcl_node[2] = GJP.ZnodeCoor ();
    lcl_node[3] = GJP.TnodeCoor ();

    int lcl2glb_offset[LORENTZs];
    for (int i = 0; i < LORENTZs; ++i)
      lcl2glb_offset[i] = lcl_sites[i] * lcl_node[i];



    //-----------------------------------------------------------------------
    //allocate space for the 4d field (?)
    //-----------------------------------------------------------------------
    int d_size_4d = GJP.VolNodeSites () * SPINORs;

    Float *d_data = (IFloat *) smalloc (d_size_4d * sizeof (IFloat));
    if (d_data == 0)
      ERR.Pointer (cname, fname, "d_data");
    VRB.Smalloc (cname, fname, "d_data", d_data, d_size_4d * sizeof (IFloat));
    //-----------------------------------------------------------------------

    //-----------------------------------------------------------------------
    //define the 4d field from the 5d field
    //-----------------------------------------------------------------------
    Lattice & lat = AlgLattice ();
    lat.Ffive2four ((Vector *) d_data, sol_5d, ls_glb / 2 - 1, ls_glb / 2);

    //-----------------------------------------------------------------------
    //Calculate the correlator
    //-----------------------------------------------------------------------
    int lcl_walls = lcl_sites[prop_dir];

    for (int lclW = 0; lclW < lcl_walls; ++lclW) {
      int lcl[LORENTZs];

      // Define hyperplane    
      lclMin[prop_dir] = lclMax[prop_dir] = lclW;

      //Loop over sites in this hyperplane in the local lattice
      for (lcl[0] = lclMin[0]; lcl[0] <= lclMax[0]; lcl[0]++)
	for (lcl[1] = lclMin[1]; lcl[1] <= lclMax[1]; lcl[1]++)
	  for (lcl[2] = lclMin[2]; lcl[2] <= lclMax[2]; lcl[2]++)
	    for (lcl[3] = lclMin[3]; lcl[3] <= lclMax[3]; lcl[3]++) {

	      int lcl_offset = siteOffset (lcl, lcl_sites) * SPINORs;

	      Float *v1 = (Float *) d_data + lcl_offset;

	      for (int vec_ind = 0; vec_ind < SPINORs; vec_ind++)
		*(j5q_pion + lclW + lcl2glb_offset[prop_dir]) +=
		  (*(v1 + vec_ind)) * (*(v1 + vec_ind));

	    }			//lcl
    }				//lclW
    sfree (d_data);

  }
//End M. Lightman
//-----------------------------------------------------------------

// YA, this does smearing of the link in the kernel of the Gaussian src smear
  void QPropW::DoLinkSmear (const QPropWGaussArg & gauss_arg)
  {
    char *fname = "DoLinkSmear()";
    VRB.Func (cname, fname);

    Lattice & lattice = AlgLattice ();
    CommonArg ca;		// if output of smearing needed, set ca.filename

    if (link_status_smeared || gauss_arg.gauss_link_smear_type == GKLS_NONE) {
      return;
    }
    if (lat_back != NULL) {
      // and as (! link_status_smeared), then lat_back is smeared link
      // which has been previouly calculated.
      // rotate lattice <-> lat_back

      Matrix *lat_tmp = new Matrix[GJP.VolNodeSites () * 4];
      if (lat_tmp == NULL) {
	ERR.Pointer (cname, cname, "lat_tmp");
      }
      // make the temporal copy of orginal link
      lattice.CopyGaugeField (lat_tmp);
      // set the smeared link
      lattice.GaugeField (lat_back);

      for (int j = 0; j < GJP.VolNodeSites () * 4; j++) {
	lat_back[j] = lat_tmp[j];
      }
      // now lat_back is orginal & lattice is smeared link

      delete[]lat_tmp;

    } else {			// then calculate smeared link

      // print plaq before smearing
      NoArg no_arg;
      AlgPlaq ap (lattice, common_arg, &no_arg);
      if (common_arg->results != 0) {
	FILE *fp;
	if ((fp = Fopen ((char *) common_arg->results, "a")) == NULL) {
	  ERR.FileA (cname, fname, (char *) common_arg->results);
	}
	Fprintf (fp, "QPropW::DoLinkSmear():Plaq_before ");
	Fclose (fp);
      }
      ap.run ();

      lat_back = new Matrix[GJP.VolNodeSites () * 4];
      if (lat_back == NULL) {
	ERR.Pointer (cname, cname, "lat_back");
      }
      // make the copy of orginal link
      lattice.CopyGaugeField (lat_back);

      // AlgSmear* as(NULL);
      // we could use pointer and then run after selection of smearing scheme,
      // if AlgSmear::run() were a virtual function.

      switch (gauss_arg.gauss_link_smear_type) {
	/*
	   case GKLS_NONE:
	   // actually this already returned
	   return;
	   break;
	 */
      case GKLS_APE:
	{
	  ApeSmearArg asa;
	  asa.coef = gauss_arg.gauss_link_smear_coeff;
	  asa.orthog = 3;	// set the smear orthogonal direction to temporal
	  AlgApeSmear as (lattice, &ca, &asa, 1);
	  for (int i = 0; i < gauss_arg.gauss_link_smear_N; i++)
	    as.run ();
	}
	break;
      case GKLS_STOUT:
	ERR.NotImplemented (cname, fname, "GKLS_STOUT yet to implement");
	/*
	   StoutSmearArg asa;
	   asa.rho = action_link_smear_coeff;
	   as = new AlgStoutSmear(lattice,&ca,&asa);
	   break;
	 */
	break;
      default:
	ERR.General (cname, fname, "unknown qp_arg.gauss_link_smear_type=%d",
		     gauss_arg.gauss_link_smear_type);
      }

      // print plaq after smear
      if (common_arg->results != 0) {
	FILE *fp;
	if ((fp = Fopen ((char *) common_arg->results, "a")) == NULL) {
	  ERR.FileA (cname, fname, (char *) common_arg->results);
	}
	Fprintf (fp, "QPropW::DoLinkSmear():Plaq_after  ");
	Fclose (fp);
      }
      ap.run ();
    }

    link_status_smeared = true;
  }

// YA, use this to set back the original link which has been replaced with
//    the smeared link for the Gaissian src smearing.
//    The smeared link is kept in lat_back.
  void QPropW::UndoLinkSmear (const QPropWGaussArg & gauss_arg)
  {
    char *fname = "UndoLinkSmear()";
    VRB.Func (cname, fname);

    if ((!link_status_smeared) || gauss_arg.gauss_link_smear_type == GKLS_NONE) {
      return;
    }
    // then lattice is smeared link & lat_back is original link
    // now rotate lattice <-> lat_back
    Lattice & lattice = AlgLattice ();
    if (lat_back == NULL)
      ERR.General (cname, fname, "lat_back=NULL, DoLinkSmear never called");
    Matrix *lat_tmp = new Matrix[GJP.VolNodeSites () * 4];
    if (lat_tmp == NULL) {
      ERR.Pointer (cname, cname, "lat_tmp");
    }
    // make the temporal copy of smeared link
    lattice.CopyGaugeField (lat_tmp);
    // set the original link back
    lattice.GaugeField (lat_back);

    for (int j = 0; j < GJP.VolNodeSites () * 4; j++) {
      lat_back[j] = lat_tmp[j];
    }
    // now lat_back is smeared & lattice is original link

    delete[]lat_tmp;

    link_status_smeared = false;
  }



//HueyWen and Peter
  QPropWGFLfuncSrc::QPropWGFLfuncSrc (Lattice & lat,
				      CgArg * arg,
				      CommonArg * common_arg,
				      Float (*fn) (int, int, int, int)
): QPropW (lat, common_arg) {
    func = fn;

    char *fname = "QPropWGFLfuncSrc(L&, CgArg*, blah)";
    char *cname = "QPropWGFLfuncSrc";
    VRB.Func (cname, fname);

    // Set the node size of the full (non-checkerboarded) fermion field
    //----------------------------------------------------------------
    size_t f_size = GJP.VolNodeSites () * lat.FsiteSize () / GJP.SnodeSites ();
    int iter = 0;
    Float true_res = 0.;

    // allocate space for the quark propagator
    //----------------------------------------------------------------
    prop =
      (WilsonMatrix *) smalloc (GJP.VolNodeSites () * sizeof (WilsonMatrix));
    if (prop == 0)
      ERR.Pointer (cname, fname, "prop");
    VRB.Smalloc (cname, fname, "prop", prop, 12 * f_size * sizeof (Float));
    FermionVectorTp src;
    FermionVectorTp sol;

    for (int spin = 0; spin < 4; spin++)
      for (int color = 0; color < GJP.Colors (); color++) {

	// initial guess
	sol.SetVolSource (color, spin);
	// set the source

	//src.SetGFLfuncSource(lat, color, spin, fn);

	// Get the prop
	//CG(lat, arg, src, sol, iter, true_res);

	// HueyWen 
	//sol.LandauGaugeFixSink(lat, 3);
	sol.LandauGaugeFixSink (lat);

	// Collect solutions in propagator.
	int j;
	for (int i = 0; i < f_size; i += SPINOR_SIZE) {
	  j = i / SPINOR_SIZE;	// lattice site
	  prop[j].load_vec (spin, color, (wilson_vector &) sol[i]);
	}

	if (common_arg->results != 0) {
	  FILE *fp;
	  if ((fp = Fopen ((char *) common_arg->results, "a")) == NULL) {
	    ERR.FileA (cname, fname, (char *) common_arg->results);
	  }
	  Fprintf (fp, "Cg iters = %d true residual = %e\n", iter, true_res);
	  Fclose (fp);
	}


      }				// End spin-color loop
  }

  QPropWGFLfuncSrc::~QPropWGFLfuncSrc () {
    char *fname = "~QPropWGFLfuncSrc()";
    char *cname = "QPropWGFLfuncSrc";
    VRB.Func (cname, fname);
    VRB.Sfree (cname, fname, "prop", prop);
    sfree (prop);
  }


//------------------------------------------------------------------
// Quark Propagator (Wilson type) with Wall Source
//------------------------------------------------------------------
QPropWWallSrc::QPropWWallSrc (Lattice & lat, CommonArg * c_arg):QPropW (lat,
	  c_arg)
  {

    char *fname = "QPropWWallSrc(L&, ComArg*)";
    cname = "QPropWWallSrc";
    VRB.Func (cname, fname);
  }
QPropWWallSrc::QPropWWallSrc (Lattice & lat, QPropWArg * arg, CommonArg * c_arg):
  QPropW (lat, arg,
	  c_arg) {

    char *fname = "QPropWWallSrc(L&, QPropWArg*, ComArg*)";
    cname = "QPropWWallSrc";
    VRB.Func (cname, fname);
  VRB.Result(cname,fname,"mass=%g\n",arg->cg.mass);
  lat.SetMassArg(arg->cg.mass);
    // get the propagator
    Run ();
  }
QPropWWallSrc::QPropWWallSrc (QPropWWallSrc & prop1, QPropWWallSrc & prop2):
  QPropW (prop1, prop2) {

    char *fname = "QPropWWallSrc(prop&, prop&)";
    cname = "QPropWWallSrc";
    VRB.Func (cname, fname);
  }

//set wall source
  void QPropWWallSrc::SetSource (FermionVectorTp & src, int spin, int color)
  {

    char *fname = "SetSource()";
    VRB.Func (cname, fname);

    src.ZeroSource ();
    src.SetWallSource (color, spin, qp_arg.t, qp_arg.flavor);
    if (GFixedSrc ())
      if (AlgLattice ().FixGaugeKind () == FIX_GAUGE_COULOMB_T)
	src.GFWallSource (AlgLattice (), spin, 3, qp_arg.t, qp_arg.flavor);
      else
	src.GaugeFixVector (AlgLattice (), spin);
  }

// === for gaussian smeared source ===================================

QPropWGaussSrc::QPropWGaussSrc (Lattice & lat, CommonArg * c_arg):QPropW (lat, c_arg)
  {
    char *fname = "QPropWGaussSrc(L&, ComArg*)";
    cname = "QPropWGaussSrc";

    VRB.Func (cname, fname);

  }

QPropWGaussSrc::QPropWGaussSrc (Lattice & lat, QPropWArg * arg, QPropWGaussArg * g_arg, CommonArg * c_arg):
  QPropW (lat, arg, c_arg), gauss_arg (*g_arg)
  {
    char *fname = "QPropWGaussSrc(L&, QPropWArg*, ComArg*)";
    cname = "QPropWGaussSrc";

    VRB.Func (cname, fname);

    // get the propagator
    Run ();
  }


// This routine should be eliminated
QPropWGaussSrc::QPropWGaussSrc (QPropWGaussSrc * prop1, QPropWGaussSrc * prop2):
  QPropW (*prop1, *prop2)
  {
    //char *fname = "QPropWGaussSrc(prop*, prop*)";
    //cname = "QPropWGaussSrc";
    //VRB.Func(cname, fname);
    gauss_arg = prop1->GaussArg ();

  }

QPropWGaussSrc::QPropWGaussSrc (QPropWGaussSrc & prop1, QPropWGaussSrc & prop2):
  QPropW (prop1, prop2)
  {
    //char *fname = "QPropWGaussSrc(prop&, prop&)";
    //cname = "QPropWGaussSrc";
    //VRB.Func(cname, fname);
    gauss_arg = prop1.GaussArg ();
  }

QPropWGaussSrc::QPropWGaussSrc (QPropW & prop1):
  QPropW (prop1) {
    // char *fname = "QPropW(prop&)";
    // cname = "QPropWGaussSrc";
    // VRB.Func(cname, fname);
    gauss_arg = prop1.GaussArg ();
  }


  QPropWGaussSrc::QPropWGaussSrc (Lattice & lat, QPropWArg * arg,
				  QPropWGaussArg * g_arg, CommonArg * c_arg,
				  char *dummy):QPropW (lat, arg, c_arg),
    gauss_arg (*g_arg)
  {
    //char *fname = "QPropWGaussSrc(L&, QPropWArg*, ComArg*)";
    //cname = "QPropWGaussSrc";
    //VRB.Func(cname, fname);
  }


//Set gaussian source
  void QPropWGaussSrc::SetSource (FermionVectorTp & src, int spin, int color)
  {
    char *fname = "SetSource()";
    VRB.Func (cname, fname);

    src.ZeroSource ();
    src.SetPointSource (color, spin, qp_arg.x, qp_arg.y, qp_arg.z, qp_arg.t);
    DoLinkSmear (gauss_arg);	// YA: smear link if needed
    src.GaussianSmearVector (AlgLattice (), spin, gauss_arg.gauss_N,
			     gauss_arg.gauss_W, qp_arg.t);
    UndoLinkSmear (gauss_arg);	// YA: get back the original link
  }

// Smear the sink of a propagator
  void QPropW::GaussSmearSinkProp (int t_sink, const QPropWGaussArg & gauss_arg)
  {
    Site site;
    FermionVectorTp tmp;
    WilsonMatrix wm;

    DoLinkSmear (gauss_arg);	// YA: smear link if needed
    for (int spin (0); spin < 4; spin++)
      for (int color (0); color < 3; color++) {
	//Copy to FermionVector
	for (site.Begin (); site.End (); site.nextSite ())
	  tmp.CopyWilsonMatSink (site.Index (), spin, color,
				 prop[site.Index ()]);
	//smear sink
	for (int s (0); s < 4; s++)
	  tmp.GaussianSmearVector (AlgLattice (), s, gauss_arg.gauss_N,
				   gauss_arg.gauss_W, t_sink);
	//Copy back to propagator
	for (site.Begin (); site.End (); site.nextSite ()) {
	  int i (site.Index () * SPINOR_SIZE);
	  prop[site.Index ()].load_row (spin, color, (wilson_vector &) tmp[i]);
	}
      }
    UndoLinkSmear (gauss_arg);	// YA: get original link back

    qp_arg.SeqSmearSink = GAUSS_GAUGE_INV;
    // This is a little dangerous. It's up to the user to know which
    // sink time slice is smeared. For the moment it's OK
  }

// Smear the sink of a propagator
  void QPropW::GaussSmearSinkProp (const QPropWGaussArg & gauss_arg)
  {
    Site site;
    FermionVectorTp tmp;
    WilsonMatrix wm;

    DoLinkSmear (gauss_arg);	// YA: smear link if needed
    for (int spin (0); spin < 4; spin++)
      for (int color (0); color < 3; color++) {
	//Copy to FermionVector
	for (site.Begin (); site.End (); site.nextSite ())
	  tmp.CopyWilsonMatSink (site.Index (), spin, color,
				 prop[site.Index ()]);
	//smear sink
	for (int s (0); s < 4; s++)
	  tmp.GaussianSmearVector (AlgLattice (), s, gauss_arg.gauss_N,
				   gauss_arg.gauss_W);
	//Copy back to propagator
	for (site.Begin (); site.End (); site.nextSite ()) {
	  int i (site.Index () * SPINOR_SIZE);
	  prop[site.Index ()].load_row (spin, color, (wilson_vector &) tmp[i]);
	}
      }
    UndoLinkSmear (gauss_arg);	// YA: get original link back

    qp_arg.SeqSmearSink = GAUSS_GAUGE_INV;
    // This is a little dangerous. It's up to the user to know which
    // sink time slice is smeared. For the moment it's OK
  }

// === for exponetial smeared source ===================================


// === for multi gaussian smeared source ===================================

QPropWMultGaussSrc::QPropWMultGaussSrc (Lattice & lat, CommonArg * c_arg):QPropW (lat, c_arg)
  {
    char *fname = "QPropWGaussSrc(L&, ComArg*)";
    cname = "QPropWGaussSrc";

    VRB.Func (cname, fname);

  }

QPropWMultGaussSrc::QPropWMultGaussSrc (Lattice & lat, QPropWArg * arg, QPropWGaussArg * g_arg, CommonArg * c_arg):
  QPropW (lat, arg, c_arg), gauss_arg (*g_arg)
  {
    char *fname = "QPropWGaussSrc(L&, QPropWArg*, ComArg*)";
    cname = "QPropWGaussSrc";

    VRB.Func (cname, fname);

    // get the propagator
    Run ();
  }


// This routine should be eliminated
QPropWMultGaussSrc::QPropWMultGaussSrc (QPropWGaussSrc * prop1, QPropWGaussSrc * prop2):
  QPropW (*prop1, *prop2)
  {
    //char *fname = "QPropWGaussSrc(prop*, prop*)";
    //cname = "QPropWGaussSrc";
    //VRB.Func(cname, fname);
    gauss_arg = prop1->GaussArg ();

  }

QPropWMultGaussSrc::QPropWMultGaussSrc (QPropWGaussSrc & prop1, QPropWGaussSrc & prop2):
  QPropW (prop1, prop2)
  {
    //char *fname = "QPropWGaussSrc(prop&, prop&)";
    //cname = "QPropWGaussSrc";
    //VRB.Func(cname, fname);
    gauss_arg = prop1.GaussArg ();
  }

QPropWMultGaussSrc::QPropWMultGaussSrc (QPropW & prop1):
  QPropW (prop1) {
    const char *fname = "QPropW(prop&)";
    cname = "QPropWGaussSrc";
    // VRB.Func(cname, fname);
    ERR.NotImplemented (cname, fname);
//  gauss_arg = prop1.GaussArg();
  }
//Set gaussian source
  void QPropWMultGaussSrc::SetSource (FermionVectorTp & src, int spin,
				      int color)
  {
    char *fname = "SetSource()";
    VRB.Func (cname, fname);

    src.ZeroSource ();

    for (int nt (0); nt < gauss_arg.nt; nt++) {
      src.SetPointSource (color, spin, qp_arg.x, qp_arg.y, qp_arg.z,
			  gauss_arg.mt[nt]);
      DoLinkSmear (gauss_arg);	// YA: smear link if needed
      src.GaussianSmearVector (AlgLattice (), spin, gauss_arg.gauss_N,
			       gauss_arg.gauss_W, gauss_arg.mt[nt]);
      UndoLinkSmear (gauss_arg);	// YA: get back the original link
    }
  }

// === for exponetial smeared source ===================================


// exponential smeared source
QPropWExpSrc::QPropWExpSrc (Lattice & lat, CommonArg * c_arg):QPropW (lat, c_arg)
  {
    char *fname = "QPropWExpSrc(L&, ComArg*)";
    cname = "QPropWExpSrc";

    VRB.Func (cname, fname);
  }

QPropWExpSrc::QPropWExpSrc (Lattice & lat, QPropWArg * arg, QPropWExpArg * e_arg, CommonArg * c_arg):QPropW (lat, arg, c_arg), exp_arg (*e_arg)
  {
    char *fname = "QPropWExpSrc(L&, QPropWArg*, ComArg*)";
    cname = "QPropWExpSrc";

    VRB.Func (cname, fname);

    // get the propagator
    Run ();
  }


QPropWExpSrc::QPropWExpSrc (QPropWExpSrc * prop1, QPropWExpSrc * prop2):
  QPropW (*prop1, *prop2) {
    //char *fname = "QPropWExpSrc(prop*, prop*)";
    //cname = "QPropWExpSrc";
    //VRB.Func(cname, fname);

  }

QPropWExpSrc::QPropWExpSrc (QPropWExpSrc & prop1, QPropWExpSrc & prop2):
  QPropW (prop1, prop2) {
    //char *fname = "QPropWExpSrc(prop&, prop&)";
    //cname = "QPropWExpSrc";
    //VRB.Func(cname, fname);

  }

QPropWExpSrc::QPropWExpSrc (QPropW & prop1):
  QPropW (prop1) {
    //CJ: Have to set exp_arg: how?
    //char *fname = "QPropW(prop&)";
    //cname = "QPropWExpSrc";
    //VRB.Func(cname, fname);
    exp_arg.exp_A = 1.2;
    exp_arg.exp_B = 0.1;
    exp_arg.exp_C = 8;
  }

// Set exponential smeared source
  void QPropWExpSrc::SetSource (FermionVectorTp & src, int spin, int color)
  {
    char *fname = "SetSource()";
    VRB.Func (cname, fname);

    src.ZeroSource ();
    //src.SetExpSource(color, spin, qp_arg.x, qp_arg.y, qp_arg.z, qp_arg.t,
    //                 exp_arg.exp_A, exp_arg.exp_B, exp_arg.exp_C);
    if (qp_arg.gauge_fix_src)
      src.GFWallSource (AlgLattice (), spin, 3, qp_arg.t);
  }

//------------------------------------------------------------------
// Quark Propagator (Wilson type) with Momentum Source
//------------------------------------------------------------------
QPropWMomSrc::QPropWMomSrc (Lattice & lat, CommonArg * c_arg):
  QPropWWallSrc (lat, c_arg) {

    char *fname = "QPropWMomSrc(L&, ComArg*)";
    cname = "QPropWMomSrc";
    VRB.Func (cname, fname);
  }
  QPropWMomSrc::QPropWMomSrc (Lattice & lat, QPropWArg * arg,
			      int *p, CommonArg * c_arg):QPropWWallSrc (lat,
									c_arg),
    mom (p)
  {

    char *fname = "QPropWMomSrc(L&, QPropWArg*, ComArg*)";
    cname = "QPropWMomSrc";
    VRB.Func (cname, fname);

//------------------------------------------------------------------
// T.Y Add
    qp_arg = *arg;
// T.Y Add
//------------------------------------------------------------------

    Run ();
  }
// copy constructor
  QPropWMomSrc::QPropWMomSrc (const QPropWMomSrc & rhs):QPropWWallSrc (rhs),
    mom (rhs.mom)
  {

    char *fname = "QPropW(const QPropW&)";
    cname = "QPropW";
    VRB.Func (cname, fname);
  }

  void QPropWMomSrc::SetSource (FermionVectorTp & src, int spin, int color)
  {
    char *fname = "SetSource()";
    VRB.Func (cname, fname);

    src.ZeroSource ();
    src.SetMomSource (color, spin, qp_arg.t, mom, qp_arg.flavor);
    if (GFixedSrc ())
      src.GFWallSource (AlgLattice (), spin, 3, qp_arg.t, qp_arg.flavor);
  }

//------------------------------------------------------------------
// Quark Propagator (Wilson type) with Momentum Cosine Source
//------------------------------------------------------------------
QPropWMomCosSrc::QPropWMomCosSrc (Lattice & lat, CommonArg * c_arg):
  QPropWWallSrc (lat, c_arg) {

    char *fname = "QPropWMomCosSrc(L&, ComArg*)";
    cname = "QPropWMomCosSrc";
    VRB.Func (cname, fname);
  }
  QPropWMomCosSrc::QPropWMomCosSrc (Lattice & lat, QPropWArg * arg,
				    const int *p,
				    CommonArg * c_arg):QPropWWallSrc (lat,
								      c_arg),
    mom (p)
  {

    char *fname = "QPropWMomCosSrc(L&, QPropWArg*, ComArg*)";
    cname = "QPropWMomCosSrc";
    VRB.Func (cname, fname);

//------------------------------------------------------------------
// T.Y Add
    qp_arg = *arg;
// T.Y Add
//------------------------------------------------------------------

    Run ();
  }

// copy constructor
  QPropWMomCosSrc::
    QPropWMomCosSrc (const QPropWMomCosSrc & rhs):QPropWWallSrc (rhs),
    mom (rhs.mom)
  {

    char *fname = "QPropW(const QPropW&)";
    cname = "QPropW";
    VRB.Func (cname, fname);
  }

  void QPropWMomCosSrc::SetSource (FermionVectorTp & src, int spin, int color)
  {
    char *fname = "SetSource()";
    VRB.Func (cname, fname);

    src.ZeroSource ();
    src.SetMomCosSource (color, spin, qp_arg.t, mom, qp_arg.flavor);
    if (GFixedSrc ())
      src.GFWallSource (AlgLattice (), spin, 3, qp_arg.t, qp_arg.flavor);
  }

//------------------------------------------------------------------
// Quark Propagator (Wilson type) with Momentum Cosine Source
// Twisted Boundary Conditions
//------------------------------------------------------------------
QPropWMomCosTwistSrc::QPropWMomCosTwistSrc (Lattice & lat, CommonArg * c_arg):
  QPropWWallSrc (lat, c_arg) {

    char *fname = "QPropWMomCosTwistSrc(L&, ComArg*)";
    cname = "QPropWMomCosTwistSrc";
    VRB.Func (cname, fname);
  }
  QPropWMomCosTwistSrc::QPropWMomCosTwistSrc (Lattice & lat, QPropWArg * arg,
					      const int *p,
					      CommonArg *
					      c_arg):QPropWWallSrc (lat, c_arg),
    mom (p)
  {

    char *fname = "QPropWMomCosTwistSrc(L&, QPropWArg*, ComArg*)";
    cname = "QPropWMomCosTwistSrc";
    VRB.Func (cname, fname);

//------------------------------------------------------------------
// T.Y Add
    qp_arg = *arg;
// T.Y Add
//------------------------------------------------------------------

    Run ();
  }
// copy constructor
  QPropWMomCosTwistSrc::
    QPropWMomCosTwistSrc (const QPropWMomCosTwistSrc & rhs):QPropWWallSrc (rhs),
    mom (rhs.mom)
  {

    char *fname = "QPropW(const QPropW&)";
    cname = "QPropW";
    VRB.Func (cname, fname);
  }

  void QPropWMomCosTwistSrc::SetSource (FermionVectorTp & src, int spin,
					int color)
  {
    char *fname = "SetSource()";
    VRB.Func (cname, fname);

    src.ZeroSource ();
    src.SetMomCosTwistSource (color, spin, qp_arg.t, mom);
    if (GFixedSrc ())
      src.GFWallSource (AlgLattice (), spin, 3, qp_arg.t);
  }


//------------------------------------------------------------------
// Quark Propagator (Wilson type) with Volume Source
//------------------------------------------------------------------
QPropWVolSrc::QPropWVolSrc (Lattice & lat, CommonArg * c_arg):QPropW (lat,
	  c_arg)
  {

    char *fname = "QPropWVolSrc(L&, ComArg*)";
    cname = "QPropWVolSrc";
    VRB.Func (cname, fname);
  }
QPropWVolSrc::QPropWVolSrc (Lattice & lat, QPropWArg * arg, CommonArg * c_arg):
  QPropW (lat, arg, c_arg) {

    char *fname = "QPropWVolSrc(L&, QPropWArg*, ComArg*)";
    cname = "QPropWVolSrc";
    VRB.Func (cname, fname);

    Run ();
  }

  void QPropWVolSrc::SetSource (FermionVectorTp & src, int spin, int color)
  {

    char *fname = "SetSource()";
    VRB.Func (cname, fname);
    src.ZeroSource ();
    src.SetVolSource (color, spin, qp_arg.flavor);
    if (GFixedSrc ())
      if (AlgLattice ().FixGaugeKind () == FIX_GAUGE_COULOMB_T)
	//   for (int t=0; t<GJP.Tnodes()*GJP.TnodeSites(); t++)
	//  src.GFWallSource(AlgLattice(), spin, 3, t, qp_arg.flavor); //assumes Coulomb gauge in T-direction!
	// else
	src.GaugeFixVector (AlgLattice (), spin);	//works for all gauges

  }

//------------------------------------------------------------------
// Quark Propagator (Wilson type) with Point Source
//------------------------------------------------------------------
QPropWPointSrc::QPropWPointSrc (Lattice & lat, CommonArg * c_arg):
  QPropW (lat, c_arg) {

    char *fname = "QPropWPointSrc(L&, ComArg*)";
    cname = "QPropWPointSrc";
    VRB.Func (cname, fname);
  }
QPropWPointSrc::QPropWPointSrc (Lattice & lat, QPropWArg * arg, CommonArg * c_arg):
  QPropW (lat, arg,
	  c_arg) {

    char *fname = "QPropWPointSrc(L&, ComArg*)";
    cname = "QPropWPointSrc";
    VRB.Func (cname, fname);

    Run ();
  }

// used for reading in a prop
  QPropWPointSrc::QPropWPointSrc (Lattice & lat, QPropWArg * arg,
				  CommonArg * c_arg, char *dummy):QPropW (lat,
									  arg,
									  c_arg)
  {

    char *fname = "QPropWPointSrc(L&, qarg*, ComArg*, char*)";
    cname = "QPropWPointSrc";
    VRB.Func (cname, fname);
  }

  void QPropWPointSrc::SetSource (FermionVectorTp & src, int spin, int color)
  {

    char *fname = "SetSource()";
    VRB.Func (cname, fname);

    src.ZeroSource ();
    src.SetPointSource (color, spin, qp_arg.x, qp_arg.y, qp_arg.z, qp_arg.t,
			qp_arg.flavor);
  }

//------------------------------------------------------------------
// Quark Propagator (Wilson type) with Box Source
//------------------------------------------------------------------
QPropWBoxSrc::QPropWBoxSrc (Lattice & lat, CommonArg * c_arg):
  QPropW (lat, c_arg) {

    char *fname = "QPropWBoxSrc(L&, ComArg*)";
    cname = "QPropWBoxSrc";
    VRB.Func (cname, fname);
  }

QPropWBoxSrc::QPropWBoxSrc (Lattice & lat, QPropWArg * arg, QPropWBoxArg * b_arg, CommonArg * c_arg):
  QPropW (lat, arg, c_arg),
    box_arg (*b_arg) {

    char *fname = "QPropWBoxSrc(L&, QPropWArg*, ComArg*)";
    cname = "QPropWBoxSrc";
    VRB.Func (cname, fname);

    //printf("TIZB gfix source %d\n", arg->gauge_fix_src);
    //TIZB, why only this was  NOT commented out > Tom || Meifeng || Chulwoo ?
    //  Run();

  }

  void QPropWBoxSrc::SetSource (FermionVectorTp & src, int spin, int color)
  {

    char *fname = "SetSource()";
    VRB.Func (cname, fname);

    if (box_arg.use_xyz_offset) {
      int src_offset[3] = { qp_arg.x, qp_arg.y, qp_arg.z };
      src.SetBoxSource (color, spin, box_arg.box_start, box_arg.box_end,
			qp_arg.t, src_offset);
    } else {
      src.SetBoxSource (color, spin, box_arg.box_start, box_arg.box_end,
			qp_arg.t);
    }
    if (GFixedSrc ())
      src.GFWallSource (AlgLattice (), spin, 3, qp_arg.t);
    else
      VRB.Warn (cname, fname, "Warning: box src not gauge fixed");
  }

// ------------------------------------------------------------------
// Quark Propagator with 4D box source
//
// Added by Hantao to handle 4D boxes, also suitable for any uniform
// point/wall/box sources.
// ------------------------------------------------------------------
  QPropW4DBoxSrc::QPropW4DBoxSrc (Lattice & lat, QPropWArg * arg,
				  QPropW4DBoxArg * b_arg, CommonArg * c_arg)
:  QPropW (lat, arg, c_arg) {
    cname = "QPropW4DBoxSrc";

    for (int mu = 0; mu < 4; ++mu) {
      box_arg.box_start[mu] = b_arg->box_start[mu];
      box_arg.box_size[mu] = b_arg->box_size[mu];
      box_arg.mom[mu] = b_arg->mom[mu];
    }

    Run();
  }

  void QPropW4DBoxSrc::SetSource (FermionVectorTp & src, int spin, int color)
  {
    const char *fname = "SetSource()";
    VRB.Func (cname, fname);

    src.Set4DBoxSource(color, spin, box_arg.box_start, box_arg.box_size,
			box_arg.mom);

    if (GFixedSrc ()) {
      // only works for t direction!
      int t_size_glb = GJP.Tnodes () * GJP.TnodeSites ();
      for (int t = 0; t < t_size_glb; t++) {
	if ((t + t_size_glb - box_arg.box_start[3]) % t_size_glb >=
	    box_arg.box_size[3])
	  continue;
	src.GFWallSource (AlgLattice (), spin, 3, t);
      }
    } else {
      VRB.Warn (cname, fname, "Warning: 4D box src not gauge fixed");
    }
  }

// ------------------------------------------------------------------
// Quark Propagator, wall source filled with Z3 boxes
//
// Added by Hantao
// ------------------------------------------------------------------
  QPropWZ3BWallSrc::QPropWZ3BWallSrc (Lattice & lat, QPropWArg * arg,
				      QPropW4DBoxArg * b_arg, CommonArg * c_arg)
:  QPropW (lat, arg, c_arg) {
    cname = "QPropWZ3BWallSrc";
    const char *fname = "QPropWZ3BWallSrc()";

    for (int mu = 0; mu < 4; ++mu) {
      box_arg.box_start[mu] = b_arg->box_start[mu];
      box_arg.box_size[mu] = b_arg->box_size[mu];
      box_arg.mom[mu] = b_arg->mom[mu];
      if (box_arg.box_size[mu] <= 0) {
	ERR.General (cname, fname, "Invalid size in %d direction: %d.\n",
		     mu, box_arg.box_size[mu]);
      }
    }
    if (box_arg.box_size[3] != 1) {
      ERR.NotImplemented (cname, fname);
    }
    if (box_arg.box_start[3] != arg->t) {
      ERR.General (cname, fname,
		   "BoxArg and QPropWArg starting time does not match.\n");
    }

    const int glb[4] = {
      GJP.XnodeSites () * GJP.Xnodes (), GJP.YnodeSites () * GJP.Ynodes (),
      GJP.ZnodeSites () * GJP.Znodes (), GJP.TnodeSites () * GJP.Tnodes (),
    };

    for (int i = 0; i < 3; ++i) {
      rand_grid[i] = (glb[i]) / box_arg.box_size[i];
    }
    rand_size = rand_grid[0] * rand_grid[1] * rand_grid[2];
    rand_num.assign (rand_size, 0);

    const Rcomplex Z3consts[3] = { 1,
      Rcomplex (-0.5, 0.5 * sqrt (3.0)),
      Rcomplex (-0.5, -0.5 * sqrt (3.0)),
    };

    VRB.Result (cname, fname, "rand_size = %d %d %d = %d\n",
		rand_grid[0], rand_grid[1], rand_grid[2], rand_size);

    for (int i = 0; i < rand_size; ++i) {
#ifndef C11
      unsigned rnd = drand48 () * 3;
#else
      LRG.AssignGenerator (0);
      unsigned rnd = LRG.Urand (3, 0);
#endif
      assert (rnd < 3);
      rand_num[i] = Z3consts[rnd];
    VRB.Result (cname, fname, "rand_num[%d]=%g %g\n", i, rand_num[i].real(),rand_num[i].imag());
    }
#ifdef USE_QMP			//should be OK as QMP is needed for any multi-node build now
    QMP_broadcast (rand_num.data (), sizeof (Rcomplex) * rand_size);
#endif

    Run ();
  }

  void QPropWZ3BWallSrc::SetSource (FermionVectorTp & src, int spin, int color)
  {
    const char *fname = "SetSource()";
    VRB.Func (cname, fname);
    VRB.Result (cname, fname, "Set Z3 boxed wall source at t=%d.\n", qp_arg.t);
    src.SetZ3BWall (color, spin, qp_arg.t, box_arg.box_size, rand_num);

    if (GFixedSrc ()) {
      src.GFWallSource (AlgLattice (), spin, 3, qp_arg.t);
    } else {
      ERR.General (cname, fname, "Warning: 4D box src not gauge fixed");
    }
  }

//------------------------------------------------------------------
// Quark Propagator (Wilson type) with Random Source
//------------------------------------------------------------------
QPropWRand::QPropWRand (Lattice & lat, CommonArg * c_arg):QPropW (lat, c_arg) {

    char *fname = "QPropWRand(L&, ComArg*)";
    cname = "QPropWRand";
    VRB.Func (cname, fname);

    rsrc = NULL;
  }

// Allocate and Delete space for the random source
  void QPropWRand::AllocateRsrc ()
  {

    char *fname = "AllocateRsrc()";
    VRB.Func (cname, fname);

    if (rsrc == NULL) {
      int rsrc_size = 2 * GJP.VolNodeSites ();
      rsrc = (Float *) smalloc (rsrc_size * sizeof (Float));
      if (rsrc == 0)
	ERR.Pointer (cname, fname, "rsrc");
      VRB.Smalloc (cname, fname, "rsrc", rsrc, rsrc_size * sizeof (Float));
    }
  }
  void QPropWRand::DeleteRsrc ()
  {

    char *fname = "DeleteRsrc()";
    VRB.Func (cname, fname);

    if (rsrc != NULL) {
      VRB.Sfree (cname, fname, "rsrc", rsrc);
      sfree (rsrc);
      rsrc = NULL;
    }
  }

QPropWRand::QPropWRand (Lattice & lat, QPropWArg * arg, QPropWRandArg * r_arg, CommonArg * c_arg):QPropW (lat, arg, c_arg), rand_arg (*r_arg)
  {
    char *fname = "QPropWRand(L&, QPropWArg*, ComArg*)";
    cname = "QPropWRand";
    VRB.Func (cname, fname);

    rsrc = NULL;
    AllocateRsrc ();

    int rsrc_size = 2 * GJP.VolNodeSites ();

    if (rand_arg.rng == GAUSS) {
      // MGE 06/10/2008
      LRG.SetSigma (0.5);
      for (int i = 0; i < rsrc_size / 2; i++) {
	LRG.AssignGenerator (i);
	rsrc[2 * i] = LRG.Grand (FOUR_D);
	rsrc[2 * i + 1] = LRG.Grand (FOUR_D);
      }
      // END MGE 06/10/2008

    }
    if (rand_arg.rng == UONE) {
      // MGE 06/10/2008
      LRG.SetInterval (6.283185307179586, 0);
      for (int i = 0; i < rsrc_size / 2; i++) {
	LRG.AssignGenerator (i);
	Float theta (LRG.Urand (FOUR_D));
	rsrc[2 * i] = cos (theta);	// real part
	rsrc[2 * i + 1] = sin (theta);	// imaginary part
      }

      // END MGE 06/10/2008
    }
    if (rand_arg.rng == ZTWO) {
      // MGE 06/10/2008  
      LRG.SetInterval (1, -1);
      for (int i = 0; i < rsrc_size / 2; i++) {
	LRG.AssignGenerator (i);
	if (LRG.Urand (FOUR_D) > 0.) {
	  rsrc[2 * i] = 1.;
	} else {
	  rsrc[2 * i] = -1.;
	}
	rsrc[2 * i + 1] = 0.0;	// source is purely real
      }
      // END MGE 06/10/2008
    }
    if (rand_arg.rng == TEST) {	// deterministic source for testing
      for (int i = 0; i < rsrc_size / 2; i++) {
	Site s (i);
	int x = s.physX (), y = s.physY (), z = s.physZ (), t = s.physT ();
	int rnum =
	  ((937 * x * x + 3826 * x * y + 7034 * z * t * t) & 0x0100) ? 1 : -1;
	rsrc[2 * i] = rnum;
	rsrc[2 * i + 1] = 0.0;
      }
    }
  }

  QPropWRand::QPropWRand (const QPropWRand & rhs):QPropW (rhs), rsrc (NULL)
  {

    char *fname = "QPropW(const QPropW&)";
    cname = "QPropW";
    VRB.Func (cname, fname);

    AllocateRsrc ();
    for (int i = 0; i < 2 * GJP.VolNodeSites (); i++)
      rsrc[i] = rhs.rsrc[i];
  }

  Complex & QPropWRand::rand_src (int i) const
  {
    return ((Complex *) rsrc)[i];
  }

  QPropWRand & QPropWRand::operator= (const QPropWRand & rhs)
  {

    char *fname = "operator=(const QPropWRand& rhs)";
    VRB.Func (cname, fname);

    if (this != &rhs) {
      QPropW::operator= (rhs);	// This copies the QPropW stuff...

      AllocateRsrc ();
      for (int i = 0; i < 2 * GJP.VolNodeSites (); i++)
	rsrc[i] = rhs.rsrc[i];
    }

    return *this;
  }

  QPropWRand::~QPropWRand () {
    char *fname = "~QPropWRand()";
    VRB.Func (cname, fname);

    DeleteRsrc ();
  }

  void QPropWRand::ShiftPropForward (int n)
  {

    char *fname = "ShiftPropForward()";
    VRB.Func (cname, fname);

    QPropW::ShiftPropForward (n);

    Float *recv_buf;
    Float *send_buf;
    int len = 12 * 12 * 2;
    int len2 = 4;
    // size of transfers in words
    recv_buf = (Float *) smalloc (len * sizeof (Float));
    if (recv_buf == 0)
      ERR.Pointer (cname, fname, "recv_buf");
    VRB.Smalloc (cname, fname, "recv_buf", recv_buf, len * sizeof (Float));

    for (int j = 0; j < n; j++) {
      // shift 1 node in t-dir.  prop -> prop
      for (int i = 0; i < 2 * GJP.VolNodeSites (); i += len2) {
	send_buf = &rsrc[i];
	getMinusData ((IFloat *) recv_buf, (IFloat *) send_buf, len2, 3);
	moveMem ((IFloat *) & rsrc[i], (IFloat *) recv_buf,
		 len2 * sizeof (IFloat));
      }

    }

    VRB.Sfree (cname, fname, "recv_buf", recv_buf);
    sfree (recv_buf);
  }
  void QPropWRand::ShiftPropBackward (int n)
  {

    char *fname = "ShiftPropBackward()";
    VRB.Func (cname, fname);

    QPropW::ShiftPropBackward (n);

    Float *recv_buf;
    Float *send_buf;
    int len = 12 * 12 * 2;
    int len2 = 4;
    // size of transfer in words
    recv_buf = (Float *) smalloc (len * sizeof (Float));
    if (recv_buf == 0)
      ERR.Pointer (cname, fname, "recv_buf");
    VRB.Smalloc (cname, fname, "recv_buf", recv_buf, len * sizeof (Float));

    for (int j = 0; j < n; j++) {
      // shift 1 node in t-dir.  prop -> prop
      for (int i = 0; i < 2 * GJP.VolNodeSites (); i += len2) {
	send_buf = &rsrc[i];
	getPlusData ((IFloat *) recv_buf, (IFloat *) send_buf, len2, 3);
	moveMem ((IFloat *) & rsrc[i], (IFloat *) recv_buf,
		 len2 * sizeof (IFloat));
      }

    }

    VRB.Sfree (cname, fname, "recv_buf", recv_buf);
    sfree (recv_buf);
  }

// Restore prop
  void QPropWRand::RestoreQProp (char *name, int mid)
  {

    char *fname = "RestoreQProp()";
    VRB.Func (cname, fname);
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
   VRB.Flow(cname,fname,"Read rsrc from file %s\n", name);
  -------------------- Quarantine ends ---------------------------*/
  }
// Save prop
  void QPropWRand::SaveQProp (char *name, int mid)
  {
    char *fname = "SaveQProp()";
    VRB.Func (cname, fname);
   /****************************************************************
     The code below is temporarily isolated for purposes of merging
     with CPS main branch,  12/09/04, Oleg Loktik 
   -------------------- Quarantine starts --------------------------
  

  QPropW::SaveQProp(name,mid) ;
  
  unsigned int* data;
  data = (unsigned int*)rsrc;
  append_data(name, data, GJP.VolNodeSites() * sizeof(Complex));
  VRB.Flow(cname,fname,"Appended rsrc to file %s\n", name);
  DeleteRsrc();
  -------------------- Quarantine ends ---------------------------*/
  }

//------------------------------------------------------------------
// Quark Propagator (Wilson type) with Random Wall Source
//------------------------------------------------------------------
QPropWRandWallSrc::QPropWRandWallSrc (Lattice & lat, CommonArg * c_arg):
  QPropWRand (lat, c_arg) {

    char *fname = "QPropWRandWallSrc(L&, ComArg*)";
    cname = "QPropWRandWallSrc";
    VRB.Func (cname, fname);
  }
  QPropWRandWallSrc::QPropWRandWallSrc (Lattice & lat, QPropWArg * arg,
					QPropWRandArg * r_arg,
					CommonArg * c_arg)
:  QPropWRand (lat, arg, r_arg, c_arg) {

    char *fname = "QPropWRandWallSrc(L&, ComArg*)";
    cname = "QPropWRandWallSrc";
    VRB.Func (cname, fname);

    Run ();
  }

  void QPropWRandWallSrc::SetSource (FermionVectorTp & src, int spin, int color)
  {

    char *fname = "SetSource()";
    VRB.Func (cname, fname);

    if (rsrc == NULL) {		//need random numbers may implemented later
      ERR.General (cname, fname, "No randrom numbers found!\n");
    }

    src.ZeroSource ();
    src.SetWallSource (color, spin, qp_arg.t, rsrc);
    if (GFixedSrc ())
      src.GFWallSource (AlgLattice (), spin, 3, qp_arg.t);
  }


//------------------------------------------------------------------
// Quark Propagator (Wilson type) with Random Volume Source
//------------------------------------------------------------------
QPropWRandVolSrc::QPropWRandVolSrc (Lattice & lat, CommonArg * c_arg):
  QPropWRand (lat, c_arg) {

    char *fname = "QPropWRandVolSrc(L&, ComArg*)";
    cname = "QPropWRandVolSrc";
    VRB.Func (cname, fname);
  }
  QPropWRandVolSrc::QPropWRandVolSrc (Lattice & lat, QPropWArg * arg,
				      QPropWRandArg * r_arg, CommonArg * c_arg)
:  QPropWRand (lat, arg, r_arg, c_arg) {

    char *fname = "QPropWRandVolSrc(L&, ComArg*)";
    cname = "QPropWRandVolSrc";
    VRB.Func (cname, fname);

    Run ();
  }

//set the random source
/*void QPropWRandVolSrc::SetSource(FermionVectorTp& src, int spin, int color) {

  char *fname = "SetSource()"; 
  VRB.Func(cname, fname);

  if (rsrc==NULL) {//need random numbers may implemented later
    ERR.General(cname,fname,"No randrom numbers found!\n") ;
  }
  
  src.SetVolSource(color, spin, rsrc);
  if (GFixedSrc()) 
    for (int t=0;t<GJP.Tnodes()*GJP.TnodeSites(); t++)
      src.GFWallSource(AlgLattice(), spin, 3, t);
}*/

//set the random source
  void QPropWRandVolSrc::SetSource (FermionVectorTp & src, int spin, int color)
  {

    char *fname = "SetSource()";
    VRB.Func (cname, fname);

    if (rsrc == NULL) {		//need random numbers may implemented later
      ERR.General (cname, fname, "No randrom numbers found!\n");
    }

    Lattice & lat = AlgLattice ();
    src.SetVolSource (color, spin, rsrc);
    if (GFixedSrc ()) {
      if (lat.FixGaugeKind () == FIX_GAUGE_COULOMB_T)
	for (int t = 0; t < GJP.Tnodes () * GJP.TnodeSites (); t++)
	  src.GFWallSource (lat, spin, 3, t);
      else if (lat.FixGaugeKind () == FIX_GAUGE_LANDAU)
	src.LandauGaugeFixSrc (lat, spin);
      else
	ERR.General (cname, fname,
		     "gauge fixing method does not work for qpropwrandvolSrc\n");
    }
  }



//------------------------------------------------------------------
// Quark Propagator (Wilson type) with Random Slab Source
//------------------------------------------------------------------
QPropWRandSlabSrc::QPropWRandSlabSrc (Lattice & lat, CommonArg * c_arg):
  QPropWRand (lat, c_arg) {

    char *fname = "QPropWRandSlabSrc(L&, ComArg*)";
    cname = "QPropWRandSlabSrc";
    VRB.Func (cname, fname);
  }

  QPropWRandSlabSrc::QPropWRandSlabSrc (Lattice & lat, QPropWArg * arg,
					QPropWSlabArg * s_arg,
					CommonArg * c_arg)
:  QPropWRand (lat, arg, &(s_arg->rand_arg), c_arg),
    slab_width (s_arg->slab_width) {

    char *fname = "QPropWRandSlabSrc(L&, ComArg*)";
    cname = "QPropWRandSlabSrc";
    VRB.Func (cname, fname);

    Run ();
  }

QPropWRandSlabSrc::QPropWRandSlabSrc (Lattice & lat, QPropWArg * arg, Float * src, CommonArg * c_arg):
  QPropWRand (lat,
	      c_arg) {

    char *fname = "QPropWRandSlabSrc(L&, QPropWArg*, Float*, CommonArg*)";
    cname = "QPropWRandSlabSrc";
    VRB.Func (cname, fname);

    for (int i = 0; i < 2 * GJP.VolNodeSites (); i++)
      rsrc[i] = src[i];

    Run ();
  }

  void QPropWRandSlabSrc::SetSource (FermionVectorTp & src, int spin, int color)
  {

    char *fname = "SetSource()";
    VRB.Func (cname, fname);

    if (rsrc == NULL) {		//need random numbers may implemented later
      ERR.General (cname, fname, "No randrom numbers found!\n");
    }
    src.ZeroSource ();

    for (int t = qp_arg.t; t < qp_arg.t + slab_width; t++) {
      src.SetWallSource (color, spin, t, rsrc);
      if (GFixedSrc ())
	src.GFWallSource (AlgLattice (), spin, 3, t);
    }
  }

//------------------------------------------------------------------
// Quark Propagator (Wilson type) with Sequential Sources
//------------------------------------------------------------------
  QPropWSeq::QPropWSeq (Lattice & lat, QPropW & q, int *p,
			QPropWArg * q_arg, CommonArg * c_arg):QPropW (lat,
								      q_arg,
								      c_arg),
    quark (q), mom (p)
  {

    char *fname = "QPropWSeq(L&, ComArg*)";
    cname = "QPropWSeq";
    VRB.Func (cname, fname);

    // Stores the quark mass of the source propagator
    // Needed by Yasumichi's not degenerate mass runs
    quark_mass = quark.Mass ();

    //if the QPropW used for constructing the sequential source
    //propagator is done using HalfFerion set the DoHalfFermion 
    //in case the user forgot to do so.
    if (q.DoHalfFermion ())
      qp_arg.do_half_fermion = 1;
  }


  QPropWSeqMesSrc::QPropWSeqMesSrc (Lattice & lat, QPropW & quark, int *p,
				    int g, QPropWArg * q_arg,
				    CommonArg * c_arg):QPropWSeq (lat, quark, p,
								  q_arg, c_arg),
    gamma (g)
  {

    char *fname = "QPropWSeq(L&,...)";
    cname = "QPropWSeq";
    VRB.Func (cname, fname);

    Run ();
  }

  void QPropWSeqMesSrc::SetSource (FermionVectorTp & src, int spin, int color)
  {

    char *fname = "SetSource()";
    cname = "QPropWSeqMesSrc";
    VRB.Func (cname, fname);

    if (color < 0 || color >= GJP.Colors ())
      ERR.General (cname, fname,
		   "Color index out of range: color = %d\n", color);

    if (spin < 0 || spin > 3)
      ERR.General (cname, fname, "Spin index out of range: spin = %d\n", spin);

    src.ZeroSource ();
    Site s;
    WilsonMatrix tmp;
    for (s.Begin (); s.End (); s.nextSite ())
      if (qp_arg.t == s.physT ()) {
	tmp = quark[s.Index ()];
	tmp.gl (gamma);		// Multiply by the meson operator
	// multiply by the mommentum factor exp(ipx)
	Complex tt (conj (mom.Fact (s)));
	tmp *= tt;
	src.CopyWilsonMatSink (s.Index (), spin, color, tmp);
      }
    // Gauge fix the source. If QPropW sink is gauge fixed
    // this has to be done!
    if (quark.GFixedSnk ()) {
      for (int ss = 0; ss < 4; ss++)
	src.GFWallSource (AlgLattice (), ss, 3, qp_arg.t);
    }
  }

  QPropWSeqBar::QPropWSeqBar (Lattice & lat, QPropW & quark, int *p,
			      ProjectType pp, QPropWArg * q_arg,
			      CommonArg * c_arg):QPropWSeq (lat, quark, p,
							    q_arg, c_arg),
    proj (pp)
  {

    char *fname = "QPropWSeqBar(L&,...)";
    cname = "QPropWSeq";
    VRB.Func (cname, fname);
  }


  QPropWSeqProtDSrc::QPropWSeqProtDSrc (Lattice & lat, QPropW & quark, int *p,
					ProjectType pp, QPropWArg * q_arg,
					QPropWGaussArg * g_arg,
					CommonArg * c_arg):QPropWSeqBar (lat,
									 quark,
									 p, pp,
									 q_arg,
									 c_arg),
    gauss_arg (*g_arg)
  {

    char *fname = "QPropWSeqProtDSrc(L&,...)";
    cname = "QPropWSeq";
    VRB.Func (cname, fname);

    Run ();

    //Multiply by gamma5 and take the dagger to make it in to quark.
    Site s;
    for (s.Begin (); s.End (); s.nextSite ()) {
      QPropW::operator[](s.Index ()).gl (-5);
      QPropW::operator[](s.Index ()).hconj ();
    }
  }

  QPropWSeqProtDSrc::QPropWSeqProtDSrc (Lattice & lat, QPropW & quark, int *p,
					ProjectType pp, QPropWArg * q_arg,
					QPropWGaussArg * g_arg,
					CommonArg * c_arg,
					char *dummy):QPropWSeqBar (lat, quark,
								   p, pp, q_arg,
								   c_arg),
    gauss_arg (*g_arg)
  {
    // char *fname = "QPropW(prop&)";
    // cname = "QPropWGaussSrc";
    // VRB.Func(cname, fname);
  }

  void QPropWSeqProtDSrc::SetSource (FermionVectorTp & src, int spin, int color)
  {

    char *fname = "SetSource()";
    VRB.Func (cname, fname);

    if (color < 0 || color >= GJP.Colors ())
      ERR.General (cname, fname,
		   "Color index out of range: color = %d\n", color);

    if (spin < 0 || spin > 3)
      ERR.General (cname, fname, "Spin index out of range: spin = %d\n", spin);

    src.ZeroSource ();
    Site s;
    Diquark diq;
    WilsonVector S;

    WilsonMatrix q;
    WilsonMatrix OqO;
    WilsonMatrix qO;
    WilsonMatrix Oq;
    for (s.Begin (); s.End (); s.nextSite ())
      for (int nt = 0; nt < gauss_arg.nt; nt++) {
	if (gauss_arg.mt[nt] == s.physT ()) {
	  int i = s.Index ();
	  q = quark[i];
	  // If DoHalfFermion is on we have non-relativistic sources
	  // Multiply the sink by 1/2(1+gamma_t) when we do Half Spinors
	  // By doing this we implement the non-relativistic sink
	  if (DoHalfFermion ())
	    q.PParProjectSink ();
	  Oq = qO = q;
	  // multiply C*gamma_5 left  C is the charge conjugation
	  Oq.ccl (5);
	  // multiply C*gamma_5 right C is the charge conjugation
	  qO.ccr (5);
	  OqO = Oq;
	  // multiply C*gamma_5 right C is the charge conjugation
	  // Shoichi's code misses a minus sign here (or in ccl)
	  OqO.ccr (5);
	  // spin is denoted as delta in notes
	  // color is denoted as d in notes
	  diq.D_diquark (OqO, q, Oq, qO, spin, color);
	  diq.Project (S, proj);

	  //multiply by the momentum factor exp(ipx)
	  Complex tt (conj (mom.Fact (s)));
	  S *= tt;

	  S.conj ();

	  //if non-relativistic sink is needed multiply  by 1/2(1+gamma_t)
	  if (DoHalfFermion ())
	    S.PParProject ();
	  //Note the order first project then multiply by gamma5 
	  //multily by gamma5
	  S.gamma (-5);

	  src.CopyWilsonVec (i, S);
	}
      }				// nt loop

    // Gauge fix the source. If quark sink is gauge fixed
    // this has to be done!
    if (quark.GFixedSnk ()) {
      for (int ss = 0; ss < 4; ss++)
	src.GFWallSource (AlgLattice (), ss, 3, qp_arg.t);
    }
    if (quark.SeqSmearSink () == GAUSS_GAUGE_INV) {
      DoLinkSmear (gauss_arg);	// YA: smear link if needed
      for (int nt = 0; nt < gauss_arg.nt; nt++) {
	for (int ss = 0; ss < 4; ss++)
	  src.GaussianSmearVector (AlgLattice (), ss,
				   quark.GaussArg ().gauss_N,
				   quark.GaussArg ().gauss_W, gauss_arg.mt[nt]);
	//source sink smearing is the same
      }
      UndoLinkSmear (gauss_arg);	// YA: get back the original link
    }
  }

  QPropWSeqProtUSrc::QPropWSeqProtUSrc (Lattice & lat, QPropW & quark, int *p,
					ProjectType pp, QPropWArg * q_arg,
					QPropWGaussArg * g_arg,
					CommonArg * c_arg):QPropWSeqBar (lat,
									 quark,
									 p, pp,
									 q_arg,
									 c_arg),
    gauss_arg (*g_arg)
  {

    char *fname = "QPropWSeqProtUSrc(L&,...)";
    cname = "QPropWSeq";
    VRB.Func (cname, fname);

    Run ();

    //Multiply by gamma5 and take the dagger to make it in to quark.
    Site s;
    for (s.Begin (); s.End (); s.nextSite ()) {
      QPropW::operator[](s.Index ()).gl (-5);
      QPropW::operator[](s.Index ()).hconj ();
    }
  }

  QPropWSeqProtUSrc::QPropWSeqProtUSrc (Lattice & lat, QPropW & quark, int *p,
					ProjectType pp, QPropWArg * q_arg,
					QPropWGaussArg * g_arg,
					CommonArg * c_arg,
					char *dummy):QPropWSeqBar (lat, quark,
								   p, pp, q_arg,
								   c_arg),
    gauss_arg (*g_arg)
  {
    // char *fname = "QPropW(prop&)";
    // cname = "QPropWGaussSrc";
    // VRB.Func(cname, fname);
  }

  void QPropWSeqProtUSrc::SetSource (FermionVectorTp & src, int spin, int color)
  {

    char *fname = "SetSource()";
    VRB.Func (cname, fname);

    if (color < 0 || color >= GJP.Colors ())
      ERR.General (cname, fname,
		   "Color index out of range: color = %d\n", color);

    if (spin < 0 || spin > 3)
      ERR.General (cname, fname, "Spin index out of range: spin = %d\n", spin);

    src.ZeroSource ();
    Site s;
    Diquark diq;
    WilsonVector S;
    WilsonMatrix q;
    WilsonMatrix OqO;

    for (s.Begin (); s.End (); s.nextSite ())
      for (int nt = 0; nt < gauss_arg.nt; nt++) {
	if (gauss_arg.mt[nt] == s.physT ()) {
	  int i (s.Index ());
	  q = quark[i];
	  // If DoHalfFermion is on we have non-relativistic sources
	  // Multiply the sink by 1/2(1+gamma_t) when we do Half Spinors
	  // By doing this we implement the non-relativistic sink
	  if (DoHalfFermion ())
	    q.PParProjectSink ();
	  OqO = q;

	  // multiply C*gamma_5 left  C is the charge conjugation
	  OqO.ccl (5);
	  // multiply C*gamma_5 right C is the charge conjugation
	  OqO.ccr (5);
	  // spin is denoted as delta in notes
	  // color is denoted as d in notes

	  diq.U_diquark (OqO, q, spin, color);
	  diq.Project (S, proj);

	  //multiply by the momentum factor exp(ipx)
	  Complex tt (conj (mom.Fact (s)));
	  S *= tt;

	  S.conj ();

	  //if non-relativistic sink is needed multiply  by 1/2(1+gamma_t)
	  if (DoHalfFermion ())
	    S.PParProject ();
	  //Note the order first project then multiply by gamma5 
	  //multily by gamma5
	  S.gamma (-5);

	  src.CopyWilsonVec (i, S);
	}			// Loop over sites
      }				// nt loop

    // Gauge fix the source. If quark sink is gauge fixed
    // this has to be done!
    if (quark.GFixedSnk ()) {
      for (int ss = 0; ss < 4; ss++)
	src.GFWallSource (AlgLattice (), ss, 3, qp_arg.t);
    }
    if (quark.SeqSmearSink () == GAUSS_GAUGE_INV) {
      DoLinkSmear (gauss_arg);	// YA: smear link if needed
      for (int nt = 0; nt < gauss_arg.nt; nt++) {
	for (int ss = 0; ss < 4; ss++)
	  src.GaussianSmearVector (AlgLattice (), ss,
				   quark.GaussArg ().gauss_N,
				   quark.GaussArg ().gauss_W, gauss_arg.mt[nt]);
	//source sink smearing is the same
      }
      UndoLinkSmear (gauss_arg);	// YA: get back the original link
    }
  }



//------------------------------------------------------------------
// Quark Propagator (Wilson type) with Volume Momentum Source
//------------------------------------------------------------------
QPropWVolMomSrc::QPropWVolMomSrc (Lattice & lat, CommonArg * c_arg):QPropW (lat,
	  c_arg)
  {

    char *fname = "QPropWVolMomSrc(L&, ComArg*)";
    cname = "QPropWVolMomSrc";
    VRB.Func (cname, fname);
  }
  QPropWVolMomSrc::QPropWVolMomSrc (Lattice & lat, QPropWArg * arg,
				    int *p, CommonArg * c_arg):QPropW (lat,
								       c_arg),
    mom (p)
  {

    char *fname = "QPropWVolMomSrc(L&, QPropWArg*, ComArg*)";
    cname = "QPropWVolMomSrc";
    VRB.Func (cname, fname);
    qp_arg = *arg;
    Run ();
  }

  QPropWVolMomSrc::QPropWVolMomSrc (const QPropWVolMomSrc & rhs):QPropW (rhs),
    mom (rhs.mom)
  {

    char *fname = "QPropWVolMomSrc(const QPropW&)";
    cname = "QPropWVolMomSrc";
    VRB.Func (cname, fname);
  }

  void QPropWVolMomSrc::SetSource (FermionVectorTp & src, int spin, int color)
  {
    char *fname = "SetSource()";
    VRB.Func (cname, fname);

    src.ZeroSource ();
    src.SetVolMomSource (color, spin, mom, qp_arg.flavor);
    if (GFixedSrc ())
      src.GaugeFixVector (AlgLattice (), spin);
  }




  //-----------------------------------------------------------------
  // TY Add Start
  int QPropW::siteOffset (const int lcl[], const int lcl_sites[]) const
  {
    int l = 4 - 1;
    int offset = lcl[l];
    while (l-- > 0)
    {
      offset *= lcl_sites[l];
      offset += lcl[l];
    }
    return offset;
  }
  // TY Add End
  //-----------------------------------------------------------------

  int QPropW::BoxSrcStart () const
  {
    ERR.NotImplemented (cname, "BoxSrcStart()");
    return 0;
  }
  int QPropW::BoxSrcEnd () const
  {
    ERR.NotImplemented (cname, "BoxSrcEnd()");
    return 0;
  }

  static QPropWGaussArg dummy_arg;
  const QPropWGaussArg & QPropW::GaussArg (void)
  {
    ERR.NotImplemented (cname, "GaussArg()");
    //dummy code to get rid or warnings
    return dummy_arg;
  }
  int QPropW::Gauss_N () const
  {
    ERR.NotImplemented (cname, "Gauss_N()");
    return 0;
  }
  Float QPropW::Gauss_W () const
  {
    ERR.NotImplemented (cname, "Gauss_W()");
    return 0;
  }

//================================================================
// Added by Meifeng Lin to do multiple sequential propagators
//================================================================

//------------------------------------------------------------------
// Quark Propagator (Wilson type) with Multiple Sequential Sources
// Added by Meifeng Lin, 10/18/2010
//------------------------------------------------------------------
  QPropWMultSeq::QPropWMultSeq (Lattice & lat, int N, QPropW ** q, int *p,
				QPropWArg * q_arg,
				CommonArg * c_arg):QPropW (lat, q_arg, c_arg),
    n_mult (N), quark (q), mom (p)
  {

    char *fname = "QPropWMultSeq(L&, ComArg*)";
    cname = "QPropWMultSeq";
    VRB.Func (cname, fname);

    // Stores the quark mass of the source propagator
    // Needed by Yasumichi's not degenerate mass runs
    quark_mass = quark[0]->Mass ();

    //if the QPropW used for constructing the sequential source
    //propagator is done using HalfFerion set the DoHalfFermion 
    //in case the user forgot to do so.
    if (quark[0]->DoHalfFermion ())
      qp_arg.do_half_fermion = 1;
  }



//-----------------------------------------------------------------------------
// Multi-sequential source propagator. -- Meifeng Lin, 10/22/2010
//-----------------------------------------------------------------------------
  QPropWMultSeqBar::QPropWMultSeqBar (Lattice & lat, int N, QPropW ** quark,
				      int *p, ProjectType pp, QPropWArg * q_arg,
				      CommonArg * c_arg):QPropWMultSeq (lat, N,
									quark,
									p,
									q_arg,
									c_arg),
    proj (pp)
  {

    char *fname = "QPropWMultSeqBar(L&,...)";
    cname = "QPropWMultSeq";
    VRB.Func (cname, fname);
  }

//-----------------------------------------------------------------------------
// Multi-sequential source propagator for the d quark. -- Meifeng Lin, 10/22/2010
//-----------------------------------------------------------------------------
  QPropWMultSeqProtDSrc::QPropWMultSeqProtDSrc (Lattice & lat, int N,
						QPropW ** quark, int *p,
						ProjectType pp,
						QPropWArg * q_arg,
						QPropWGaussArg * g_arg,
						CommonArg * c_arg,
						int *t):QPropWMultSeqBar (lat,
									  N,
									  quark,
									  p, pp,
									  q_arg,
									  c_arg),
    gauss_arg (*g_arg), time (t)
  {

    char *fname = "QPropWMultSeqProtDSrc(L&,...)";
    cname = "QPropWMultSeq";
    VRB.Func (cname, fname);

    Run ();

    //Multiply by gamma5 and take the dagger to make it in to quark.
    Site s;
    for (s.Begin (); s.End (); s.nextSite ()) {
      QPropW::operator[](s.Index ()).gl (-5);
      QPropW::operator[](s.Index ()).hconj ();
    }
  }

//-----------------------------------------------------------------------------
// Multi-sequential source propagator for the d quark. -- Meifeng Lin, 10/22/2010
// A dummy function used to load quark propagators. 
//-----------------------------------------------------------------------------
  QPropWMultSeqProtDSrc::QPropWMultSeqProtDSrc (Lattice & lat, int N,
						QPropW ** quark, int *p,
						ProjectType pp,
						QPropWArg * q_arg,
						QPropWGaussArg * g_arg,
						CommonArg * c_arg,
						char
						*dummy):QPropWMultSeqBar (lat,
									  N,
									  quark,
									  p, pp,
									  q_arg,
									  c_arg),
    gauss_arg (*g_arg)
  {
    //   char *fname = "QPropWMultSeqProtDSrc(Lattice&, int,  QPropW**, int *, ProjectType, QPropWArg*, QPropWGaussArg *, CommonArg*, char*)";
    //   VRB.Func(cname, fname);
  }

//--------------------------------------------------------------------------
// Source for the multiple-sequential propagator. 
// The source is the superposition of individual sequential d-quark sources.
// The only difference between this and QPropWSeqProtDSrc::SetSource is that
// the input is multiple quark propagators, hence
// q = quark[i] is replaced by 
// q= quark[nt]->operator[](i)
// Perhaps there is an easier way to implement this? 
// -- Meifeng Lin, 10/22/2010
//--------------------------------------------------------------------------
  void QPropWMultSeqProtDSrc::SetSource (FermionVectorTp & src, int spin,
					 int color)
  {

    char *fname = "SetSource()";
    VRB.Func (cname, fname);

    if (color < 0 || color >= GJP.Colors ())
      ERR.General (cname, fname,
		   "Color index out of range: color = %d\n", color);

    if (spin < 0 || spin > 3)
      ERR.General (cname, fname, "Spin index out of range: spin = %d\n", spin);

    src.ZeroSource ();
    Site s;
    Diquark diq;
    WilsonVector S;

    WilsonMatrix q;
    WilsonMatrix OqO;
    WilsonMatrix qO;
    WilsonMatrix Oq;
    for (int n = 0; n < n_mult; n++) {
      for (s.Begin (); s.End (); s.nextSite ())
	for (int nt = 0; nt < gauss_arg.nt; nt++) {
	  if (gauss_arg.mt[nt] == s.physT () || time[n] == s.physT ()) {
	    int i = s.Index ();
	    q = quark[n]->operator[](i);

	    // If DoHalfFermion is on we have non-relativistic sources
	    // Multiply the sink by 1/2(1+gamma_t) when we do Half Spinors
	    // By doing this we implement the non-relativistic sink
	    if (DoHalfFermion ())
	      q.PParProjectSink ();
	    Oq = qO = q;
	    // multiply C*gamma_5 left  C is the charge conjugation
	    Oq.ccl (5);
	    // multiply C*gamma_5 right C is the charge conjugation
	    qO.ccr (5);
	    OqO = Oq;
	    // multiply C*gamma_5 right C is the charge conjugation
	    // Shoichi's code misses a minus sign here (or in ccl)
	    OqO.ccr (5);
	    // spin is denoted as delta in notes
	    // color is denoted as d in notes
	    diq.D_diquark (OqO, q, Oq, qO, spin, color);
	    diq.Project (S, proj);

	    //multiply by the momentum factor exp(ipx)
	    Complex tt (conj (mom.Fact (s)));
	    S *= tt;

	    S.conj ();

	    //if non-relativistic sink is needed multiply  by 1/2(1+gamma_t)
	    if (DoHalfFermion ())
	      S.PParProject ();
	    //Note the order first project then multiply by gamma5 
	    //multily by gamma5
	    S.gamma (-5);

	    src.CopyWilsonVec (i, S);
	  }
	}			// nt loop

      // Gauge fix the source. If quark sink is gauge fixed
      // this has to be done!
      if (quark[n]->GFixedSnk ()) {
	for (int ss = 0; ss < 4; ss++)
	  src.GFWallSource (AlgLattice (), ss, 3, time[n]);
      }
      if (quark[n]->SeqSmearSink () == GAUSS_GAUGE_INV) {
	DoLinkSmear (gauss_arg);	// YA: smear link if needed
	for (int nt = 0; nt < gauss_arg.nt; nt++) {
	  for (int ss = 0; ss < 4; ss++)
	    src.GaussianSmearVector (AlgLattice (), ss,
				     quark[n]->GaussArg ().gauss_N,
				     quark[n]->GaussArg ().gauss_W, time[n]);
	  //source sink smearing is the same
	}
	UndoLinkSmear (gauss_arg);	// YA: get back the original link
      }
    }				// loop over forward propagators
  }

//--------------------------------------------------------------------------------
// Multi-sequential source propagator for the u quark. -- Meifeng Lin, 10/22/2010
//--------------------------------------------------------------------------------
  QPropWMultSeqProtUSrc::QPropWMultSeqProtUSrc (Lattice & lat, int N,
						QPropW ** quark, int *p,
						ProjectType pp,
						QPropWArg * q_arg,
						QPropWGaussArg * g_arg,
						CommonArg * c_arg,
						int *t):QPropWMultSeqBar (lat,
									  N,
									  quark,
									  p, pp,
									  q_arg,
									  c_arg),
    gauss_arg (*g_arg), time (t)
  {

    char *fname = "QPropWSeqProtUSrc(L&,...)";
    cname = "QPropWSeq";
    VRB.Func (cname, fname);

    Run ();

    //Multiply by gamma5 and take the dagger to make it in to quark.
    Site s;
    for (s.Begin (); s.End (); s.nextSite ()) {
      QPropW::operator[](s.Index ()).gl (-5);
      QPropW::operator[](s.Index ()).hconj ();
    }
  }

  QPropWMultSeqProtUSrc::QPropWMultSeqProtUSrc (Lattice & lat, int N,
						QPropW ** quark, int *p,
						ProjectType pp,
						QPropWArg * q_arg,
						QPropWGaussArg * g_arg,
						CommonArg * c_arg,
						char
						*dummy):QPropWMultSeqBar (lat,
									  N,
									  quark,
									  p, pp,
									  q_arg,
									  c_arg),
    gauss_arg (*g_arg)
  {
    //   char *fname = "QPropWMultSeqProtUSrc(Lattice&, int,  QPropW**, int *, ProjectType, QPropWArg*, QPropWGaussArg *, CommonArg*, char*)";
    //   VRB.Func(cname, fname);
  }


//--------------------------------------------------------------------------
// Source for the multiple-sequential propagator. 
// The source is the superposition of individual sequential d-quark sources.
// The only difference between this and QPropWSeqProtUSrc::SetSource is that
// the input is multiple quark propagators, hence
// q = quark[i] is replaced by 
// q= quark[nt]->operator[](i)
// Perhaps there is an easier way to implement this? 
// -- Meifeng Lin, 10/22/2010
//--------------------------------------------------------------------------
  void QPropWMultSeqProtUSrc::SetSource (FermionVectorTp & src, int spin,
					 int color)
  {

    char *fname = "SetSource()";
    VRB.Func (cname, fname);

    if (color < 0 || color >= GJP.Colors ())
      ERR.General (cname, fname,
		   "Color index out of range: color = %d\n", color);

    if (spin < 0 || spin > 3)
      ERR.General (cname, fname, "Spin index out of range: spin = %d\n", spin);

    src.ZeroSource ();
    Site s;
    Diquark diq;
    WilsonVector S;
    WilsonMatrix q;
    WilsonMatrix OqO;

    for (int n = 0; n < n_mult; n++) {
      for (s.Begin (); s.End (); s.nextSite ())
	for (int nt = 0; nt < gauss_arg.nt; nt++) {
	  if (gauss_arg.mt[nt] == s.physT () || time[n] == s.physT ()) {
	    int i (s.Index ());
	    q = quark[n]->operator[](i);
	    // If DoHalfFermion is on we have non-relativistic sources
	    // Multiply the sink by 1/2(1+gamma_t) when we do Half Spinors
	    // By doing this we implement the non-relativistic sink
	    if (DoHalfFermion ())
	      q.PParProjectSink ();
	    OqO = q;

	    // multiply C*gamma_5 left  C is the charge conjugation
	    OqO.ccl (5);
	    // multiply C*gamma_5 right C is the charge conjugation
	    OqO.ccr (5);
	    // spin is denoted as delta in notes
	    // color is denoted as d in notes

	    diq.U_diquark (OqO, q, spin, color);
	    diq.Project (S, proj);

	    //multiply by the momentum factor exp(ipx)
	    Complex tt (conj (mom.Fact (s)));
	    S *= tt;

	    S.conj ();

	    //if non-relativistic sink is needed multiply  by 1/2(1+gamma_t)
	    if (DoHalfFermion ())
	      S.PParProject ();
	    //Note the order first project then multiply by gamma5 
	    //multily by gamma5
	    S.gamma (-5);

	    src.CopyWilsonVec (i, S);
	  }			// Loop over sites
	}			// nt loop

      // Gauge fix the source. If quark sink is gauge fixed
      // this has to be done!
      if (quark[n]->GFixedSnk ()) {
	for (int ss = 0; ss < 4; ss++)
	  src.GFWallSource (AlgLattice (), ss, 3, time[n]);
      }
      if (quark[n]->SeqSmearSink () == GAUSS_GAUGE_INV) {
	DoLinkSmear (gauss_arg);	// YA: smear link if needed
	for (int nt = 0; nt < gauss_arg.nt; nt++) {
	  for (int ss = 0; ss < 4; ss++)
	    src.GaussianSmearVector (AlgLattice (), ss,
				     quark[n]->GaussArg ().gauss_N,
				     quark[n]->GaussArg ().gauss_W, time[n]);
	  //source sink smearing is the same
	}
	UndoLinkSmear (gauss_arg);	// YA: get back the original link
      }
    }				//loop over forward propagators
  }


CPS_END_NAMESPACE
