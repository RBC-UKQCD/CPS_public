#include <config.h>
#include <stdlib.h>		// exit()
#include <stdio.h>
#include <string.h>
#include <alg/common_arg.h>
#include <comms/glb.h>
#include <comms/scu.h>
// #include <util/data_io.h>
#include <comms/sysfunc_cps.h>

#include <fcntl.h>		// read and write control flags,
#include <unistd.h>		// close(). These are needed for io parts to
			// compile on PCs
#include <alg/qpropw.h>
#include <util/qcdio.h>

#include <util/qioarg.h>
#include <util/sproj_tr.h>
#include <util/time_cps.h>

//YA
#include <alg/alg_plaq.h>
#include <alg/alg_smear.h>
#include <alg/no_arg.h>

#define VOLFMT QIO_VOLFMT

CPS_START_NAMESPACE
// Restore prop
  void
QPropW::RestoreQProp (char *name, int mid)
{

  char *fname = "RestoreQProp()";

  VRB.Func (cname, fname);
  VRB.Result (cname, fname, "name=%s mid=%d\n", name, mid);
  int float_size = sizeof (WilsonMatrix) * GJP.VolNodeSites () / sizeof (Float);

  if (prop == NULL)
    Allocate (PROP);
  if (mid)
    Allocate (MIDPROP);

  int sc_part_file_exist = CheckSCfile (name);

  // we need to store the source
  Float *read_source = (Float *) smalloc (cname, fname, "read_source",
					  GJP.VolNodeSites () * 288 *
					  sizeof (Float));
	src_hyper=0;

#ifdef USE_QIO

  if (sc_part_file_exist) {
    //        for(int spin=qpropw_arg.StartSrcSpin;spin<qpropw_arg.EndSrcSpin;spin++){
    //          for(int color=qpropw_arg.StartSrcColor;color<qpropw_arg.EndSrcColor;color++){
    for (int spin = 0; spin < 4; spin++) {
      for (int color = 0; color < 3; color++) {
	char file[256];

	snprintf (file, 256, "%ss%dc%d", name, spin, color);
//        ReLoadSC(file,spin,color);
	qio_readPropagator readPropQio (file, spin, color, &prop[0],
					read_source, float_size,
					float_size, VOLFMT);
		src_hyper=readPropQio.hyper_n;
		if (src_hyper)
		for(int i = 0;i<src_hyper;i++){
			hyper_lower[i] = readPropQio.hyper_lower[i];
			hyper_upper[i] = readPropQio.hyper_upper[i];
		}
      }
    }
  } else {
    qio_readPropagator readPropQio (qp_arg.file, QIO_FULL_SOURCE,
				    &prop[0], read_source, GJP.argc (),
				    GJP.argv (), VOLFMT);
		src_hyper=readPropQio.hyper_n;
		if (src_hyper)
		for(int i = 0;i<src_hyper;i++){
			hyper_lower[i] = readPropQio.hyper_lower[i];
			hyper_upper[i] = readPropQio.hyper_upper[i];
		}
  }
#else
  ERR.General(cname,fname,"Needs QIO to read in propagators");
#endif // USE_QIO


  //Get the lattice form the Alg base class
  Lattice & Lat = this->AlgLattice ();

  if (Lat.Fclass () == F_CLASS_DWF) {
    Float renFac = 5. - GJP.DwfHeight ();

    for (int ii (0); ii < GJP.VolNodeSites (); ++ii)
      *(prop + ii) *= renFac;
  }

  int glb_walls = GJP.TnodeSites () * GJP.Tnodes ();

  if ((Lat.Fclass () == F_CLASS_DWF) && mid) {

    // Set the node size of the full (non-checkerboarded) fermion field
    //----------------------------------------------------------------
    int ls = GJP.SnodeSites ();
    int ls_glb = GJP.SnodeSites () * GJP.Snodes ();
    int f_size = GJP.VolNodeSites () * Lat.FsiteSize () / GJP.SnodeSites ();
    int f_size_5d = f_size * ls;

    FermionVectorTp sol;
    FermionVectorTp src;
    FermionVectorTp midsol;
    Vector *src_4d = (Vector *) src.data ();
    Vector *sol_4d = (Vector *) sol.data ();
    Vector *midsol_4d = (Vector *) midsol.data ();
    Vector *src_5d = (Vector *) smalloc (cname, fname, "src_5d",
					 f_size_5d * sizeof (IFloat));
    Vector *sol_5d = (Vector *) smalloc (cname, fname, "sol_5d",
					 f_size_5d * sizeof (IFloat));

    j5q_pion =
      (Float *) smalloc (cname, fname, "d_j5q_pion_p",
			 glb_walls * sizeof (Float));
    conserved =
      (Float *) smalloc (cname, fname, "d_conserved_p",
			 glb_walls * sizeof (Float));

    Float *flt_p = (Float *) j5q_pion;

    for (int i = 0; i < glb_walls; i++)
      *flt_p++ = 0.0;
    flt_p = (Float *) conserved;
    for (int i = 0; i < glb_walls; i++)
      *flt_p++ = 0.0;

    for (int spn = 0; spn < 4; spn++)
      for (int col = 0; col < 3; col++) {
	Float *src_4d_f = (Float *) src_4d;
	double src_norm = 0., src_temp;

	// store the source
	for (int index (0); index < GJP.VolNodeSites (); ++index)
	  for (int mm (0); mm < 4; ++mm)
	    for (int cc (0); cc < GJP.Colors (); ++cc) {
	      // now same ordering as propagator [volume][spin][color][solution_spin][solution_color][ReIm]
//                 *(save_source + 288*index + 72*mm + 24*cc + 6*spn + 2*col)       = src[24*index + 6*mm + 2*cc];
//                 *(save_source + 288*index + 72*mm + 24*cc + 6*spn + 2*col + 1)   = src[24*index + 6*mm + 2*cc+1];
	      src_temp =
		src_4d_f[24 * index + 6 * mm + 2 * cc] =
		*(read_source + 288 * index + 72 * mm +
		  24 * cc + 6 * spn + 2 * col);
	      src_norm += src_temp * src_temp;
	      src_temp =
		src_4d_f[24 * index + 6 * mm + 2 * cc + 1] =
		*(read_source + 288 * index + 72 * mm +
		  24 * cc + 6 * spn + 2 * col + 1);
	      src_norm += src_temp * src_temp;
	    }
	VRB.Result (cname, fname, "spn=%d col=%d src_norm=%e\n", spn,
		    col, src_temp);
	Lat.Ffour2five (src_5d, src_4d, 0, ls_glb - 1);
	SaveRow (spn, col, sol, midsol);
	UnfixSol (sol);
	Lat.Fsolfour2five (sol_5d, sol_4d, src_5d, &(qp_arg.cg));


//      iter = Lat.FmatInv(sol_5d, src_5d, &(qp_arg.cg), &true_res, CNV_FRM_YES, PRESERVE_NO);


	//-----------------------------------------------------------------
	// TY Add Start
//    if(qp_arg.save_ls_prop) 
//       for(int nls(0);nls<GJP.SnodeSites();nls++) SaveQPropLs(sol_5d, qp_arg.file, nls);
//    spnclr_cnt++;

//    MeasConAxialOld(sol_5d);
	// TY Add End
	//-----------------------------------------------------------------

	//-----------------------------------------------------------------
	// M. Lightman
	{

	  MeasConAxialOld (sol_5d);

	  MeasJ5qPion (sol_5d);
	}
	// End M. Lightman
	//-----------------------------------------------------------------

	// prop on walls
	Lat.Ffive2four (sol_4d, sol_5d, ls_glb - 1, 0);
	FixSol (sol);
	// midpoint prop
	if (StoreMidprop ()) {
	  Lat.Ffive2four (midsol_4d, sol_5d, ls_glb / 2 - 1, ls_glb / 2);
	  FixSol (midsol);
	}
      }
    //-----------------------------------------------------------------
    // TY Add Start
    // Print out conserved axial results
//   int time_size = GJP.TnodeSites()*GJP.Tnodes();
    for (int t (0); t < glb_walls; t++)
      slice_sum ((Float *) & conserved[t], 1, 99);
    if (common_arg->results != 0) {
      FILE *fp;

      if ((fp = Fopen ((char *) common_arg->results, "a")) == NULL) {
	ERR.FileA (cname, fname, (char *) common_arg->results);
      }
      Fprintf (fp, "Conserved Axial w_spect\n");
      for (int t = 0; t < glb_walls; t++) {
	Fprintf (fp, "%d = %.16e\n", t, conserved[t]);
      }
      Fclose (fp);
    }
    // TY Add End
    //-----------------------------------------------------------------

    //-----------------------------------------------------------------
    // M. Lightman
    // Print out J5q Pion contraction
    for (int t (0); t < glb_walls; t++)
      slice_sum ((Float *) & j5q_pion[t], 1, 99);
    if (common_arg->results != 0) {
      FILE *fp1;

      if ((fp1 = Fopen ((char *) common_arg->results, "a")) == NULL) {
	ERR.FileA (cname, fname, (char *) common_arg->results);
      }
      Fprintf (fp1, "J5q Pion Contraction\n");
      for (int t = 0; t < glb_walls; t++) {
	Fprintf (fp1, "%d = %.16e\n", t, j5q_pion[t]);
      }
      Fclose (fp1);
    }
    // End M. Lightman
    //-----------------------------------------------------------------

    sfree (cname, fname, "d_j5q_pion_p", j5q_pion);
    sfree (cname, fname, "d_conserved_p", conserved);

    sfree (cname, fname, "sol_5d", sol_5d);
    sfree (cname, fname, "src_5d", src_5d);

  } else {
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
}

#if 1
  //-----------------------------------------------------------------
  // TY Add Start
// Save 5d prop at each ls
void
QPropW::SaveQPropLs (Vector * sol_5d, char *name, int ls)
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
					 (wilson_vector &) sol_5d[i
								  + skip_buf]);
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
#endif

// Restore 5d prop at ls
void
QPropW::RestoreQPropLs (char *name, int ls)
{

  char *fname = "RestoreQPropLs()";

  VRB.Func (cname, fname);

  // Flag set if sequential propagator 
  int seq_src = ((SrcType () == PROT_U_SEQ) ||
		 (SrcType () == PROT_D_SEQ) || (SrcType () == MESSEQ));

  if (qp_arg.save_ls_prop == 1) {
    FermionVectorTp sol;
    FermionVectorTp midsol;

    int fv_size = GJP.Colors () * 4 * 2 * sizeof (Float) * GJP.VolNodeSites ();
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
    int f_size = sizeof (WilsonMatrix) * GJP.VolNodeSites () / sizeof (Float);

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
void
QPropW::RestoreQPropLs_ftom (char *name)
{

  char *fname = "RestoreQPropLs()";

  VRB.Func (cname, fname);

  // Flag set if sequential propagator 
  int seq_src = ((SrcType () == PROT_U_SEQ) ||
		 (SrcType () == PROT_D_SEQ) || (SrcType () == MESSEQ));

  if (qp_arg.save_ls_prop == 1) {
    FermionVectorTp sol;
    FermionVectorTp midsol;

    int fv_size = GJP.Colors () * 4 * 2 * sizeof (Float) * GJP.VolNodeSites ();
    char sname[100];

    // Allocate 5d memory
    if (propls == NULL) {
      propls =
	(WilsonMatrix *) smalloc (GJP.VolNodeSites () *
				  GJP.SnodeSites () * sizeof (WilsonMatrix));
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
	      propls[s + shft_buf].load_row (spn, col, (wilson_vector &)
					     sol[i]);
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

void
QPropW::SaveQProp (char *name, int mid)
{

  char *fname = "SaveQProp()";

  VRB.Func (cname, fname);
  // save prop
//   if (do_cg && qp_arg.save_prop) {
  if (mid)
    ERR.General (cname, fname, "midpoint store not implemented");

  //Start timing how long it takes to save
  Float dtime_last = dclock ();
  WilsonMatrix *save_prop;

  char propType[256], sourceType[256];

//propOutfile[256];
  char *propOutfile = name;

  char gfixInfo[256];

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


  char fermionInfo[256];

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



  sprintf (propType,
	   "4D propagator, mass=%0.4f, StpCond=%0.0E,\nBC=%s%s%s%s,\n%s,\n%s",
	   qp_arg.cg.mass, qp_arg.cg.stop_rsd,
	   ((GJP.Xbc () == BND_CND_PRD) ? "P" : "A"),
	   ((GJP.Ybc () == BND_CND_PRD) ? "P" : "A"),
	   ((GJP.Zbc () == BND_CND_PRD) ? "P" : "A"),
	   ((GJP.Tbc () == BND_CND_PRD) ? "P" : "A"), gfixInfo, fermionInfo);

  //sprintf(sourceType, "fullSource");
  // be a bit more sophisticated
  sprintf (sourceType, "%s-source at t=%i", SourceType_map[SrcType ()].name,
	   SourceTime ());

//       sprintf(propOutfile,qp_arg.file);

  //in case of DWF, renormalize first
  if (AlgLattice ().Fclass () == F_CLASS_DWF) {

    Float renFac = 1. / (5. - GJP.DwfHeight ());

    save_prop =
      (WilsonMatrix *) smalloc (cname, fname, "save_prop",
				GJP.VolNodeSites () * sizeof (WilsonMatrix));

    for (int ii (0); ii < GJP.VolNodeSites (); ++ii)
      *(save_prop + ii) = renFac * prop[ii];


  } else
    save_prop = &prop[0];

  Float *save_source = (Float *) smalloc (cname, fname, "save_source",
					  GJP.VolNodeSites () * 288 *
					  sizeof (Float));
  FermionVectorTp src;

#ifdef USE_QIO
  Float qio_time = -dclock ();

  for (int spn = 0; spn < 4; spn++)
    for (int col = 0; col < 3; col++) {
//              Float *src_4d_f=(Float *)src_4d;
      SetSource (src, spn, col);
      // store the source
      for (int index (0); index < GJP.VolNodeSites (); ++index)
	for (int mm (0); mm < 4; ++mm)
	  for (int cc (0); cc < GJP.Colors (); ++cc) {
	    // now same ordering as propagator [volume][spin][color][solution_spin][solution_color][ReIm]
	    *(save_source + 288 * index + 72 * mm + 24 * cc +
	      6 * spn + 2 * col) = src[24 * index + 6 * mm + 2 * cc];
	    *(save_source + 288 * index + 72 * mm + 24 * cc +
	      6 * spn + 2 * col + 1) = src[24 * index + 6 * mm + 2 * cc + 1];
	  }
    }

  // always writes the full 4D source
  //qio_writePropagator writePropQio(propOutfile, QIO_FULL_SOURCE, save_prop, save_source,
  //                       qp_arg.ensemble_id, qp_arg.ensemble_label, qp_arg.seqNum, propType, sourceType,
  //                       GJP.argc(), GJP.argv(), VOLFMT);

  // write a t-slice/slices or hypercube in some cases for the source
  qio_writePropagator writePropQio;

  writePropQio.setHeader (qp_arg.ensemble_id, qp_arg.ensemble_label,
			  qp_arg.seqNum, propType, sourceType);

#if 1
  // just writing hypercube? SrcType doesn't work
  if (src_hyper) {
    VRB.Result (cname, fname,
	      " source: only write (%d %d %d %d) to (%d %d %d %d) to file\n",
	hyper_lower[0], hyper_lower[1], hyper_lower[2], hyper_lower[3],
	hyper_upper[0], hyper_upper[1], hyper_upper[2], hyper_upper[3]);
    writePropQio.setSourceHypercube (hyper_lower, hyper_upper);
  }
#endif

  switch (qp_arg.save_prop) {
  case 1:
    writePropQio.write_12pairs (propOutfile, QIO_FULL_SOURCE, save_prop,
				save_source, VOLFMT);
    break;
  case 2:
    for (int spn = 0; spn < 4; spn++) {
      for (int col = 0; col < 3; col++) {
#if 1
	char file[256];

	sprintf (file, "%ss%dc%d", propOutfile, spn, col);
	writePropQio.write_pair (file, QIO_FULL_SOURCE, save_prop,
				 save_source, spn, col, VOLFMT);
#else
	string file (propOutfile);

	file += '';
	writePropQio.write_pair (file.c_str (), QIO_FULL_SOURCE,
				 save_prop, save_source, spn, col, VOLFMT);
#endif
      }
    }
    break;
  default:
    ERR.General (cname, fname, "invalid save_prop in qp arge\n");
  }
  qio_time += dclock ();
  print_time ("QPropW::Run", "qio_writePropagator", qio_time);


#endif // USE_QIO

  if (AlgLattice ().Fclass () == F_CLASS_DWF)
    sfree (save_prop);

  sfree (save_source);
  // the old storage function
  //SaveQProp(qp_arg.file,PROP); 

  //Print out time taken to save
  Float dtime_this = dclock ();
  Float time_tmp = dtime_this - dtime_last;
  Float hr_tmp = time_tmp / 3600.0;
  Float min_tmp = (time_tmp - 3600.0 * hr_tmp) / 60.0;
  Float sec_tmp = time_tmp - 3600.0 * hr_tmp - 60.0 * min_tmp;

  VRB.Result (cname, fname,
	      "Time taken to save: %d hours %d minutes %f seconds.\n",
	      hr_tmp, min_tmp, sec_tmp);

}



CPS_END_NAMESPACE
