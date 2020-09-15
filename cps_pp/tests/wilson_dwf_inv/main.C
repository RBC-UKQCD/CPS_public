#include<omp.h>
#include<config.h>
#include <util/qcdio.h>
#include <math.h>
#include<util/lattice.h>
#include<util/gjp.h>
#include<util/verbose.h>
#include<util/dirac_op.h>
#include<util/wilson.h>
#include<util/error.h>
#include<util/time_cps.h>
#include<comms/scu.h>
#include<alg/alg_hmd.h>
#include<alg/do_arg.h>
#include <util/omp_wrapper.h>
#include <util/command_line.h>

//#include<cps.h>

//#define OMP(A) #pragma omp A


USING_NAMESPACE_CPS DoArg do_arg;
DoArgExt doext_arg;
//MobiusArg mobius_arg;


//const char *f_wilson_test_filename = CWDPREFIX("f_wilson_test");
//const char *psi_filename = CWDPREFIX("psi");

static int nx, ny, nz, nt, ns;
static CgArg cg_arg;


static void SetZmobiusPC(int flag)
{
    switch (flag){
        case 0:
        GJP.ZMobius_PC_Type (ZMOB_PC_ORIG);
        break;
        case 1:
        GJP.ZMobius_PC_Type (ZMOB_PC_SYM1);
        break;
        case 2:
        GJP.ZMobius_PC_Type (ZMOB_PC_SYM2);
        break;
#if 0
        case 5:
        GJP.ZMobius_PC_Type (ZMOB_PC_SYM3);
        break;
        case 3:
        GJP.ZMobius_PC_Type (ZMOB_PC_SYM1_MIT);
        break;
        case 4:
        GJP.ZMobius_PC_Type (ZMOB_PC_SYM2_MIT);
        break;
#endif
    default:
        ERR.General("","SetZmobiusPC()","%d Not implemented\n",flag);
    }
}

//void run_inv(Lattice &lat, DiracOp &dirac, StrOrdType str_ord, char *out_name, int DO_CHECK);
void run_inv (Lattice & lat, StrOrdType str_ord, char *out_name, int DO_CHECK);

int main (int argc, char *argv[])
{

  Start (&argc, &argv);
  char *fname = "main()";
#pragma omp parallel default(shared)
  {
    int tnum = omp_get_num_threads ();


#pragma omp for
    for (int i = 0; i < 100; i++) {
      if (!UniqueID ())
	printf ("thread %d of %d i=%d\n", omp_get_thread_num (), tnum, i);
    }
  }


  //----------------------------------------------------------------
  // Initializes all Global Job Parameters
  //----------------------------------------------------------------
  DoArg do_arg;
  char *out_file = NULL;
  CommandLine::is(argc,argv);
//  if (argc > 1)
    out_file = CommandLine::arg();
  if (!do_arg.Decode (CommandLine::arg(), "do_arg")) {
    ERR.General ("", fname, "Decoding of do_arg failed\n");
  }
  do_arg.Encode ("do_arg.dat", "do_arg");
  if (!doext_arg.Decode (CommandLine::arg(), "doext_arg")) {
    ERR.General ("", fname, "Decoding of doext_arg failed\n");
  }
  doext_arg.Encode ("doext_arg.dat", "doext_arg");
  if (!cg_arg.Decode (CommandLine::arg(), "cg_arg")) {
    ERR.General ("", fname, "Decoding of cg_arg failed\n");
  }
#ifdef USE_QUDA
  if (!QudaParam.Decode (CommandLine::arg(), "QudaParam")) {
    printf ("Bum quda_arg\n");
    exit (-1);
  }
#endif
  SetZmobiusPC(CommandLine::arg_as_int());

  GJP.Initialize (do_arg);
  GJP.InitializeExt (doext_arg);


  {
    GwilsonFmobius lat;
//    DiracOpMobius dirac(lat,NULL,NULL,&cg_arg,CNV_FRM_NO);
    run_inv (lat, DWF_4D_EOPREC_EE, out_file, 1);
  }

  End ();
}

//void run_inv(Lattice &lat, DiracOp &dirac, StrOrdType str_ord, char *out_name, int DO_CHECK){
void run_inv (Lattice & lat, StrOrdType str_ord, char *out_name, int DO_CHECK)
{
  FILE *fp;
  const char *fname = "run_inv()";
  double dtime;
  int DO_IO = 1;
  if (out_name == NULL)
    DO_IO = 0;
  if (DO_IO)
    fp = Fopen (ADD_ID, out_name, "w");
  else
    fp = stdout;
  VRB.Result ("", "", "DO_CHECK=%d DO_IO=%d\n", DO_CHECK, DO_IO);

  Vector *result =
    (Vector *) smalloc ("", "", "result",
			GJP.VolNodeSites () * lat.FsiteSize () * sizeof (IFloat));
  Vector *X_out =
    (Vector *) smalloc ("", "", "X_out",
			GJP.VolNodeSites () * lat.FsiteSize () * sizeof (IFloat));
  Vector *X_out2 = (Vector *) smalloc ("", "", "X_out2",
			GJP.VolNodeSites () * lat.FsiteSize () * sizeof (IFloat));
  Vector *tmp = (Vector *) smalloc ("", "", "tmp",
			GJP.VolNodeSites () * lat.FsiteSize () * sizeof (IFloat));

  int s_size = 1;
  if (lat.F5D ())
    s_size = GJP.SnodeSites ();

  int s[5];
  Vector *X_in =
    (Vector *) smalloc (GJP.VolNodeSites () * lat.FsiteSize () *
			sizeof (IFloat));
  if (!X_in)
    ERR.Pointer ("", "", "X_in");
  memset (X_in, 0, GJP.VolNodeSites () * lat.FsiteSize () * sizeof (IFloat));
#if 1
  lat.RandGaussVector (X_in, 1.0);
#else

//    lat.RandGaussVector(X_in,1.0);
  Matrix *gf = lat.GaugeField ();
  IFloat *gf_p = (IFloat *) lat.GaugeField ();
  int fsize = lat.FsiteSize () / s_size;
  VRB.Result ("", "main()", "fsize=%d", fsize);

  for (s[4] = 0; s[4] < s_size; s[4]++)
    for (s[3] = 0; s[3] < GJP.NodeSites (3); s[3]++)
      for (s[2] = 0; s[2] < GJP.NodeSites (2); s[2]++)
	for (s[1] = 0; s[1] < GJP.NodeSites (1); s[1]++)
	  for (s[0] = 0; s[0] < GJP.NodeSites (0); s[0]++) {

	    int n = lat.FsiteOffset (s) + s[4] * GJP.VolNodeSites ();

	    int crd = 1.;
	    if (CoorX () == 0 && CoorY () == 0 && CoorZ () == 0 && CoorT () == 0)
	      crd = 1.0;
	    else
	      crd = 0.0;
	    if (s[0] != 0) crd = 0.;
	    if (s[1] != 0) crd = 0.;
	    if (s[2] != 0) crd = 0.;
	    if (s[3] != 0) crd = 0.;
	    if (s[4] != 0) crd = 0.;

	    IFloat *X_f = (IFloat *) (X_in) + (n * fsize);
	    for (int v = 0; v < fsize; v += 1) {
	      if ((v%6)  == 0)
		*(X_f + v) = crd;
	      else
		*(X_f + v) = 0;
	    }
	  }
#endif

  Vector *out;
  Float true_res;

  for (int k = 0; k < 1; k++) {
    double maxdiff = 0.,maxdiff2=0.;
    VRB.Result ("", fname, "k=%d\n", k);
    if (k == 0)
      out = result;
    else
      out = X_out;
    memset ((char *) out, 0,
	    GJP.VolNodeSites () * lat.FsiteSize () * sizeof (IFloat));
    lat.Fconvert (out, str_ord, CANONICAL);
    lat.Fconvert (X_in, str_ord, CANONICAL);
    int offset = GJP.VolNodeSites () * lat.FsiteSize () / (2 * 6);
#if 1
    dtime = -dclock ();
    int iter =
      lat.FmatInv (out, X_in, &cg_arg, &true_res, CNV_FRM_NO, PRESERVE_YES,
		   0);
// testing if the guess works
      lat.FmatInv (out, X_in, &cg_arg, &true_res, CNV_FRM_NO, PRESERVE_YES,
		   0);
    dtime += dclock ();
#else
    lat.Fdslash (result, X_in, &cg_arg, CNV_FRM_NO, 0);
#endif

    if (1) {
//      if (k == 0) {
	memset ((char *) X_out2, 0,
		GJP.VolNodeSites () * lat.FsiteSize () * sizeof (IFloat));
	lat.Fdslash (X_out2, out, &cg_arg, CNV_FRM_NO, 0);
//      }
    }
    lat.Fconvert (out, CANONICAL, str_ord);
    lat.Fconvert (X_in, CANONICAL, str_ord);
    lat.Fconvert (X_out2, CANONICAL, str_ord);
//      if (lat.F5D())
    if (lat.Fclass () == F_CLASS_DWF)
      X_out2->FTimesV1PlusV2 (-0.5 / (5.0 - GJP.DwfHeight ()), X_out2, out,
			      GJP.VolNodeSites () * lat.FsiteSize ());
    else if (lat.Fclass () == F_CLASS_WILSON)
      X_out2->FTimesV1PlusV2 (-0.5 / (cg_arg.mass + 4.0), X_out2, out,
			      GJP.VolNodeSites () * lat.FsiteSize ());

    Float dummy;
    Float dt = 2;
    if (DO_CHECK) {
      for (s[4] = 0; s[4] < s_size; s[4]++)
	for (s[3] = 0; s[3] < GJP.NodeSites (3); s[3]++)
	  for (s[2] = 0; s[2] < GJP.NodeSites (2); s[2]++)
	    for (s[1] = 0; s[1] < GJP.NodeSites (1); s[1]++)
	      for (s[0] = 0; s[0] < GJP.NodeSites (0); s[0]++) {

//                  int n = lat.FsiteOffset(s)*lat.SpinComponents()*GJP.SnodeSites();
		int n =
		  (lat.FsiteOffset (s) +
		   GJP.VolNodeSites () * s[4]) * lat.SpinComponents ();
		for (int i = 0; i < (3 * lat.SpinComponents ()); i++) {
		  double re_re = *((IFloat *) & out[n] + i * 2);
		  double in_re = *((IFloat *) & X_in[n] + i * 2);
		  double re_im = *((IFloat *) & out[n] + i * 2 + 1);
		  double in_im = *((IFloat *) & X_in[n] + i * 2 + 1);
		  if ((re_re * re_re + re_im * re_im + in_re * in_re +
		       in_im * in_im) > 1e-8)
		    if ( (!k) && DO_IO) {
			Fprintf (ADD_ID, fp, " %d %d %d %d %d (%d) ",
				 CoorX () * GJP.NodeSites (0) + s[0],
				 CoorY () * GJP.NodeSites (1) + s[1],
				 CoorZ () * GJP.NodeSites (2) + s[2],
				 CoorT () * GJP.NodeSites (3) + s[3],
				 CoorS () * GJP.NodeSites (4) + s[4], i, n);
			Fprintf (ADD_ID, fp, " ( %0.7e %0.7e ) (%0.7e %0.7e)",
				 *((IFloat *) & out[n] + i * 2),
				 *((IFloat *) & out[n] + i * 2 + 1),
				 *((IFloat *) & X_in[n] + i * 2),
				 *((IFloat *) & X_in[n] + i * 2 + 1));
#if 1
		      Fprintf (ADD_ID, fp, " ( %0.2e %0.2e )\n",
#if 1
			       *((IFloat *) & X_out2[n] + i * 2) - *((IFloat *) & X_in[n] + i * 2),
			       *((IFloat *) & X_out2[n] + i * 2 + 1) - *((IFloat *) & X_in[n] + i * 2 + 1));
#else
			       *((IFloat *) & X_out2[n] + i * 2),
			       *((IFloat *) & X_out2[n] + i * 2 + 1));
#endif
#else
		      Fprintf (ADD_ID, fp, "\n");
#endif
		    }		//DO_IO
		  double diff =
		    *((IFloat *) & X_out2[n] + i * 2) - *((IFloat *) & X_in[n] +
							  i * 2);
		  if (fabs (diff) > maxdiff)
		    maxdiff = fabs (diff);
		  diff =
		    *((IFloat *) & X_out2[n] + i * 2 + 1) -
		    *((IFloat *) & X_in[n] + i * 2 + 1);
		  if (fabs (diff) > maxdiff)
		    maxdiff = fabs (diff);

		  diff =
		    *((IFloat *) & out[n] + i * 2) - *((IFloat *) & result[n] +
							  i * 2);
		  if (fabs (diff) > maxdiff2)
		    maxdiff2 = fabs (diff);
		  diff = *((IFloat *) & out[n] + i * 2 + 1) - *((IFloat *) & result[n] + i * 2 + 1);
		  if (fabs (diff) > maxdiff2)
		    maxdiff2 = fabs (diff);
		}
	      }
      VRB.Result ("", "run_inv()",
		  "Max diff between X_in and M*X_out = %0.2e result and out = %0.2e\n", maxdiff,maxdiff2);
    }
  }
  if (DO_IO)
    Fclose (fp);

  sfree (X_in);
  sfree (result);
  sfree (X_out);
  sfree (X_out2);
}
