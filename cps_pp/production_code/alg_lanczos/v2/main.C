#include <conf.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <alg/alg_int.h>
#include <alg/int_arg.h>
#include <alg/alg_hmc.h>
#include <alg/hmc_arg.h>
#include <alg/alg_plaq.h>
#include <alg/alg_lanczos.h>
#include <alg/alg_meas.h>
#include <alg/alg_remez.h>
#include <alg/alg_eig.h>
#include <alg/alg_fix_gauge.h>
#include <alg/qpropw.h>
#include <alg/meson.h>
#include <alg/nuc2pt.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/error.h>
#include <util/command_line.h>
#include <util/ReadLatticePar.h>
#include <util/ReadU1LatticePar.h>
#include <util/qio_readLattice.h>
#include <util/qcdio.h>
#include <util/eigen_container.h>



USING_NAMESPACE_CPS 
using namespace std;


DoArg do_arg;
DoArgExt doext_arg;

//HmcArg hmc_arg;
//ActionGaugeArg gauge_arg;

// It is very imporatnt that mobius_arg is global and exits for ever
// as in GJP, we use zmobius_c_coeff and zmobius_b_coeff is pointing the contetnts of mobius_arg
MobiusArg mobius_arg;
//MobiusArg mobius_arg2;
//MdwfArg mdwf_arg;
//RemezArg remez_arg;

NoArg no_arg;
CommonArg common_arg;
LanczosArg lanczos_arg;
QPropWArg qp_arg;


void movefloattoFloat (Float * out, float *in, int f_size);


static void SetZmobiusPC(int flag)
{
    switch (flag){
	case 5:
    	GJP.ZMobius_PC_Type (ZMOB_PC_SYM3);
	break;
	case 0:
    	GJP.ZMobius_PC_Type (ZMOB_PC_ORIG);
	break;
	case 1:
    	GJP.ZMobius_PC_Type (ZMOB_PC_SYM1);
	break;
	case 2:
    	GJP.ZMobius_PC_Type (ZMOB_PC_SYM2);
	break;
	case 3:
    	GJP.ZMobius_PC_Type (ZMOB_PC_SYM1_MIT);
	break;
	case 4:
    	GJP.ZMobius_PC_Type (ZMOB_PC_SYM2_MIT);
	break;
    default:
	ERR.General("","SetZmobiusPC()","%d Not recognized\n",flag);
    }
}

EigenCache *ecache;
char *cname;

void CalMesons (QPropW & q, char *out_dir, int traj, char *src_str);
void CalNucleons (QPropW & q, char *out_dir, int traj, char *src_str);

char evecname_bc[1024];


void meas_bound (Lattice & lattice, Float *meas_high, Float *meas_low,Float mass)
{
  char *fname = "meas_bound(lattice&)";
  Float etime = time_elapse ();
  EigArg eig_arg;
  if (!eig_arg.Decode ("eig_arg.vml", "eig_arg")) {
    ERR.General (fname, fname, "Decoding of eig_arg failed\n");
  }
  eig_arg.Mass.Mass_val[0]=mass;
  eig_arg.mass=mass;
  eig_arg.Encode ("eig_arg.dat", "eig_arg");
if (meas_high)
{
  eig_arg.RitzMatOper =  NEG_MATPCDAG_MATPC;
  AlgEig ritz(lattice,&common_arg,&eig_arg);
  Float eigen[1][1];
// Segfault. Can't be bother to fix it at the moment
//  ritz.run((Float**)eigen,(Vector**)NULL);
//  VRB.Result("",fname,"high_value =%0.16g\n",eigen[0][0]);
//  *meas_high = eigen[0][0];
  ritz.run();
}
if (meas_low)
{
  eig_arg.RitzMatOper =  MATPCDAG_MATPC;
  AlgEig ritz(lattice,&common_arg,&eig_arg);
//  Float eigen[1][1];
//  ritz.run((Float**)eigen,(Vector**)NULL);
//  VRB.Result("",fname,"low =%0.16g\n",eigen[0][0]);
  ritz.run();
}
  etime = time_elapse ();
  if (!UniqueID ())
    printf ("Time for Ritz %g\n", etime);
}

void comp_read_eigenvectors (Lattice & lattice, int eig_start=0,int num_eig=0)
{
  char *fname = "comp_read_eigenvectors(lattice&, int)";
  VRB.Func("",fname);
  /*
   *  Compute eigen- vectors and values
   *
   */
  Float etime = time_elapse ();
  char cache_name[1024];
  snprintf (cache_name, 1024, "cache_0_mass%g", lanczos_arg.mass);
  ecache = new EigenCache (cache_name);
  char evecname[1024];
  snprintf (evecname, 1024, "%s/eig4dee.mass%g.ls%d.pc%d",
	    lanczos_arg.file, lanczos_arg.mass, GJP.SnodeSites (), GJP.ZMobius_PC_Type ());
  lanczos_arg.file = evecname;



  snprintf (evecname_bc, 1024, "%s.bc%d%d%d%d", lanczos_arg.file, GJP.Bc (0),
	    GJP.Bc (1), GJP.Bc (2), GJP.Bc (3));


  AlgLanczos eig (lattice, &common_arg, &lanczos_arg, ecache);
  int Ncb = eig.NumChkb (lanczos_arg.RitzMat_lanczos);
  int fsize = GJP.VolNodeSites () * lattice.FsiteSize () * Ncb / 2 / 2;	//last 2 for single prec.;
  EigenCacheList.push_back (ecache);
  int neig;


#if 1
  if (lanczos_arg.nt_lanczos_vectors > 0) {

    int init_flag = 0;		// 0 : wall start, 1 :  compression, 2 : compress & decompress, 3: refinement using Ritz
    int ncompress = 1;		// number of evecs in the compressed linear combination

    char *comp_file;

    neig = lanczos_arg.nk_lanczos_vectors + lanczos_arg.np_lanczos_vectors;
    ecache->alloc (evecname_bc, neig, fsize);
    eig.run (init_flag, ncompress, comp_file);
    neig = lanczos_arg.nt_lanczos_vectors;
    ecache->set_neig (neig);
  } else {
    if (num_eig>0) neig=num_eig;
    else
    neig = mobius_arg.cg.neig;
    if (neig > 0) {
      ecache->alloc (evecname_bc, neig, fsize);
      {				//read in only
	const int n_fields = GJP.SnodeSites ();
	const int f_size_per_site = lattice.FsiteSize () / n_fields / 2;
	EigenContainer eigcon (lattice, evecname_bc, neig, f_size_per_site / 2,
			       n_fields, ecache);
	// factor of 2 for single-prec.
	// have to do this if stride != 1 FIX!
	for (int iev = 0; iev < neig; iev++) {
	  Vector *evec = eigcon.nev_load (iev+eig_start);
	  ecache->set_index (iev);
	}
      }
    }
  }

  etime = time_elapse ();
  if (!UniqueID ())
    printf ("Time for Lanczos %g\n", etime);
#endif



  // Read eigenvectors for small Ls mobius in mobius_arg
  //-------------------------------------------------------
  VRB.FuncEnd("",fname);

}

void truncate_it(CommonArg *common_arg, const char stem[], int traj)
{
    char fnbuf[1024];
    sprintf(fnbuf, "%s.%d", stem, traj);
    FILE *truncate_it = Fopen(fnbuf, "w");
    Fclose(truncate_it);
    common_arg->set_filename(fnbuf);
}

int main (int argc, char *argv[])
{

  cname = argv[0];
  char *fname = "main()";

  Start (&argc, &argv);

  if (argc < 3) {
    if (!UniqueID ())
      printf
	("(exe) working_dir ZmobPc meas_arg (io_sparse) but argc=%d\n",
	 argc);
    exit (-1);
  }

  chdir (argv[1]);
//  int do_zmob_sm = atoi (argv[5]);
  int do_meas = atoi(argv[3]);
  int ntraj = atoi(argv[4]);
  int io_sparse=16;
  if (argc>5) io_sparse = atoi(argv[5]);
  setQioSparseNum(io_sparse);

  if (!do_arg.Decode ("do_arg.vml", "do_arg")) {
    ERR.General (fname, fname, "Decoding of do_arg failed\n");
  }
  // make a record of what was run
  do_arg.Encode ("do_arg.dat", "do_arg");

  if (!doext_arg.Decode ("doext_arg.vml", "doext_arg")) {
    ERR.General (fname, fname, "Decoding of doext_arg failed\n");
  }
  // make a record of what was run
  doext_arg.Encode ("doext_arg.dat", "doext_arg");

  if (!mobius_arg.Decode ("mobius_arg.vml", "mobius_arg")) {
    mobius_arg.Encode ("mobius_arg.dat", "mobius_arg");
    ERR.General (fname, fname, "Decoding of mobius_arg failed\n");
  }
  mobius_arg.Encode ("mboius_arg.dat", "mobius_arg");

  if (!qp_arg.Decode ("qpropw_arg.vml", "qpropw_arg")) {
    ERR.General (fname, fname, "Decoding of qpropw_arg failed\n");
  }
  qp_arg.Encode ("qpropw_arg.dat", "qpropw_arg");

  if (!lanczos_arg.Decode ("lanczos_arg.vml", "lanczos_arg")) {
    lanczos_arg.Encode ("lanczos_arg.dat", "lanczos_arg");
    if (!UniqueID ())
      printf ("Decoding of lanczos_arg failed\n");
    exit (-1);
  }
  lanczos_arg.Encode ("lanczos_arg.dat", "lanczos_arg");


  GJP.Initialize (do_arg);
  GJP.InitializeExt (doext_arg);
  VRB.Result("","main()","GJP.SaveStride()=%d\n",GJP.SaveStride());
  VRB.Level (do_arg.verbose_level);

  // Solve  Small Ls with Zmobius
  //--------------------------------------------
 // if (do_zmob_sm) 
  {
    GJP.SnodeSites (mobius_arg.ls);
#if 0
    GnoneFzmobius lattice;
    GJP.ZMobius_b (mobius_arg.zmobius_b_coeff.zmobius_b_coeff_val,
		   mobius_arg.ls);
    GJP.ZMobius_c (mobius_arg.zmobius_c_coeff.zmobius_c_coeff_val,
		   mobius_arg.ls);
#endif
#if 0
    GnoneFmobius lattice;
    GJP.Mobius_b (mobius_arg.mobius_b_coeff);
    GJP.Mobius_c (mobius_arg.mobius_c_coeff);
#endif
    GnoneFstag lattice;
//    GnoneFdwf lattice;


    SetZmobiusPC(atoi(argv[2]));

    VRB.Result(cname,"main()", "GJP.ZMobius_PC_Type() = %d\n",GJP.ZMobius_PC_Type());

if (do_meas){
    Float high = 0.,low=0.;
    meas_bound (lattice,&high,&low,mobius_arg.cg.mass);
//  cps::sync();exit(-1);
}
    comp_read_eigenvectors (lattice);

    CgArg cg_save = qp_arg.cg;
    qp_arg.cg = mobius_arg.cg;

    CommonArg carg;

    QPropWPointSrc qp (lattice, &qp_arg, &carg);
    // 
    system ("mkdir -p zmob_sm");
    CalMesons (qp, "zmob_sm", 0, "point");
    CalNucleons (qp, "zmob_sm", 0, "point");

    qp_arg.cg = cg_save;
  }




  EigenCacheListCleanup ();

  //End();
  return 0;
}


void movefloattoFloat (Float * out, float *in, int f_size)
{

  float flt;
  for (int i = 0; i < f_size; i++) {
    flt = in[i];
    out[i] = (Float) flt;
  }
};


void CalNucleons (QPropW & q, char *out_dir, int traj, char *src_str)
{
  Nuc2pt nuc (NUC_G5C, POINT);
  nuc.calcNucleon (q);

  char file[256];
  sprintf (file, "%s/nucleon_%s.dat.%d", out_dir, src_str, traj);
  FILE *fp = Fopen (&file[0], "w");

  {
    Nuc2pt nuc (NUC_G5C, POINT);
    nuc.calcNucleon (q);
    nuc.Print (fp);
  }

  Fclose (fp);


}

void CalMesons (QPropW & q, char *out_dir, int traj, char *src_str)
{

  char file[256];
  sprintf (file, "%s/meson_%s.dat.%d", out_dir, src_str, traj);
  FILE *fp = Fopen (&file[0], "w");

  int mu, nu;
  //vector and scalar
  for (mu = 0; mu < 4; mu++) {
    char str[256];
    sprintf (str, "GAM_%d", mu + 1);
    Meson mes (mu, str);
    mes.Zero ();
    mes.setMass (q.Mass (), q.Mass ());
    mes.calcMeson (q, q);
    if (!UniqueID ())
      mes.Print (fp);
  }

  //pseudoscalar
  {
    mu = -5;
    char str[256];
    sprintf (str, "GAM_5");
    Meson mes (mu, str);
    mes.Zero ();
    mes.setMass (q.Mass (), q.Mass ());
    mes.calcMeson (q, q);
    if (!UniqueID ())
      mes.Print (fp);
  }
  //tensor
  for (mu = 0; mu < 4; mu++) {
    for (nu = mu + 1; nu < 4; nu++) {
      char str[256];
      sprintf (str, "GAM_%d%d", mu + 1, nu + 1);
      Meson mes1 (mu, nu, str);
      mes1.Zero ();
      mes1.setMass (q.Mass (), q.Mass ());
      mes1.calcMeson (q, q);
      if (!UniqueID ())
	mes1.Print (fp);

      sprintf (str, "GAM_%d%d", nu + 1, mu + 1);
      Meson mes2 (nu, mu, str);
      mes2.Zero ();
      mes2.setMass (q.Mass (), q.Mass ());
      mes2.calcMeson (q, q);
      if (!UniqueID ())
	mes2.Print (fp);
    }
  }
  //axialvector
  nu = -5;
  for (mu = 0; mu < 4; mu++) {
    char str[256];
    sprintf (str, "GAM_%d5", mu + 1);
    Meson mes (mu, nu, str);
    mes.Zero ();
    mes.setMass (q.Mass (), q.Mass ());
    mes.calcMeson (q, q);
    if (!UniqueID ())
      mes.Print (fp);
  }

  // FIXME:  PUT IN THE  A0  SCALAR MESON, or WE WILL MISS A LOT !!

  Fclose (fp);

}

