#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2003-07-24 16:53:53 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_s_spect/alg_s_spect.C,v 1.2 2003-07-24 16:53:53 zs Exp $
//  $Id: alg_s_spect.C,v 1.2 2003-07-24 16:53:53 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.9  2002/03/11 22:25:42  anj
//  This should now be the correct, fully merged code from our two versions. Anj
//
//  Revision 1.6.2.1  2002/03/08 16:34:59  anj
//  Checking in the Columbia code branch on tag Columbia4_1_1_test-branch, to be
//  merged with the UKQCD head branch shortly.  Anj
//
//  Revision 1.6  2001/08/16 10:49:40  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.4  2001/07/03 17:00:45  anj
//
//  Multiple minor alterations to change some #include's from referring to
//  files relative to the top-level source directory to referring to files
//  relative to the source-file positions.  This alteration makes the code
//  backwards compatable with the make structure of QCDSP, although this
//  may have to be changed to a more usual form in the future. Anj.
//
//  Revision 1.3  2001/06/21 15:40:09  anj
//  Updated the _TARTAN ifdefs, using PARALLEL instead (where appropriate).Anj
//
//  Revision 1.2  2001/06/19 18:11:29  anj
//  Serious ANSIfication.  Plus, degenerate double64.h files removed.
//  Next version will contain the new nga/include/double64.h.  Also,
//  Makefile.gnutests has been modified to work properly, propagating the
//  choice of C++ compiler and flags all the way down the directory tree.
//  The mpi_scu code has been added under phys/nga, and partially
//  plumbed in.
//
//  Everything has newer dates, due to the way in which this first alteration was handled.
//
//  Anj.
//
//  Revision 1.2  2001/05/25 06:16:00  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: alg_s_spect.C,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_s_spect/alg_s_spect.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// alg_s_spect.C
//
// Code for all alg classes relevant to staggered fermion
// spectroscopy. The type of glue or fermion is given as
// an argument of type Lattice& to the constructor. If the 
// fermion type is not F_CLASS_STAG the constructors exit with
// a general error.
//
//------------------------------------------------------------------


CPS_END_NAMESPACE
#include <stdlib.h>	// exit()
#include <stdio.h>
#include <alg/quark_prop_s.h>
#include <alg/meson_prop_s.h>
#include <alg/mom_meson_p_s.h>
#include <alg/nucl_prop_s.h>
#include <alg/nlocal_prop_s.h>
#include <alg/aots_s.h>
#include <alg/alg_s_spect.h>
#include <alg/common_arg.h>
#include <alg/s_spect_arg.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/smalloc.h>
#include <util/vector.h>
#include <util/verbose.h>
#include <util/error.h>
//#include <util/mom.h>
#include <alg/myenum.h>
CPS_START_NAMESPACE
#ifdef PARALLEL
CPS_END_NAMESPACE
#include <sysfunc.h>
CPS_START_NAMESPACE
#endif

QuarkPropSMng AlgStagQuark::sqpm;

static char *class_str = "Class type was %d instead of F_CLASS_STAG";

static void zero_buffer(Float *buf, int len)
{ 
  for(int i = 0; i < len; i++)
    *buf++ = 0;
}

static void write_to_file(const Float *data_p, FILE *fp, 
			  int num_IFloats, int num_slices, int mom, 
			  HadronType type, BndCndType bc)
{
  // printf("write_to_file: no_of_momenta = %d \n",mom);
  int unit;
  switch(type) {
    case SMESON:
	unit = 4; break;
    case SMOMMESON:
        unit = 4*mom; break;
    case SNUCLEON: 
	unit = 1; break;
    case SNONLOCAL:
	unit = 4; break;
    default: 
	unit = 0; break;
  }

  /*
  printf("write_to_file: mom=%d unit =%d num_IFloats=%d num_slices=%d \n",mom,uni
t,num_IFloats,num_slices);

  printf("data_p[0] = %e %e %e \n",data_p[0],data_p[1],data_p[2]);

  FILE *fp1=fopen("total_correlator.dat", "a");
  for (int ii=0; ii < num_IFloats; ii++) fprintf(fp1,"%e\n",(IFloat)data_p[ii]);
  fclose(fp1);
  */
  
  int i, k;

  //-------------------------------------------------------------
  // * write C(t = 0)
  //-------------------------------------------------------------
  for (i = 0; i < unit; ++i) {
    Float tmp = data_p[2*i] / num_slices;
    fprintf(fp, "%e\n", (IFloat)tmp);
  }
  
  //-------------------------------------------------------------
  // * write (C(t) + C(N_t-t))/2 for t = 1 to (N_t/2)-1
  //-------------------------------------------------------------
  int half_length = num_IFloats / 4;
  int N_t = num_IFloats / (2*unit);
  
  if (type == SMESON || type == SMOMMESON) {
    for (i = unit; i < half_length; ++i) {
      k = unit*(N_t - i/unit) + i%unit;
      Float tmp = (data_p[2*i] + data_p[2*k])*0.5/num_slices;
      fprintf(fp, "%e\n", (IFloat)tmp);
    }
  }
  else if (bc == BND_CND_APRD)  {
    for (i = unit; i < half_length; ++i) {
      k = unit*(N_t - i/unit) + i%unit;
      Float tmp;

      if (((i/unit) % 2) == 0)
  	tmp = (data_p[2*i] - data_p[2*k])*0.5/num_slices;
      else
  	tmp = (data_p[2*i] + data_p[2*k])*0.5/num_slices;

      fprintf(fp, "%e\n", (IFloat)tmp);
    }
  }

  else {
    for (i = unit; i < half_length; ++i) {
      k = unit*(N_t - i/unit) + i%unit;
      Float tmp;

      if (((i/unit) % 2) == 0)
  	tmp = (data_p[2*i] + data_p[2*k])*0.5/num_slices;
      else
  	tmp = (data_p[2*i] - data_p[2*k])*0.5/num_slices;

      fprintf(fp, "%e\n", (IFloat)tmp);
    }
  }

  //-------------------------------------------------------------
  // * write C(t=N_t/2)
  //-------------------------------------------------------------
  for (i = half_length; i < half_length + unit; ++i) {
    Float tmp = data_p[2*i] / num_slices;
    fprintf(fp, "%e\n", (IFloat)tmp);
  }
}
  
//------------------------------------------------------------------
//------------------------------------------------------------------
//
// AlgStagQuark is derived from Alg and is relevant to  
// the staggered quark propagator.
// The type of fermion is determined by the argument to the 
// constructor. If the fermion type is not F_CLASS_STAG the 
// constructor will exit with a general error.
//
//------------------------------------------------------------------
//------------------------------------------------------------------

//------------------------------------------------------------------
// Constructor 
//------------------------------------------------------------------
AlgStagQuark::AlgStagQuark(Lattice& latt, 
			   CommonArg *c_arg,
			   StagQuarkArg *arg, Aots &aots) 
: Alg(latt, c_arg)
{
  cname = "AlgStagQuark";
  char *fname = "AlgStagQuark(L&,CommonArg*,StagQuarkArg*)";
  VRB.Func(cname,fname);

  // Check fermion type
  //----------------------------------------------------------------
  if(latt.Fclass() != F_CLASS_STAG)
    ERR.General(cname,fname, class_str, int(latt.Fclass()));


  // Initialize the argument pointer
  //----------------------------------------------------------------
  if(arg == 0)
    ERR.Pointer(cname,fname, "arg");

  arg->src.origin[arg->src.dir] = arg->src.end[arg->src.dir] = aots.slice();
  alg_stag_quark_arg = arg;
}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
AlgStagQuark::~AlgStagQuark() {
  char *fname = "~AlgStagQuark()";
  VRB.Func(cname,fname);
}


//------------------------------------------------------------------
// Calculate propagator
//------------------------------------------------------------------
void AlgStagQuark::run()
{
  char *fname = "run()";
  VRB.Func(cname,fname);

  // Set the Lattice pointer
  //----------------------------------------------------------------
  Lattice& lat = AlgLattice();

  QuarkPropS q1(lat, *alg_stag_quark_arg);
  q1.setupQuarkPropS();
  q1.getQuarkPropS((char *)common_arg->results);
}


//------------------------------------------------------------------
// Destroy propagator
//------------------------------------------------------------------
void AlgStagQuark::free()
{
  char *fname = "destroy()";
  VRB.Func(cname,fname);

  // Set the Lattice pointer
  //----------------------------------------------------------------
  Lattice& lat = AlgLattice();

  sqpm.destroyQuarkPropS(alg_stag_quark_arg->qid);
}



//------------------------------------------------------------------
//------------------------------------------------------------------
//
// AlgStagMeson is derived from Alg and is relevant to  
// the staggered meson propagator.
// The type of fermion is determined by the argument to the 
// constructor. If the fermion type is not F_CLASS_STAG the 
// constructor will exit with a general error.
//
//------------------------------------------------------------------
//------------------------------------------------------------------

//------------------------------------------------------------------
// Constructor 
//------------------------------------------------------------------
AlgStagMeson::AlgStagMeson(Lattice& latt, 
			   CommonArg *c_arg,
			   StagMesonArg *arg, Aots& a) 
:   Alg(latt, c_arg), aots(a) 
{
  cname = "AlgStagMeson";
  char *fname = "AlgStagMeson(L&,CommonArg*,StagMesonArg*)";
  VRB.Func(cname,fname);

  // Check fermion type
  //----------------------------------------------------------------
  if(latt.Fclass() != F_CLASS_STAG)
    ERR.General(cname,fname, class_str, int(latt.Fclass()));

  
  // Initialize the argument pointer
  //----------------------------------------------------------------
  if(arg == 0)
    ERR.Pointer(cname,fname, "arg");
  alg_stag_meson_arg = arg;
}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
AlgStagMeson::~AlgStagMeson() {
  char *fname = "~AlgStagMeson()";
  VRB.Func(cname,fname);
}


//------------------------------------------------------------------
// Calculate propagator
//------------------------------------------------------------------
void AlgStagMeson::run()
{
  char *fname = "run()";
  VRB.Func(cname,fname);

  // Set the Lattice pointer
  //----------------------------------------------------------------
  Lattice& lat = AlgLattice();

  MesonPropS m1(lat, *alg_stag_meson_arg);
  m1.getHadronPropS();

  // download propagator and AOTS 
  //----------------------------------------------------------------
  if(aots.begin()) {
    alg_stag_meson_arg->meson_buf = 
		(Float *)smalloc(m1.propLenTotal()*sizeof(Float));
    if(alg_stag_meson_arg->meson_buf == 0)
          ERR.Pointer(cname,fname, "alg_stag_meson_arg->meson_buf");
    VRB.Smalloc(cname,fname, "alg_stag_meson_arg->meson_buf", 
		alg_stag_meson_arg->meson_buf, m1.propLenTotal()*sizeof(Float));
    zero_buffer(alg_stag_meson_arg->meson_buf, m1.propLenTotal());
  }

  m1.download_prop(SMESON, alg_stag_meson_arg->meson_buf);

  if (aots.last()) {
    char *data_file = common_arg->results ? 
		      (char *)common_arg->results : CAST_AWAY_CONST("meson.def"); 
    FILE *fp;
    if( NULL == (fp = fopen(data_file, "a")) ) {
      ERR.FileA(cname,fname, data_file);
    }

    write_to_file(alg_stag_meson_arg->meson_buf, fp, 
   	          m1.propLenTotal(), aots.numSlices(), 1, SMESON, m1.bcd());
    fclose(fp);

    VRB.Sfree(cname,fname, "alg_stag_meson_arg->meson_buf", 
			    alg_stag_meson_arg->meson_buf);
    sfree(alg_stag_meson_arg->meson_buf);
    alg_stag_meson_arg->meson_buf = 0;
  }
}



//------------------------------------------------------------------
//------------------------------------------------------------------
//
// AlgStagMomMeson is derived from Alg and is relevant to
// the staggered meson propagator.
// The type of fermion is determined by the argument to the
// constructor. If the fermion type is not F_CLASS_STAG the
// constructor will exit with a general error.
//
//------------------------------------------------------------------
//------------------------------------------------------------------

//------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------
AlgStagMomMeson::AlgStagMomMeson(Lattice& latt,
                           CommonArg *c_arg,
                           StagMomMesonArg *arg, Aots& a)
:   Alg(latt, c_arg), aots(a)
{
  cname = "AlgStagMomMeson";
  char *fname = "AlgStagMomMeson(L&,CommonArg*,StagMomMesonArg*, Aots& )";
  VRB.Func(cname,fname);

  // Check fermion type
  //----------------------------------------------------------------
  if(latt.Fclass() != F_CLASS_STAG)
    ERR.General(cname,fname, class_str, int(latt.Fclass()));


  // Initialize the argument pointer
  //----------------------------------------------------------------
  if(arg == 0)
    ERR.Pointer(cname,fname, "arg");
  alg_stag_mom_meson_arg = arg;
}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
AlgStagMomMeson::~AlgStagMomMeson() {
  char *fname = "~AlgStagMomMeson()";
  VRB.Func(cname,fname);
}


//------------------------------------------------------------------
// Calculate propagator
//------------------------------------------------------------------
void AlgStagMomMeson::run()
{
  char *fname = "run()";
  VRB.Func(cname,fname);

  // Set the Lattice pointer
  //----------------------------------------------------------------
  Lattice& lat = AlgLattice();

  MomMesonPropS m1(lat, *alg_stag_mom_meson_arg);
  m1.getHadronPropS();

  // download propagator and AOTS
  //----------------------------------------------------------------
  if(aots.begin()) {
    alg_stag_mom_meson_arg->meson_buf =
      (Float *)smalloc(m1.propLenTotal()*sizeof(Float));
    if(alg_stag_mom_meson_arg->meson_buf == 0)
      ERR.Pointer(cname,fname, "alg_stag_mom_meson_arg->meson_buf");
    VRB.Smalloc(cname,fname, "alg_stag_mom_meson_arg->meson_buf",
         alg_stag_mom_meson_arg->meson_buf, m1.propLenTotal()*sizeof(Float));
    zero_buffer(alg_stag_mom_meson_arg->meson_buf, m1.propLenTotal());
  }


  m1.download_prop(SMOMMESON, alg_stag_mom_meson_arg->meson_buf);

  if (aots.last()) {
    char *data_file = common_arg->results ?
                      (char *)common_arg->results : (char *) "mom_meson.def";
    FILE *fp;
    if( NULL == (fp = fopen(data_file, "a")) ) {
      ERR.FileA(cname,fname, data_file);
    }

    write_to_file(alg_stag_mom_meson_arg->meson_buf, fp,
                  m1.propLenTotal(), aots.numSlices(), 
	alg_stag_mom_meson_arg->no_of_momenta,SMOMMESON, m1.bcd());
    fclose(fp);

    VRB.Sfree(cname,fname, "alg_stag_mom_meson_arg->meson_buf",
                            alg_stag_mom_meson_arg->meson_buf);
    sfree(alg_stag_mom_meson_arg->meson_buf);
    alg_stag_mom_meson_arg->meson_buf = 0;
  }
}



//------------------------------------------------------------------
//------------------------------------------------------------------
//
// AlgStagNucleon is derived from Alg and is relevant to  
// the staggered nucleon propagator.
// The type of fermion is determined by the argument to the 
// constructor. If the fermion type is not F_CLASS_STAG the 
// constructor will exit with a general error.
//
//------------------------------------------------------------------
//------------------------------------------------------------------

//------------------------------------------------------------------
// Constructor 
//------------------------------------------------------------------
AlgStagNucleon::AlgStagNucleon(Lattice& latt, 
			       CommonArg *c_arg,
			       StagNucleonArg *arg, Aots& a) :
			       Alg(latt, c_arg), aots(a) 
{
  cname = "AlgStagNucleon";
  char *fname = "AlgStagNucleon(L&,CommonArg*,StagNucleonArg*)";
  VRB.Func(cname,fname);

  // Check fermion type
  //----------------------------------------------------------------
  if(latt.Fclass() != F_CLASS_STAG)
    ERR.General(cname,fname, class_str, int(latt.Fclass()));


  // Initialize the argument pointer
  //----------------------------------------------------------------
  if(arg == 0)
    ERR.Pointer(cname,fname, "arg");
  alg_stag_nucleon_arg = arg;
}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
AlgStagNucleon::~AlgStagNucleon() {
  char *fname = "~AlgStagNucleon()";
  VRB.Func(cname,fname);
}


//------------------------------------------------------------------
// Calculate propagator
//------------------------------------------------------------------
void AlgStagNucleon::run()
{
  char *fname = "run()";
  VRB.Func(cname,fname);

  // Set the Lattice pointer
  //----------------------------------------------------------------
  Lattice& lat = AlgLattice();

  NucleonPropS n1(lat, *alg_stag_nucleon_arg);
  n1.getHadronPropS();

  // download propagator and AOTS 
  //----------------------------------------------------------------
  if(aots.begin()) {
    alg_stag_nucleon_arg->nucleon_buf = 
		(Float *)smalloc(n1.propLenTotal()*sizeof(Float));
    if(alg_stag_nucleon_arg->nucleon_buf == 0)
          ERR.Pointer(cname,fname, "alg_stag_nucleon_arg->nucleon_buf");
    VRB.Smalloc(cname,fname, "alg_stag_nucleon_arg->nucleon_buf", 
		alg_stag_nucleon_arg->nucleon_buf, n1.propLenTotal()*sizeof(Float));
    zero_buffer(alg_stag_nucleon_arg->nucleon_buf, n1.propLenTotal());
  }

  n1.download_prop(SNUCLEON, alg_stag_nucleon_arg->nucleon_buf);

  if (aots.last()) {
    char *data_file = common_arg->results ? 
		      (char *)common_arg->results : CAST_AWAY_CONST("nucleon.def"); 
    FILE *fp;
    if( NULL == (fp = fopen(data_file, "a")) ) {
      ERR.FileA(cname,fname, data_file);
    }
   
    write_to_file(alg_stag_nucleon_arg->nucleon_buf, fp, 
		  n1.propLenTotal(), aots.numSlices(), 1, SNUCLEON, n1.bcd());
    fclose(fp);

    VRB.Sfree(cname,fname, "alg_stag_nucleon_arg->nucleon_buf", 
			    alg_stag_nucleon_arg->nucleon_buf);
    sfree(alg_stag_nucleon_arg->nucleon_buf);
    alg_stag_nucleon_arg->nucleon_buf = 0;
  }
}



//------------------------------------------------------------------
//------------------------------------------------------------------
//
// AlgStagNonLocal is derived from Alg and is relevant to  
// the staggered non-local hadron propagator.
// The type of fermion is determined by the argument to the 
// constructor. If the fermion type is not F_CLASS_STAG the 
// constructor will exit with a general error.
//
//------------------------------------------------------------------
//------------------------------------------------------------------

//------------------------------------------------------------------
// Constructor 
//------------------------------------------------------------------
AlgStagNonLocal::AlgStagNonLocal(Lattice& latt, 
				 CommonArg *c_arg,
				 StagNonLocalArg *arg, Aots& a) :
				 Alg(latt, c_arg), aots(a) 
{
  cname = "AlgStagNonLocal";
  char *fname = "AlgStagNonLocal(L&,CommonArg*,StagNonLocalArg*)";
  VRB.Func(cname,fname);

  // Check fermion type
  //----------------------------------------------------------------
  if(latt.Fclass() != F_CLASS_STAG)
    ERR.General(cname,fname, class_str, int(latt.Fclass()));


  // Initialize the argument pointer
  //----------------------------------------------------------------
  if(arg == 0)
    ERR.Pointer(cname,fname, "arg");
  alg_stag_non_local_arg = arg;
}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
AlgStagNonLocal::~AlgStagNonLocal() {
  char *fname = "~AlgStagNonLocal()";
  VRB.Func(cname,fname);
}


//------------------------------------------------------------------
// Calculate propagator
//------------------------------------------------------------------
void AlgStagNonLocal::run()
{
  char *fname = "run()";
  VRB.Func(cname,fname);

  // Set the Lattice pointer
  //----------------------------------------------------------------
  Lattice& lat = AlgLattice();

  NLocalPropS nlc(lat, *alg_stag_non_local_arg);
  nlc.getHadronPropS();

  // download propagator and AOTS 
  //----------------------------------------------------------------
  if(aots.begin()) {
    alg_stag_non_local_arg->nlocal_buf = 
		(Float *)smalloc(nlc.propLenTotal()*sizeof(Float));

    if(alg_stag_non_local_arg->nlocal_buf == 0)
          ERR.Pointer(cname,fname, "alg_stag_non_local_arg->nlocal_buf");
    VRB.Smalloc(cname,fname, "alg_stag_non_local_arg->nlocal_buf", 
		alg_stag_non_local_arg->nlocal_buf, 
		nlc.propLenTotal()*sizeof(Float));

    zero_buffer(alg_stag_non_local_arg->nlocal_buf, nlc.propLenTotal());
  }

  nlc.download_prop(SNONLOCAL, alg_stag_non_local_arg->nlocal_buf);

  if (aots.last()) {
    char *data_file = common_arg->results ? 
		      (char *)common_arg->results : CAST_AWAY_CONST("nonlocal.def"); 
    FILE *fp;
    if( NULL == (fp = fopen(data_file, "a")) ) {
      ERR.FileA(cname,fname, data_file);
    }

    write_to_file(alg_stag_non_local_arg->nlocal_buf, fp, 
		nlc.propLenTotal(), aots.numSlices(), 1, SNONLOCAL, nlc.bcd());
    fclose(fp);

    VRB.Sfree(cname,fname, "alg_stag_non_local_arg->nlocal_buf", 
			    alg_stag_non_local_arg->nlocal_buf);
    sfree(alg_stag_non_local_arg->nlocal_buf);
    alg_stag_non_local_arg->nlocal_buf = 0;
  }
}

CPS_END_NAMESPACE
