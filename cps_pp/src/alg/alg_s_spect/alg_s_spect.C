#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2008/02/08 18:35:05 $
//  $Header: /space/cvs/cps/cps++/src/alg/alg_s_spect/alg_s_spect.C,v 1.13 2008/02/08 18:35:05 chulwoo Exp $
//  $Id: alg_s_spect.C,v 1.13 2008/02/08 18:35:05 chulwoo Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $RCSfile: alg_s_spect.C,v $
//  $Revision: 1.13 $
//  $Source: /space/cvs/cps/cps++/src/alg/alg_s_spect/alg_s_spect.C,v $
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
#include <util/qcdio.h>
#include <alg/quark_prop_s.h>
#include <alg/meson_prop_s.h>
#include <alg/mom_meson_p_s.h>
#include <alg/nucl_prop_s.h>
#include <alg/nlocal_prop_s.h>
#include <alg/nlocal_propmes_s.h>
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
#include <alg/myenum.h>

#ifdef PARALLEL
#include <comms/sysfunc_cps.h>
#endif

CPS_START_NAMESPACE

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

#if TARGET==cpsMPI
    using MPISCU::fprintf;
#endif

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

  
  int i, k;

  //-------------------------------------------------------------
  // for the Non Local staggered meson, just print out every result!
  // other modification is done by anlysis script 
  //-------------------------------------------------------------
  if (type == NLSTAG) {
   Float tmp;    
    for(i = 0 ; i < num_IFloats ; i+=2) {
      //printf("i = %d\n",i);
      tmp = data_p[i] / num_slices;
      Fprintf(fp,"%e    ",(IFloat)tmp);
      tmp = data_p[i+1] / num_slices;
      Fprintf(fp,"%e\n",(IFloat)tmp);
      //printf("data_p[%d] = %f data_p[%d] = %f\n",i,data_p[i]/num_slices,i+1,data_p[i+1]/num_slices);
    }
    return ;
  }
  //-------------------------------------------------------------
  // * write C(t = 0)
  //-------------------------------------------------------------
  for (i = 0; i < unit; ++i) {
    Float tmp = data_p[2*i] / num_slices;
    Fprintf(fp, "%e\n", (IFloat)tmp);
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
      Fprintf(fp, "%e\n", (IFloat)tmp);
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

      Fprintf(fp, "%e\n", (IFloat)tmp);
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

      Fprintf(fp, "%e\n", (IFloat)tmp);
    }
  }

  //-------------------------------------------------------------
  // * write C(t=N_t/2)
  //-------------------------------------------------------------
  for (i = half_length; i < half_length + unit; ++i) {
    Float tmp = data_p[2*i] / num_slices;
    Fprintf(fp, "%e\n", (IFloat)tmp);
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
  if(!latt.FstagType())
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
//  Lattice& lat = AlgLattice();

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
  if(!latt.FstagType())
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
		(unsigned long)smalloc(m1.propLenTotal()*sizeof(Float));
    if(alg_stag_meson_arg->meson_buf == 0)
          ERR.Pointer(cname,fname, "alg_stag_meson_arg->meson_buf");
    VRB.Smalloc(cname,fname, "alg_stag_meson_arg->meson_buf", 
		(void *)alg_stag_meson_arg->meson_buf, m1.propLenTotal()*sizeof(Float));
    zero_buffer((Float *)alg_stag_meson_arg->meson_buf, m1.propLenTotal());
  }

  m1.download_prop(SMESON, (Float *)alg_stag_meson_arg->meson_buf);

  if (aots.last()) {
    char *data_file = common_arg->results ? 
		      (char *)common_arg->results : CAST_AWAY_CONST("meson.def"); 
    FILE *fp;
    if( NULL == (fp = Fopen(data_file, "a")) ) {
      ERR.FileA(cname,fname, data_file);
    }

    write_to_file((Float *)alg_stag_meson_arg->meson_buf, fp, 
   	          m1.propLenTotal(), aots.numSlices(), 1, SMESON, m1.bcd());
    Fclose(fp);

    VRB.Sfree(cname,fname, "alg_stag_meson_arg->meson_buf", 
			    (void *)alg_stag_meson_arg->meson_buf);
    sfree((Float *)alg_stag_meson_arg->meson_buf);
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
  if(!latt.FstagType())
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
      (unsigned long)smalloc(m1.propLenTotal()*sizeof(Float));
    if(alg_stag_mom_meson_arg->meson_buf == 0)
      ERR.Pointer(cname,fname, "alg_stag_mom_meson_arg->meson_buf");
    VRB.Smalloc(cname,fname, "alg_stag_mom_meson_arg->meson_buf",
         (void *)alg_stag_mom_meson_arg->meson_buf, m1.propLenTotal()*sizeof(Float));
    zero_buffer((Float *)alg_stag_mom_meson_arg->meson_buf, m1.propLenTotal());
  }


  m1.download_prop(SMOMMESON, (Float *)alg_stag_mom_meson_arg->meson_buf);

  if (aots.last()) {
    char *data_file = common_arg->results ?
                      (char *)common_arg->results : (char *) "mom_meson.def";
    FILE *fp;
    if( NULL == (fp = Fopen(data_file, "a")) ) {
      ERR.FileA(cname,fname, data_file);
    }

    write_to_file((Float *)alg_stag_mom_meson_arg->meson_buf, fp,
                  m1.propLenTotal(), aots.numSlices(), 
	alg_stag_mom_meson_arg->no_of_momenta,SMOMMESON, m1.bcd());
    Fclose(fp);

    VRB.Sfree(cname,fname, "alg_stag_mom_meson_arg->meson_buf",
                            (void *)alg_stag_mom_meson_arg->meson_buf);
    sfree((void *)alg_stag_mom_meson_arg->meson_buf);
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
  if(!latt.FstagType())
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
		(unsigned long)smalloc(n1.propLenTotal()*sizeof(Float));
    if(alg_stag_nucleon_arg->nucleon_buf == 0)
          ERR.Pointer(cname,fname, "alg_stag_nucleon_arg->nucleon_buf");
    VRB.Smalloc(cname,fname, "alg_stag_nucleon_arg->nucleon_buf", 
		(void *)alg_stag_nucleon_arg->nucleon_buf, n1.propLenTotal()*sizeof(Float));
    zero_buffer((Float *)alg_stag_nucleon_arg->nucleon_buf, n1.propLenTotal());
  }

  n1.download_prop(SNUCLEON, (Float *)alg_stag_nucleon_arg->nucleon_buf);

  if (aots.last()) {
    char *data_file = common_arg->results ? 
		      (char *)common_arg->results : CAST_AWAY_CONST("nucleon.def"); 
    FILE *fp;
    if( NULL == (fp = Fopen(data_file, "a")) ) {
      ERR.FileA(cname,fname, data_file);
    }
   
    write_to_file((Float *)alg_stag_nucleon_arg->nucleon_buf, fp, 
		  n1.propLenTotal(), aots.numSlices(), 1, SNUCLEON, n1.bcd());
    Fclose(fp);

    VRB.Sfree(cname,fname, "alg_stag_nucleon_arg->nucleon_buf", 
			    (void *)alg_stag_nucleon_arg->nucleon_buf);
    sfree((void *)alg_stag_nucleon_arg->nucleon_buf);
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
  if(!latt.FstagType())
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
		(unsigned long)smalloc(nlc.propLenTotal()*sizeof(Float));

    if(alg_stag_non_local_arg->nlocal_buf == 0)
          ERR.Pointer(cname,fname, "alg_stag_non_local_arg->nlocal_buf");
    VRB.Smalloc(cname,fname, "alg_stag_non_local_arg->nlocal_buf", 
		(void *)alg_stag_non_local_arg->nlocal_buf, 
		nlc.propLenTotal()*sizeof(Float));

    zero_buffer((Float *)alg_stag_non_local_arg->nlocal_buf, nlc.propLenTotal());
  }

  nlc.download_prop(SNONLOCAL, (Float *)alg_stag_non_local_arg->nlocal_buf);

  if (aots.last()) {
    char *data_file = common_arg->results ? 
		      (char *)common_arg->results : CAST_AWAY_CONST("nonlocal.def"); 
    FILE *fp;
    if( NULL == (fp = Fopen(data_file, "a")) ) {
      ERR.FileA(cname,fname, data_file);
    }

    write_to_file((Float *)alg_stag_non_local_arg->nlocal_buf, fp, 
		nlc.propLenTotal(), aots.numSlices(), 1, SNONLOCAL, nlc.bcd());
    Fclose(fp);

    VRB.Sfree(cname,fname, "alg_stag_non_local_arg->nlocal_buf", 
			    (void *)alg_stag_non_local_arg->nlocal_buf);
    sfree((void *)alg_stag_non_local_arg->nlocal_buf);
    alg_stag_non_local_arg->nlocal_buf = 0;
  }
}

//------------------------------------------------------------------
//------------------------------------------------------------------
//
// AlgNLStagMeson is derived from Alg and is relevant to  
// the staggered non-local hadron propagator.
// The type of fermion is determined by the argument to the 
// constructor. If the fermion type is not F_CLASS_STAG the 
// constructor will exit with a general error.
// This class use NLMesonPropS class writted by chateau
//
//------------------------------------------------------------------
//------------------------------------------------------------------

//------------------------------------------------------------------
// Constructor 
//------------------------------------------------------------------
AlgNLStagMeson::AlgNLStagMeson(Lattice& latt, 
				 CommonArg *c_arg,
				 NLStagMesonArg *arg, Aots& a) :
				 Alg(latt, c_arg), aots(a) 
{
  cname = "AlgNLStagMeson";
  char *fname = "AlgNLStagMeson(L&,CommonArg*,NLStagMesonArg*)";
  VRB.Func(cname,fname);

  // Check fermion type
  //----------------------------------------------------------------
  /*if(latt.Fclass() != F_CLASS_STAG)
    ERR.General(cname,fname, class_str, int(latt.Fclass()));*/


  // Initialize the argument pointer
  //----------------------------------------------------------------
  if(arg == 0)
    ERR.Pointer(cname,fname, "arg");
  alg_stag_non_local_arg = arg;
}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
AlgNLStagMeson::~AlgNLStagMeson() {
  char *fname = "~AlgNLStagMeson()";
  VRB.Func(cname,fname);
}


//------------------------------------------------------------------
// Calculate propagator
//------------------------------------------------------------------
void AlgNLStagMeson::run()
{
  char *fname = "run()";
  VRB.Func(cname,fname);

  // Set the Lattice pointer
  //----------------------------------------------------------------
  Lattice& lat = AlgLattice();

  NLSMesonPropS nlsm(lat, *alg_stag_non_local_arg);
  nlsm.getHadronPropS();

  // download propagator and AOTS 
  //----------------------------------------------------------------
  if(aots.begin()) {
    alg_stag_non_local_arg->nlocal_buf = 
		(unsigned long)smalloc(nlsm.propLenTotal()*sizeof(Float));

    if(alg_stag_non_local_arg->nlocal_buf == 0)
          ERR.Pointer(cname,fname, "alg_stag_non_local_arg->nlocal_buf");
    VRB.Smalloc(cname,fname, "alg_stag_non_local_arg->nlocal_buf", 
		(void *)alg_stag_non_local_arg->nlocal_buf, 
		nlsm.propLenTotal()*sizeof(Float));

    zero_buffer((Float *)alg_stag_non_local_arg->nlocal_buf, nlsm.propLenTotal());
  }

//  nlsm.download_prop(SNONLOCAL, alg_stag_non_local_arg->nlocal_buf);
  nlsm.download_prop(NLSTAG, (Float *)alg_stag_non_local_arg->nlocal_buf);

  if (aots.last()) {
    char *data_file = common_arg->results ? 
		      (char *)common_arg->results : CAST_AWAY_CONST("nonlocal.def"); 
    FILE *fp;
    if( NULL == (fp = Fopen(data_file, "a")) ) {
      ERR.FileA(cname,fname, data_file);
    }

    write_to_file((Float *)alg_stag_non_local_arg->nlocal_buf, fp, 
		  nlsm.propLenTotal(), aots.numSlices(), 1, NLSTAG, nlsm.bcd());
    Fclose(fp);

    VRB.Sfree(cname,fname, "alg_stag_non_local_arg->nlocal_buf", 
			    (void *)alg_stag_non_local_arg->nlocal_buf);
    sfree((void *)alg_stag_non_local_arg->nlocal_buf);
    alg_stag_non_local_arg->nlocal_buf = 0;
  }
}

CPS_END_NAMESPACE
