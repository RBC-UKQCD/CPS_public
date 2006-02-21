#include<config.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------
//
// alg_s_spect.h
//
// Header file for all alg classes relevant to staggered fermion
// spectroscopy. The type of glue or fermion is given as
// an argument of type Lattice& to the constructor. If the 
// fermion type is not F_CLASS_STAG the constructors exit with
// a general error.
//
//------------------------------------------------------------------


#ifndef INCLUDED_ALG_S_SPECT_H
#define INCLUDED_ALG_S_SPECT_H

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/smalloc.h>
#include <util/pmalloc.h>
#include <alg/alg_base.h>
#include <alg/quark_prop_s.h>
#include <alg/common_arg.h>
#include <alg/s_spect_arg.h>
CPS_START_NAMESPACE

// Forward declaration
class Aots;

//------------------------------------------------------------------
//
// AlgStagQuark is derived from Alg and is relevant to  
// the staggered quark propagator.
// The type of fermion is determined by the argument to the 
// constructor. If the fermion type is not F_CLASS_STAG the 
// constructor will exit with a general error.
//
//------------------------------------------------------------------
class AlgStagQuark : public Alg
{
 private:
    char *cname;

    StagQuarkArg *alg_stag_quark_arg;
        // The argument structure for
        // the quark propagator

    static QuarkPropSMng sqpm; 
	// The quark propagator manager is initilaized
	// before entering the main function, and responsible for
	// registering and destroying quark propagators 

 public:
    AlgStagQuark(Lattice & latt, CommonArg *c_arg, 
		 StagQuarkArg *arg, Aots &aots);

    virtual ~AlgStagQuark();

    void run(void);

    void free(void);
};


//------------------------------------------------------------------
//
// AlgStagMeson is derived from Alg and is relevant to  
// the staggered meson propagator.
// The type of fermion is determined by the argument to the 
// constructor. If the fermion type is not F_CLASS_STAG the 
// constructor will exit with a general error.
//
//------------------------------------------------------------------
class AlgStagMeson : public Alg
{
 private:
    char *cname;

    StagMesonArg *alg_stag_meson_arg;
        // The argument structure for
        // the meson propagator

    Aots & aots;
	// Control info to Average Over Time Slices

 public:
    AlgStagMeson(Lattice & latt, CommonArg *c_arg, 
	         StagMesonArg *arg, Aots& a);

    virtual ~AlgStagMeson();

    void run(void);

};

//------------------------------------------------------------------
// added by manke: 05/04/2001
// AlgStagMomMeson is derived from Alg and is relevant to  
// the staggered meson propagator.
// The type of fermion is determined by the argument to the 
// constructor. If the fermion type is not F_CLASS_STAG the 
// constructor will exit with a general error.
//
//------------------------------------------------------------------
class AlgStagMomMeson : public Alg
{
 private:
    char *cname;

    StagMomMesonArg *alg_stag_mom_meson_arg;
        // The argument structure for
        // the meson propagator

    Aots & aots;
	// Control info to Average Over Time Slices

 public:
    AlgStagMomMeson(Lattice & latt, CommonArg *c_arg, 
	         StagMomMesonArg *arg, Aots& a);

    virtual ~AlgStagMomMeson();

    void run(void);

};


//------------------------------------------------------------------
//
// AlgStagNucleon is derived from Alg and is relevant to  
// the staggered nucleon propagator.
// The type of fermion is determined by the argument to the 
// constructor. If the fermion type is not F_CLASS_STAG the 
// constructor will exit with a general error.
//
//------------------------------------------------------------------
class AlgStagNucleon : public Alg
{
 private:
    char *cname;

    StagNucleonArg *alg_stag_nucleon_arg;
        // The argument structure for
        // the nucleon propagator
 
    Aots & aots;
	// Control info to Average Over Time Slices

 public:
    AlgStagNucleon(Lattice & latt, CommonArg *c_arg, 
		   StagNucleonArg *arg, Aots& a);

    virtual ~AlgStagNucleon();

    void run(void);

};



//------------------------------------------------------------------
//
// AlgStagNonLocal is derived from Alg and is relevant to  
// the staggered non-local hadron propagator.
// The type of fermion is determined by the argument to the 
// constructor. If the fermion type is not F_CLASS_STAG the 
// constructor will exit with a general error.
//
//------------------------------------------------------------------
class AlgStagNonLocal : public Alg
{
 private:
    char *cname;

    StagNonLocalArg *alg_stag_non_local_arg;
        // The argument structure for
        // the non-local hadron propagator
 
    Aots & aots;
	// Control info to Average Over Time Slices

 public:
    AlgStagNonLocal(Lattice & latt, CommonArg *c_arg, 
		    StagNonLocalArg *arg, Aots& a);

    virtual ~AlgStagNonLocal();

    void run(void);

};


class AlgNLStagMeson : public Alg
{
 private:
    char *cname;

    NLStagMesonArg *alg_stag_non_local_arg;
        // The argument structure for
        // the non-local hadron propagator
 
    Aots & aots;
        // Control info to Average Over Time Slices

 public:
    AlgNLStagMeson(Lattice & latt, CommonArg *c_arg, 
                    NLStagMesonArg *arg, Aots& a);

    virtual ~AlgNLStagMeson();

    void run(void);

};




#endif





CPS_END_NAMESPACE
