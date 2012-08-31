#include <config.h>
/*!\file
  \brief  Definitions of the AlgThreePt class.

  $Id: alg_threept.h,v 1.7 2012-08-31 04:55:08 chulwoo Exp $
*/
//------------------------------------------------------------------

CPS_START_NAMESPACE
#ifndef INCLUDED_ALG_3PT_H
#define INCLUDED_ALG_3PT_H
CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/smalloc.h>
#include <util/pmalloc.h>
#include <alg/alg_base.h>
#include <alg/common_arg.h>

//#include <alg/threept_arg.h>
#include "threept_arg.h"
#include "threept_prop_arg.h"
//#include <alg/qpropw.h>
#include "qpropw.h"

CPS_START_NAMESPACE

//------------------------------------------------------------------
//! Class implementing the calculation of meson three point functions.
/*!
  This uses Wilson, clover or domain wall fermions.

  \ingroup alg
*/
// PLEASE NOTE that the propagators passed to this method will
// be altered to change propagators with periodic (P) and 
// antiperiodic (A) boundary conditions into P+A and P-A if
// necessary.  This is done so that extra memory isn't taken
// up storing separate P+A and P-A propagators.
//------------------------------------------------------------------
class AlgThreePt : public Alg
{
 private:
    char* cname;

    ThreePtArg* alg_threept_arg;
        // arguments for the three point calculation
    ThreePtPropArg* alg_threept_prop_arg;
        // pointers to propagators used in three point calculation
    int f_size;
        // Node checkerboard size of the fermion field
	enum GammaType {
	  SCALAR = 0, //  no gamma indices
	  VECTOR,     // one gamma index
	  TENSOR      // two gamma indices
	};
	typedef enum GammaType GammaType;
	enum TraceType {
	  TR = 0,    // one full trace  (spin and color)
	  TRTR,      // two full traces 
	  TR_MX,     // one spin trace  (two color traces)
	  TRTR_MX    // two spin traces (one color trace)
	};
	typedef enum TraceType TraceType;

	char gam [3][5];  // strings describing gamma structure
	char gam2[3][5];  
	char tra [4][10];

	int do_susy;
	FILE *fp;
	FILE *fp_mres_ZA;
	FILE *fp_pipi;

 public:
    AlgThreePt(Lattice & latt, CommonArg* c_arg, ThreePtArg* arg, ThreePtPropArg* p_arg);

    virtual ~AlgThreePt();

    void run(void);

    void pipi(QPropW& q_src);
    void pipi(QPropW& q_u, QPropW& q_d_mom, int mom_num, int mom_dir);
 
    void     spectrum(QPropW& q1_src, QPropW& q2_src, int bc=2);
    void     spectrum(QPropW& q_s, QPropW& q_u_mom, int mom_num, int mom_dir);
    void box_spectrum(QPropW& q1_src, QPropW& q2_src, int bc=2);
    void box_spectrum(QPropW& q_s, QPropW& q_u_mom, int mom_num, int mom_dir);

    void k_to_vac(QPropW& q_str, QPropW& q_src,
				  QPropW& q_ppl);

    void figure8(QPropW& q_str, QPropW& q_src,
				 QPropW& q_snk1, QPropW& q_snk2,
				 int is_light=0);

    void eye(QPropW& q_str, QPropW& q_spc,
			 QPropW& q_snk,
			 QPropW& q_ppl);
// copied from 5_0_3-wme to ensure backward compatibillity, bu CJ
    void figure8_spectator_old(QPropW& q_str, QPropW& q_spc,
						   QPropW& q_snk1, QPropW& q_snk2, QPropW& q_snk3);

    void figure8_spectator(QPropW& q_str, QPropW& q_spc,
						   QPropW& q_snk1, QPropW& q_snk2, QPropW& q_snk3, int mom_num=0, int mom_dir=0);

    void figure8_vacuum(QPropW& q_str, QPropW& q_src,
						QPropW& q_snk1, QPropW& q_vac, QPropW& q_snk2);

	void eye_vacuum(QPropW& q_str, QPropW& q_spc,
					QPropW& q_vac, QPropW& q_snk,
					QPropW& q_ppl);

};


#endif
CPS_END_NAMESPACE
