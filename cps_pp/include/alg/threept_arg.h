#include <config.h>
/*!\file
  \brief Definition of the ThreePtArg structure.

  $Id: threept_arg.h,v 1.6 2010-04-05 19:55:48 chulwoo Exp $ 
*/
//---------------------------------------------------------------------------

CPS_START_NAMESPACE
#ifndef INCLUDED_3PT_ARG_H
#define INCLUDED_3PT_ARG_H          //!< Prevent multiple inclusion.
CPS_END_NAMESPACE
#include <alg/cg_arg.h>
#include <util/vector.h>
#include <alg/qpropw_arg.h>
#include <util/vml/types.h>
#include <util/vml/vml.h>
#include <util/enum.h>
#include <util/defines.h>
CPS_START_NAMESPACE

//! The structure holds parameters relevant to the three-point functions measurement.
/*!  \ingroup algargs */
class VML;
class ThreePtArg {
 public:

  bool Encode(char *filename,char *instance);
  bool Decode(char *filename,char *instance);
  bool Vml(VML *vmls,char *instance);

  char *results; /*!< filename for output */
  char *results_mres_ZA; /*!< filename for mres/ZA output */
  char *results_pipi; /*!< filename for Pi Pi correlator output */

  CgArg cg;		/*!< Parameters for fermion matrix inversion to
				  compute the quark propagator. */

  RandomType rng;      /*!< Type of RNG for random source */

  int gauge_fix; /*!< Type of gauge fixing to apply */

  int t_src;	/*!< Source timeslice.  We have been setting this to 0 and shifting the lattice instead */
  int t_snk;    /*!< Sink timeslice.  This would be 64 if the source timeslice is 0 and you are doing P+A and P-A */
  int t_shift;  /*!< Amount to shift lattice.  Must be a multiple of the number of t sites on a node. Only does anything if you're generating propagators */
  int box_len;  /*!< Length of box source cube */

  int t_op;		/*!< The location of the random propagator */
  int width;    /*!< Width of random propagator slab (starting at t_op) */
  int num_hits; /*!< Number of random propagators to do */
  int do_susy;  /*!< Compute non-vector matrix elements? */
  
  int num_light;	/*!< The number of light masses */
  Float l_mass[20];	/*!< The list of light masses */
  int num_strange;      /*!< The number of strange masses */
  Float s_mass[20];     /*!< The list of strange masses */
  int num_heavy;    /*!< The number of heavy masses */
  Float h_mass[20]; /*!< The list of heavy masses */
  
  int num_tK;   /*!< The number of different kaon source times */
  int tK[20];   /*!< The list of kaon source times */

  int do_zero_mom;  /*!< Whether or not to do contractions with zero momentum pions. */
  int do_first_mom;  /*!< Whether or not to do the first non-zero lattice momentum. */
  int do_second_mom; /*!< Whether or not to do the second non-zero lattice momentum (i.e. sqrt(2) times the first non-zero lattice momentum). */
  int do_third_mom; /*!< Whether or not to do the third non-zero lattice momentum (i.e. sqrt(3) times the first non-zero lattice momentum). */
  int do_pipi_non_zero_tot_mom; /*!< Whether or not to do two pion correlators where the two pions have non-zero total momentum. */

  int do_p_plus_a_kaon; /*!< Whether or not to use P+A (set this flag to 1) or just P (set this flag to 0) for the strange quark propagator in the kaon. */
  int do_kaon_at_walls; /*!< Whether or not to do three point functions with kaons with sources at 0 and time_size. (If not, just do kaons with sources at tK[j] in the three point functions). */
  int do_kaons_tK; /*!< Whether or not to generate light propagators with sources at tK[j] in order that kaon correlators with sources at tK[j] can be computed. */

  int chkpoints; /*!< Whether or not to do timing checkpoints in the run function. */

  char* ensemble_label;
  char* ensemble_id;
  int seqNum; /*!< Trajectory number. */

  ThreePtArg ( ) ;
 
};

/* the xdr functions */

#ifdef __cplusplus
extern "C" {
#endif

#if defined(__STDC__) || defined(__cplusplus)
  extern  bool_t vml_ThreePtArg (VML *, char *instance, ThreePtArg*);

#else /* K&R C */
  extern  bool_t vml_ThreePtArg (VML *, char *instance, ThreePtArg*);

#endif /* K&R C */

#ifdef __cplusplus
}
#endif


#endif /* !INCLUDED_3PT_ARG_H */
CPS_END_NAMESPACE
