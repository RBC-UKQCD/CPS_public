#include<config.h>
CPS_START_NAMESPACE
/*  pbp_arg.h */

/*  The structure type PbpArg holds the parameters relevant
    to the PsiBar Psi measurement. */

#ifndef INCLUDED_PBP_ARG_H
#define INCLUDED_PBP_ARG_H


#define MAX_PBP_MASSES  100     /* The maximum number of masses
                                for which PsiBarPsi is measured
                                (for the same source vector). 
                                Used when PatternType = ARRAY */


enum PatternType {ARRAY = 0,	/* Use masses from mass[] table. */
		  LIN   = 1,	/* Use step pattern using +=     */
		  LOG   = 2 };	/* Use step pattern using *=     */


struct PbpArg {

  PatternType	pattern_kind;	/* Specifies the pattern used
                                to scan the mass values. For
                                each mass value PsiBarPsi is measured
                                using the same source vector. */
 
  int n_masses;                 /* The number of masses for which
				PsiBarPsi is measured */

  Float	mass_start;		/* Relevant to the LIN and LOG patterns.
                                This is the starting mass for these 
                                patterns. */                   

  Float	mass_step;		/* Relevant to the LIN and LOG patterns.
                                This is the mass step size for these 
                                patterns:
                                          += mass_step	(LIN)
				          *= mass_step	(LOG)   */
 
  Float mass[MAX_PBP_MASSES];   /* Relevant to the ARRAY pattern.
				The array of mass values for which 
                                PsiBarPsi is measured. */

  int max_num_iter;             /* The maximum number of conjugate
			        gradient iterations to do. */


  Float stop_rsd;               /*  The residual for the stopping 
				condition. */

  int src_u_s;                  /* Relevant to Domain Wall Fermions only.
                                It is the location s on the wall of the 
                                upper two components of the source.
                                For spread out DWF this is the global
                                coordinate, i.e. its range is from 0 to ls-1. 
                                ls = GJP.SnodeSites() * GJP.Snodes()
                                For maximum localization use 0. */

  int src_l_s;                  /* Relevant to Domain Wall Fermions only.
                                It is the location s on the wall of the 
                                lower two components of the source. 
                                For spread out DWF this is the global
                                coordinate, i.e. its range is from 0 to ls-1. 
                                ls = GJP.SnodeSites() * GJP.Snodes()
                                For maximum localization use ls-1. */

  int snk_u_s;                  /* Relevant to Domain Wall Fermions only.
                                It is the location s on the wall of the 
                                upper two components of the sink.
                                For spread out DWF this is the global
                                coordinate, i.e. its range is from 0 to ls-1. 
                                ls = GJP.SnodeSites() * GJP.Snodes()
                                For maximum localization use ls-1.
                                Not used if snk_loop set. */

  int snk_l_s;                  /* Relevant to Domain Wall Fermions only.
                                It is the location s on the wall of the 
                                lower two components of the sink. 
                                For spread out DWF this is the global
                                coordinate, i.e. its range is from 0 to ls-1. 
                                ls = GJP.SnodeSites() * GJP.Snodes()
                                For maximum localization use 0.
                                Not used if snk_loop set. */

  int snk_loop;                /* Relevant to Domain Wall Fermions only.
                                If this flag is set ( i.e. != 0 ) then 
                                the sink of the upper spin components
                                takes values s and the sink of the
                                lower spin components takes values ls-s-1
                                for s = 0, 1, ..., ls-1.
                                For spread out DWF this is the global
                                coordinate.
                                ls = GJP.SnodeSites() * GJP.Snodes() */

};

#endif /* !INCLUDED_PBP_ARG_H */
CPS_END_NAMESPACE
