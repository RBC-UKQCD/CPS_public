/*!\file
  \brief  Definition of the PbpArg structure.

  $Id: pbp_arg.x,v 1.3 2005-05-18 06:17:35 chulwoo Exp $
*/
/*---------------------------------------------------------------------------*/

/*! The maximum number of masses for which the condensate can be measured.*/
enum PBP_LIM { MAX_PBP_MASSES = 100     };
                                /* Used when PatternType = ARRAY */

/*! A structure holding the parameters relevant to the condensate measurement.*/
/*!  \ingroup algargs */
class PbpArg {

  PatternType	pattern_kind;	/*!< Specifies the pattern used
                                to obtain the mass values. For
                                each mass value the condensate is measured
                                using the same source vector. */
 
  int n_masses;                 /*!< The number of masses for which
				the condensate is measured */

  Float	mass_start;		/*!< 
				  This is the starting mass of the progressions
				  defining the masses in the ::LIN and ::LOG
				  patterns:
				*/                   

  Float	mass_step;		/*!<
				  This is the increment in the mass for the
				  arithmetic progression, or the factor of the
				  mass in the geometric progression
				  defining the masses in the ::LIN and ::LOG
				  patterns:
				  
				*/

  Float mass[MAX_PBP_MASSES];   /*!< Relevant to the ::ARRAY pattern:
				The array of mass values for which 
                                the condensate is measured. */

  int max_num_iter;             /*!< The maximum number of conjugate
			        gradient iterations to do. */


  Float stop_rsd;               /*!<  The target residual for the stopping 
				condition. */

  int src_u_s;                  /*!< Relevant to Domain Wall Fermions only.
				   It is the global 5th direction coordinate s
				   of the location on the wall of the
                                 location on the wall of the 
                                upper two components of the source.
                                For maximum localization use 0. */

    int src_l_s;                  /*!< Relevant to Domain Wall Fermions only.
				   It is the global 5th direction coordinate s
				   of the location on the wall of the
                                location  on the wall of the  
                                lower two components of the source. 
                                For maximum localization use the maximum value Ls-1. */

  int snk_u_s;                  /*!< Relevant to Domain Wall Fermions only.
				   It is the global 5th direction coordinate s
				   of the location on the wall of the
                                upper two components of the sink.
                                For maximum localization use the maximum value
				Ls-1.
                                Not used if snk_loop set. */

  int snk_l_s;                  /*!< Relevant to Domain Wall Fermions only.
				   It is the global 5th direction coordinate s
				   of the location on the wall of the
                                lower two components of the sink. 
                                For spread out DWF this is the global
                                coordinate, i.e. its range is from 0 to ls-1. 
                                ls = GJP.SnodeSites() * GJP.Snodes()
                                For maximum localization use 0.
                                Not used if snk_loop set. */

  int snk_loop;                /*!< Relevant to Domain Wall Fermions only.
                                If this flag is set ( i.e. != 0 ) then 
                                the sink of the upper spin components
                                takes values s and the sink of the
                                lower spin components takes values Ls-s-1
                                for s = 0, 1, ..., Ls-1, where s is the
				global lattice 5th direction coordinate,
			       */
};

