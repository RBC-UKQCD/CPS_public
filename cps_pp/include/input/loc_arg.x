class LocArg {

  PatternType	pattern_kind;	/*!< Specifies the pattern used
                                to obtain the locations . For
                                each mass value the condensate is measured
                                using the same source vector. */
 
  int n_loc; 
  int	start;
  int	step;
  int locs<>;
  int max_src; /* number of maximum measurement per job */
};

