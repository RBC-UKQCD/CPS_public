
class FixGaugeArg {
  FixGaugeType fix_gauge_kind;   /* The kind of gauge fixing*/
  int hyperplane_start;          /* The full lattice coordinate of the first */
                                 /* hyperplane of gauge fixing matrices.*/
  int hyperplane_step;           /* The coordinate step between hyperplanes.*/
  int hyperplane_num;            /* The number of hyperplanes.*/
  Float stop_cond;               /* The stopping condition.*/
  int max_iter_num;              /* Maximum number of iterations.*/
                                 /* It exits if reached. If is set to 0 then */
                                 /* the number of interations is not checked*/
                                 /* and will continue untill the stopping */
                                 /* condition is satisfied.*/
};

