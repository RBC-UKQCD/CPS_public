/*  mom_arg.h*/

class MomArg {
  int no_of_momenta; /* number of different momenta    */
  int deg;           /* control flag: average over degenerate momenta on/off*/
  int dir;           /* propagation direction*/
  int src_begin[4];  /* source */
  int src_end[4];

};
