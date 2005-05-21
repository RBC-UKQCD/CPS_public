/*  pot_arg.h*/

class PotArg {
  int prop_dir;    
  /* "temporal" direction along which the potential*/
  /* is measure as exponential decay*/

  int therm_steps; /* additional thermalisation steps for each new xi0*/
  int therm_meas;  /* measure wline and plaq every therm_meas steps*/
  int xi; /* renormalised anisotropy*/
  int check; /* check W(z,x=check) against W(z,xi*t)*/
  int sweep; /* sweep number which is needed in output*/
  int max_X; /* maximal extent of Q-Q pair in X-direction*/
  int max_Y; /* maximal extent of Q-Q pair in Y-direction*/
  int max_Z; /* maximal extent of Q-Q pair in Z-direction*/
  int max_T; /* maximal extent of Q-Q pair in T-direction*/

};


