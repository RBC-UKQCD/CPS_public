/*
  Class used to return information regarding the MD force.
  L1 - The L1 norm of the force
  L2 - The L2 norm of the force
  Linf - The Linfinity norm of the force

  Placed print routine in this class to minimise rewriting if code evolves.
*/

#ifndef FORCE_ARG_H
#define FORCE_ARG_H

#include <config.h>
#include <alg/enum.h>

CPS_START_NAMESPACE

class ForceArg {
public:
  char *cname;
  Float L1;
  Float L2;
  Float Linf;

  ForceArg();
  ForceArg(Float L1, Float L2, Float Linf);
  ~ForceArg();

  void print(Float dt, char *label);
};

CPS_END_NAMESPACE

#endif /* !FORCE_ARG_H */
