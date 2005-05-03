/*!
  Arguments for Hyp smearing; tolerance is the
  stopping condition on the su3 projection, if orthog
  is 0 to 3 then the smearing will only be in the
  direction orthogonal to that direction.
  c1 to c3 are the smearing coefficients.
*/
class HypSmearArg
{
  Float tolerance;
  int   orthog;
  Float c1;
  Float c2;
  Float c3;

  memfun HypSmearArg();
};

