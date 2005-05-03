/*!
  Coefficients for an asqtad style smearing of the
  lattice. If orthog is 0 to 3 then the smearing will 
  only be in the direction orthogonal to that direction.
*/
class KineticSmearArg
{
  int   orthog;
  Float single_link;
  Float three_link;
  Float five_link;
  Float seven_link;
  Float lepage;

  memfun KineticSmearArg();
};
