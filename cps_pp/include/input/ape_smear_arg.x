/*! Input paramters for AlgApeSmear */
class ApeSmearArg
{     
  //! tolerance for the SU(3) projection
  Float tolerance;

  //! smear in hyper-plane orthoganal to this direction
  int   orthog;

  //! ape smearing coefficient 
  Float coef;
  
  memfun ApeSmearArg();
};
