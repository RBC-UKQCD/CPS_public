/*
  File:  threept_prop_arg.h

  Defines a structure threept_prop_arg which holds pointers to QPropW
  objects to pass to an alg_threept object.

  M. Lightman 08/27/2009
*/

#ifndef INCLUDED_3PT_PROP_ARG_H
#define INCLUDED_3PT_PROP_ARG_H

#include "qpropw.h"

CPS_START_NAMESPACE

struct ThreePtPropArg {

  QPropW* q_light_tpi[20][2][4][3];
  QPropW* q_light_tK[20][2][20];
  QPropW* q_strange_tpi[20][2];
  QPropW* q_strange_tK[20][2][20];
  QPropW* h_mass[20]; /*!< Pointers to heavy quark propagators.*/

};

CPS_END_NAMESPACE

#endif /* !INCLUDED_3PT_PROP_ARG_H */
