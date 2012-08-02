#ifndef GAUGE_AGG_H
#define GAUGE_AGG_H
#ifdef SCIDAC
#include <asq_data_types.h>
#else
#include <util/asq_data_types.h>
#endif

struct gauge_agg{
  unsigned long src;
  unsigned long dest;
  IFloat mat[18];
};

typedef struct gauge_agg_cb{
  //Index for initial position of field
  unsigned long src;
  //Index for "transported field"
  unsigned long dest;
  //Index for the gauge link
  unsigned long gauge;
} ind_agg;

struct hop_pointer {
  unsigned long src;
  unsigned long dest;
};
//CPS_END_NAMESPACE
#endif
