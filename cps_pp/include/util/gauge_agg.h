#ifndef GAUGE_AGG_H
#define GAUGE_AGG_H
#ifdef SCIDAC
#include <asq_data_types.h>
#else
#include <util/asq_data_types.h>
#endif
struct gauge_agg{
  int src;
  int dest;
  IFloat mat[18];
};

typedef struct gauge_agg_cb{
  //Index for initial position of field
  int src;
  //Index for "transported field"
  int dest;
  //Index for the gauge link
  int gauge;
} ind_agg;

struct hop_pointer {
  int src;
  int dest;
};
//CPS_END_NAMESPACE
#endif
