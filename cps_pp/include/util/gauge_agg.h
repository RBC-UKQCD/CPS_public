#ifndef GAUGE_AGG_H
#define GAUGE_AGG_H
#include <util/asq_data_types.h>
struct gauge_agg{
  int src;
  int dest;
  IFloat mat[18];
};

struct gauge_agg_cb{
  //Index for initial position of field
  int src;
  //Index for "transported field"
  int dest;
  //Index for the padded field.  This assumes that the field will be
  //transported in all 8 possible directions.  Instead of storing
  //each direction in a single block, fields transported to each
  //lattice site from different directions are stored together
  //This allows for faster linear combination in the p4 action.
  int dest2;
  //Index for the gauge link
  int gauge_index;
  //Determines if the gauge link needs to be complex conjugated
  int dagger;
};

struct hop_pointer {
  int src;
  int dest;
};
//CPS_END_NAMESPACE
#endif
