#ifndef _QLA_RANDOM_H
#define _QLA_RANDOM_H

#include <qla_types.h>

/* random number structures */

typedef struct {
  /* We assume long is at least 32 bits */
  unsigned long r0,r1,r2,r3,r4,r5,r6;
  unsigned long multiplier,addend,ic_state;
  float scale;
} QLA_RandomState;

/* Generic random number generator returning a uniformly distributed
   random value on [0,1] */

QLA_F_Real QLA_random(QLA_RandomState *prn_pt);
QLA_F_Real QLA_gaussian(QLA_RandomState *prn_pt);
void QLA_seed_random(QLA_RandomState *prn_pt, int seed, QLA_Int index);

#endif /* _QLA_RANDOM_H */
