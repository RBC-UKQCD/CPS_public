#ifndef INCLUDED_MY_UTIL_KL3_H
#define INCLUDED_MY_UTIL_KL3_H

// compute the xyzt index from a single number, useful for OpenMP.
static void compute_coord(int x[4], const int size[4], int i)
{
    for(int j = 0; j < 4; ++j) {
        x[j] = i % size[j];
        i /= size[j];
    }
}

// compute the xyzt index from a single number, useful for OpenMP.
//
// Note 1: this version is for P+A/P-A propagators.
//
// Note 2: this version also computes correct coordinates for P (or A)
// only propagators if i is smaller than the 4D volume.
static void compute_coord_ap(int x[4], const int size[4], int i, int glb_t)
{
    for(int j = 0; j < 4; ++j) {
        x[j] = i % size[j];
        i /= size[j];
    }
    if(i != 0) {
        x[3] += glb_t;
    }
}

// compute the xyzt index from a single number, useful for OpenMP.
static void compute_coord(int x[4],
                          const int delta[4], const int low[4],
                          int i)
{
    for(int j = 0; j < 4; ++j) {
        x[j] = i % delta[j] + low[j];
        i /= delta[j];
    }
}

static int compute_id(const int x[4], const int size[4])
{
    int ret = 0;
    for(int j = 3; j >= 0; --j) {
        ret = ret * size[j] + x[j];
    }
    return ret;
}

// type of propagators to be use in contractions.
// make sure PROP_TYPE and ptype_str[] is consistant.
enum PROP_TYPE {
    PROP_P, // periodic 
    PROP_A, // antiperiodic 
    PROP_PA // P+A
};

static const char *ptype_str[] = {"P", "A", ""};

#endif
