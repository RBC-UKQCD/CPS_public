//
//   Flconv.h
//               convert   IEEE32BIG <->  TIDSP


#ifndef __FCONVERT__
#define __FCONVERT__

//typedef unsigned long int type32 ;
typedef uint32_t type32 ;

void ti2ieee (type32 *,int);
//void ieee2ti (type32 *,int);


void byterevn(type32 w[], int n);

#endif 
