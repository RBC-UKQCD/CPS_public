/******************************************************************************/
/*                                                                            */
/*                          Global Barriers using SPI                         */
/*                                                         30.06.07 S.Krieg   */
/******************************************************************************/

#ifndef _GLB_BARRIER_H_
#define _GLB_BARRIER_H_

/******************************************************************************/
/* Barrier using GLOBINT_CHANNEL_APP_BARRIER hardcoded.                       */
/******************************************************************************/

#define GLOB_BARRIER GlobBarrier()

int GlobBarrier(void);

#endif
