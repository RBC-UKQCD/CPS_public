/* Abstracted threading model
 * int bagel_thread_create(func, arg)
 * void bagel_thread_join (int)
 */
#include "wfm.h"
#include <stdio.h>
#include <stdlib.h>

void wfm::CoreCount(int num)
{
#ifdef USE_THREADS_NONE
 if ( num != 1 ) {
  if ( 1 ) printf("Bagel: CoreCount wrong disabled thread configuration\n");
   exit(-1);
 }
#endif
#ifdef USE_THREADS_BGL
 if ( ( num < 1 ) || num > 2 ) { 
  if ( 1 ) printf("Bagel: CoreCount wrong for BG/L threading\n");
   exit(-1);
 }
#endif
 if ( num < 1 ) { 
  if ( 1 ) printf("Bagel: Bad CoreCount\n");
   exit(-1);
 }
 nthread = num;
}

#if defined( USE_THREADS_FAKE) || defined (USE_THREADS_NONE)
void wfm::thread_create(void *(*func)(void *), void *arg,int id)
{
if ( 1 ) printf("Fake threads\n"); fflush(stdout);
 func(arg);
}
void wfm::thread_join  (int id) 
{
}
#endif

#ifdef USE_THREADS_POSIX

#include <vector>
#include <pthread.h>

std::vector<pthread_t> threads;

void wfm::thread_create(void *(*func)(void *), void *arg,int id)
{
  if ( id+1 > threads.size()) 
  threads.resize(id+1);

  if ( id < 0  ) { 
   if ( 1 ) printf("Bagel threading : core number out of range\n");
    exit(-1);
  }
  pthread_create(&threads[id],NULL,func,arg);
}
void wfm::thread_join  (int id) 
{
  pthread_join(threads[id],NULL);
}
#endif
#ifdef USE_THREADS_BGL

#include <vector>
#include <rts.h>

std::vector<BGL_co_handle_t> threads;

void wfm::thread_create(void *(*func)(void *), void *arg,int id)
{
if ( 1 ) printf("BGL threads\n"); fflush(stdout);
  if ( id+1 > threads.size()) 
  threads.resize(id+1);

  if ( id < 0  ) { 
   if ( 1 ) printf("Bagel threading : core number out of range\n");
    exit(-1);
  }

  rts_dcache_evict();
  threads[id] = co_start(func,arg);
}
void wfm::thread_join  (int id) 
{
  co_join(threads[id]);
}

#endif
