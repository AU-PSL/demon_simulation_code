/*===- Parallel.h - libSimulation -=============================================
*
*                                  DEMON
*
* This file is distributed under the BSD Open Source License. See LICENSE.TXT 
* for details.
*
*===-----------------------------------------------------------------------===*/

#ifndef PARALLEL_H
#define PARALLEL_H

#define BLOCK_VALUE_TIME currentTimeStep
#define BLOCK_VALUE_DIST currentDist

#ifdef __APPLE__

#include <dispatch/dispatch.h>
#define DISPATCH_QUEUES

typedef size_t cloud_index;
#define BEGIN_PARALLEL_FOR(i,e,num,step) \
dispatch_apply(num/step, dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0), ^(cloud_index i) { \
    i *= step;
#define END_PARALLEL_FOR });

#define SEMAPHORES dispatch_semaphore_t *semaphores;
#define SEMAPHORES_MALLOC(num) semaphores(new dispatch_semaphore_t[num])
#define SEMAPHORES_INIT(num) \
dispatch_apply(num, dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0), ^(cloud_index i) { \
    semaphores[i] = dispatch_semaphore_create(1); \
});
#define SEMAPHORES_FREE(num) \
dispatch_apply(cloud->n/2, dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0), ^(cloud_index i) { \
    dispatch_release(semaphores[i]); \
});
#define SEMAPHORE_WAIT(i) dispatch_semaphore_wait(semaphores[i], DISPATCH_TIME_FOREVER);
#define SEMAPHORE_SIGNAL(i) dispatch_semaphore_signal(semaphores[i]);

#undef BLOCK_VALUE_TIME
#define BLOCK_VALUE_TIME currTimeStep
#undef BLOCK_VALUE_DIST
#define BLOCK_VALUE_DIST currDist

#else

typedef unsigned int cloud_index;
#define BEGIN_PARALLEL_FOR(i,e,num,step) for (cloud_index i = 0, e = num; i < e; i += step) {
#define END_PARALLEL_FOR }

#endif

#endif // PARALLEL_H
