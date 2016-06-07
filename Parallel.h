/**
* @file  Parallel.h
* @brief Gives the code parallelization support
*
* @license This file is distributed under the BSD Open Source License. 
*          See LICENSE.TXT for details. 
**/

#ifndef PARALLEL_H
#define PARALLEL_H

/*===- OpenMP -------------------------------------------------------------===*/
// Implements OpenMP parallelization if the compiler contains support.
#ifdef _OPENMP
#include <omp.h>

#define STRING(s) #s

typedef int cloud_index;

// Parallelize for loops. kind determines the scheduling of iterations to the 
// threads.
#define BEGIN_PARALLEL_FOR(i,e,num,step,kind)\
_Pragma(STRING(omp parallel for schedule(kind))) \
for (cloud_index i = 0; i < num; i += step) {

#define END_PARALLEL_FOR }

// Thread synronization routines.
#define SEMAPHORES omp_lock_t *locks;

#define SEMAPHORES_MALLOC(num) , locks(new omp_lock_t[num])

#define SEMAPHORES_INIT(num) \
_Pragma("omp parallel for schedule(static)") \
for (cloud_index i = 0; i < (num); i++) \
omp_init_lock(locks + i);

#define SEMAPHORES_FREE(num) \
_Pragma("omp parallel for schedule(static)") \
for (cloud_index i = 0; i < (num); i++) \
omp_destroy_lock(locks + i);

#define SEMAPHORE_WAIT(i) omp_set_lock(locks + i);

#define SEMAPHORE_SIGNAL(i) omp_unset_lock(locks + i);

/*===- libDispatch --------------------------------------------------------===*/
// Implements libDispatch parallelization. This is used as the default on all 
// Apple targets. This can be used on BSD, Linux and Windows targets if a port 
// of libDispatch is avalible for those patforms with block support in the 
// compiler. libDispatch is avalible here http://libdispatch.macosforge.org/
#elif defined (__APPLE__)
#include <dispatch/dispatch.h>

#define DISPATCH_QUEUES

typedef size_t cloud_index;

// Parallelize for loops. kind is unused.
#define BEGIN_PARALLEL_FOR(i,e,num,step,kind) \
dispatch_apply((num)/step, dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0), ^(cloud_index i) { \
    i *= step;

#define END_PARALLEL_FOR });

// Thread synronization routines.
#define SEMAPHORES dispatch_semaphore_t *semaphores;

#define SEMAPHORES_MALLOC(num) , semaphores(new dispatch_semaphore_t[num])

#define SEMAPHORES_INIT(num) \
dispatch_apply(num, dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0), ^(cloud_index i) { \
    semaphores[i] = dispatch_semaphore_create(1); \
});

#define SEMAPHORES_FREE(num) \
dispatch_apply(num, dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0), ^(cloud_index i) { \
    dispatch_release(semaphores[i]); \
});

#define SEMAPHORE_WAIT(i) dispatch_semaphore_wait(semaphores[i], DISPATCH_TIME_FOREVER);

#define SEMAPHORE_SIGNAL(i) dispatch_semaphore_signal(semaphores[i]);

/*===- scalar -------------------------------------------------------------===*/
// If no parallelization is availible fallback to single threaded code.
#else

typedef unsigned int cloud_index;

// Loop over in a serial fasion. kind is unused.
#define BEGIN_PARALLEL_FOR(i,e,num,step,kind) for (cloud_index i = 0, e = num; i < e; i += step) {

#define END_PARALLEL_FOR }

// Thread synronization routines. Since there is only one thread these expand to
// nothing.
#define SEMAPHORES

#define SEMAPHORES_MALLOC(num)

#define SEMAPHORES_INIT(num)

#define SEMAPHORES_FREE(num)

#define SEMAPHORE_WAIT(i)

#define SEMAPHORE_SIGNAL(i)

#endif

#endif // PARALLEL_H
