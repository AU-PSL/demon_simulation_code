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

typedef unsigned int cloud_index;
#define begin_parallel_for(i,e,num,step) for (cloud_index i = 0, e = num; i < e; i += step) {
#define end_parallel_for }

#endif // PARALLEL_H
