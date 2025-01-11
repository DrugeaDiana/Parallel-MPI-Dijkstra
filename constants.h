#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

#ifndef NUM_NODES
#define NUM_NODES 10000
#endif
#ifndef DENSITY
#define DENSITY 0.6
#endif
#define NUM_EDGES NUM_NODES * (NUM_NODES - 1) * DENSITY
#ifndef SEED
#define SEED 42
#endif
#ifndef NUM_THREADS
#define NUM_THREADS 4
#endif
#define INF_GEN USHRT_MAX
#define INF 999999999
#define MAX_COST 10000

#endif // __CONSTANTS_H__