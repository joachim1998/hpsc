#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "omp.h"

static int32_t s_randtbl[32] =
  {
    3,

    -1726662223, 379960547, 1735697613, 1040273694, 1313901226,
    1627687941, -179304937, -2073333483, 1780058412, -1989503057,
    -615974602, 344556628, 939512070, -1249116260, 1507946756,
    -812545463, 154635395, 1388815473, -1926676823, 525320961,
    -1009028674, 968117788, -123449607, 1284210865, 435012392,
    -2017506339, -911064859, -370259173, 1132637927, 1398500161,
    -205601318,
  };

typedef struct Seed{
  int32_t randtbl[32];
  int32_t *fptr;
  int32_t *rptr;
  int32_t *end_ptr;
  int32_t *state;
} Seed; 

// The random function you could used in your project instead of the C 'rand' function. 
unsigned int my_rand(Seed *seed) 
{
  int32_t *fptr = seed->fptr;
  int32_t *rptr = seed->rptr;
  int32_t *end_ptr = seed->end_ptr;
  unsigned int val;

  val = *fptr += (unsigned int) *rptr;
  ++fptr;
  if (fptr >= end_ptr) {
    fptr = seed->state;
    ++rptr;
  }
  else {
    ++rptr;
    if (rptr >= end_ptr)
      rptr = seed->state;
  }
  seed->fptr = fptr;
  seed->rptr = rptr;

  return val;
}

void init_seed(Seed *seed)
{
  for(int i = 0; i < 32; ++i) {
    seed->randtbl[i] = s_randtbl[i];
  }
  seed->fptr = &(seed->randtbl[2]);
  seed->rptr = &(seed->randtbl[1]);
  seed->end_ptr = &(seed->randtbl[sizeof (seed->randtbl) / sizeof (seed->randtbl[0])]);
  seed->state = &(seed->randtbl[1]);

  unsigned int init = (time(NULL) << omp_get_thread_num());
  seed->state[0] = init;
  int32_t *dst = seed->state;
  int32_t word = init;
  int kc = 32;
  for(int i = 1; i < kc; ++i) {
    long int hi = word / 127773;
    long int lo = word % 127773;
    word = 16807 * lo - 2836 * hi;
    if (word < 0)
      word += 2147483647;
    *++dst = word;
  }
  seed->fptr = &(seed->state[3]);
  seed->rptr = &(seed->state[0]);
  kc *= 10;
  while (--kc >= 0) {
    my_rand(seed);
  }
}

int main(int argc, char **argv)
{
#pragma omp parallel
  {
    // This variable MUST be local (each thread must have its own seed variable). You can use init_seed function to initialize it.
    Seed seed;
    init_seed(&seed);
      
#pragma omp for
    for(int i = 0; i < 10; ++i) {
      // The 'pragma omp critical' could be removed. It is used only to avoid making a mess in printf. In your project, remove it!
#pragma omp critical
      printf("Thread %u, rand = %u\n", omp_get_thread_num(), my_rand(&seed));
    }
  }
}
