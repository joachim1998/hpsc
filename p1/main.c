#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include "omp.h"
#include <stdint.h>

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

/*
*INPUT: int x, int y, int N
*OUTPUT: int x*N+y
*
*This function will transfrom the x and y coordinates from a matrix
*into the position in a single row array and will return the result
*/
int getPosition(int x, int y, int N){
	return x * N + y; 
}

/*
*INPUT: int nb
*OUTPUT: boolean
*
*This function will tell us if a cell is conducting or not
*/
bool isConductingCell(int nb){
	if(nb == 1)
		return true;
	return false;
}

/*
*INPUT: int* array, int N
*OUTPUT: boolean
*
*This function will check if the grid is conductible
*/
bool isConductingGrid(int* array, int N){
	bool val = false;

	for(int i=0; i<N; i++){
		if(array[getPosition(i, N-1, N)] == 2)
			val = true;
	}
	return val;
}

/*
*INPUT: int* array, int x, int y, int N
*OUTPUT: /
*
*see statement in the recursive function part
*/
void recurs_func(int* array, int x, int y, int N){
	if(isConductingCell(array[getPosition(x, y, N)])){
		array[getPosition(x, y, N)] = 2;

		int left_neighbor = y - 1;
		int right_neighbor = y + 1;
		int up_neighbor = x - 1;
		int down_neighbor = x + 1;

		if(left_neighbor >= 0)
			recurs_func(array, x, left_neighbor, N);
		if(right_neighbor < N)
			recurs_func(array, x, right_neighbor, N);
		if(up_neighbor >= 0)
			recurs_func(array, up_neighbor, y, N);
		if(down_neighbor < N)
			recurs_func(array, down_neighbor, y, N);
	}
}

/*
*INPUT: int N, int* array
*OUTPUT: /
*
*This function will get a ppm image from the grid
*/
void getPPM(int N, int* array){
	FILE* fp = fopen("grid.ppm", "wb");
	fprintf(fp, "P6\n%d %d\n255\n", N, N);

	for(int i=0; i<N; i++){
		for(int j=0; j<N; j++){
			static unsigned char color[3];

			if(array[getPosition(i, j, N)] == 2){ //connected cell
				color[0] = 255;
				color[1] = 0;
				color[2] = 0;
			}
			else if(array[getPosition(i, j, N)] == 1){ //conducting cell
				color[0] = 0;
				color[1] = 0;
				color[2] = 0;
			}
			else{
				color[0] = 255;
				color[1] = 255;
				color[2] = 255;
			}

			fwrite(color, 1, 3, fp);
		}
	}

	fclose(fp);
}

/*
*INPUT: int deadline, Seed* seed, int N, int nb_fibers
*OUTPUT: int nb_conduc
*
*see statement in the percolation algorithm part
*It will return 1 if the grid is conducting or 0 if not.
*/
int core(int deadline, Seed *seed, int N, int nb_fibers){
	int x, y, dir;

	//generate the array of N² integers fill with 0.
	int* array = calloc (N*N, sizeof(int));
	if(!array)
		exit(-1);

	for(int i=0; i<nb_fibers; i++){

		x = my_rand(seed) % N;
		y = my_rand(seed) % N;
		dir = my_rand(seed) % 2;

		array[getPosition(x, y, N)] = 1;

		if(dir == 0){
		    int left_neighbor = y - 1;
		    int right_neighbor = y + 1;

		    if(left_neighbor >= 0)
		    	array[getPosition(x, left_neighbor, N)] = 1;

		    if(right_neighbor < N)
		    	array[getPosition(x, right_neighbor, N)] = 1;
		}
		else{
		    int up_neighbor = x - 1;
		    int down_neighbor = x + 1;

		    if(up_neighbor >= 0)
		    	array[getPosition(up_neighbor, y, N)] = 1;

		    if(down_neighbor < N)
		    	array[getPosition(down_neighbor, y, N)] = 1;
		}
	}

	//Loop over all the leftmost cells
	for(int i=0; i<N; i++){
		recurs_func(array, i, 0, N);
	}

	//Check if connected
	int nb_conduc = 0;
	if(isConductingGrid(array, N))
		nb_conduc = 1;

	if(deadline == 0)
		getPPM(N, array);

	free(array);

	return nb_conduc;
}

int main(int argc, char** argv){
    if (argc < 4)
		exit(-1);

	int flag, nb_samples, N, nb_conduc = 0;
	double d;

    sscanf(argv[1], "%d", &flag);
    sscanf(argv[2], "%d", &N);
    sscanf(argv[3], "%lf", &d);

    if(flag == 1)
    	sscanf(argv[4], "%d", &nb_samples);

    //computing the number of fibers
    int nb_fibers = N * N * d;

    if(flag == 1){
		#pragma omp parallel default(none) shared(nb_fibers, N, nb_samples, nb_conduc, flag)
		{
			Seed seed;
			init_seed(&seed);

			#pragma omp for schedule(static)
			for(int i=0; i<nb_samples; i++){
				int conduc = core(flag, &seed, N, nb_fibers);

				#pragma omp atomic
				nb_conduc += conduc;
			}
		}
		double prob_cond = ((double) nb_conduc / (double) nb_samples) * 100;
		printf("The probability of conduction is : %.2lf%%\n", prob_cond);
	}
	else{
		Seed seed;
		init_seed(&seed);
		if(core(flag, &seed, N, nb_fibers) == 1)
			printf("The grid is conducting\n");
		else printf("The grid is NOT conducting\n");
	}
	
    return 0;
}