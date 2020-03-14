#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "project2_Paquay_Bonhomme_main.h"
#include "project2_Paquay_Bonhomme_explicit.h"
#include "project2_Paquay_Bonhomme_implicit.h"

int main(int argc, char *argv[]){
  
  if(argc != 4){
    printf("3 arguments needed\n");
    exit(-1);
  }

  int flag;
  sscanf(argv[3], "%d", &flag);

  //read the .dat file
  FILE* map_file = fopen(argv[2], "rb");
  Map map;

  //read a,b,X and Y
  fread(&map, sizeof(Map), 1, map_file);
  double delta_xh = map.a/(map.X-1);
  double delta_yh = map.b/(map.Y-1);

  //read h
  double* h = malloc(map.X*map.Y*sizeof(double));
  fread(h, sizeof(double), map.X*map.Y, map_file);
  fclose(map_file);

  //read the .txt file
  int nb_args = 11;
  FILE* param = fopen(argv[1], "r");
  Parameters parameters;

  int nb = fscanf(param, "%lf\n" "%lf\n" "%lf\n"
                         "%lf\n" "%lf\n" "%lf\n"
                         "%lf\n" "%lf\n" "%u\n"
                         "%u\n" "%lf\n",
                         &parameters.g, &parameters.gamma, &parameters.delta_x,
                         &parameters.delta_y, &parameters.delta_t, &parameters.T_max,
                         &parameters.A, &parameters.f, &parameters.S,
                         &parameters.s, &parameters.r_tresh);
  fclose(param);

  if(nb != nb_args){
    printf("Wrong number or type of parameters in the parameters file.\n");
    exit(-1);
  }

  int myrank, np;
  MPI_Status status;

  MPI_Init(&argc , &argv);
  MPI_Comm_size(MPI_COMM_WORLD , &np); //total number of MPI process
  MPI_Comm_rank(MPI_COMM_WORLD , &myrank); //which process

  if(flag == 0){
    explicitEuler(myrank, np, parameters, map, h, delta_xh, delta_yh, status);
  }
  if(flag == 1){
    implicitEuler(myrank, np, parameters, map, h, delta_xh, delta_yh, status);
  }

  MPI_Finalize();

  return 0 ; 
}