#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <mpi.h>
#include <time.h>
#include <stdbool.h>
//#include "euler_explicit.h" //p-e faire un fichier header??

typedef struct{
  double g;
  double gamma;
  double delta_x;
  double delta_y;
  double delta_t;
  double T_max;
  double A;
  double f;
  unsigned int S;
  unsigned int s;
  double r_tresh;
}Parameters;

typedef struct{
  double a;
  double b;
  unsigned int X;
  unsigned int Y;
}Map;

int main(int argc, char *argv[]){
  if(argc != 4)
    exit(-1);

  //READ the arguments and parse them => A FAIRE
  // on hardcode pour le moment
  string map_str = "map.dat";
  string param_str = "param.txt";

  //read the .dat file
  FILE* map_file = fopen(&map_str, "rb");
  Map map;

  //read a,b,X and Y
  fread(&map, sizeof(Map), 1, map_file);
  double delta_xh = map.a/(map.X-1);
  double delta_yh = map.b/(map.Y-1);

  //read h
  double* h = malloc(X*Y*sizeof(double));
  fread(h, sizeof(double), X*Y, map_file);
  fclose(map_file);

  //read the .txt file
  FILE* param = fopen(&param_str, "r");
  Parameters parameters;

  fread(&parameters, sizeof(Parameters), 1, param);
  fclose(param);

  int myrank, np;

  MPI_Init(&argc , &argv);
  MPI_Comm_size(MPI_COMM_WORLD , &np); //le nombre total de process MPI
  MPI_Comm_rank(MPI_COMM_WORLD , &myrank); //dans quel process on est

  if(flag == 0){
    explicitEuler(myrank, np, parameters, map, h, delta_xh, delta_yh);
  }
  else{
    //implicitEuler();
  }



  	
  return 0 ; 
}