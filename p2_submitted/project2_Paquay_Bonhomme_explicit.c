#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>

#include "project2_Paquay_Bonhomme_main.h"

#define PI 3.14159265358979323846

//Write the .dat file
void toWrite_exp(double* eta, double* u, double* v, unsigned int nb_lines_tot, unsigned int nb_columns, int saveNb){
  char* eta_name = malloc(25*sizeof(char));
  char* u_name = malloc(25*sizeof(char));
  char* v_name = malloc(25*sizeof(char));
  sprintf(eta_name, "results/eta_exp_%d.dat", saveNb);
  sprintf(u_name, "results/u_exp_%d.dat", saveNb);
  sprintf(v_name, "results/v_exp_%d.dat", saveNb);

  FILE* eta_file = fopen(eta_name, "wb");
  FILE* u_file = fopen(u_name, "wb");
  FILE* v_file = fopen(v_name, "wb");

  //eta
  fwrite(&nb_columns, sizeof(unsigned int), 1, eta_file);
  fwrite(&nb_lines_tot, sizeof(unsigned int), 1, eta_file);
  fwrite(eta, sizeof(double), nb_lines_tot*nb_columns, eta_file);
  fclose(eta_file);
  free(eta_name);

  //u
  unsigned int nb_columns_u = nb_columns + 1;
  fwrite(&nb_columns_u, sizeof(unsigned int), 1, u_file);
  fwrite(&nb_lines_tot, sizeof(unsigned int), 1, u_file);
  fwrite(u, sizeof(double), nb_lines_tot*nb_columns_u, u_file);
  fclose(u_file);
  free(u_name);

  //v
  unsigned int nb_lines_v = nb_lines_tot + 1;
  fwrite(&nb_columns, sizeof(unsigned int), 1, v_file);
  fwrite(&nb_lines_v, sizeof(unsigned int), 1, v_file);
  fwrite(v, sizeof(double), nb_lines_v*nb_columns, v_file);
  fclose(v_file);
  free(v_name);
}

//Compute the value of the interpolated H on a given position on the grid
double compute_h_exp(Parameters parameters, double delta_xh, double delta_yh, double *h, double x, double y,
                unsigned int Y,unsigned int X) {
  if (x<0)
    x=0;
  if (y<0)
    y=0;

  int h_k = (int)floor(x/delta_xh);
  int h_l = (int)floor(y/delta_yh);

  if (h_k <0)
    h_k = 0;

  if (h_l <0)
    h_l = 0;
  
  int h_k1 = h_k+1;
  int h_l1 = h_l+1;
    
  if (h_k +1 >=X)
    h_k1 = X-1;

  if (h_l+1>=Y)
    h_l1 = Y-1;

  double x_k = h_k*delta_xh;
  double y_l = h_l*delta_yh;
  double x_k1 = (h_k + 1) * delta_xh;
  double y_l1 = (h_l + 1) * delta_yh;
  
  double h_xy = ((x_k1 - x)*(y_l1 - y) * h[h_l*X + h_k]
                + (x_k1 - x)*(y -y_l) * h[h_l1*X + h_k]
                + (x - x_k)*(y_l1 - y) *h[h_l*X + h_k1] 
                + (x - x_k)*(y - y_l) * h[h_l1*X + h_k1])/(delta_xh*delta_yh);
 
  return h_xy;
}

//Check if there is a left boundary condition for a given position
bool checkCL_left(int x){
  if(x == 0)
    return 1;
  return 0;
}

//Check if there is a right boundary condition for a given position
bool checkCL_right(int x, int nb_columns){
  if(x == nb_columns)
    return 1;
  return 0;
}

//Check if there is a top boundary condition for a given position
bool checkCL_top(int y, int nb_lines){
  if(y == nb_lines - 1)
    return 1;
  return 0;
}

//Check if there is a bottom boundary condition for a given position
bool checkCL_bottom(int y){
  if(y == 0)
    return 1;
  return 0;
}

//Compute u at time t+1
void computeNextU(double* u_t, double* u_t1, double* eta_t, Parameters parameters, int nb_lines, int nb_col){
  #pragma omp parallel default(shared) 
  {
    #pragma omp for schedule(static)
    for(int it = 0; it < nb_lines * (nb_col + 1); it++){ 
      int j = ((int)floor(it/(nb_col + 1))) * (nb_col + 1);
      int i = it%(nb_col + 1);
      int j_eta = ((int)floor(it/(nb_col+1))) * nb_col;

      if(checkCL_left(i) || checkCL_right(i, nb_col))
          u_t1[i + j] = 0;
      else{
          u_t1[i + j] = (-parameters.g*(eta_t[i + j_eta] - eta_t[(i - 1) + j_eta])/(parameters.delta_x)- 
                          parameters.gamma*u_t[i + j])*parameters.delta_t + u_t[i + j];
        }
    }
  }
}

//Compute v at time t+1
void computeNextV(double* v_t, double* v_t1, double* eta_t, Parameters parameters, int nb_lines,
                  int nb_col, int myrank, int np, double* eta_under, double t){ 
  #pragma omp parallel default(shared)
  {
    #pragma omp for schedule(static)
    for(int it = 0; it < nb_lines * nb_col; it++){ 
      int j = ((int)floor(it/nb_col)) * nb_col;
      int i = it%nb_col;

      if(checkCL_top(j/nb_col, nb_lines) && myrank == np-1){
        if(parameters.s == 0){ 
          v_t1[i + j + nb_col] = parameters.A * sin(2*PI*parameters.f*t);
        }
        else
          v_t1[i + j + nb_col] = parameters.A * sin(2*PI*parameters.f*t) * exp(-t/500);

        if(checkCL_bottom(j) && np != 1) //if we only have one line, we are both in the top and in the bottom
          v_t1[i + j] = (-parameters.g*(eta_t[i + j] - eta_under[i])/(parameters.delta_y)- 
                        parameters.gamma*v_t[i + j])*parameters.delta_t + v_t[i + j];
        else if(checkCL_bottom(j)) //if we only have one process
          v_t1[i + j] = 0;
        else
          v_t1[i + j] = (-parameters.g*(eta_t[i + j]  - eta_t[i + (j-nb_col)])/(parameters.delta_y)- 
                        parameters.gamma*v_t[i + j])*parameters.delta_t + v_t[i + j];
      }
      else if(checkCL_bottom(j)){
        if(myrank == 0)
          v_t1[i + j] = 0;
        else
          v_t1[i + j] = (-parameters.g*(eta_t[i + j] - eta_under[i])/(parameters.delta_y)- 
                        parameters.gamma*v_t[i + j])*parameters.delta_t + v_t[i + j];
      }
      else
        v_t1[i + j] = (-parameters.g*(eta_t[i + j] - eta_t[i + (j-nb_col)])/(parameters.delta_y)- 
                      parameters.gamma*v_t[i + j])*parameters.delta_t + v_t[i + j];
    }
  }
}

//Compute eta at time t+1
void computeNextEta(double* eta_t, double* eta_t1, double* u_t, double* v_t, Parameters parameters,
                   int nb_lines, int nb_col, int myrank, int np, double* v_above, double* H_PDi_j,
                   double* H_MDi_j, double* H_i_PDj, double* H_i_MDj){
  #pragma omp parallel default(shared)
  {
    #pragma omp for schedule(static)
    for(int it = 0; it < nb_lines * nb_col; it++){ 
      int j = ((int)floor(it/nb_col)) * nb_col;
      int i = it%nb_col;
      int j_u = ((int)floor(it/nb_col)) * (nb_col+1);

      if(checkCL_top(j/nb_col, nb_lines) && myrank != np-1){
        eta_t1[i + j] = eta_t[i + j] - ((H_PDi_j[i + j]*u_t[i + 1 + j_u] - H_MDi_j[i + j]*u_t[i + j_u])/parameters.delta_x 
                        + (H_i_PDj[i + j]*v_above[i] - H_i_MDj[i + j]*v_t[i + j])/parameters.delta_y)*parameters.delta_t;
      }
      else{
        eta_t1[i + j] = eta_t[i + j] - ((H_PDi_j[i + j]*u_t[i + 1 + j_u] - H_MDi_j[i + j]*u_t[i + j_u])/parameters.delta_x 
                        + (H_i_PDj[i + j]*v_t[i + j + nb_col] - H_i_MDj[i + j]*v_t[i + j])/parameters.delta_y)*parameters.delta_t;
      }
    }
  }
}

//Store in the vectors H_PDi_j, H_MDi_j, H_i_PDj and H_i_MDj the interpolated h
void getH_exp(int nb_lines, int nb_columns, double* H_PDi_j, double* H_MDi_j,
         double* H_i_PDj, double* H_i_MDj, Parameters parameters, double delta_xh, 
         double delta_yh, double* h, unsigned int Y, unsigned int X, int nb_lines_above){
  #pragma omp parallel default(shared)
  {
    #pragma omp for schedule(static)
    for(int it = 0; it < nb_lines * nb_columns; it++){ 
      int j = ((int)floor(it/nb_columns)) * nb_columns;
      int i = it%nb_columns;

      H_PDi_j[i + j] = compute_h_exp(parameters, delta_xh, delta_yh, h, (i+0.5)*parameters.delta_x ,  //PDi= i+0.5
                                (j/nb_columns + nb_lines_above)*parameters.delta_y , Y, X);

      if(i == 0)
        H_MDi_j[i + j] = compute_h_exp(parameters, delta_xh, delta_yh, h, (i-0.5)*parameters.delta_x ,  //PMi= i-0.5
                                  (j/nb_columns + nb_lines_above)*parameters.delta_y , Y, X);
      else
        H_MDi_j[i + j] = H_PDi_j[i + j - 1];

      H_i_PDj[i + j] = compute_h_exp(parameters, delta_xh, delta_yh, h, i*parameters.delta_x ,        //PDj= j+0.5
                                (j/nb_columns + 0.5 + nb_lines_above)*parameters.delta_y , Y, X);

      if(j == 0)
        H_i_MDj[i + j] = compute_h_exp(parameters, delta_xh, delta_yh, h, i*parameters.delta_x ,        //PMj= j-0.5
                                  (j/nb_columns - 0.5 + nb_lines_above)*parameters.delta_y , Y, X);
      else
        H_i_MDj[i + j] = H_i_PDj[i + j - nb_columns];
    }
  }
}

//Return the first line of matrix v
void getFirstLineV(double* first_line_v, double* v, int nb_columns){
  #pragma omp parallel default(shared)
  {
    #pragma omp for schedule(static)
    for(int i = 0; i < nb_columns; i++){
          first_line_v[i] = v[i];
    }
  }
}

//Return the last line of the matrix eta
void getLastLineEta(double* last_line_eta, double* eta, int nb_columns, int nb_lines){
  #pragma omp parallel default(shared)
  {
    #pragma omp for schedule(static)
    for(int i = 0; i < nb_columns; i++){
      last_line_eta[i] = eta[(nb_lines-1)*nb_columns + i];
    }
  }
}

//Gather the vector eta, u and v in a single process in order to write them
void gather_exp(int length_eta, int length_u, int length_v, int myrank, double* eta,
            double* u, double* v, int nb_lines_tot, int nb_columns,
            int* arg_length_eta, int* arg_pos_eta, int* arg_length_u, int* arg_pos_u,
            int* arg_length_v, int* arg_pos_v, int saveNb){
  double* eta_final;
  double* u_final;
  double* v_final;

  if(myrank == 0){
    eta_final = malloc(nb_lines_tot * nb_columns * sizeof(double));
    u_final = malloc(nb_lines_tot * (nb_columns + 1) * sizeof(double));
    v_final = malloc((nb_lines_tot + 1) * nb_columns * sizeof(double));
  }

  MPI_Gatherv(eta, length_eta, MPI_DOUBLE, eta_final, arg_length_eta, arg_pos_eta, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gatherv(u, length_u, MPI_DOUBLE, u_final, arg_length_u, arg_pos_u, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gatherv(v, length_v, MPI_DOUBLE, v_final, arg_length_v, arg_pos_v, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if(myrank == 0){
    toWrite_exp(eta_final, u_final, v_final, (unsigned int) nb_lines_tot, (unsigned int) nb_columns, saveNb);
    free(eta_final);
    free(u_final);
    free(v_final);
  }
}

//Do the explicit resolution
void explicitEuler(int myrank, int np, Parameters parameters, Map map, double* h, double delta_xh,
                  double delta_yh, MPI_Status status){
  //Get the number of lines for the current process
  int nb_columns = (int) map.a/parameters.delta_x + 1;
  int nb_lines_tot = (int) map.b/parameters.delta_y + 1;

  int nb_lines = floor(nb_lines_tot/np);
  int mod = nb_lines_tot - np*nb_lines;

  if(myrank < mod)
    nb_lines += 1;

  //get the length
  int length_eta = nb_lines * nb_columns;
  int length_u = nb_lines * (nb_columns + 1);
  int length_v;

  //get the vectors at time 0
  double* eta_t = calloc(nb_lines * nb_columns, sizeof(double));
  double* u_t = calloc(nb_lines * (nb_columns + 1), sizeof(double));
  double* v_t;
  if(myrank == np-1){
    v_t = calloc((nb_lines + 1) * nb_columns, sizeof(double)); //add the above line
    length_v = (nb_lines + 1) * nb_columns;
  }
  else{
    v_t = calloc(nb_lines * nb_columns, sizeof(double));
    length_v = nb_lines* nb_columns;
  }

  //Get the number of lines above
  int nb_lines_above;
  if(myrank < mod)
    nb_lines_above = myrank * nb_lines;
  else
    nb_lines_above = mod * (nb_lines + 1) + (myrank - mod) * nb_lines;
  
  //Get the h interpolated
  double* H_PDi_j = malloc(nb_lines * nb_columns * sizeof(double));
  double* H_MDi_j = malloc(nb_lines * nb_columns * sizeof(double));
  double* H_i_PDj = malloc(nb_lines * nb_columns * sizeof(double));
  double* H_i_MDj = malloc(nb_lines * nb_columns * sizeof(double));

  getH_exp(nb_lines, nb_columns, H_PDi_j, H_MDi_j, H_i_PDj, H_i_MDj, parameters, delta_xh, 
       delta_yh, h, map.Y, map.X, nb_lines_above);

  //Get the vectors in order to gather
  int* arg_length_eta; 
  int* arg_length_u;
  int* arg_length_v;
  int* arg_pos_eta;
  int* arg_pos_u;
  int* arg_pos_v;

  if(np != 1 && parameters.S != 0){
    if(myrank == 0){
      arg_length_eta = malloc(np * sizeof(int));
      arg_length_u = malloc(np * sizeof(int));
      arg_length_v = malloc(np * sizeof(int));
    }

    MPI_Gather(&length_eta, 1, MPI_INT, arg_length_eta, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&length_u, 1, MPI_INT, arg_length_u, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&length_v, 1, MPI_INT, arg_length_v, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if(myrank == 0){
      arg_pos_eta = malloc(np * sizeof(int));
      arg_pos_u = malloc(np * sizeof(int));
      arg_pos_v = malloc(np * sizeof(int));

      arg_pos_eta[0] = 0;
      arg_pos_u[0] = 0;
      arg_pos_v[0] = 0;

      for(int i = 1; i < np; i++){
        arg_pos_eta[i] = arg_pos_eta[i - 1] + arg_length_eta[i - 1];
        arg_pos_u[i] = arg_pos_u[i - 1] + arg_length_u[i - 1];
        arg_pos_v[i] = arg_pos_v[i - 1] + arg_length_v[i - 1];
      }
    }
  }

  //Prepare the vectors needed for the time loop
  double* eta_t1 = malloc(nb_lines * nb_columns * sizeof(double));
  double* u_t1 = malloc(nb_lines * (nb_columns + 1) * sizeof(double));
  double* v_t1;
  if(myrank == np - 1)
    v_t1 = malloc((nb_lines + 1) * nb_columns * sizeof(double));
  else
    v_t1 = malloc(nb_lines * nb_columns * sizeof(double));

  double* first_line_v = malloc(nb_columns * sizeof(double));
  double* eta_under = malloc(nb_columns * sizeof(double));
  double* v_above = malloc(nb_columns * sizeof(double));
  double* last_line_eta = malloc(nb_columns * sizeof(double));

  //Start the time loop
  for(int t = 1; t <= parameters.T_max/parameters.delta_t; t++){
    //Send and receive all the lines needed to compute the next submatrix
    if(myrank == np-1){
      if(np > 1){
        //To send
        getFirstLineV(first_line_v, v_t, nb_columns);
        MPI_Send(first_line_v, nb_columns, MPI_DOUBLE, myrank-1, 23, MPI_COMM_WORLD);

        //To receive
        MPI_Recv(eta_under, nb_columns, MPI_DOUBLE, myrank-1, 23, MPI_COMM_WORLD, &status);
      }
    }
    else{
      if(myrank == 0){
        //To receive
        MPI_Recv(v_above, nb_columns, MPI_DOUBLE, 1, 23, MPI_COMM_WORLD, &status);

        //To send
        getLastLineEta(last_line_eta, eta_t, nb_columns, nb_lines);
        MPI_Send(last_line_eta, nb_columns, MPI_DOUBLE, 1, 23, MPI_COMM_WORLD);
      }
      else{
        //To receive UNDER
        MPI_Recv(v_above, nb_columns, MPI_DOUBLE, myrank+1, 23, MPI_COMM_WORLD, &status);

        //To send UNDER
        getLastLineEta(last_line_eta, eta_t, nb_columns, nb_lines);
        MPI_Send(last_line_eta, nb_columns, MPI_DOUBLE, myrank+1, 23, MPI_COMM_WORLD);

        //To send ABOVE
        getFirstLineV(first_line_v, v_t, nb_columns);
        MPI_Send(first_line_v, nb_columns, MPI_DOUBLE, myrank-1, 23, MPI_COMM_WORLD);

        //To receive ABOVE
        MPI_Recv(eta_under, nb_columns, MPI_DOUBLE, myrank-1, 23, MPI_COMM_WORLD, &status);
      }
    }

    //Compute the next submatrix
    computeNextU(u_t, u_t1, eta_t, parameters, nb_lines, nb_columns);
    computeNextV(v_t, v_t1, eta_t, parameters, nb_lines, nb_columns, myrank, np, eta_under, t*parameters.delta_t);
    computeNextEta(eta_t, eta_t1, u_t, v_t, parameters, nb_lines, nb_columns, myrank, np, v_above,
                  H_PDi_j, H_MDi_j, H_i_PDj, H_i_MDj);
  
    //Copy for the next loop
    memcpy(u_t, u_t1, nb_lines * (nb_columns + 1) * sizeof(double));
    memcpy(eta_t, eta_t1, nb_lines * nb_columns * sizeof(double));

    if(myrank == np-1)
      memcpy(v_t, v_t1, (nb_lines + 1) * nb_columns * sizeof(double));
    else
      memcpy(v_t, v_t1, nb_lines * nb_columns * sizeof(double));

    //Save the data
    if(parameters.S!=0 && t % parameters.S == 0 && np != 1)
      gather_exp(length_eta, length_u, length_v, myrank, eta_t, u_t, v_t, nb_lines_tot, nb_columns,
            arg_length_eta, arg_pos_eta, arg_length_u, arg_pos_u, arg_length_v, arg_pos_v, t);
  }

  free(eta_t);
  free(u_t);
  free(v_t);

  free(first_line_v);
  free(eta_under);
  free(v_above);
  free(last_line_eta);

  if(myrank == 0 && np != 1 && parameters.S!=0){
    free(arg_length_eta);
    free(arg_pos_eta);
    free(arg_pos_u);
    free(arg_length_u);
    free(arg_pos_v);
    free(arg_length_v);
  }
}