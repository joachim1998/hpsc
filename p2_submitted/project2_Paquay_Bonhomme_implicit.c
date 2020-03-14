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
void toWrite_imp(double* eta, double* u, double* v, unsigned int nb_lines_tot, unsigned int nb_columns, int saveNb){
  char* eta_name = malloc(25*sizeof(char));
  char* u_name = malloc(25*sizeof(char));
  char* v_name = malloc(25*sizeof(char));
  sprintf(eta_name, "results/eta_imp_%d.dat", saveNb);
  sprintf(u_name, "results/u_imp_%d.dat", saveNb);
  sprintf(v_name, "results/v_imp_%d.dat", saveNb);

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
  unsigned int nb_columns_u = nb_columns + 1 - 2;
  fwrite(&nb_columns_u, sizeof(unsigned int), 1, u_file);
  fwrite(&nb_lines_tot, sizeof(unsigned int), 1, u_file);
  fwrite(u, sizeof(double), nb_lines_tot*nb_columns_u, u_file);
  fclose(u_file);
  free(u_name);

  //v
  unsigned int nb_lines_v = nb_lines_tot + 1 - 2;
  fwrite(&nb_columns, sizeof(unsigned int), 1, v_file);
  fwrite(&nb_lines_v, sizeof(unsigned int), 1, v_file);
  fwrite(v, sizeof(double), nb_lines_v*nb_columns, v_file);
  fclose(v_file);
  free(v_name);
}

//Compute the value of the interpolated H on a given position on the grid
double compute_h_imp(Parameters parameters, double delta_xh, double delta_yh, double *h, double x, double y,
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

//Store in the vectors H_i and H_j the interpolated h on a given direction (i or j)
void getH_imp(int nb_lines_i, int nb_lines_j, int nb_columns_i, int nb_columns_j, double* H_i, double* H_j, Parameters parameters, double delta_xh, 
         double delta_yh, double* h, unsigned int Y, unsigned int X, int nb_lines_above){ //OK
  #pragma omp parallel default(shared)
  {
    #pragma omp for schedule(static)
    for(int it = 0; it < nb_lines_i * nb_columns_i; it++){ //for H_i
      int j = ((int)floor(it/nb_columns_i)) * nb_columns_i;
      int i = it%nb_columns_i;

      H_i[i + j] = compute_h_imp(parameters, delta_xh, delta_yh, h, (i-0.5)*parameters.delta_x,
                                (j/nb_columns_i + nb_lines_above)*parameters.delta_y , Y, X);
    }

    #pragma omp for schedule(static)
    for(int it = 0; it < nb_lines_j*nb_columns_j; it++){ //for H_j
      int j = ((int)floor(it/nb_columns_j)) * nb_columns_j;
      int i = it%nb_columns_j;

      H_j[i + j] = compute_h_imp(parameters, delta_xh, delta_yh, h, i*parameters.delta_x ,
                            (j/nb_columns_j - 0.5 + nb_lines_above)*parameters.delta_y , Y, X);
    }
  }
}

//Store in the vectors H_i_tot and H_j_tot all interpolated h map on a given direction (i and j) for all the processes
void interpolateH(int myrank, int np,int nb_columns_tot, int nb_lines_tot, Parameters parameters, double delta_xh, double delta_yh,
                 Map map, double* h, double* H_i_tot, double* H_j_tot){
  //Get the h interpolated for all the process
  //Nb of lines of h computed by the process
  int nb_lines = floor(nb_lines_tot/np);
  int mod = nb_lines_tot - np*nb_lines;

  if(myrank < mod)
    nb_lines += 1;

  int nb_lines_i = nb_lines;
  int nb_lines_j = nb_lines;
  int nb_columns_tot_i = nb_columns_tot + 1;
  int nb_columns_tot_j = nb_columns_tot;

  if(myrank == np-1)
    nb_lines_j += 1;

  //get the length in order to gather the h
  int length_i = nb_lines_i * nb_columns_tot_i;
  int length_j = nb_lines_j * nb_columns_tot_j;

  //Get the number of lines above
  int nb_lines_above;
  if(myrank < mod)
    nb_lines_above = myrank * nb_lines;
  else
    nb_lines_above = mod * (nb_lines + 1) + (myrank - mod) * nb_lines;

  //Get the h interpolated
  double* H_i = malloc(nb_lines_i * nb_columns_tot_i * sizeof(double));
  double* H_j = malloc(nb_lines_j * nb_columns_tot_j * sizeof(double));
  
  getH_imp(nb_lines_i, nb_lines_j, nb_columns_tot_i, nb_columns_tot_j, H_i, H_j, parameters, delta_xh, 
       delta_yh, h, map.Y, map.X, nb_lines_above);

  int* arg_length_i = malloc(np * sizeof(int));
  int* arg_length_j = malloc(np * sizeof(int));
  MPI_Allgather(&length_i, 1, MPI_INT, arg_length_i, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Allgather(&length_j, 1, MPI_INT, arg_length_j, 1, MPI_INT, MPI_COMM_WORLD);

  int* arg_pos_i = malloc(np * sizeof(int));
  int* arg_pos_j = malloc(np * sizeof(int));
  arg_pos_i[0] = 0;
  arg_pos_j[0] = 0;

  for(int i = 1; i < np; i++){
    arg_pos_i[i] = arg_pos_i[i - 1] + arg_length_i[i - 1];
    arg_pos_j[i] = arg_pos_j[i - 1] + arg_length_j[i - 1];
  }

  MPI_Allgatherv(H_i, length_i, MPI_DOUBLE, H_i_tot, arg_length_i, arg_pos_i, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Allgatherv(H_j, length_j, MPI_DOUBLE, H_j_tot, arg_length_j, arg_pos_j, MPI_DOUBLE, MPI_COMM_WORLD);

  free(arg_length_i);
  free(arg_length_j);
  free(arg_pos_i);
  free(arg_pos_j);
  free(H_i);
  free(H_j);
}

//Gather all the sub-vectors on every process in all the process
void getWholeVector(int myrank, int np, double* sub_vec, int size_sub, double* all_vec){
  if(np == 1)
    memcpy(all_vec, sub_vec, size_sub * sizeof(double));

  else{
    int* arg_length = malloc(np * sizeof(int));
    MPI_Allgather(&size_sub, 1, MPI_INT, arg_length, 1, MPI_INT, MPI_COMM_WORLD);

    int* arg_pos = malloc(np * sizeof(int));
    arg_pos[0] = 0;

    for(int i = 1; i < np; i++){
      arg_pos[i] = arg_pos[i - 1] + arg_length[i - 1];
    }

    MPI_Allgatherv(sub_vec, size_sub, MPI_DOUBLE, all_vec, arg_length, arg_pos, MPI_DOUBLE, MPI_COMM_WORLD);

    free(arg_length);
    free(arg_pos);
  }
}

//Add the vec1 and vec2 with a coef to vec2 and store the result in vec_sol
void addTwoVector(double* vec1, double* vec2, double* vec_sol, int size, double coef_vec2){
  #pragma omp parallel default(shared)
  {
    #pragma omp for schedule(static)
    for(int i = 0; i < size; i++){
        vec_sol[i] = vec1[i] + coef_vec2 * vec2[i];
    }
  }
}

//Add the vec1 and vec2 with a coef to vec2 and store the result in vec_sol and vec_sol1
void addTwoVectorSpe(double* vec1, double* vec2, double* vec_sol, double* vec_sol1, int size, double coef_vec2){
  #pragma omp parallel default(shared)
  {
    #pragma omp for schedule(static)
    for(int i = 0; i < size; i++){
        vec_sol[i] = vec1[i] + coef_vec2 * vec2[i];
        vec_sol1[i] = vec_sol[i];
    }
  }
}

//Return the dot product of two vector
double multTwoVector(double* vec1, double* vec2, int size){
  double sum = 0;
  #pragma omp parallel default(shared) reduction(+:sum)
  {
    #pragma omp for schedule(static)
    for(int i = 0; i < size; i++){
        sum += vec1[i] * vec2[i];
    }
  }
  return sum;
}

//Divide every element of a vector by a coef
void divideVec(double* vec, int size, double coef){
  #pragma omp parallel default(shared)
  {
    #pragma omp for schedule(static)
    for(int i = 0; i < size; i++){
        vec[i] = vec[i] / coef;
    }
  }
}

//Add the boundary condition to a vector
void addCL(double* vec, int nb_elem_eta, int nb_col_eta, Parameters parameters, double t, double* H_j){
  #pragma omp parallel default(shared)
  {
    #pragma omp for schedule(static)
    for(int it = nb_elem_eta - nb_col_eta; it < nb_elem_eta; it++){
      int j = it - (nb_elem_eta - nb_col_eta);

      double cl;
      if(parameters.s == 0)
        cl = parameters.A * sin(2*PI*parameters.f*t);
      else
        cl = parameters.A * sin(2*PI*parameters.f*t) * exp(-t/500);

      vec[it] -= (H_j[nb_elem_eta + j] * cl)/parameters.delta_y;
    }
  }
}

//Check if the position is a boundary condition
bool checkCl(int pos, int nb_elem_eta, int nb_elem_u, int nb_elem_v, int nb_col_eta){
  if(pos >= nb_elem_eta && pos < nb_elem_eta + nb_elem_u){ //range u
    int newj = pos-nb_elem_eta;

    if(newj - (newj/(nb_col_eta+1))*(nb_col_eta+1) == 0 || newj - (newj/(nb_col_eta+1))*(nb_col_eta+1) == nb_col_eta){
      return 1;
    }
  }
  if(pos >= nb_elem_eta + nb_elem_u){ //range v
    int newj = pos-nb_elem_eta-nb_elem_u;

    if((newj >= 0 && newj < nb_col_eta) || (newj >= nb_elem_v-nb_col_eta && newj < nb_elem_v)){
      return 1;
    }
  }
  return 0;
}

//Return the number of boundary conditions before a coordinate (in i or j), this coordinate is included
int getNbCl(int pos, int nb_col_eta, int nb_lines_eta, int nb_elem_eta, int nb_elem_u, int nb_elem_v){
  if(pos < nb_elem_eta) //range eta
    return 0;

  if(pos >= nb_elem_eta && pos < nb_elem_eta + nb_elem_u) //range u
    return 1 + (pos-nb_elem_eta)/(nb_col_eta+1) + (pos+1-nb_elem_eta)/(nb_col_eta+1);

  if(pos >= nb_elem_eta + nb_elem_u && pos < nb_elem_eta + nb_elem_u + nb_col_eta)
    return 2*nb_lines_eta + (pos - (nb_elem_eta + nb_elem_u)) % nb_col_eta + 1;

  if(pos >= nb_elem_eta + nb_elem_u + nb_col_eta && pos < nb_elem_eta + nb_elem_u + nb_elem_v - nb_col_eta)
    return 2*nb_lines_eta + nb_col_eta;

  else
    return 2*nb_lines_eta + nb_col_eta + (pos - (nb_elem_eta + nb_elem_u + nb_elem_v - nb_col_eta)) % nb_col_eta + 1;
 
}

//Update the sol vector with the index of the columns that are not a 0 (in the matrix A) from a given line
void getNon0Pos(int j, int nb_col_eta, int nb_elem_eta, int nb_elem_u, int* sol){
  if(j < nb_elem_eta){
    sol[0] = j;
    sol[1] = j+nb_elem_eta + j/nb_col_eta;
    sol[2] = sol[1] + 1;
    sol[3] = j+nb_elem_eta + nb_elem_u;
    sol[4] = sol[3] + nb_col_eta;
  }
  else if(j >= nb_elem_eta && j < nb_elem_eta + nb_elem_u){
    int newj = j - nb_elem_eta -1;
    sol[0] = newj - newj/(nb_col_eta+1);
    sol[1] = sol[0] + 1;
    sol[2] = j;
  }
  else{
    int newj = j - nb_elem_eta - nb_elem_u - nb_col_eta;
    sol[0] = newj;
    sol[1] = newj + nb_col_eta;
    sol[2] = j;
  }
}

//Return the coef of the matrix A on a given (j,i) coordinates
double getCoeffA(int j, int i, int nb_columns_A, int nb_elem_eta, int nb_elem_u, int nb_col_eta,
                 double* H_i, double* H_j, Parameters parameters){
  if(j < nb_elem_eta){ //range eta
    if(i == j){
      return 1/parameters.delta_t;
    }
    if(i == j + nb_elem_eta + j/nb_col_eta){
      return -H_i[i - nb_elem_eta]/parameters.delta_x;
    }
    if(i == j + nb_elem_eta + j/nb_col_eta + 1){
      return H_i[i - nb_elem_eta]/parameters.delta_x;
    }
    if(i == j + nb_elem_eta + nb_elem_u){
      return -H_j[i - nb_elem_eta - nb_elem_u]/parameters.delta_y;
    }
    if(i == j + nb_elem_eta + nb_elem_u + nb_col_eta){
      return H_j[i - nb_elem_eta - nb_elem_u]/parameters.delta_y;
    }
  }
  else if(j >= nb_elem_eta + 1 && j < nb_elem_eta + nb_elem_u - 1){ //range u
    int new_j = j - nb_elem_eta;
    int decal = new_j/(nb_col_eta + 1);

    if(i == new_j - decal){
      return parameters.g / parameters.delta_x;
    }
    if(i == new_j - decal - 1){
      return -parameters.g / parameters.delta_x;
    }
    if(i == new_j + nb_elem_eta){
      return 1/parameters.delta_t + parameters.gamma;
    }
  }
  else if(j >= nb_elem_eta + nb_elem_u){ //range v
    int new_j = j - (nb_elem_eta + nb_elem_u);

    if(i == new_j - nb_col_eta){
      return -parameters.g/parameters.delta_y;
    }
    if(i == new_j){
      return parameters.g/parameters.delta_y;
    }
    if(i == new_j + nb_elem_eta + nb_elem_u){
      return 1/parameters.delta_t + parameters.gamma;
    }
  }
  return 0;
}

//Store in the vector sol_vec the results of the product A*vec where A is a matrix
void matVecProduct(double* vec, double* sol_vec, int nb_lines, int nb_lines_above, int N, bool trans,
                   int nb_elem_eta, int nb_elem_u, int nb_elem_v, int nb_col_eta, int nb_lines_eta,
                   double* H_i, double* H_j, Parameters parameters){

  int nb_CL_before = getNbCl(nb_lines_above-1, nb_col_eta, nb_lines_eta, nb_elem_eta, nb_elem_u, nb_elem_v);

  #pragma omp parallel default(shared)
  {
    #pragma omp for schedule(static)
    for(int j = nb_lines_above; j < nb_lines_above + nb_lines; j++){
      if(checkCl(j, nb_elem_eta, nb_elem_u, nb_elem_v, nb_col_eta)){
        continue;
      }

      int nb_CL = getNbCl(j, nb_col_eta, nb_lines_eta, nb_elem_eta, nb_elem_u, nb_elem_v);
      int pos_vec = j - nb_lines_above - nb_CL + nb_CL_before;
   
      int size_non_0;
      if(nb_CL == 0) //range eta
        size_non_0 = 5;
      else //range u and v
        size_non_0 = 3;

      int* non_0_pos = malloc(size_non_0 * sizeof(int));
      getNon0Pos(j, nb_col_eta, nb_elem_eta, nb_elem_u, non_0_pos);

      //Do the line of A * vec
      for(int it = 0; it < size_non_0; it++){
        int i = non_0_pos[it];

        if(checkCl(i, nb_elem_eta, nb_elem_u, nb_elem_v, nb_col_eta)){
          continue;
        }

        double prod;
        if(!trans){
          prod = getCoeffA(j,i, N, nb_elem_eta, nb_elem_u, nb_col_eta, H_i, H_j, parameters) *
                 vec[i - getNbCl(i, nb_col_eta, nb_lines_eta, nb_elem_eta, nb_elem_u, nb_elem_v)];
        }
        else{
          prod = getCoeffA(i,j, N, nb_elem_eta, nb_elem_u, nb_col_eta, H_i, H_j, parameters) *
                 vec[i - getNbCl(i, nb_col_eta, nb_lines_eta, nb_elem_eta, nb_elem_u, nb_elem_v)];
        }
        if(prod == 0)
          continue;

        #pragma omp atomic
        sol_vec[pos_vec] += prod;
      }
      free(non_0_pos);
    }
  }
}

//Set all the elements of a vector to 0
void setVecToZero(double* vec, int size){
  #pragma omp parallel default(shared)
  {
    #pragma omp for schedule(static)
    for(int i = 0; i < size; i++){
      vec[i] = 0;
    }
  }
}

//Separated the vector x_tot into eta, u and v in order to do the writing
void gather_imp(int myrank, int np, double* x, int nb_lines_A, int nb_columns_A, int nb_lines_H, int nb_columns_H, double* x_tot, int saveNb){
  if(np > 1)
    getWholeVector(myrank, np, x, nb_lines_A, x_tot);
  else
    memcpy(x, x_tot, nb_lines_A * sizeof(double));

  if(myrank == 0){
    int size_eta = nb_lines_H * nb_columns_H;
    int size_u = nb_lines_H * (nb_columns_H + 1) - 2*nb_lines_H;
    int size_v = (nb_lines_H + 1) * nb_columns_H - 2*nb_columns_H;

    //Get eta, u and v
    double* eta = malloc(size_eta * sizeof(double));
    double* u = malloc(size_u * sizeof(double));
    double* v = malloc(size_v * sizeof(double));

    #pragma omp parallel default(shared)
    {
      #pragma omp for schedule(static)
      for(int i = 0; i < size_eta + size_u + size_v; i++){
        if(i < size_eta)
          eta[i] = x_tot[i];

        else if(i >= size_eta && i < size_eta + size_u)
          u[i-size_eta] = x_tot[i];

        else
          v[i-size_eta - size_u] = x_tot[i];
      }
    }
    toWrite_imp(eta, u, v, (unsigned int) nb_lines_H, (unsigned int) nb_columns_H, saveNb);
    free(u);
    free(v);
    free(eta);
  }
}

//Do the implicit resolution
void implicitEuler(int myrank, int np, Parameters parameters, Map map, double* h,
                  double delta_xh, double delta_yh, MPI_Status status){
  //Get the number of lines for the current process
  int nb_columns_H = (int) map.a/parameters.delta_x + 1;
  int nb_lines_H = (int) map.b/parameters.delta_y + 1;

  //Get the H in all the process
  double* H_i_tot = malloc(nb_lines_H * (nb_columns_H + 1) * sizeof(double));
  double* H_j_tot = malloc((nb_lines_H + 1) * nb_columns_H * sizeof(double));

  interpolateH(myrank, np, nb_columns_H, nb_lines_H, parameters, delta_xh, delta_yh,
               map, h, H_i_tot, H_j_tot);

  //Nb elements in b, x0, p,...
  int nb_elem_eta = nb_columns_H * nb_lines_H;
  int nb_elem_u = (nb_columns_H + 1) * nb_lines_H;
  int nb_elem_v = nb_columns_H * (nb_lines_H + 1);

  int nb_columns_A = nb_elem_eta + nb_elem_u + nb_elem_v;

  int nb_lines_tot_A = nb_columns_A; //square matrix

  int nb_lines_A = (int)floor(nb_columns_A/np); //per process
  int mod = nb_lines_tot_A - np*nb_lines_A;

  if(myrank < mod)
    nb_lines_A += 1;

  //Get the number of lines of A above
  int nb_lines_above;
  if(myrank < mod)
    nb_lines_above = myrank * nb_lines_A;
  else
    nb_lines_above = mod * (nb_lines_A + 1) + (myrank - mod) * nb_lines_A;

  //get the number of boundary conditions
  int nb_CL_tot = 2*nb_lines_H + 2*nb_columns_H;
  int nb_CL_process = getNbCl(nb_lines_A - 1 + nb_lines_above, nb_columns_H, nb_lines_H, nb_elem_eta, nb_elem_u, nb_elem_v);
  int nb_CL_before = getNbCl(nb_lines_above-1, nb_columns_H, nb_lines_H, nb_elem_eta, nb_elem_u, nb_elem_v);

  //get the vector size without the boundary conditions
  int whole_vector_size = nb_columns_A - nb_CL_tot;
  int sub_vector_size = nb_lines_A - nb_CL_process + nb_CL_before;

  //Set the initial condition
  double* x = calloc(sub_vector_size, sizeof(double));
  double* r = malloc(sub_vector_size * sizeof(double));
  double* p = malloc(sub_vector_size * sizeof(double));
  double* vec_tot = malloc(whole_vector_size * sizeof(double));
  double* vec_sub = malloc(sub_vector_size * sizeof(double));
  double* vec_sub2 = malloc(sub_vector_size * sizeof(double));

  //Start the time loop
  for(int t = 1; t <= parameters.T_max/parameters.delta_t; t++){
    //get the whole vector x 
    getWholeVector(myrank, np, x, sub_vector_size, vec_tot);

    //get the vector Ax
    setVecToZero(vec_sub, sub_vector_size);
    matVecProduct(vec_tot, vec_sub, nb_lines_A, nb_lines_above, nb_columns_A, 0, nb_elem_eta, nb_elem_u, nb_elem_v, nb_columns_H, nb_lines_H, H_i_tot, H_j_tot, parameters);

    //get the whole vector Ax
    getWholeVector(myrank, np, vec_sub, sub_vector_size, vec_tot);

    //get AtAx
    setVecToZero(vec_sub, sub_vector_size);
    matVecProduct(vec_tot, vec_sub, nb_lines_A, nb_lines_above, nb_columns_A, 1, nb_elem_eta, nb_elem_u, nb_elem_v, nb_columns_H, nb_lines_H, H_i_tot, H_j_tot, parameters);
    
    //get the whole vector b
    getWholeVector(myrank, np, x, sub_vector_size, vec_tot);
    divideVec(vec_tot, whole_vector_size, parameters.delta_t);

    //add the CL
    addCL(vec_tot, nb_elem_eta, nb_columns_H, parameters, t*parameters.delta_t, H_j_tot);

    //get Atb
    setVecToZero(vec_sub2, sub_vector_size);
    matVecProduct(vec_tot, vec_sub2, nb_lines_A, nb_lines_above, nb_columns_A, 1, nb_elem_eta, nb_elem_u, nb_elem_v, nb_columns_H, nb_lines_H, H_i_tot, H_j_tot, parameters);

    //r = Atb - AtAx && p = r
    addTwoVectorSpe(vec_sub2, vec_sub, r, p, sub_vector_size, -1);

    double r0_ = multTwoVector(r, r, sub_vector_size);

    double r0, norm_ri;
    MPI_Allreduce(&r0_, &r0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    r0 = sqrt(r0);
    norm_ri = r0;

    //**************** START OF CONJUGATE GRADIENT *******************
    while(norm_ri / r0 >= parameters.r_tresh){
      //get the whole vector p
      getWholeVector(myrank, np, p, sub_vector_size, vec_tot);

      //get the vector Ap
      setVecToZero(vec_sub, sub_vector_size);
      matVecProduct(vec_tot, vec_sub, nb_lines_A, nb_lines_above, nb_columns_A, 0, nb_elem_eta, nb_elem_u, nb_elem_v, nb_columns_H, nb_lines_H, H_i_tot, H_j_tot, parameters);

      //get the whole vector Ap
      getWholeVector(myrank, np, vec_sub, sub_vector_size, vec_tot);

      //get AtAp
      setVecToZero(vec_sub, sub_vector_size);
      matVecProduct(vec_tot, vec_sub, nb_lines_A, nb_lines_above, nb_columns_A, 1, nb_elem_eta, nb_elem_u, nb_elem_v, nb_columns_H, nb_lines_H, H_i_tot, H_j_tot, parameters);

      double r_all_ = multTwoVector(r, r, sub_vector_size);
      double r_all;
      MPI_Allreduce(&r_all_, &r_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); //ri * ri

      double p_all_ = multTwoVector(p, vec_sub, sub_vector_size);
      double p_all;
      MPI_Allreduce(&p_all_, &p_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); //pi * AtApi

      double alpha = r_all / p_all; //alpha

      addTwoVector(x, p, x, sub_vector_size, alpha); //xi+1
      addTwoVector(r, vec_sub, r, sub_vector_size, -alpha); //ri+1

      double r_all1_ = multTwoVector(r,r, sub_vector_size);
      double r_all1;
      MPI_Allreduce(&r_all1_, &r_all1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); //ri1 * ri1

      double beta = r_all1 / r_all; //beta

      addTwoVector(r, p, p, sub_vector_size, beta); //pi+1

      norm_ri = sqrt(r_all1);
    }
    //**************** END OF CONJUGATE GRADIENT *******************

    //writing part
    if(parameters.S!=0 && t % parameters.S == 0 ){
      double* x_tot = malloc(whole_vector_size * sizeof(double));
      gather_imp(myrank, np, x, sub_vector_size, whole_vector_size, nb_lines_H, nb_columns_H, x_tot, t);
      free(x_tot);
    }
  }

  free(x);
  free(r);
  free(p);
  free(H_i_tot);
  free(H_j_tot);
  free(vec_sub);
  free(vec_sub2);
  free(vec_tot);
}