#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <mpi.h>

//
// allocate a a flattened matrix of "n" elements
//
double *allocMatrix( size_t n)
{
   double *m;
   m = (double *)malloc( n*sizeof(double));
   if( m==NULL) {
      perror( "failed to allocate matrix; ");
   }
   return m;
}

//
// initialise the values of the given matrix "out" of size "nxn" with 0s
//
void init( double *out, size_t n)
{
   size_t i,j;

   for( i=0; i<n; i++) {
      for( j=0; j<n; j++) {
         out[i*n+j] = 0;
      }
   }

}

//
// print the values of a given matrix "out" of size "nxn"
//
void print( double *out, size_t n)
{
   size_t i,j,maxn;

   maxn = (n < 20 ? n : 20);
   for( i=0; i<maxn; i++) {
      printf( "|");
      for( j=0; j<maxn; j++) {
         printf( " %7.2f", out[i*n+j]);
      }
      if( maxn < n) {
         printf( "...|\n");
      } else {
         printf( "|\n");
      }
   }
   if( maxn < n) {
         printf( "...\n");
      }
}

//
// individual step of the 5-point stencil
// computes values in matrix "out" from those in matrix "in"
// assuming both are of size "nxn"
//
void relax( double *in, double *out, size_t n, int rank, int start, int stop)
{

   // size_t i,j;
   // for( i=1; i<n-1; i++) {
   //    for( j=1; j<n-1; j++) {
   //       out[i*n+j] = 0.25*in[(i-1)*n+j]      // upper neighbour
   //                    + 0.25*in[i*n+j]        // center
   //                    + 0.125*in[(i+1)*n+j]   // lower neighbour
   //                    + 0.175*in[i*n+(j-1)]   // left neighbour
   //                    + 0.2*in[i*n+(j+1)];    // right neighbour
   //    }
   // }
   size_t i;

   for( i = start; i < stop; i++){
      //printf("rank %d is at i = %d\n", rank, i);
      if(i >= n && i <= n*n-n && i%n != 0 && i%n != n-1){
         out[i] = 0.25*in[i-n]                // upper
                  + 0.25*in[i]                // center
                  + 0.125*in[i+n]             // lower
                  + 0.175*in[i-1]             // left
                  + 0.2*in[i+1];              // right    
      printf("rank %d is at i = %d. VALUE = %f\n", rank, i, out[i]);

      }else{
        printf("rank %d dit not compute i = %d\n", rank, i);
      }
   }
}

int main (int argc, char *argv[])
{
   double *a,*b;
   size_t n=0;
   int i;
   int max_iter;

   if( argc < 3) {
      printf("call should have two arguments \"%s <n> <iter>\"\n", argv[0]);
      exit(1);
   }
   if( sscanf( argv[1], "%zu", &n) != 1) {
      printf("non size_t value for matrix size\n");
      exit(1);
   }

   if( sscanf( argv[2], "%d", &max_iter) != 1) {
      printf("non int value for # iterations\n");
      exit(1);
   }

   a = allocMatrix(n*n);
   b = allocMatrix(n*n);

   init( a, n);
   init( b, n);

   a[n/4] = 100.0;;
   b[n/4] = 100.0;;

   a[(n*3)/4] = 1000.0;;
   b[(n*3)/4] = 1000.0;;
                                                                                                                                            
//printf( "size   : n = %zu => %d M elements (%d MB)\n",
           //n, (int)(n*n/1000000), (int)(n*n*sizeof(double) / (1024*1024)));
   //printf( "iter   : %d\n", max_iter);

   //print(a, n);

   int my_rank, p;
   int local_start, local_stop, local_n;
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
   MPI_Comm_size(MPI_COMM_WORLD, &p);
   local_n = floor(n*n/p);
   local_start = my_rank * local_n;
   if (my_rank == p -1){
      local_n = local_n + n*n%p;
   }
   //local_start = my_rank * local_n;
   local_stop = local_start + local_n;
   double *local_matrix;
   local_matrix = allocMatrix(local_n);
   printf("Local_n: %d for rank: %d\n", local_n, my_rank);
   for( i=0; i<max_iter; i++) {
      //relax( a, local_matrix, n, my_rank, local_start, local_stop);
      //MPI_Allgather(&local_matrix, 5, MPI_DOUBLE, &a, 5, MPI_DOUBLE, MPI_COMM_WORLD);
      //if(my_rank == 0){
        //printf("matrix a of rank %d\n", my_rank); 
        //print(a, n);
      //}else{
        //printf("matrix a of rank %d\n", my_rank);
        //print(a,n);

      //}
      if(my_rank == 0){
              MPI_Bcast(a, n*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
              relax(a, local_matrix, n, my_rank, local_start, local_stop);
              int i;
              for(i = 0; i < local_stop; i++){
                      a[i] = local_matrix[i];
              }
              for(i = 1; i < p; i++){
                      if(i == p-1){
                        double *temp;
                        temp = allocMatrix(local_n + n*n%p);
                        MPI_Recv(temp, local_n + n*n%p, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        for(int j = 0; j< local_n + n*n%p; j ++){
                           a[local_n*i + j] = temp[j];
                        }
                      }else{
			 double *temp;
                        temp = allocMatrix(local_n);
                        MPI_Recv(temp, local_n, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        for(int j = 0; j < local_n; j ++){
                           a[local_n*i + j] = temp[j];
                        }
                      }
              }
             print(a, n);
      }
      else{
              MPI_Bcast(b, n*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
              relax(b, local_matrix, n, my_rank, local_start, local_stop);
              MPI_Send(local_matrix, local_n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
      }

   }
   printf("HERE\n");
   //if(my_rank == 0){
      printf("rank %d called finalize!\n", my_rank);
      MPI_Finalize();
   //}
   //printf( "Matrix after %d iterations:\n", i);
   //print( b, n);
   printf("OR HERE\n");

  // print(a,n);
   return 0;
}

