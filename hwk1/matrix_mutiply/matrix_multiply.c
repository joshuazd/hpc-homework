//*****************************************************************************
//
// This program calculates the product of a square matrix with itself:
//
// B = A * A
//
// Please keep all code in this single file.
//
//
//*****************************************************************************

#include <stdio.h>
#include <stdlib.h>


void multiply(int *a, int *b,int n) {
    int r,c,i;
    for(r=0;r<n;++r) {
        for(c=0;c<n;++c) {
            for(i=0;i<n;++i) {
                b[n*r+c] += a[n*r+i]*a[c+n*i]; 
            }
        }
    }
}

void print(int *b,int n) {
    printf("\n");
    int r,c;
    for(r=0;r<n;++r) {
        for(c=0;c<n;++c) {
            printf("%d ",b[n*r+c]);
        }
        printf("\n");
    }
}

void read_matrix(FILE *fp_in,int *m,int n) {
    int r,c;
    for(r=0;r<n;++r) {
        for(c=0;c<n;++c) {
            int k;
            fscanf(fp_in,"%d",&k);
            m[n*r+c] = k;
        }
    }
}


int main(int argc, char ** argv)
{
   
   // check command line arguments
    if ( argc != 3 ) {
        printf("This program computes the product of n x n matrix with itself\n");
        printf("Usage: ./matrix_multiply filename n\n");
        exit(0);
    }

   // TODO: parse input arguments
    int n = atoi(argv[2]);
    FILE *fp_in;
    fp_in = fopen(argv[1],"r");
    if(fp_in == NULL) {
        printf("%s not opened, exiting.\n",argv[1]);
        exit(0);
    }

   // TODO: dynamically allocate space for matrix_A (input matrix) in 1d array
    int * matrix_A;  // declare input matrix
    matrix_A = malloc(n*n * sizeof(int));

   // TODO: call function to read data from file and copy into matrix_A
    read_matrix(fp_in,matrix_A,n);
    fclose(fp_in);
    
    // TODO: dynamically allocate space for matrix_B (output matrix) in 1d array
    int * matrix_B;  // declare output matrix
    matrix_B = malloc(n*n * sizeof(int));
    int r,c;
    for(r=0;r<n;++r)
        for(c=0;c<n;++c)
            matrix_B[n*r+c] = 0;
   // TODO: call function to perform matrix multiplication ( matrix_B = matrix_A * matrix_A )
    multiply(matrix_A,matrix_B,n);
   
    // TODO: call function to write results (matrix_B) to stdout
    print(matrix_B,n);
   // TODO: free space allocated for matrix_A and matrix_B
    free(matrix_A);
    free(matrix_B);
    return 0;

}
