#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>
#include <sys/time.h>
#include <sched.h>

void populateMatrix(int** mat, int row, int col){
    for (int j = 0; j<(row); j++){
		for (int k = 0; k<col; k++){
            int item = (rand() % 50) + 1;
            mat[j][k] = item;
		}
	}
}

void printMatrix(int** mat,int row, int col){
    for (int j = 0; j<(row); j++){
		for (int k = 0; k<col; k++){
            printf("%d\t",mat[j][k]);
		}
        printf("\n");
	}
}

void populateVector(int* vec, int n){
    for (int k = 0; k<n; k++){
        int item = (rand() % 50) + 1;
        vec[k] = item;
    }
}

void printVector(double* vec, int n){
    for (int k = 0; k<n; k++) printf("%f\t",vec[k]);
    printf("\n");
}

void sample_pearson_cor(float* mat, float* vecY, int size){
    float sumX = 0;
    float sumY = 0;
    float sumX2 = 0;
    float sumY2 = 0;
    float sumXY = 0;
    float numerator = 0;
    double denominator = 0;
    double ans = 0;
    
    // USE THIS for testing (same algorithm for computing the pearson correlation)
    for(int i = 0; i<size; i++){
        sumX += mat[i];
        sumY += vecY[i];
        sumX2 += mat[i] * mat[i];
        sumY2 += vecY[i] * vecY[i];
        sumXY += mat[i] * vecY[i];
    }

    numerator = (size*sumXY) - (sumX*sumY);
    denominator = (((size*sumX2)-(sumX*sumX))*((size*sumY2)-(sumY*sumY)));
    denominator = sqrt(denominator);
    // printf("%f/%f\n",numerator,denominator);
    ans = numerator/denominator;
    printf("r = %f\n",ans);
}

// Function that will compute the pearson correlation coefficient of a given matrix and a vector
double* pearson_cor(int** mat, int* vecY, int row, int col){
    // Creation of vector x;
    double *vecX = (double *)malloc(col * sizeof(double));

    // Outer loop that iterates from 0-size, to compute the pearson_cor of each column in the matrix
    for (int i = 0; i<col; i++){
        float sumX = 0;
        float sumY = 0;
        float sumX2 = 0;
        float sumY2 = 0;
        float sumXY = 0;
        float numerator = 0;
        double denominator = 0;
        double ans = 0;
        // Inner loop that iterates from 0-size to get the summation of values inside the vectors that will be needed to compute
        // the pearson_cor
        for(int j = 0; j<row; j++){
            sumX += mat[j][i];
            sumY += vecY[j];
            sumX2 += mat[j][i] * mat[j][i];
            sumY2 += vecY[j] * vecY[j];
            sumXY += mat[j][i] * vecY[j];
        }
        // Using the summations to get the final value
        numerator = (row*sumXY) - (sumX*sumY);
        denominator = (((row*sumX2)-(sumX*sumX))*((row*sumY2)-(sumY*sumY)));
        denominator = sqrt(denominator);
        ans = numerator/denominator;
        // stores the final value to the appropriate position of vector X
        vecX[i] = ans;
    }
    // printVector(vecX,col);
    return vecX;
    
}

// void setCPU(int core_id, pthread_attr_t * attr){
//     // pthread_attr_t attr;
//     cpu_set_t cpuset;

//     // Initialize thread attributes
//     pthread_attr_init(attr);

//     // Clear the CPU set
//     CPU_ZERO(&cpuset);
//     // Set the CPU core
//     CPU_SET(core_id, &cpuset);

//     // Set CPU affinity attribute
//     pthread_setaffinity(0, sizeof(cpu_set_t), &cpuset);
// }

typedef struct ARG{
    int** mat;
    int* vecY; 
    int row;
    int col;
    int cpu;
    double* result;
    // int tnum;
}args;

void* pearson_cor_thread(void* arg){
    args * temp;
    temp = (args *) arg;
    
    // assign thread to a single CPU
    cpu_set_t cpuset;
    // Clear the CPU set
    CPU_ZERO(&cpuset);
    // Set the CPU core
    CPU_SET(temp->cpu, &cpuset);
    sched_setaffinity(0,sizeof(cpuset),&cpuset);

    temp->result = pearson_cor(temp->mat,temp->vecY,temp->row,temp->col);
    // printf("Thread %d Done!\n",temp->tnum);
    pthread_exit(NULL);
}


int main(){

    time_t t;
    int n,num_threads;
    int cols = 0;
    struct timeval start_time, end_time;
    double elapsed;
    int cpu_to_use = 4;
    int curr_cpu = 1;


    printf("Input n: ");
    scanf("%d",&n);
    printf("Input t: ");
    scanf("%d",&num_threads);

    args* arguments = (args*) malloc(num_threads * sizeof(args));   
    pthread_t tid[num_threads];
    // pthread_attr_t attr[num_threads];
    
    // Matrix Creation
    int **matrix = (int **)malloc(n * sizeof(int*));
    for(int i = 0; i < n; i++){
        matrix[i] = (int*)malloc(n * sizeof(int));
    } 

    // Populate the matrix
    populateMatrix(matrix,n,n);

    // Creation of vector y;
    int *vecY = (int *)malloc(n * sizeof(int));

    // Populate vector y
    populateVector(vecY,n);

    // Split matrix into submatirces
    cols = n/num_threads;
    int*** submatrices = (int***) malloc(num_threads*sizeof(int**));
    printf("Submatrices will be %d x %d\n",n,cols);

    int start = 0;
    // Calculating the excess column and adding it to the last submatrix
    int excess_column = n - (cols*num_threads);
    printf("Excess column: %d\n",excess_column);
    for(int i=0; i<num_threads; i++){
        submatrices[i] = (int**) malloc(n*sizeof(int*));
        for (int j = 0; j<n; j++){
            int end=cols;
            if(i==(num_threads-1) && num_threads!=1){
                submatrices[i][j] = (int*) malloc(((cols+excess_column)*sizeof(int)));
                end = cols+excess_column;
            }else{
                submatrices[i][j] = (int*) malloc(cols*sizeof(int));
            }
            for (int k = 0; k<end; k++){
                submatrices[i][j][k] = matrix[j][start+k];
            }
	}
        start += cols;
    }

    double ** vecX = (double**) malloc(num_threads*sizeof(double*));

    for(int j = 0; j<num_threads; j++){
        if(j == (num_threads-1) && num_threads!=1){
            vecX[j] = (double*) malloc((cols+excess_column)*sizeof(double));
        }else{
            vecX[j] = (double*) malloc((cols)*sizeof(double));
        }
        if (vecX[j] == NULL) {
        printf("Memory allocation failed for vecX[%d]\n", j);
        }
    }

    double* finalR = (double*) malloc(n*sizeof(double));

    gettimeofday(&start_time,NULL);
    for(int j = 0; j<num_threads; j++){

        arguments[j].mat = submatrices[j];
        arguments[j].vecY = vecY;
        arguments[j].row = n;
        arguments[j].cpu = curr_cpu;

        if (j == (num_threads - 1) && num_threads != 1) {
            arguments[j].col = cols + excess_column;
        } else {
            arguments[j].col = cols;
        }
        // Assign thread to a CPU core
        if(curr_cpu == 3) curr_cpu=0;

         // Create thread
        pthread_create(&tid[j],NULL,pearson_cor_thread,(void *) &arguments[j]);
        curr_cpu++;
    }

    // Wait for the threads to finish
    for(int k = 0; k<num_threads; k++){
		pthread_join(tid[k], NULL);
        // pthread_attr_destroy(&attr[k]);
	}

    // store the output of threaded pearson_cor into a single vector
    int index = 0;
    for(int i = 0; i<num_threads; i++){
        vecX[i] = arguments[i].result;
        int column = cols;
        if(i== (num_threads - 1) && num_threads != 1) column+=excess_column;
        for(int j = 0; j<column; j++){
            finalR[index] = vecX[i][j];
            index++;
        }
    }
    gettimeofday(&end_time,NULL);

    double elapsed_sec = (end_time.tv_sec - start_time.tv_sec) + (end_time.tv_usec - start_time.tv_usec) / 1000000.0;
    printf("Run Time: %f\n", elapsed_sec);
    
    free(arguments);
    free(vecY);
    free(finalR);

    for(int i = 0; i<num_threads; i++){
        for(int j = 0; j<n; j++){
            free(submatrices[i][j]);
        }
        free(submatrices[i]);
    }
    free(submatrices);

    for(int i = 0; i<num_threads; i++){
        free(vecX[i]); 
    } 
    free(vecX); 

    for(int j = 0; j<n; j++){
        free(matrix[j]);
    }
    free(matrix);
}

// add -lm -lpthread when compiling