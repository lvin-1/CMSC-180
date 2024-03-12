#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>
#include <sys/time.h>

void populateMatrix(int** mat, int row, int col){
    for (int j = 0; j<(row); j++){
		for (int k = 0; k<col; k++){
            int item = (rand() % 50) + 1;
            mat[j][k] = item;
            // mat[row][col]++;
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
    double *vecX = (double *)malloc(row * sizeof(double));

    // Outer loop that iterates from 0-size, to compute the pearson_cor of each column in the matrix
    for (int i = 0; i<row; i++){
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
        for(int j = 0; j<col; j++){
            sumX += mat[i][j];
            sumY += vecY[j];
            sumX2 += mat[i][j] * mat[i][j];
            sumY2 += vecY[j] * vecY[j];
            sumXY += mat[i][j] * vecY[j];
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

typedef struct ARG{
    int** mat;
    int* vecY; 
    int row;
    int col;
    double* result;
    // int tnum;
}args;

void* pearson_cor_thread(void* arg){
    args * temp;
    temp = (args *) arg;
    temp->result = pearson_cor(temp->mat,temp->vecY,temp->row,temp->col);
    // printf("Thread %d Done!\n",temp->tnum);
    pthread_exit(NULL);
}

int** transposeMatrix(int ** mat, int n){
    int ** transposedMat = (int**) malloc(n*sizeof(int*));
    for (int i = 0; i<n; i++){
        transposedMat[i] = (int*) malloc(n*sizeof(int));
    }

    for(int j = 0; j<n; j++){
        for(int k = 0; k<n; k++){
            transposedMat[k][j] = mat[j][k];
        }
    }

    return transposedMat;
}

/**/

int main(){

    // SAMPLE VALUES OF CHECKING

    // float sampleX[10] = {3.63 , 3.02 , 3.82 , 3.42 , 3.59 , 2.87 , 3.03 , 3.46 , 3.36 , 3.30};
    // float sampleY[10] = {53.1 , 49.7 , 48.4 , 54.2 , 54.9 , 43.7 , 47.2 , 45.2 , 54.4 , 50.4};

    // sample_pearson_cor(sampleX,sampleY,10);

    time_t t;
    int n,num_threads;
    int rows = 0;
    int ***submatrices;
    printf("Input n: ");
    scanf("%d",&n);
    printf("Input t: ");
    scanf("%d",&num_threads);

    args* arguments = (args*) malloc(sizeof(args)*num_threads);   
    pthread_t tid[num_threads];
    
    // Matrix Creation
    int **matrix = (int **)malloc(n * sizeof(int*));
    for(int i = 0; i < n; i++){
        matrix[i] = (int*)malloc(n * sizeof(int));
    } 

    // Populate the matrix
    populateMatrix(matrix,n,n);

    matrix = transposeMatrix(matrix,n); //Transposing the matrix

    // Creation of vector y;
    int *vecY = (int *)malloc(n * sizeof(int));

    // Populate vector y
    populateVector(vecY,n);

    // Split matrix into submatirces
    rows = n/num_threads;
    submatrices = (int***) malloc(num_threads*sizeof(int**));
    printf("Submatrices will be %d x %d\n",rows,n);

    int start = 0;
    // Calculating the excess row and adding it to the last submatrix
    int excess_row = n - (rows*num_threads);
    printf("Excess row: %d\n",excess_row);
    for(int i=0; i<num_threads; i++){
        int end=rows;
        if(i==(num_threads-1) && num_threads!=1){
            submatrices[i] = (int**) malloc(((rows+excess_row)*sizeof(int*)));
            end = rows+excess_row;
        }else{
            submatrices[i] = (int**) malloc(rows*sizeof(int*));
        }
        for (int j = 0; j<end; j++){
            submatrices[i][j] = (int*) malloc(n*sizeof(int));
            for (int k = 0; k<n; k++){
                submatrices[i][j][k] = matrix[start+j][k];
            }
	}
        start += rows;
    }

    double ** vecX = (double**) malloc(num_threads*sizeof(double*));

    for(int j = 0; j<num_threads; j++){
        if(j == (num_threads-1) && num_threads!=1){
            vecX[j] = (double*) malloc((rows+excess_row)*sizeof(double));
        }else{
            vecX[j] = (double*) malloc((rows)*sizeof(double));
        }
        if (vecX[j] == NULL) {
        printf("Memory allocation failed for vecX[%d]\n", j);
        }
    }

    double* finalR = (double*) malloc(n*sizeof(double));

    struct timeval start_time, end_time;
    double elapsed;
    gettimeofday(&start_time,NULL);
    for(int j = 0; j<num_threads; j++){

        arguments[j].mat = submatrices[j];
        arguments[j].vecY = vecY;
        arguments[j].col = n;

        if (j == (num_threads - 1) && num_threads != 1) {
            arguments[j].row = rows + excess_row;
        } else {
            arguments[j].row = rows;
        }
         // Create thread
        pthread_create(&tid[j],NULL,pearson_cor_thread,(void *) &arguments[j]);
    }

    // Wait for the threads to finish
    for(int k = 0; k<num_threads; k++){
		pthread_join(tid[k], NULL);
	}

    // store the output of threaded pearson_cor into a single vector
    int index = 0;
    for(int i = 0; i<num_threads; i++){
        vecX[i] = arguments[i].result;
        int r = rows;
        if(i== (num_threads - 1) && num_threads != 1) r+=excess_row;
        for(int j = 0; j<r; j++){
            finalR[index] = vecX[i][j];
            index++;
        }
    }
    gettimeofday(&end_time,NULL);


    // elapsed = (end_time.tv_sec - start_time.tv_sec);
    // elapsed += (end_time.tv_nsec - start_time.tv_nsec) / 1000000000.0;
    double elapsed_sec = (end_time.tv_sec - start_time.tv_sec) + (end_time.tv_usec - start_time.tv_usec) / 1000000.0;
    printf("Run Time: %f\n", elapsed_sec);
    
    // printVector(finalR,n);

    // Free allocated memory
    // printf("Freeing arguments\n");
    free(arguments);
    // printf("Freed arguments\n");

    // printf("Freeing vecX\n");
    for(int i = 0; i<num_threads; i++){
        free(vecX[i]); 
    } 
    free(vecX); 
    // printf("Freed vecX\n");

    // printf("Freeing submatrices\n");
    for(int i = 0; i<num_threads; i++){
        int r = rows;
        if(i== (num_threads - 1) && num_threads != 1) r+=excess_row;
        for(int j = 0; j<r; j++){
            free(submatrices[i][j]);
        }
        free(submatrices[i]);
    }
    free(submatrices);
    // printf("Successfully freed submatrices\n");

    // printf("Freeing matrix\n");
    // for(int i = 0; i < n; i++){
    //     free(matrix[i]);
    // }
    // free(matrix);
    // printf("Successfully freed matrix\n");

    // printf("Freeing vecY\n");
    free(vecY);
    free(finalR);
    // printf("Freed vecY\n");
}

// add -lm -lpthread when compiling