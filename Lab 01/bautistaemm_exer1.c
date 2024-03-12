#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

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

void printVector(float* vec, int n){
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
void pearson_cor(int** mat, int* vecY, float* vecX, int size){
    // Outer loop that iterates from 0-size, to compute the pearson_cor of each column in the matrix
    for (int i = 0; i<size; i++){
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
        for(int j = 0; j<size; j++){
        sumX += mat[j][i];
        sumY += vecY[j];
        sumX2 += mat[j][i] * mat[j][i];
        sumY2 += vecY[j] * vecY[j];
        sumXY += mat[j][i] * vecY[j];
        }
        // Using the summations to get the final value
        numerator = (size*sumXY) - (sumX*sumY);
        denominator = (((size*sumX2)-(sumX*sumX))*((size*sumY2)-(sumY*sumY)));
        denominator = sqrt(denominator);
        ans = numerator/denominator;
        // stores the final value to the appropriate position of vector X
        vecX[i] = ans;
    }
    
}

/**/

int main(){

    // SAMPLE VALUES OF CHECKING

    float sampleX[10] = {3.63 , 3.02 , 3.82 , 3.42 , 3.59 , 2.87 , 3.03 , 3.46 , 3.36 , 3.30};
    float sampleY[10] = {53.1 , 49.7 , 48.4 , 54.2 , 54.9 , 43.7 , 47.2 , 45.2 , 54.4 , 50.4};

    sample_pearson_cor(sampleX,sampleY,10);

    time_t t;
    int n;
    printf("Input n: ");
    scanf("%d",&n);

    // Matrix Creation
    int **matrix = (int **)malloc(n * sizeof(int*));
    for(int i = 0; i < n; i++){
        matrix[i] = (int*)malloc(n * sizeof(int));
    } 

    // Populate the matrix
    populateMatrix(matrix,n,n);

    // Print Matrix
    // printMatrix(matrix,n,n);
    printf("\n");
    // Creation of vector y;
    int *vecY = (int *)malloc(n * sizeof(int));

    // Populate vector y
    populateVector(vecY,n);
    // printVector(vecY,n);

    // Creation of vector x;
    float *vecX = (float *)malloc(n * sizeof(float));
    t = clock();
    pearson_cor(matrix,vecY,vecX,n);
    t = clock() - t;

    double runTime = ((double)t)/CLOCKS_PER_SEC;
    printf("Run Time: %f\n", runTime);
    
    // printVector(vecX,n);

    // Free allocated memory
    for(int i = 0; i < n; i++){
       free(matrix[i]);
    }
    free(matrix);
    free(vecY); 
    free(vecX); 


}

// add -lm when compiling