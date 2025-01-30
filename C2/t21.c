#include <stdio.h>
#include <stdlib.h>

// Function to compute the scalar product of two arrays
double scalar_product(double *array1, double *array2, int length) {
    double product = 0.0;
    for (int i = 0; i < length; i++) {
        product += array1[i] * array2[i];
    }
    return product;
}

int main() {
    int N;

    // Prompt for length of  vectors
    printf("Enter length of vector: ");
    scanf("%d", &N);

    // Dynamically allocate memory for two arrays of length N
    double *array1 = (double *)calloc(N, sizeof(double));
    double *array2 = (double *)calloc(N, sizeof(double));

    // Check if memory allocation was successful
    if (array1 == NULL || array2 == NULL) {
        perror("Memory allocation failed.\n");
        exit(1);
    }

    // Prompt the user to enter values for each array
    printf("Enter values for the first vector:\n");
    for (int i = 0; i < N; i++) {
        printf("array1[%d]: ", i);
        scanf("%lf", &array1[i]);
    }

    printf("Enter values for the second vector:\n");
    for (int i = 0; i < N; i++) {
        printf("array2[%d]: ", i);
        scanf("%lf", &array2[i]);
    }

    // Calculate the scalar product
    double product = scalar_product(array1, array2, N);

    // Print the result
    printf("The scalar product of the two vectors is: %.2f\n", product);

    // Free the allocated memory
    free(array1);
    free(array2);

    return 0;
}
