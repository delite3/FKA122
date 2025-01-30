#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Function to dynamically allocate an nx3 2D array
double **allocate_2d_array(int n) {
    double **array = (double **)malloc(n * sizeof(double *));
    if (array == NULL) {
        printf("Memory allocation failed.\n");
        exit(1);
    }
    for (int i = 0; i < n; i++) {
        array[i] = (double *)malloc(3 * sizeof(double));
        if (array[i] == NULL) {
            printf("Memory allocation failed.\n");
            exit(1);
        }
    }
    return array;
}

// Function to calculate the Euclidean distance between two points
double calculate_distance(double *point1, double *point2) {
    return sqrt(pow(point2[0] - point1[0], 2) +
                pow(point2[1] - point1[1], 2) +
                pow(point2[2] - point1[2], 2));
}

// Function to write the 2D array to a CSV file
void write_to_csv(double **points, int n, const char *filename) {
    FILE *file = fopen(filename, "w");  // Open the file for writing
    if (file == NULL) {
        printf("Error opening file.\n");
        return;
    }

    // Write each point to the file in CSV format
    for (int i = 0; i < n; i++) {
        fprintf(file, "%.2f,%.2f,%.2f\n", points[i][0], points[i][1], points[i][2]);
    }

    fclose(file);  // Close the file
    printf("Data written to %s\n", filename);
}

int main() {
    int n;

    // Prompt the user to enter the number of points
    printf("Enter the number of points (rows): ");
    scanf("%d", &n);

    // Allocate the nx3 array
    double **points = allocate_2d_array(n);

    // Fill the array with values of your choice (example values here)
    printf("Enter values for each point (3 coordinates per point):\n");
    for (int i = 0; i < n; i++) {
        printf("Point %d:\n", i + 1);
        for (int j = 0; j < 3; j++) {
            printf("Coordinate %d: ", j + 1);
            scanf("%lf", &points[i][j]);
        }
    }

    // Calculate the distance between the first and second point as an example
    if (n >= 2) { // Ensure there are at least two points
        double distance = calculate_distance(points[0], points[1]);
        printf("Distance between Point 1 and Point 2: %.2f\n", distance);
    } else {
        printf("Not enough points to calculate a distance.\n");
    }

    // Write the matrix to a CSV file
    write_to_csv(points, n, "points.csv");

    // Free the allocated memory
    for (int i = 0; i < n; i++) {
        free(points[i]);
    }
    free(points);

    return 0;
}
