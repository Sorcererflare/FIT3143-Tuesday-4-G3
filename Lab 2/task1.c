#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

// is_prime: uses sqrt() to check if a number is prime
bool is_prime(int k) {
    if (k <= 1) return false;
    if (k == 2) return true;
    if (k % 2 == 0) return false; 

    int limit = (int)sqrt(k);
    for (int i = 3; i <= limit; i += 2) {
        if (k % i == 0) {
            return false;
        }
    }
    return true;
}

int main() {
	
    int n;
    printf("Enter an integer n: ");
    if (scanf("%d", &n) != 1) {
        return 1;
    }

    FILE *output_file = fopen("task1.txt", "w");
    if (output_file == NULL) {
        fprintf(stderr, "Error opening file.\n");
        return 1;
    }

    
    // include 2
    if (n >= 2) {
        fprintf(output_file, "2\n");
    }
    
	clock_t start_time = clock();
	
    // check numbers from 3 to n
    for (int i = 3; i < n; i += 2) {
        if (is_prime(i)) {
            fprintf(output_file, "%d\n", i);
        }
    }

    clock_t end_time = clock();
    double time_taken = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;

    fprintf(output_file, "Serial execution time: %f seconds\n", time_taken);
    fclose(output_file);
    
    return 0;
}