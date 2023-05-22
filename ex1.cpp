#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

void print_matrix(int **matrix, int n, int m) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            printf("%d ", matrix[i][j]);
        }
        printf("\n");
    }
}

int **allocate_matrix(int n, int m) {
    int **matrix = new int*[n];
    for (int i = 0; i < n; i++) {
        matrix[i] = new int[m];
        for (int j = 0; j < m; j++)
            matrix[i][j] = 0;
    }
    return matrix;
}

void free_matrix(int **matrix, int n, int m) {
    for (int i = 0; i < n; i++) {
        delete[] matrix[i];
    }
    delete[] matrix;
}

bool binary_matrix_inverse(int **_matrix, int **inverse, int n) {
    int **matrix = allocate_matrix(n, n);

    // Create I matrix, and clone matrix
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            matrix[i][j] = _matrix[i][j];
            inverse[i][j] = (i == j) ? 1 : 0;
        }
    }

    // Gaussian
    for (int i = 0; i < n; i++) {
        if (matrix[i][i] == 0) {
            int swap_row = -1;
            for (int j = i + 1; j < n; j++) {
                if (matrix[j][i] == 1) {
                    swap_row = j;
                    break;
                }
            }

            if (swap_row == -1) {
                return false;
            }

            // Rows swap
            for (int j = 0; j < n; j++) {
                int temp = matrix[i][j];
                matrix[i][j] = matrix[swap_row][j];
                matrix[swap_row][j] = temp;

                temp = inverse[i][j];
                inverse[i][j] = inverse[swap_row][j];
                inverse[swap_row][j] = temp;
            }
        }

        // Remove element below main diagonal
        for (int j = i + 1; j < n; j++) {
            if (matrix[j][i] == 1) {
                for (int k = 0; k < n; k++) {
                    matrix[j][k] ^= matrix[i][k];
                    inverse[j][k] ^= inverse[i][k];
                }
            }
        }
    }

    // Remove the element above the main diagonal
    for (int i = n - 1; i > 0; i--) {
        for (int j = i - 1; j >= 0; j--) {
            if (matrix[j][i] == 1) {
                for (int k = 0; k < n; k++) {
                    matrix[j][k] ^= matrix[i][k];
                    inverse[j][k] ^= inverse[i][k];
                }
            }
        }
    }
    free_matrix(matrix, n, n);
    return true;
}

// n,k x k,m = n,m
void multiply(int** mat1, int** mat2, int** res, int n, int k, int m, int modulo = 2) {
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < m; ++j) {
            res[i][j] = 0;
        }
    }

    // Perform matrix multiplication
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < m; ++j) {
            for(int z = 0; z < k; ++z) {
                res[i][j] += mat1[i][z] * mat2[z][j];
                if(modulo > 0) { 
                    res[i][j] %= modulo;
                }
            }
        }
    }
}

bool checkm(int** m, int** c, int **G, int n, int k, int t)
{
    int count = 0;
    int** cp = allocate_matrix(1, n);
    multiply(m, G, cp, 1, k, n);
    for (int i = 0; i < n; i++)
    {
        count += (cp[0][i] != c[0][i]);
    }
    free_matrix(cp, 1, n);
    return count == t;
}

void generate_random_vector(int* a, int n, int t) {
    for (int i = 0; i < t; i++)
        a[i] = 1;
    for (int i = t; i < n; i++)
        a[i] = 0;

    // Now shuffle the array to randomize the locations of 1s
    for (int i = n - 1; i > 0; i--) {
        // Pick a random index from 0 to i
        int j = rand() % (i + 1);

        // Swap e[i] with the element at random index
        int temp = a[i];
        a[i] = a[j];
        a[j] = temp;
    }
}

void ISD_encode(int** m, int** c, int **G, int n, int k, int t)
{
    int **e = allocate_matrix(1, n);
    generate_random_vector(e[0], n, t);

    multiply(m, G, c, 1, k, n);
    for (int i = 0; i < n; i++)
    {
        c[0][i] = (c[0][i] + e[0][i])%2;
    }
    free_matrix(e, 1, n);
}

void generate_random_information_set(int** Gs, int** G, int** cs, int** c, int n, int k) {
    int** info_set = allocate_matrix(1, n);
    for (int i = 0; i < n; i++) {
        info_set[0][i] = 0;
    }
    for (int i = 0; i < k; i++) {
        int index;
        do {
            index = rand() % n;
        } while (info_set[0][index] == 1);
        info_set[0][index] = 1;
        for (int j = 0; j < k; j++)
        {
            Gs[j][i] = G[j][index];
            cs[0][i] = c[0][index];
        }
    }
    free_matrix(info_set, 1, n);
}

bool ISD_decode(int** m, int** c, int** G, int n, int k, int t, int max_iterations = -1)
{
    int **cs = allocate_matrix(1, k);
    int **Gs = allocate_matrix(k, k);
    int **inv_Gs = allocate_matrix(k, k);

    int count = 0;
    for (int i = 0; i < max_iterations || max_iterations == -1; i++)
    {
        generate_random_information_set(Gs, G, cs, c, n, k);
        if (binary_matrix_inverse(Gs, inv_Gs, k))
        {
            multiply(cs, inv_Gs, m, 1, k, k);
            if (checkm(m, c, G, n, k, t))
            {
                free_matrix(cs, 1, k);
                free_matrix(Gs, k, k);
                free_matrix(inv_Gs, k, k);
                return true;
            }
        }
    }
    
    free_matrix(cs, 1, k);
    free_matrix(Gs, k, k);
    free_matrix(inv_Gs, k, k);

    return false;
}

void test_ISD(int n, int k, int t)
{
    int **G = allocate_matrix(k, n);
    int **m = allocate_matrix(1, k);
    int **c = allocate_matrix(1, n);
    int **mdecode = allocate_matrix(1, k);

    // random values
    for (int i = 0; i < k; i++)
    {
        generate_random_vector(G[i], n, 60*n/100);
    }
    generate_random_vector(m[0], k, k/2);

    ISD_encode(m, c, G, n, k, t);

    ISD_decode(mdecode, c, G, n, k, t);

    free_matrix(G, k, n);
    free_matrix(m, 1, k);
    free_matrix(c, 1, n);
    free_matrix(mdecode, 1, k);
}

int main(int argc, char *argv[]) {
    int n = 1000, k = 500, t = 10;
    if (argc > 1)
        n = atoi(argv[1]);
    if (argc > 2)
        k = atoi(argv[2]);
    if (argc > 3)
        t = atoi(argv[3]);
    clock_t tStart = clock();
    srand(time(NULL));
    test_ISD(n, k, t);
    printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
    return 0;
}