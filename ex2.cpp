#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>

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

void multiply(int** mat1, int** mat2, int** res, int n, int k, int m, int modulo = 2) {
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < m; ++j) {
            res[i][j] = 0;
        }
    }

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

void transpose(int** a, int** res, int n, int m) {
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < m; ++j) {
            res[j][i] = a[i][j];
        }
    }
}

void add(int** mat1, int** mat2, int** res, int n, int m, int modulo = 2)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            res[i][j] = mat1[i][j] + mat2[i][j];
            if (modulo > 0)
            {
                res[i][j] = res[i][j]%modulo;
            }
        }
    }
}

void copy(int **dst, int **src, int n, int m)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            dst[i][j] = src[i][j];
        }
    }
}

void hstack(int **mat1, int **mat2, int **res, int n, int m)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            res[i][j] = mat1[i][j];
            res[i][j+m] = mat2[i][j];
        }
    }
}

int sum(int** mat, int n, int m)
{
    int result = 0;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            result += mat[i][j];
        }
    }
    return result;
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

void circulant(int* a, int ** res, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            res[i][j] = a[(n - j + i)%n];
        }
    }
}

bool BitFlipping(int** h0, int** h1, int** s, int** resu, int **resv, int n, int T, int w, int max_iter = 1000, bool smooth = true)
{
    bool res = true;
    int iter = 0;
    int **H = allocate_matrix(n, 2*n);
    int **uv = allocate_matrix(2*n, 1);
    int **flipu = allocate_matrix(1, n);
    int **flipv = allocate_matrix(1, n);
    int **total = allocate_matrix(1, 2*n);
    int **syndrome = allocate_matrix(1, n);
    int **syndrome_temp = allocate_matrix(n, 1);
    int **flipped_position = allocate_matrix(2*n, 1);

    hstack(h0, h1, H, n, n);
    copy(syndrome, s, 1, n);
    int lastlasts = 0, lasts = sum(syndrome, 1, n);

    while ((sum(resu, 1, n) != w || sum(resv, 1, n) != w && lasts != 0) && (iter < max_iter))
    {
        iter += 1;
        multiply(syndrome, H, total, 1, n, 2*n, -1);
        // std::cout << sum(resu, 1, n) << " " << sum(resv, 1, n) << " " << lasts << " " << T << " " << sum(flipped_position, 2*n, 1) << std::endl;
        for (int i = 0 ; i < n; i++)
        {
            flipu[0][i] = int(total[0][i] >= T);
            flipv[0][i] = int(total[0][i+n] >= T);
            flipped_position[i][0] = int(total[0][i] >= T);
            flipped_position[i+n][0] = int(total[0][i+n] >= T);
        }
        add(resu, flipu, resu, 1, n);
        add(resv, flipv, resv, 1, n);

        multiply(H, flipped_position, syndrome_temp, n, 2*n, 1);
        for (int i = 0; i < n; i++)
        {
            syndrome[0][i] = (syndrome[0][i] + syndrome_temp[i][0]) % 2;
        }

        if (smooth){
            if (sum(flipped_position, 2*n, 1) == 0)
            {
                T = T-1;
            }
            else if (sum(syndrome, 1, n) == lastlasts || sum(syndrome, 1, n) == lasts)
            {
                T = std::max(T+1, 2)*2;
                for (int i = 0; i < n; i++){
                    resu[0][i] = 0;
                    resv[0][i] = 0;
                    syndrome[0][i] = s[0][i];
                }
            }
            else if (sum(flipped_position, 2*n, 1) > 50*w)
            {
                T += sum(flipped_position, 2*n, 1) / w / 20;
                for (int i = 0; i < n; i++){
                    resu[0][i] = 0;
                    resv[0][i] = 0;
                    syndrome[0][i] = s[0][i];
                }
            }
        }

        lastlasts = lasts;
        lasts = sum(syndrome, 1, n);
    }
    
    for (int i = 0; i < n; i++)
    {
        uv[i][0] = resu[0][i];
        uv[i+n][0] = resv[0][i];
    }
    multiply(H, uv, syndrome_temp, n, 2*n, 1);
    transpose(syndrome_temp, syndrome, n, 1);
    for (int i = 0; i < n; i++)
    {
        if (syndrome[0][i] != s[0][i])
            res =  false;
    }

    free_matrix(H, n, 2*n);
    free_matrix(uv, 2*n, 1);
    free_matrix(flipu, 1, n);
    free_matrix(flipv, 1, n);
    free_matrix(total, 1, 2*n);
    free_matrix(syndrome, 1, n);
    free_matrix(syndrome_temp, n, 1);
    free_matrix(flipped_position, 2*n, 1);

    return res;
}

void test_BitFlipping(int n, int w, int T, bool smooth)
{
    int **e0 = allocate_matrix(1, n);
    int **e1 = allocate_matrix(1, n);
    int **h0 = allocate_matrix(1, n);
    int **mh0 = allocate_matrix(n, n);
    int **mh0t = allocate_matrix(n, n);
    int **h1 = allocate_matrix(1, n);
    int **mh1 = allocate_matrix(n, n);
    int **mh1t = allocate_matrix(n, n);
    int **s = allocate_matrix(1, n);
    int **s0 = allocate_matrix(1, n);
    int **s1 = allocate_matrix(1, n);
    int **resu = allocate_matrix(1, n);
    int **resv = allocate_matrix(1, n);

    generate_random_vector(e0[0], n, w);
    generate_random_vector(e1[0], n, w);
    generate_random_vector(h0[0], n, w);
    generate_random_vector(h1[0], n, w);
    circulant(h0[0], mh0, n);
    circulant(h1[0], mh1, n);
    transpose(mh0, mh0t, n, n);
    transpose(mh1, mh1t, n, n);


    multiply(e0, mh0t, s0, 1, n, n);
    multiply(e1, mh1t, s1, 1, n, n);
    add(s0, s1, s, 1, n);

    BitFlipping(mh0, mh1, s, resu, resv, n, T, w, 1000, smooth);

    free_matrix(e0, 1, n);
    free_matrix(e1, 1, n);
    free_matrix(h0, 1, n);
    free_matrix(h1, 1, n);
    free_matrix(mh0, n, n);
    free_matrix(mh1, n, n);
    free_matrix(mh0t, n, n);
    free_matrix(mh1t, n, n);
    free_matrix(s, 1, n);
    free_matrix(s0, 1, n);
    free_matrix(s1, 1, n);
    free_matrix(resu, 1, n);
    free_matrix(resv, 1, n);
}

int main(int argc, char *argv[])
{
    int k = 26, w = 39, n = 4813;
    bool smooth = true;
    if (argc > 1)
        k = atoi(argv[1]);
    if (argc > 2)
        w = atoi(argv[2]);
    if (argc > 3)
        n = atoi(argv[3]);
    if (argc > 4)
        smooth = atoi(argv[4]);
    clock_t tStart = clock();
    srand(time(NULL));
    test_BitFlipping(n, w, k, smooth);
    printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
}