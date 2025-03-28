#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

double cpuSecond()
{
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    return ((double)ts.tv_sec + (double)ts.tv_nsec * 1.e-9);
}

void matrix_vector_product(double *a, double *b, double *c, size_t m, size_t n)
{
    for (int i = 0; i < m; i++)
    {
        c[i] = 0.0;
        for (int j = 0; j < n; j++)
            c[i] += a[i * n + j] * b[j];
    }
}

void matrix_vector_product_omp(double *a, double *b, double *c, size_t m, size_t n, int threads)
{
#pragma omp parallel num_threads(threads)
    {
        int nthreads = omp_get_num_threads();
        int threadid = omp_get_thread_num();
        int items_per_thread = m / nthreads;
        int lb = threadid * items_per_thread;
        int ub = (threadid == nthreads - 1) ? (m - 1) : (lb + items_per_thread - 1);
        for (int i = lb; i <= ub; i++)
        {
            c[i] = 0.0;
            for (int j = 0; j < n; j++)
                c[i] += a[i * n + j] * b[j];
        }
    }
}

double run_serial(size_t n, size_t m)
{
    double *a, *b, *c;
    a = (double*)malloc(sizeof(*a) * m * n);
    b = (double*)malloc(sizeof(*b) * n);
    c = (double*)malloc(sizeof(*c) * m);

    if (a == NULL || b == NULL || c == NULL)
    {
        free(a);
        free(b);
        free(c);
        printf("Error allocate memory!\n");
        exit(1);
    }

    for (size_t i = 0; i < m; i++)
    {
        for (size_t j = 0; j < n; j++)
            a[i * n + j] = i + j;
    }

    for (size_t j = 0; j < n; j++)
        b[j] = j;

    double t = cpuSecond();
    matrix_vector_product(a, b, c, m, n);
    t = cpuSecond() - t;

    printf("Elapsed time (serial): %.6f sec.\n", t);
    free(a);
    free(b);
    free(c);
    return t;
}

double run_parallel(size_t n, size_t m, int threads, double t1)
{
    double *a, *b, *c;

    a = (double*)malloc(sizeof(*a) * m * n);
    b = (double*)malloc(sizeof(*b) * n);
    c = (double*)malloc(sizeof(*c) * m);

    if (a == NULL || b == NULL || c == NULL)
    {
        free(a);
        free(b);
        free(c);
        printf("Error allocate memory!\n");
        exit(1);
    }

    for (size_t i = 0; i < m; i++)
    {
        for (size_t j = 0; j < n; j++)
            a[i * n + j] = i + j;
    }

    for (size_t j = 0; j < n; j++)
        b[j] = j;
    

    double t = cpuSecond();
    matrix_vector_product_omp(a, b, c, m, n, threads);
    t = cpuSecond() - t;

    printf("Elapsed time (%d threads): %.6lf sec. S_p - %.6lf\n", threads, t, t1/t);
    free(a);
    free(b);
    free(c);
    return t;
}

int main(int argc, char *argv[])
{
    size_t M = 20000;
    size_t N = 20000;
    double t1;
    int threads[] = {1, 2, 4, 7, 8, 16, 20, 40};

    
    for(int z=0; z<2; z++){
        printf("size: %lux%lu\n", M, N);
        t1 = run_serial(M, N);
        for(int i=0; i<8; i++)
            run_parallel(M, N, threads[i], t1);
        M += M;
        N += N;
    }

    return 0;
}