#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <time.h>
#include <omp.h>

using namespace std;

int n_size = 2000;
double t = 0.00001;
double e = 0.000003;
vector<int> n_threads = {1,2,4,7,8,16,20,40};

double cpuSecond()
{
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    return ((double)ts.tv_sec + (double)ts.tv_nsec * 1.e-9);
}

int main(){
    vector<double> b(n_size, n_size + 1);
    vector<vector<double>> A(n_size, vector<double>(n_size));
    for (int i=0; i<n_size; i++)
        for (int j=0; j<n_size; j++){
            if (i != j)
                A[i][j] = 1;
            else
                A[i][j] = 2;
        }

    for(int z=0; z<n_threads.size(); z++){
        vector<double> x(n_size, 0);
        vector<double> result(n_size, 0);

        double counter = 0;
        double res1 = 0, res2 = 0;
        double critery = 0.1;
        
        double time = cpuSecond();
        #pragma omp parallel num_threads(n_threads[z])
        {
            int nthreads = omp_get_num_threads();
            int threadid = omp_get_thread_num();
            int items_per_thread = n_size / nthreads;
            int lb = threadid * items_per_thread;
            int ub = (threadid == nthreads) ? n_size : (lb + items_per_thread);

            while (critery >= e)
            {
                for (int i=lb; i<ub; i++){
                    result[i] = 0;
                    for (int j=0; j<n_size; j++)
                        result[i] += A[i][j] * x[j];
                    result[i] -= b[i];
                    result[i] *= t;
                    x[i] -= result[i];
                }            

                #pragma omp barrier
                for (int i=lb; i<ub; i++){
                    result[i] = 0;
                    for (int j=0; j<n_size; j++)
                        result[i] += A[i][j] * x[j];
                    result[i] -= b[i];
                }

                double thread_res = 0;
                for (int i=lb; i<ub; i++)
                    thread_res += pow(result[i], 2);
                
                #pragma omp atomic
                res1 += thread_res;

                thread_res = 0;
                for (int i=lb; i<ub; i++)
                    thread_res += pow(b[i], 2);

                #pragma omp atomic
                res2 += thread_res;

                #pragma omp barrier
                #pragma omp single
                {
                    // cout << "Итерация - " << counter << "    \tКритерий - " << critery << endl;
                    critery = sqrt(res1) / sqrt(res2);
                    res1 = 0, res2 = 0;
                    counter++;
                }
                #pragma omp barrier
            }
        }
        time = cpuSecond() - time;
        cout << "Время работы на " << n_threads[z] << " потоках - " << time << " секунд\n";
    }
}