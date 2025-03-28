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

double cpuSecond(){
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    return ((double)ts.tv_sec + (double)ts.tv_nsec * 1.e-9);
}

int main(){
    vector<double> b(n_size, n_size+1);
    vector<vector<double>> A(n_size, vector<double>(n_size));
    for(int i=0; i<n_size; i++)
        for(int j=0; j<n_size; j++){
            if(i != j)
                A[i][j] = 1;
            else
                A[i][j] = 2;
        }

    for(int z=0; z<n_threads.size(); z++){
        vector<double> x(n_size, 0);

        double counter = 0;
        double res1 = 0, res2 = 0;
        double critery = 0.1;
        vector<double> result(n_size, 0);

        double time = cpuSecond();
        while (critery >= e){
            #pragma omp parallel for num_threads(n_threads[z])
            for(int i=0; i<n_size; i++){
                result[i] = 0;
                for (int j=0; j<n_size; j++)
                    result[i] += A[i][j] * x[j];
                result[i] -= b[i];
                result[i] *= t;
                x[i] -= result[i];
            }

            #pragma omp parallel for num_threads(n_threads[z])
            for(int i=0; i<n_size; i++){
                result[i] = 0;
                for (int j = 0; j<n_size; j++)
                    result[i] += A[i][j] * x[j];
                result[i] -= b[i];
            }

            #pragma omp parallel for reduction(+:res1) num_threads(n_threads[z])
            for(int i=0; i<n_size; i++)
                res1 += pow(result[i], 2);
            
            #pragma omp parallel for reduction(+:res2) num_threads(n_threads[z])
            for(int i=0; i<n_size; i++)
                res2 += pow(b[i], 2);

            // cout << "Итерация - " << counter << "    \tКритерий - " << critery << endl;
            critery = sqrt(res1) / sqrt(res2);
            res1 = 0, res2 = 0;
            counter++;
        }

        time = cpuSecond() - time;
        cout << "Время работы на " << n_threads[z] << " потоках - " << time << " секунд\n";
    }
}