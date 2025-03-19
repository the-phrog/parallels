#include <iostream>
#include <cmath>
#include <vector>

#ifndef USE_DOUBLE
using num = float;
#else
using num = double;
#endif

using namespace std;

const int n = pow(10, 7);
const num pi = 3.14159265358979323846;

int main() {
    vector<num> arr(n);
    num sum = 0;
    
    for(int i=0; i<n; i++){
        arr[i] = sin(2*pi * i/n);
        sum += arr[i];
    }

    cout << sum << endl;
    return 0;
}
