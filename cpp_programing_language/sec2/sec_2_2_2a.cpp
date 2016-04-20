#include<iostream>
#include<cmath>

#define PI 3.14159265

using namespace std;

void print_sin_value(double x) {
    cout << "Sin(" << x << " deg) = " << sin(x*PI/180.00e0) << "\n";
}

int main () {
    const int m=1;
    
    int n=10;
    double x {1.0e0};

    print_sin_value(x);
    cout << ++x << "\n";

}
