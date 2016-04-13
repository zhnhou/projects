#include <iostream>
#include <cmath>

using namespace std;

double square_root(double x) {
    return sqrt(x);
}

void print_square_root(double x) {
    cout << "The square root of " << x << " is " << square_root(x) << "\n";
}

int main() {

    print_square_root(1.234e0);

}
