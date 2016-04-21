#include<iostream>

int count_x(char* p, char x) {
    if (p == nullptr) return 0;

    int ct = 0;

    // here it seems the end of a char is identified as 0

    for (; *p!=0; ++p) {
        std::cout << *p << '\n';
        if (*p == x) {
            ++ct;
        }
    }
    return ct;
}

void main() {

    char* p = "eventually";
    char  x = 'e';

    std::cout << count_x(p, x) << "\n";
} 
