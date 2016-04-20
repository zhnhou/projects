#include<iostream>

bool accept() {

    char    answer;
    
    std::cout << "Do you want to proceed (y or n)?\n";
    std::cin  >> answer;

    if (answer == 'y') return true;
    return false;
}

bool accept3() {
    
    char    answer;
    int     itry = 0;

    while (itry < 3) {
        std::cout << "Do you want to proceed (y or n)? \n";
        std::cin  >> answer;

        switch (answer) {
        case 'y':
            return true;
        case 'n':
            return false;
        default:
            std::cout << "Sorry, I don't understand that. \n";
            ++ itry;
        }
    }

    std::cout << "I will take that as no. \n";
    return false;
}

void main() {

    accept3();
    
}
