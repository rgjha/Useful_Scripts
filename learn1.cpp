#include <iostream> // Basic I/O functions
#include <fstream>
using namespace std;

int main()
{
    int x;
    int *ptr_p;
    
    x = 5;
    ptr_p = &x;
    
    cout << ptr_p << endl ;
    return 0;
}
