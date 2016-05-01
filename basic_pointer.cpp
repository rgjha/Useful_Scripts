#include <iostream> // Basic I/O functions
#include <fstream>
using namespace std;

int main()
{
    
    int data[2] = {3,4} ;
    struct Point {
        int x;
        int y;
    };
    
    Point* p;      // declare pointer to a Point struct
    
    p = new Point; // dynamically allocate a Point
    p->x = data[0];  // set the field values.
    cout << " X is " << (*p).x << endl ;
    p->y = 34;
    (*p).x = data[1];
    //cout << " X is " << p->x << endl ;
    cout << " X is " << (*p).x << endl ;
    return 0; 
}
