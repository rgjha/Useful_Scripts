#include <iostream>
#include <ctime>
#include <locale>
using namespace std ;
class Solution {
public:
    int reverse(int x) {
        if (x == INT_MIN)
            return 0;
        if (x < 0)
            return -reverse(-x);
        
        int num_reverse = 0; // store reversed integer
        while (x != 0) {
            if (num_reverse > INT_MAX / 10 || 10 * num_reverse > INT_MAX - x % 10)
            {
                cout << " Overflow Detected" << endl ;
                return 0; // OR is double Pipe // * :) //
            }
            num_reverse = num_reverse * 10 + x % 10;
            
            x = x / 10; // Move to a digit place on left //
        }
        
        return num_reverse;
    }
};


const std::string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d and %X", &tstruct);
    
    return buf;
}

using namespace std ;
int main () {
    Solution sol;
    int i;
    cout << "Enter an integer value :: " << endl ;
    cin >> i ;
    cout << "Reversed Number is :: " << sol.reverse(i) << endl ;

        /*time_t t = time(0);   // get time now
        struct tm * now = localtime( & t );
        cout << (now->tm_year + 1900) << '-'
        << (now->tm_mon + 1) << '-'
        <<  now->tm_mday
        << endl;*/
    
    
        std::cout << "Current Date and Time :: " << currentDateTime() << std::endl;
        //getchar();  // wait for keyboard input

    

return 0;
    
}
