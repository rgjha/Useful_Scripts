#include <iostream>
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

using namespace std ;
int main () {
    Solution sol;
    int i;
    cout << "Please enter an integer value :: " << endl ;
    cin >> i ;
    cout << "Reversed Number is :: " << sol.reverse(i) << endl ;
    return 0;
    
}