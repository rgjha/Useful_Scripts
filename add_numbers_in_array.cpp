/* Finding largest possible sum of any pair in a given array */
/* Find the largest and 2nd largest element of the array */

/* Aug 14, 2015 */

#include<iostream>
using namespace std;

/* Function which returns largest pair sum. Assumes that
 there are at-least  two elements in arr[] */


int LSP(int arr[], int n) /* Largest Sum Pair (LSP) */
{
    // Initialize First & Second largest element
    int first, second;
    if (arr[0] > arr[1])
    {
        first = arr[0];
        second = arr[1];
    }
    else
    {
        first = arr[1];
        second = arr[0];
    }
    
    // Remaining elements of the array  //
    
    for (int i = 2; i<n; i ++)
    {
        
        if (arr[i] > first)
        {
            second = first;
            first = arr[i];
        }
        
        
        else if (arr[i] > second && arr[i] != first)
            second = arr[i];
    }
    return (first + second);
}


/* Main program to test function */
int main()
{
    int arr[] = {12, 34, 10, 6, 45, 34, 43, 11, 10, 40};
    int n = sizeof(arr)/sizeof(arr[0]);
    cout << "N is " << n << endl ;
    cout << "Max Pair Sum is " << LSP(arr, n) << endl ;
    return 0;
}