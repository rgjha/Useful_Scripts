/* Simple sort program */
/* R.Jha, Aug 2015 */ 

#include <iostream>
using namespace std;

// Declaration of swap function
void swap(int &m, int &n);

int main()
{
    int N = 6;
    int a[6] = {52, 1002,
        2001, 0, 79, 99};
    
    // Selection Sort
    for (int i = 0; i < (N - 1); i++)
    {
        int minIndex = i;
        
        // Find the index of the minimum element
        for (int j = i + 1; j < N; j++)
        {
            if (a[j] < a[minIndex])
            {
                minIndex = j;
            }
        }
        
        // Swap if i-th element not already smallest
        if (minIndex > i)
        {
            swap(a[i], a[minIndex]);
        }
    }
    
    
    // Print sorted results
    for (int i = 0; i < N; i++)
    {
        cout << i << " " << a[i] << endl;
    }
    
    return 0;
}

void swap(int &x, int &y)
{
    int temp;
    
    temp = x;
    x = y;
    y = temp;
}

