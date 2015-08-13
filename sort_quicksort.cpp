/* Thu Aug 13 19:07:59 2015 /* Quick sort */
#include<stdio.h>
#include <iostream>
#include <fstream>
#include <string>

// Swapping two integers //
void swap(int* a, int* b)
{
    int t = *a;
    *a = *b;
    *b = t;
}

/* This function takes last element as pivot, places the pivot element at its
 correct position in sorted array, and places all smaller (smaller than pivot)
 to left of pivot and all greater elements to right of pivot */

int partition (int arr[], int l, int h)
{
    int x = arr[h]; // pivot
    int i = (l - 1); // Index of smaller element
    
    for (int j = l; j <= h- 1; j++)
    {
        // If current element is smaller than or equal to pivot
        if (arr[j] <= x)
        {
            i++; // increment index of smaller element
            swap(&arr[i], &arr[j]); // Swap current element with index
        }
    }
    swap(&arr[i + 1], &arr[h]);
    return (i + 1);
}

/* arr[] --> Array to be sorted, l --> Starting index, h --> Ending index */
void quickSort(int arr[], int l, int h)
{
    if (l < h)
    {
        int p = partition(arr, l, h); /* Partitioning index */
        quickSort(arr, l, p - 1);
        quickSort(arr, p + 1, h);
    }
}

/* Function to print some array */
void printArray(int arr[], int size)
{
    int i;
    for (i=0; i < size; i++)
        printf("%d ", arr[i]);
    printf("\n");
}

// Main program to test functions
int main()
{
    
    int arr[] = {102, 54, 34, 46, 51, 78};
    int n = sizeof(arr)/sizeof(arr[0]);
    quickSort(arr, 0, n-1);
    printf("The sorted array is ~ : \n");
    printArray(arr, n);
    return 0;
}





