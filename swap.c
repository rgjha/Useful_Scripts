/*
 * C Program to swap set of two numbers *
*/

#include <stdio.h>
void swap(int *xp, int *yp)
{
    if (xp == yp) // Check if the two addresses are same
    
    return;
    *xp = *xp + *yp;
    *yp = *xp - *yp;
    *xp = *xp - *yp;
}

void swap1(int *xp, int *yp)  // Bitwise XOR operator ^ // 
{
    *xp = *xp ^ *yp;
    *yp = *xp ^ *yp;
    *xp = *xp ^ *yp;
}

int main()
{
    int x = 10;
    int y = 5;
    swap1(&x, &y);
    printf("After swap x is = %d", x);
    putchar('\n');
    fflush(stdout);
    printf("After swap y is = %d", y);
    putchar('\n');
    fflush(stdout);
    return 0;
}



