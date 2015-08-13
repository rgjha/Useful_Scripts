#include <iostream> // Basic I/O functions
#include <fstream>
using namespace std;

int main()
{
    int array_size = 100000; // define the size of character array
    char * array = new char[array_size]; // allocating an array of 1kb
    int position = 0; //Increment to fill characters in the array
    
    ifstream fin("example.dat"); //Opening an input stream for example.dat//

    if(fin.is_open())
    {
        //File opened sucessfully//
        //This loop run until end of file (eof) does not occur
        while(!fin.eof() && position < array_size)
        {
            fin.get(array[position]); // Reading file //
            position++;
        }
        array[position-1] = '\0'; //placing character array terminating character
        
        cout << "Array...is...." << endl << endl;
        for(int i = 0; array[i] != '\0'; i++)
        {
            cout << array[i];
        }
    }
    else //File could not be opened
    {
        cout << "File could not be opened." << endl;
    }
    return 0;
}
