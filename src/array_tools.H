#ifndef array_tools_H
#define array_tools_H

#include "Types.H"
#include <iostream>

template <class T, size_t N> inline
size_t size_of_array(const T(&)[N]) {return N;}

template <class T>
T **allocate_2Darray(int nRows, int nCols)
{
   T **array;

   array = new T*[nRows];
   for( int i = 0 ; i < nRows ; i++ )
      array[i] = new T [nCols];

   return array;
}

template <class T>
void free_2Darray(T** array)
{
   delete [] *array;
   delete [] array;
}

template<class T>
void print_array(const T *array, int size)
{
    for (int i = 0; i < size; i++)
        std::cout<<array[i]<<" ";
    std::cout<<endl;
}

template <class T>
void full_transpose(int ncol, int nrow, T *B, T* A)
{
    for(int i=0;i<ncol;i++) c_tcopy(nrow,&B[i*nrow],1,&A[i],ncol);
}

#endif
