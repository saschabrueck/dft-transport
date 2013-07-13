#ifndef vector_tools_H
#define vector_tools_H

#include <vector>
#include <iterator>
#include <iostream>
#include "vectorTypes.H"

template <class T>
void print_Vector1D(Vector1D<T>& Vec)
{
   typename std::vector<T>::const_iterator it;
   for (it=Vec.begin(); it!=Vec.end(); ++it)
       std::cout << " " << *it ;
   std::cout << endl;
}

template <class T>
void print_Vector2D(Vector2D<T>& Vec)
{
   typename std::vector< std::vector<T> >::const_iterator it_outer;
   for (it_outer=Vec.begin(); it_outer!=Vec.end(); ++it_outer){
       typename std::vector<T>::const_iterator it_inner=it_outer->begin(),end_inner=it_outer->end();
       std::copy(it_inner,end_inner,std::ostream_iterator<T>(std::cout, " "));
       std::cout << endl;
   };
}

template <class T>
void initialize_Vector2D(Vector2D<T>& Vec, int size)
{
   for(int i = 0; i<size; i++){
      vector<T> empty_vec;
      Vec.push_back(empty_vec);
   };
}

template <class T>
void convert_Vector2D_to_Vector1D(Vector2D<T>& Vec2D, Vector1D<T>& Vec1D)
{
   size_t reserve_size = 0;
   for(unsigned i=0; i<Vec2D.size(); ++i)
      reserve_size += Vec2D[i].size();

   Vec1D.reserve(reserve_size);

   for(unsigned i=0; i<Vec2D.size(); ++i)
      Vec1D.insert(Vec1D.end(), Vec2D[i].begin(), Vec2D[i].end());
}

template <class T>
void Vector1D_abs(Vector1D<T>& Vec1D)
{
   typename std::vector<T>::iterator i;
   for(i = Vec1D.begin(); i != Vec1D.end(); ++i)
      if(*i < 0) *i *= -1;
}

#endif
