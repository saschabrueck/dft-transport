#ifndef vectorTypes_H
#define vectorTypes_H

#include <vector>

template <class T>
struct Vector1D : public std::vector<T> {};
template<class T>
struct Vector2D : public std::vector< std::vector<T> > {};

#endif
