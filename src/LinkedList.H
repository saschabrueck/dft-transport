#ifndef __LINKEDLIST
#define __LINKEDLIST

#include <mpi.h>
#include <iostream>
#include "CSR.H"

class LinkedList {
    
  public:    
    /*Constructors*/
    LinkedList();
    
    /*Operations*/
    void insert(double,int,int,int,int);
    void comp_insert(double,double,int,int,int,int);
    void add_to_list(int,int,int,int,int*,int*,int,double,double,int);
    
    LinkedList* next;
    double r_value,i_value;
    int index_i,index_j,type,added,left_right;

    private:
    
};

#endif
