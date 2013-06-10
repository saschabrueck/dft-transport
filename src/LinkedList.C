#include "LinkedList.H"

LinkedList::LinkedList()
{
    next       = NULL;
    index_j    = INF;
    i_value    = 0.0;
    type       = 0;
    added      = 0;
    left_right = 0;
}

/************************************************************************************************/

void LinkedList::insert(double val,int i,int j,int t,int lr)
{
    r_value    = val;
    index_i    = i;
    index_j    = j;
    type       = t;
    left_right = lr;
    next       = new LinkedList();
}

/************************************************************************************************/

void LinkedList::comp_insert(double r_val,double i_val,int i,int j,int t,int lr)
{
    r_value    = r_val;
    i_value    = i_val;
    index_i    = i;
    index_j    = j;
    type       = t;
    left_right = lr;
    next       = new LinkedList();
}

/************************************************************************************************/

void LinkedList::add_to_list(int i0,int iN,int j0,int jN,int* no_element,int* ind_i,int ind_shift,
                             double rval,double ival,int added_value)
{
    
    LinkedList *node=this, *ptr, *prev;
    int i,j;

    for(i=i0;i<iN;i++){
        while(node->index_i<i){
            prev = node;
            node = node->next;
        }
        for(j=j0;j<jN;j++){
	  
            while((node->index_j<j)&(node->index_i==i)){
                prev = node;
                node = node->next;
	     }
	    
            if((j<node->index_j)|(i<node->index_i)){
                ptr = new LinkedList();
		ptr->r_value  = rval;
                ptr->i_value  = ival;
                ptr->index_i  = i;
                ptr->index_j  = j;
                ptr->added    = added_value;
                prev->next    = ptr;
                ptr->next     = node;
                node          = ptr;
 		
                  	
                ind_i[i-ind_shift]++;
                *no_element   = *no_element+1;
            }else{
                node->r_value = node->r_value+rval;
                node->i_value = node->i_value+ival;
                node->added   = added_value;
            }
        }
    }
}

/************************************************************************************************/
