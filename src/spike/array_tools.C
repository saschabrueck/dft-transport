#include "array_tools.H"

void full_conjugate_transpose(int n1, int n2, CPX *in, CPX* out)
{
    for (int j=0;j<n2;j++)
        for (int i=0;i<n1;i++)
            out[j*n1+i]=conj(in[i*n2+j]);
}
