#include <iostream>
#include "Eigenvalue.h"
using namespace luMath;
int main()
{
    setlocale(LC_ALL, "Rus");
    Eigenvalue<double> data;

    switch (data.getTask())
    {
    case Eigenvalue<double>::TASK::EIGENVALUES:
        data.GetFrobenius();
        break;
    case Eigenvalue<double>::TASK::EIGENVECTORS:
        //data.getDeterminant();
        break;
    }
    
    return 0;
}

