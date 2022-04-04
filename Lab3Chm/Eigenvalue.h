#ifndef EIGENVALUE_H
#define EIGENVALUE_H
#include <iostream>
#include "../../../2/Lab2Chm/Lab2Chm/Matrix.h"
#include "../../../2/Lab2Chm/Lab2Chm/Vector.h"

template<class T>
class Eigenvalue 
{
public:
    enum class TASK
    {
        EIGENVALUES = 1,
        EIGENVECTORS,
    };
private:
    TASK _task;
    int m; // ����������� ���������� �������
    Matrix<T>* A;
    Matrix<T>* P;         // ������� ����������
    Vector<Vector<T>>* x; // ��������� 
    T _eigenvalue;
    std::ofstream* _fout;
    
public:

    Eigenvalue() 
    {
        std::ifstream* _fin = new std::ifstream("input.txt");
        _fout = new std::ofstream("output.txt");
        int c;
        *_fin >> c; // ����������� ��� ������
        _task = static_cast<TASK>(c);

        *_fin >> m;
        std::cout << "\n\t������� �������: " << m;

        T* array = new T[m * (m + 1)];
        for (int i = 0; i < (m + 1) * m; i++)
            *_fin >> array[i];

        delete _fin;
    }

    ~InputData()
    {

        
    }


};

#endif
