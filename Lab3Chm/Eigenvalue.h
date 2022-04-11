#ifndef EIGENVALUE_H
#define EIGENVALUE_H
#include <iostream>
#include <fstream>
#include "../../../2/Lab2Chm/Lab2Chm/Matrix.h"
#include "../../../2/Lab2Chm/Lab2Chm/Vector.h"

namespace luMath
{
    double unit_matrix_initer(size_t m, size_t n, size_t r, size_t c)
    {
        return r == c;
    }

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
        int                m;           // размерность квадратной матрицы
        Matrix<T>* A;           // исходная матрица
        Matrix<T>* P;           // матрица Фробениуса
        Vector<Vector<T>*>* x;           // собственные вектора 
        Vector<T>* eigenvalues; // собственные числа
        Vector<T>* k;           // кратность собственных чисел/векторов
        std::ofstream* _fout;

    public:

        Eigenvalue()
        {
            std::ifstream _fin("input.txt");
            _fout = new std::ofstream("output.txt");
            int c;
            _fin >> c; // считывается тип задачи
            _task = static_cast<TASK>(c);

            _fin >> m;
            std::cout << "\n\tПорядок матрицы: " << m;

            T* array = new T[m * m];
            for (int i = 0; i < m * m; i++)
                _fin >> array[i];

            _fin.close();

            A = new Matrix<T>(m, array);
            P = new Matrix<T>(m);
           
            //Vector<T>* vectors = new Vector<T>[m];
            //x = new Vector<T>(m, vectors);
            //delete vectors;

            eigenvalues = new Vector<T>(m);
            k = new Vector<T>(m);
        }

        ~Eigenvalue()
        {
            delete A;
            delete P;

            delete x;
            delete eigenvalues;
            delete k;
        }
        
        TASK getTask() { return _task; }

        Matrix<T> GetFrobenius()
        {
            std::cout << "\n" << std::setw(10) << *A;
            Matrix<T> P(*A);
            Matrix<T> tempA(*A);
            Matrix<T> M(m), M_1(m);
            Matrix<T> E(m, unit_matrix_initer);
            for (int k = m - 1; k >= 0; k--)
            {
                // Вычисляем M(k)
                for (int i = 0; i < m; i++)
                    for (int j = 0; j < m; j++)
                    {
                        if (i == (k-1) && j == (k-1)) // m(k,k) = 1 / a^(n-k-1)_(k+1,k)
                            M[i][j] = 1 / tempA[k][k-1];
                        else if (i == (k-1)) // m(k,j) = -a^(n-k-1)_(k+1,j) / a^(n-k-1)_(k+1,k); j = 1,2,...,n; j != k
                            M[i][j] = -tempA[k][j] / tempA[k][k-1];
                        else //m(i,j) = e(i,j); i = 1,2,...,n; j = 1,2,...,n; i!=k
                            M[i][j] = E[i][j];
                        std::cout << "\n"<< std::setw(10) << M << "\n";
                    }
                P *= M; // ~A^k = A^(k-1) * M_(n-k)
                std::cout << "\n~A^" << k-1 << " = A^" << k-2 <<" * M_" << k << "\n" << std::setw(10) << P;
                // Вычисляем M^-1(k)
                for (int i = 0; i < m; i++)
                    for (int j = 0; j < m; j++)
                    {
                        if (i == (k - 1)) // m(k,j) = a^(n-k-1)_(k+1,j); j = 1,2,...,n
                            M_1[i][j] = tempA[k][j];
                        else //m(i,j) = e(i,j); i = 1,2,...,n; j = 1,2,...,n; i!=k
                            M_1[i][j] = E[i][j];
                        std::cout << "\n" << std::setw(10) << M_1 << "\n";
                    }
                P = M_1 * P; // A^k = M^(-1)_(n-k) * ~A^k 
                std::cout << "\nA^" << k - 1 << " = M^-1_" << k << " * ~A^" << m-k << "\n" << std::setw(10) << P;
                tempA = P;
            }
        
            return P;
        }

        
    };
}
#endif
