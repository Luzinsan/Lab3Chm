#ifndef EIGENVALUE_H
#define EIGENVALUE_H
#include <iostream>
#include <fstream>
#include "../../../2/Lab2Chm/Lab2Chm/Matrix.h"
#include "../../../2/Lab2Chm/Lab2Chm/Vector.h"

#include "Polynomial.h"
#include "PolStr.h"
#define EPS 1E-3
#define MAX 10

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

        static T getDeterminant(const Matrix<T>& _Matrix) // Через метод Гаусса
        {
            Matrix<T> A(_Matrix);
            T determinant = 1;
            // Прямой ход метода Гаусса - преобразование матрицы к треугольному виду
            for (int i = 0; i < A.getRows(); i++) // проходим по всем строкам
            {
                T coeff = A[i][i]; // запоминаем коэффициент по диагонали
                determinant *= coeff;
                for (int j = i; j < A.getRows(); j++) // проходим по всем элементам текущей строки, включая вектор коэффициентов
                    A[i][j] /= coeff;
                for (int j = i + 1; j < A.getRows(); j++)
                {
                    coeff = A[j][i]; 
                    for (int k = i; k < A.getCols(); k++) 
                        A[j][k] -= coeff * A[i][k]; 
                }
            }
            return determinant;
        }

        
        void MethodDanilevsky() // Получение собственных чисел
        {
            *_fout << "\t\tМетод Данилевского для нахождения собственных чисел.\n";
            Matrix<T> P = GetFrobenius();
            *_fout << "\n\tМатрица Фробениуса:\n"
                << P;
           
            T* array = new T[m+1];
           
            array[m] = (m % 2 == 0) ? 1 : -1;
            for (int i = m-1; i >= 0; i--)
                array[i] = P[0][m-i-1];
            
            Polynomial<T> pol(m+1, array);
            delete[] array;
            

            std::string ss = pol.to_string();
            
            const char* polStr = CreatePolStr(ss.c_str(), 0);
            Vector<T> eigenvalue(m);
            if (GetError() == ERR_OK && polStr)
            {
                int n = 0; //номер корня
                double x0 = 0; // Начальное приближение
                double x1 = x0;
                bool flag = false; // было ли на предыдущей итерации найдено приближение корня в промежутке (-EPS, EPS)
                T root1 = EvalPolStr(polStr, x0, 0);
                T root2 = root1;
                for (x1 = EPS; n < m && x1 < MAX ; x1 += EPS)
                {
                    root1 = root2;
                    root2 = EvalPolStr(polStr, x1, 0);
                    if (abs(root2) <= EPS && abs(root2) < abs(root1))
                    {
                        x0 = x1;
                        flag = true;
                    }
                    else if (flag == true)
                    {
                        int k = 0;
                        flag = false;
                        eigenvalue[n] = x0;
                        k++; n++;

                        std::cout << EvalPolStr(polStr, x0, k) << "\n";
                        double eps_b = EPS;
                        while (abs(EvalPolStr(polStr, x0, k)) < eps_b)
                        {
                            std::cout << abs(EvalPolStr(polStr, x0, k)) << "\n";
                            k++;
                            eps_b *= 10;
                        }
                        for (int i = 0; i < k - 1; i++)
                            eigenvalue[n++] = x0;
                    }
                }
                for (x1 = EPS; n < m && x1 > -MAX; x1 -= EPS)
                {
                    root1 = root2;
                    root2 = EvalPolStr(polStr, x1, 0);
                    if (abs(root2) <= EPS && abs(root2) < abs(root1))
                    {
                        x0 = x1;
                        flag = true;
                    }
                    else if (flag == true)
                    {
                        int k = 0;
                        flag = false;
                        eigenvalue[n] = x0;
                        k++; n++;

                        std::cout << EvalPolStr(polStr, x0, k) << "\n";
                        double eps_b = EPS;
                        while (abs(EvalPolStr(polStr, x0, k)) < eps_b)
                        {
                            std::cout << abs(EvalPolStr(polStr, x0, k)) << "\n";
                            k++;
                            eps_b *= 10;
                        }
                        for (int i = 0; i < k-1; i++)
                            eigenvalue[n++] = x0;
                    }
                }
                
                n = 1;
                int k;
                int i;
                Matrix<T> E(m, unit_matrix_initer);
                for (i = 1; i < m; i++) 
                {
                    k = 1; 
                    if (eigenvalue[i] == eigenvalue[i-1])
                        k++;
                    else
                    {
                        std::cout << "\n\tСобственное число #" << n++ << ": " << eigenvalue[i-1]
                            << "\n\t-> Кратность: " << k;
                        std::cout << "\n\tПроверка: \n" << (*A) << "\t-\n" << E << "\n\t*" << eigenvalue[i - 1] << "\n\t";
                        std::cout << "\n\tПроверка: " << getDeterminant(((*A) - E * eigenvalue[i - 1]));
                        k = 1;
                    }
                }
                std::cout << "\n\n\tСобственное число #" << n << ": " << eigenvalue[i - 1]
                    << "\n\t-> Кратность: " << k;
                std::cout << "\n\tПроверка: \n" << (*A) << "\t-\n" << E << "\n\t*" << eigenvalue[i - 1] << "\n\t";
                std::cout << "\n\tПроверка: " << getDeterminant(((*A) - E * eigenvalue[i - 1]));
                std::cout << "\n" << std::setw(5) << eigenvalue;
            }
        }

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
