#ifndef VECTOR_H
#define VECTOR_H

#include <vector>
#include <complex>

class Vector
{
public:
    Vector();
};

enum class MATRIX_TYPE{DENSE,SYMM,HER,UTRI,LTRI,CSR,CSR3};

template<typename T,MATRIX_TYPE type,class Enable = void>
class Matrix{
private:
    size_t row;
    size_t col;
    std::vector<T> elem;
public:
    Matrix() = default;
    Matrix(size_t m,size_t n):row(m),col(n),elem(m*n,T()){}
};

template<typename T>
class Matrix<T,MATRIX_TYPE::HER,std::enable_if_t<std::is_same_v<T,std::complex>>>{
private:
    size_t row;
    size_t col;
    std::vector<T> elem;
public:
    Matrix() = default;
    Matrix(size_t m,size_t n):row(m),col(n),elem(m*n,T()){}
};

Matrix<double,MATRIX_TYPE::HER> mat1(4,5);

#endif // VECTOR_H
