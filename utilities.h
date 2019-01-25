#ifndef UTILITIES_H
#define UTILITIES_H

#include "ndimmatrix/matrix.h"

inline std::tuple<MatDoub,double> inv_det_2x2_Matrix(const MatDoub &A){
    assert(A.rows() == 2 && A.cols() == 2);
    const double a = A(0,0);
    const double b = A(0,1);
    const double c = A(1,0);
    const double d = A(1,1);
    MatDoub K(2,2);
    K(0,0) = d; K(0,1) = -b;
    K(1,0) = -c; K(1,1) = a;
    double det = (a*d-b*c);
    K/=det;
    return std::tuple<MatDoub,double>(K,det);
}
inline std::tuple<MatDoub,double> inv_det_3x3_Matrix(const MatDoub &A){
    assert(A.rows() == 3 && A.cols() == 3);
    double a = A(0,0);double b = A(0,1);double c = A(0,2);
    double d = A(1,0);double e = A(1,1);double f = A(1,2);
    double g = A(2,0);double h = A(2,1);double i = A(2,2);

    double aa = (e*i - f*h);
    double bb = -(d*i - f*g);
    double cc = (d*h - e*g);

    double dd = -(b*i - c*h);
    double ee = (a*i - c*g);
    double ff = -(a*h - b*g);

    double gg = (b*f - c*e);
    double hh = -(a*f - c*d);
    double ii = (a*e - b*d);

    double det = a*aa + b*bb + c*cc;

    MatDoub R(3,3);

    R(0,0) = aa/det; R(0,1) = dd/det; R(0,2) = gg/det;
    R(1,0) = bb/det; R(1,1) = ee/det; R(1,2) = hh/det;
    R(2,0) = cc/det; R(2,1) = ff/det; R(2,2) = ii/det;

    return std::tuple<MatDoub,double>(R,det);
}
inline double det_3x3_Matrix(const MatDoub &A){
    assert(A.rows() == 3 && A.cols() == 3);
    double a = A(0,0);double b = A(0,1);double c = A(0,2);
    double d = A(1,0);double e = A(1,1);double f = A(1,2);
    double g = A(2,0);double h = A(2,1);double i = A(2,2);

    double aa = (e*i - f*h);
    double bb = -(d*i - f*g);
    double cc = (d*h - e*g);

    return a*aa + b*bb + c*cc;
}
inline double det_2x2_Matrix(const MatDoub &A){
    assert(A.rows() == 2 && A.cols() == 2);
    return A(0,0)*A(1,1) - A(1,0)*A(0,1);
}

#endif // UTILITIES_H
