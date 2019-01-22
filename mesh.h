#ifndef MESH_H
#define MESH_H

#include <QPointF>
#include "ndimmatrix/matrix.h"
#include "element.h"

typedef Matrix<uint32_t,2,Matrix_Type::GEN,Matrix_Storage_Scheme::FULL> MatUint32;
typedef Matrix<double,2,Matrix_Type::GEN,Matrix_Storage_Scheme::FULL> MatDoub;
typedef Matrix<double,1,Matrix_Type::GEN,Matrix_Storage_Scheme::FULL> VecDoub;

//enum class ELEMENT_TYPE{TRI3,QUAD4};

template<ELEMENT_TYPE etype>
struct rectangular_mesh
{
    static constexpr ELEMENT_TYPE element_type = etype;

    uint32_t element_number;
    double x_lenght;
    double y_lenght;
    std::vector<QPointF> nodes_coordinates;
    MatUint32 element_connect;

    //rectangular_mesh(double Lx,double Ly,uint32_t Nx,uint32_t Ny,std::enable_if_t<ELEMENT_TYPE::TRI3 == etype,bool> diag = true);
    rectangular_mesh(double Lx,double Ly,uint32_t Nx,uint32_t Ny);
    rectangular_mesh(double Lx,double Ly,uint32_t Nx,uint32_t Ny,bool diag);

    /*
     * Lx: longitud en direccion x;
     * Ly: longitud en direccion y
     * Nx: numero de elementos rectangulares en direccion x
     * Ny: numero de elementos rectangulares en direccion y
     * diag: en caso de elementos triangulares define la diagonal, true = '/', false = '\'
    */
};
extern template struct rectangular_mesh<ELEMENT_TYPE::QUAD4>;
extern template struct rectangular_mesh<ELEMENT_TYPE::TRI3>;


//TODO: cambiar esta funcion de header
//dN: matriz dN(dim,nodxelem): derivadas de funciones de forma elemento master
//X: vector de coordenadas por elemento X(nelem,COORD), COORD(nodxelem,dim): coordenadas de los nodos por elemento
inline std::vector<MatDoub> transpose_jacobian(const MatDoub &dN,const std::vector<MatDoub> &X){
    //size_t dim = dN.rows();
    //size_t nodxelem = dN.cols();
    size_t nelem = X.size();
    //std::vector<MatDoub> v_JT(nelem,MatDoub(dim,dim));
    std::vector<MatDoub> v_JT(nelem);
    std::transform(X.begin(),X.end(),v_JT.begin(),[&dN](const auto &x){return dN*x;});
    return v_JT;
}

VecDoub quad4_shape_functions(const QPointF &p){
    return {0.25*(1 - p.x() - p.y() + p.x()*p.y()),0.25*(1 + p.x() - p.y() - p.x()*p.y()),
            0.25*(1 + p.x() + p.y() + p.x()*p.y()),0.25*(1 - p.x() + p.y() - p.x()*p.y())};
}
MatDoub quad4_shape_functions_gradient(const QPointF &p){
    return {{0.25*(-1+p.y()),0.25*(1-p.y()),0.25*(1+p.y()),-0.25*(1+p.y())},
            {0.25*(-1+p.x()),0.25*(-1-p.x()),0.25*(1+p.x()),0.25*(1-p.x())}};
}

struct gaussian_cuadrature{
    ELEMENT_TYPE element_type;
    uint32_t gauss_points_number;
    std::vector<QPointF> gauss_points;
    std::vector<double> weights;

    gaussian_cuadrature(ELEMENT_TYPE etype):element_type(etype){
        switch (etype) {
        case ELEMENT_TYPE::QUAD4:{
            gauss_points_number = 4;
            std::vector<QPointF> gauss_points = {QPointF(-1/std::sqrt(3),-1/std::sqrt(3)),
                                                 QPointF(1/std::sqrt(3),-1/std::sqrt(3)),
                                                 QPointF(1/std::sqrt(3),1/std::sqrt(3)),
                                                 QPointF(-1/std::sqrt(3),1/std::sqrt(3))};
            weights = {1.0,1.0,1.0,1.0};
            break;
        }
        default:
            throw ("gaussian_cuadrature: error invalid element type");
            //break;
        }
    }
    VecDoub integration_points_function_evaluation(double (*f)(const QPointF &x)){
        VecDoub vres(gauss_points_number);
        std::transform(gauss_points.begin(),gauss_points.end(),weights.begin(),
                       vres.begin(),[&f](auto &p,auto &w){return w*f(p);});
        return vres;
    }
    double integrate(double (*f)(const QPointF &x)){
        return std::inner_product(gauss_points.begin(),gauss_points.end(),weights.begin(),0.0,std::plus<double>(),
                                  [&f](auto &p,auto &w){return w*f(p);});
    }

};

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
//inline uint8_t gauss_integration_points(ELEMENT_TYPE etype){
//    switch (etype) {
//    case ELEMENT_TYPE::TRI3:
//        return 1;
//    case ELEMENT_TYPE::QUAD4:
//        return 4;
//    //default:
//        //return 0;
//    }
//    return 0;
//}
//inline std::vector<QPointF> gauss_integration_points(ELEMENT_TYPE etype){
//    switch (etype) {
//       case ELEMENT_TYPE::QUAD4:
//        std::vector<QPointF> gauss_points = {QPointF(-1/std::sqrt(3),-1/std::sqrt(3)),
//                                             QPointF(1/std::sqrt(3),-1/std::sqrt(3)),
//                                             QPointF(1/std::sqrt(3),1/std::sqrt(3)),
//                                             QPointF(-1/std::sqrt(3),1/std::sqrt(3))};

//    }
//    //    switch (etype) {
//    //    case ELEMENT_TYPE::TRI3:
//    //        return 1;
//    //    case ELEMENT_TYPE::QUAD4:
//    //        return 4;
//    //    //default:
//    //        //return 0;
//    //    }
//    //    return 0;
//}
#endif // MESH_H
