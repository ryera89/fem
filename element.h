#ifndef ELEMENT_H
#define ELEMENT_H

#include "node.h"
#include "ndimmatrix/matrix.h"
#include <array>


typedef Matrix<double,1> VecDoub;
typedef Matrix<double,2> MatDoub;

enum class ELEMENT_TYPE{TRI3,QUAD4};

//funciones de forma del elemento master [QUAD4]
inline VecDoub quad4_master_element_shape_functions(const QPointF &p){
    return {0.25*(1 - p.x() - p.y() + p.x()*p.y()),0.25*(1 + p.x() - p.y() - p.x()*p.y()),
            0.25*(1 + p.x() + p.y() + p.x()*p.y()),0.25*(1 - p.x() + p.y() - p.x()*p.y())};
}
//Derivadas de las funciones de forma del elemento master [QUAD4]
inline MatDoub quad4_master_element_shape_functions_gradient(const QPointF &p){
    return {{0.25*(-1+p.y()),0.25*(1-p.y()),0.25*(1+p.y()),-0.25*(1+p.y())}, //Ni,x
            {0.25*(-1+p.x()),0.25*(-1-p.x()),0.25*(1+p.x()),0.25*(1-p.x())}}; //Ni,y
}

//TODO comentariar esto que no se va a usar
template<ELEMENT_TYPE etype,uint8_t nnum>
struct element{
    static constexpr ELEMENT_TYPE type = etype;
    static constexpr uint8_t nod_number = nnum;
    std::array<node,nnum> nodes;
};


#endif // ELEMENT_H
