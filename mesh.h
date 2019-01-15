#ifndef MESH_H
#define MESH_H

#include <QPointF>
#include "ndimmatrix/matrix.h"
#include "element.h"

typedef Matrix<uint32_t,2,Matrix_Type::GEN,Matrix_Storage_Scheme::FULL> Matrix_uint32;

//enum class ELEMENT_TYPE{TRI3,QUAD4};

template<ELEMENT_TYPE etype>
struct rectangular_mesh
{
    static constexpr ELEMENT_TYPE element_type = etype;

    uint32_t element_number;
    double x_lenght;
    double y_lenght;
    std::vector<QPointF> nodes_coordinates;
    Matrix_uint32 element_connect;

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

#endif // MESH_H
