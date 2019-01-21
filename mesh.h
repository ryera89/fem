#ifndef MESH_H
#define MESH_H

#include <QPointF>
#include "ndimmatrix/matrix.h"
#include "element.h"

typedef Matrix<uint32_t,2,Matrix_Type::GEN,Matrix_Storage_Scheme::FULL> MatUint32;
typedef Matrix<double,2,Matrix_Type::GEN,Matrix_Storage_Scheme::FULL> MatDoub;

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

inline uint8_t gauss_integration_points(ELEMENT_TYPE etype){
    switch (etype) {
    case ELEMENT_TYPE::TRI3:
        return 1;
    case ELEMENT_TYPE::QUAD4:
        return 4;
    //default:
        //return 0;
    }
    return 0;
}
#endif // MESH_H
