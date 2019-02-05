#ifndef MESH_H
#define MESH_H

#include <QPointF>
#include "element.h"

//TODO programar bien esta structura
template<ELEMENT_TYPE etype>
struct rectangular_mesh
{
    static constexpr ELEMENT_TYPE element_type = etype;

    uint32_t dofxnode; //grados de libertad por elemento
    uint32_t element_number;
    double x_lenght;
    double y_lenght;
    std::vector<QPointF> nodes_coordinates;
    Matrix<uint32_t,2> element_connect;

    std::vector<uint32_t> m_interior_nodes;
    std::vector<uint32_t> m_left_bottom_nodes;
    std::vector<uint32_t> m_right_top_nodes;
    std::vector<uint32_t> m_corner_nodes;

    std::vector<uint32_t> m_interior_dof;
    std::vector<uint32_t> m_left_bottom_dof;
    std::vector<uint32_t> m_right_top_dof;
    std::vector<uint32_t> m_corner_dof;

    //rectangular_mesh(double Lx,double Ly,uint32_t Nx,uint32_t Ny,std::enable_if_t<ELEMENT_TYPE::TRI3 == etype,bool> diag = true);
    rectangular_mesh(double Lx,double Ly,uint32_t Nx,uint32_t Ny);
    rectangular_mesh(double Lx,double Ly,uint32_t Nx,uint32_t Ny,bool diag);

    void split_mesh_nodes_and_dof();

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
