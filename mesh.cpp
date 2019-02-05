#include "mesh.h"

template <>
rectangular_mesh<ELEMENT_TYPE::QUAD4>::rectangular_mesh(double Lx, double Ly, uint32_t Nx, uint32_t Ny):
    element_number(Nx*Ny),x_lenght(Lx),y_lenght(Ly),nodes_coordinates((Nx+1)*(Ny+1)),element_connect(element_number,4)
{
    double dx = Lx/Nx;
    double dy = Ly/Ny;
    for (uint32_t i = 0; i < Ny; ++i){
        for (uint32_t j = 0 ; j < Nx; ++j){
            uint32_t ii = i*Nx+j;
            element_connect(ii,0) = ii+i;
            element_connect(ii,1) = ii+i+1;
            element_connect(ii,2) = ii+i+Nx+1;
            element_connect(ii,3) = ii+i+Nx+2;

        }
    }
    for (uint32_t i = 0; i < Ny+1; ++i){
        for (uint32_t j = 0 ; j < Nx+1; ++j){
            uint32_t ii = i*(Nx+1) + j;
            nodes_coordinates[ii].setX(j*dx);
            nodes_coordinates[ii].setY(i*dy);
        }
    }
}

template <>
rectangular_mesh<ELEMENT_TYPE::TRI3>::rectangular_mesh(double Lx, double Ly, uint32_t Nx, uint32_t Ny,
                                                       bool diag):element_number(2*Nx*Ny),
                                                                    x_lenght(Lx),y_lenght(Ly),
                                                                    nodes_coordinates((Nx+1)*(Ny+1)),
                                                                    element_connect(element_number,3)
{
    double dx = Lx/Nx;
    double dy = Ly/Ny;
    uint32_t tmp = 2*Nx;

    //std::vector<element<ELEMENT_TYPE::TRI,3>> elements(nelem);
    //elements.reserve(nelem);

    if (diag){
    /*  __ __ __ __
         * | /| /| /| /|
         * |/_|/_|/_|/_|
         * | /| /| /| /|
         * |/_|/_|/_|/_|
         * | /| /| /| /|
         * |/_|/_|/_|/_|
         * | /| /| /| /|
         * |/_|/_|/_|/_|
        */
    for (uint32_t i = 0; i < Ny; ++i){
        for (uint32_t j = 0 ; j < Nx; ++j){
            uint32_t ii = i*tmp+2*j;
            uint32_t iii = ii+1;
            element_connect(ii,0) = i*Nx+j+i;
            element_connect(ii,1) = (i+1)*Nx+j+i+1;
            element_connect(ii,2) = (i+1)*Nx+j+i+2;

            element_connect(iii,0) = i*Nx+j+i;
            element_connect(iii,1) = i*Nx+j+i+1;
            element_connect(iii,2) = (i+1)*Nx+j+i+2;
        }
    }
    }else{
        /*  __ __ __ __
         * |\ |\ |\ |\ |
         * |_\|_\|_\|_\|
         * |\ |\ |\ |\ |
         * |_\|_\|_\|_\|
         * |\ |\ |\ |\ |
         * |_\|_\|_\|_\|
         * |\ |\ |\ |\ |
         * |_\|_\|_\|_\|
        */
        for (uint32_t i = 0; i < Ny; ++i){
            for (uint32_t j = 0 ; j < Nx; ++j){
                uint32_t ii = i*tmp+2*j;
                uint32_t iii = ii+1;
                element_connect(ii,0) = i*Nx+j+i;
                element_connect(ii,1) = i*Nx+j+i+1;
                element_connect(ii,2) = (i+1)*Nx+j+i+1;

                element_connect(iii,0) = i*Nx+j+i+1;
                element_connect(iii,1) = (i+1)*Nx+j+i+1;
                element_connect(iii,2) = (i+1)*Nx+j+i+2;

            }
        }
    }
    for (uint32_t i = 0; i < Ny+1; ++i){
        for (uint32_t j = 0 ; j < Nx+1; ++j){
            uint32_t ii = i*(Nx+1) + j;
            nodes_coordinates[ii].setX(j*dx);
            nodes_coordinates[ii].setY(i*dy);
        }
    }

}

template<ELEMENT_TYPE etype>
void rectangular_mesh<etype>::split_mesh_nodes_and_dof()
{
    auto [ixmin,ixmax] = std::minmax_element(nodes_coordinates.begin(),nodes_coordinates.end(),
            [](const QPointF &v1,const QPointF &v2){return v1.x() < v2.x();});
    auto [iymin,iymax] = std::minmax_element(nodes_coordinates.begin(),nodes_coordinates.end(),
            [](const QPointF &v1,const QPointF &v2){return v1.y() < v2.y();});

    double xmin = ixmin->x();
    double xmax =  ixmax->x();
    double ymin = *iymin->y();
    double ymax = *iymax->y();
}

//template<ELEMENT_TYPE etype>
//rectangular_mesh<etype>::rectangular_mesh(double Lx, double Ly, uint32_t Nx, uint32_t Ny,
//                                          std::enable_if_t<ELEMENT_TYPE::TRI3 == etype, bool> diag)
//{

//}
