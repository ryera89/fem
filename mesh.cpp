#include "mesh.h"

using std::abs;
using std::isless;
using std::isgreater;

template <>
rectangular_mesh<ELEMENT_TYPE::QUAD4>::rectangular_mesh(double Lx, double Ly, uint32_t Nx, uint32_t Ny):
    element_number(Nx*Ny),nodes_number((Nx+1)*(Ny+1)),x_lenght(Lx),y_lenght(Ly),nodes_coordinates(nodes_number),
    element_connect(element_number,4)
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
                                                       bool diag):element_number(2*Nx*Ny),nodes_number((Nx+1)*(Ny+1)),
                                                                    x_lenght(Lx),y_lenght(Ly),
                                                                    nodes_coordinates(nodes_number),
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
    double xmax = ixmax->x();
    double ymin = iymin->y();
    double ymax = iymax->y();

    double xtoler = 1.0e-8*x_lenght;
    double ytoler = 1.0e-8*y_lenght;

    std::vector<uint32_t> left_edge_nodes;
    std::vector<uint32_t> right_edge_nodes;
    std::vector<uint32_t> bottom_edge_nodes;
    std::vector<uint32_t> top_edge_nodes;
    for (uint32_t i = 0; i < nodes_number; ++i){
        if (is_corner_node(nodes_coordinates[i],xtoler,ytoler,xmin,xmax,ymin,ymax)){
            m_corner_nodes.push_back(i);
            continue;
        }
        if (is_left_edge_node(nodes_coordinates[i],xtoler,ytoler,xmin,xmax,ymin,ymax)){
            left_edge_nodes.push_back(i);
            continue;
        }
        if (is_right_edge_node(nodes_coordinates[i],xtoler,ytoler,xmin,xmax,ymin,ymax)){
            right_edge_nodes.push_back(i);
            continue;
        }
        if (is_top_edge_node(nodes_coordinates[i],xtoler,ytoler,xmin,xmax,ymin,ymax)){
            top_edge_nodes.push_back(i);
            continue;
        }
        if (is_bottom_edge_node(nodes_coordinates[i],xtoler,ytoler,xmin,xmax,ymin,ymax)){
            bottom_edge_nodes.push_back(i);
            continue;
        }
        //interior node
        m_interior_nodes.push_back(i);
    }
    if (top_edge_nodes.size() != bottom_edge_nodes.size()) throw ("split_mesh_nodes_and_dof: error number "
        "of nodes of top and bottom edges do not coincide");
    if (left_edge_nodes.size() != right_edge_nodes.size()) throw ("split_mesh_nodes_and_dof: error number "
        "of nodes of left and right edges do not coincide");
    if (m_corner_nodes.size() != 4) throw ("split_mesh_nodes_and_dof: error there must be 4 corner nodes");

     //chequeo de correspondencia uno-a-uno
    std::vector<double> top_edge_node_xcoor(top_edge_nodes.size());
    std::vector<double> bot_edge_node_xcoor(bottom_edge_nodes.size());

    std::vector<double> left_edge_node_ycoor(left_edge_nodes.size());
    std::vector<double> right_edge_node_ycoor(right_edge_nodes.size());

    for (uint32_t i = 0; i < top_edge_nodes.size(); ++i) top_edge_node_xcoor[i] = nodes_coordinates[top_edge_nodes[i]].x();
    for (uint32_t i = 0; i < bottom_edge_nodes.size(); ++i) bot_edge_node_xcoor[i] = nodes_coordinates[bottom_edge_nodes[i]].x();
    for (uint32_t i = 0; i < left_edge_nodes.size(); ++i) left_edge_node_ycoor[i] = nodes_coordinates[left_edge_nodes[i]].y();
    for (uint32_t i = 0; i < right_edge_nodes.size(); ++i) right_edge_node_ycoor[i] = nodes_coordinates[right_edge_nodes[i]].y();

    std::sort(top_edge_node_xcoor.begin(),top_edge_node_xcoor.end());
    std::sort(bot_edge_node_xcoor.begin(),bot_edge_node_xcoor.end());
    std::sort(left_edge_node_ycoor.begin(),left_edge_node_ycoor.end());
    std::sort(right_edge_node_ycoor.begin(),right_edge_node_ycoor.end());

    for (uint32_t i = 0; i < top_edge_node_xcoor.size(); ++i){
        if (isgreater(abs(top_edge_node_xcoor[i]-bot_edge_node_xcoor[i]),xtoler))
            throw ("split_mesh_nodes_and_dof: error there is not one-to-one correspondence"
            " between top and bottom edge nodes");
    }
    for (uint32_t i = 0; i < left_edge_node_ycoor.size(); ++i){
        if (isgreater(abs(left_edge_node_ycoor[i]-right_edge_node_ycoor[i]),ytoler))
            throw ("split_mesh_nodes_and_dof: error there is not one-to-one correspondence"
            " between left and right edge nodes");
    }

    //ordering border nodes
    std::sort(top_edge_nodes.begin(),top_edge_nodes.end(),[this](uint32_t ii,uint32_t jj)
                                 {return nodes_coordinates[ii].x() < nodes_coordinates[jj].x();});
    std::sort(bottom_edge_nodes.begin(),bottom_edge_nodes.end(),[this](uint32_t ii,uint32_t jj)
                                 {return nodes_coordinates[ii].x() < nodes_coordinates[jj].x();});
    std::sort(left_edge_nodes.begin(),left_edge_nodes.end(),[this](uint32_t ii,uint32_t jj)
                                 {return nodes_coordinates[ii].y() < nodes_coordinates[jj].y();});
    std::sort(right_edge_nodes.begin(),right_edge_nodes.end(),[this](uint32_t ii,uint32_t jj)
                                 {return nodes_coordinates[ii].y() < nodes_coordinates[jj].y();});

    m_right_top_nodes


}

template<ELEMENT_TYPE etype>
bool rectangular_mesh<etype>::is_corner_node(const QPointF &coor,double xtoler,double ytoler,
                                             double xmin,double xmax,double ymin,double ymax)
{
    bool bot_left_corner = isless(abs(xmin - coor.x()),xtoler) && isless(abs(ymin - coor.y()),ytoler);
    bool top_left_corner = isless(abs(xmin - coor.x()),xtoler) && isless(abs(ymax - coor.y()),ytoler);
    bool bot_right_corner = isless(abs(xmax - coor.x()),xtoler) && isless(abs(ymin - coor.y()),ytoler);
    bool top_right_corner = isless(abs(xmax - coor.x()),xtoler) && isless(abs(ymax - coor.y()),ytoler);
    return bot_left_corner || top_left_corner || bot_right_corner || top_right_corner;
}

template<ELEMENT_TYPE etype>
bool rectangular_mesh<etype>::is_right_edge_node(const QPointF &coor, double xtoler, double ytoler, double xmin,
                                            double xmax, double ymin, double ymax)
{
    return isless(abs(xmax - coor.x()),xtoler) && !is_corner_node(coor,xtoler,ytoler,xmin,xmax,ymin,ymax);
}

template<ELEMENT_TYPE etype>
bool rectangular_mesh<etype>::is_top_edge_node(const QPointF &coor, double xtoler, double ytoler, double xmin, double xmax, double ymin, double ymax)
{
    return isless(abs(ymax - coor.y()),ytoler) && !is_corner_node(coor,xtoler,ytoler,xmin,xmax,ymin,ymax);
}

template<ELEMENT_TYPE etype>
bool rectangular_mesh<etype>::is_bottom_edge_node(const QPointF &coor, double xtoler, double ytoler, double xmin, double xmax, double ymin, double ymax)
{
    return isless(abs(ymin - coor.y()),ytoler) && !is_corner_node(coor,xtoler,ytoler,xmin,xmax,ymin,ymax);
}

template<ELEMENT_TYPE etype>
bool rectangular_mesh<etype>::is_left_edge_node(const QPointF &coor, double xtoler, double ytoler, double xmin, double xmax, double ymin, double ymax)
{
    return isless(abs(xmin - coor.x()),xtoler) && !is_corner_node(coor,xtoler,ytoler,xmin,xmax,ymin,ymax);
}

//template<ELEMENT_TYPE etype>
//rectangular_mesh<etype>::rectangular_mesh(double Lx, double Ly, uint32_t Nx, uint32_t Ny,
//                                          std::enable_if_t<ELEMENT_TYPE::TRI3 == etype, bool> diag)
//{

//}
