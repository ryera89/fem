#ifndef MESH_H
#define MESH_H

#include <QPointF>
#include "element.h"

using std::abs;
using std::isless;
using std::isgreater;

//TODO programar bien esta structura
template<ELEMENT_TYPE etype>
struct rectangular_mesh
{
    static constexpr ELEMENT_TYPE element_type = etype;

    uint32_t m_dofxnode; //grados de libertad por elemento
    uint32_t m_element_number;
    uint32_t m_nodes_number;
    double m_x_lenght;
    double m_y_lenght;
    std::vector<QPointF> m_nodes_coordinates;
    Matrix<uint32_t,2> m_element_connect;

    std::vector<uint32_t> m_interior_nodes;
    std::vector<uint32_t> m_left_bottom_nodes;
    std::vector<uint32_t> m_right_top_nodes;
    std::vector<uint32_t> m_corner_nodes;

    std::vector<uint32_t> m_interior_dof;
    std::vector<uint32_t> m_left_bottom_dof;
    std::vector<uint32_t> m_right_top_dof;
    std::vector<uint32_t> m_corner_dof;

    //rectangular_mesh(double Lx,double Ly,uint32_t Nx,uint32_t Ny,std::enable_if_t<ELEMENT_TYPE::TRI3 == etype,bool> diag = true);
    //template<typename = typename std::enable_if_t<etype == ELEMENT_TYPE::QUAD4>>
    rectangular_mesh(uint32_t dofxnode,double Lx,double Ly,uint32_t Nx,uint32_t Ny):m_dofxnode(dofxnode), m_element_number(Nx*Ny),
                                                                                    m_nodes_number((Nx+1)*(Ny+1)),m_x_lenght(Lx),
                                                                                    m_y_lenght(Ly),m_nodes_coordinates(m_nodes_number),
                                                                                    m_element_connect(m_element_number,4)
    {
        double dx = Lx/Nx;
        double dy = Ly/Ny;
        for (uint32_t i = 0; i < Ny; ++i){
            for (uint32_t j = 0 ; j < Nx; ++j){
                uint32_t ii = i*Nx+j;
                m_element_connect(ii,0) = ii+i;
                m_element_connect(ii,1) = ii+i+1;
                m_element_connect(ii,2) = ii+i+Nx+1;
                m_element_connect(ii,3) = ii+i+Nx+2;

            }
        }
        for (uint32_t i = 0; i < Ny+1; ++i){
            for (uint32_t j = 0 ; j < Nx+1; ++j){
                uint32_t ii = i*(Nx+1) + j;
                m_nodes_coordinates[ii].setX(j*dx);
                m_nodes_coordinates[ii].setY(i*dy);
            }
        }

        //split_mesh_nodes_and_dof();
    }
    //template<typename = typename std::enable_if_t<etype == ELEMENT_TYPE::TRI3>>
    rectangular_mesh(uint32_t dofxnode,double Lx,double Ly,uint32_t Nx,uint32_t Ny,bool diag):m_dofxnode(dofxnode),m_element_number(2*Nx*Ny),
                                                                                              m_nodes_number((Nx+1)*(Ny+1)),m_x_lenght(Lx),
                                                                                              m_y_lenght(Ly),m_nodes_coordinates(m_nodes_number),
                                                                                              m_element_connect(m_element_number,3)
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
                    m_element_connect(ii,0) = i*Nx+j+i;
                    m_element_connect(ii,1) = (i+1)*Nx+j+i+1;
                    m_element_connect(ii,2) = (i+1)*Nx+j+i+2;

                    m_element_connect(iii,0) = i*Nx+j+i;
                    m_element_connect(iii,1) = i*Nx+j+i+1;
                    m_element_connect(iii,2) = (i+1)*Nx+j+i+2;
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
                    m_element_connect(ii,0) = i*Nx+j+i;
                    m_element_connect(ii,1) = i*Nx+j+i+1;
                    m_element_connect(ii,2) = (i+1)*Nx+j+i+1;

                    m_element_connect(iii,0) = i*Nx+j+i+1;
                    m_element_connect(iii,1) = (i+1)*Nx+j+i+1;
                    m_element_connect(iii,2) = (i+1)*Nx+j+i+2;

                }
            }
        }
        for (uint32_t i = 0; i < Ny+1; ++i){
            for (uint32_t j = 0 ; j < Nx+1; ++j){
                uint32_t ii = i*(Nx+1) + j;
                m_nodes_coordinates[ii].setX(j*dx);
                m_nodes_coordinates[ii].setY(i*dy);
            }
        }
        //split_mesh_nodes_and_dof();
    }

    void split_mesh_nodes_and_dof()
    {
        auto [ixmin,ixmax] = std::minmax_element(m_nodes_coordinates.begin(),m_nodes_coordinates.end(),
                                                  [](const QPointF &v1,const QPointF &v2){return v1.x() < v2.x();});
        auto [iymin,iymax] = std::minmax_element(m_nodes_coordinates.begin(),m_nodes_coordinates.end(),
                                                  [](const QPointF &v1,const QPointF &v2){return v1.y() < v2.y();});


        double xmin = (*ixmin).x();
        double xmax = (*ixmax).x();
        double ymin = (*iymin).y();
        double ymax = (*iymax).y();

        double xtoler = 1.0e-8*m_x_lenght;
        double ytoler = 1.0e-8*m_y_lenght;

        std::vector<uint32_t> left_edge_nodes;
        std::vector<uint32_t> right_edge_nodes;
        std::vector<uint32_t> bottom_edge_nodes;
        std::vector<uint32_t> top_edge_nodes;
        for (uint32_t i = 0; i < m_nodes_number; ++i){
            if (is_corner_node(m_nodes_coordinates[i],xtoler,ytoler,xmin,xmax,ymin,ymax)){
                m_corner_nodes.push_back(i);
                continue;
            }
            if (is_left_edge_node(m_nodes_coordinates[i],xtoler,ytoler,xmin,xmax,ymin,ymax)){
                left_edge_nodes.push_back(i);
                continue;
            }
            if (is_right_edge_node(m_nodes_coordinates[i],xtoler,ytoler,xmin,xmax,ymin,ymax)){
                right_edge_nodes.push_back(i);
                continue;
            }
            if (is_top_edge_node(m_nodes_coordinates[i],xtoler,ytoler,xmin,xmax,ymin,ymax)){
                top_edge_nodes.push_back(i);
                continue;
            }
            if (is_bottom_edge_node(m_nodes_coordinates[i],xtoler,ytoler,xmin,xmax,ymin,ymax)){
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

        for (uint32_t i = 0; i < top_edge_nodes.size(); ++i) top_edge_node_xcoor[i] = m_nodes_coordinates[top_edge_nodes[i]].x();
        for (uint32_t i = 0; i < bottom_edge_nodes.size(); ++i) bot_edge_node_xcoor[i] = m_nodes_coordinates[bottom_edge_nodes[i]].x();
        for (uint32_t i = 0; i < left_edge_nodes.size(); ++i) left_edge_node_ycoor[i] = m_nodes_coordinates[left_edge_nodes[i]].y();
        for (uint32_t i = 0; i < right_edge_nodes.size(); ++i) right_edge_node_ycoor[i] = m_nodes_coordinates[right_edge_nodes[i]].y();

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
                  {return m_nodes_coordinates[ii].x() < m_nodes_coordinates[jj].x();});
        std::sort(bottom_edge_nodes.begin(),bottom_edge_nodes.end(),[this](uint32_t ii,uint32_t jj)
                  {return m_nodes_coordinates[ii].x() < m_nodes_coordinates[jj].x();});
        std::sort(left_edge_nodes.begin(),left_edge_nodes.end(),[this](uint32_t ii,uint32_t jj)
                  {return m_nodes_coordinates[ii].y() < m_nodes_coordinates[jj].y();});
        std::sort(right_edge_nodes.begin(),right_edge_nodes.end(),[this](uint32_t ii,uint32_t jj)
                  {return m_nodes_coordinates[ii].y() < m_nodes_coordinates[jj].y();});

        m_right_top_nodes.insert(m_right_top_nodes.end(),right_edge_nodes.begin(),right_edge_nodes.end());
        m_right_top_nodes.insert(m_right_top_nodes.end(),top_edge_nodes.begin(),top_edge_nodes.end());

        m_left_bottom_nodes.insert(m_left_bottom_nodes.end(),left_edge_nodes.begin(),left_edge_nodes.end());
        m_left_bottom_nodes.insert(m_left_bottom_nodes.end(),bottom_edge_nodes.begin(),bottom_edge_nodes.end());

        //DOFs dofxelem
        if (m_dofxnode == 1){
            m_right_top_dof = m_right_top_nodes;
            m_left_bottom_dof = m_left_bottom_nodes;
            m_corner_dof = m_corner_nodes;
            m_interior_dof = m_interior_nodes;
        }else{
            if (m_dofxnode == 2){
                m_right_top_dof.clear();
                m_left_bottom_dof.clear();
                m_corner_dof.clear();
                m_interior_dof.clear();
                m_right_top_dof.reserve(2*m_right_top_nodes.size());
                m_left_bottom_dof.reserve(2*m_left_bottom_nodes.size());
                m_corner_dof.reserve(2*m_corner_nodes.size());
                m_interior_dof.reserve(2*m_interior_nodes.size());
                for (uint32_t i = 0; i < m_right_top_nodes.size(); ++i){
                    m_right_top_dof.push_back(2*m_right_top_nodes[i]);
                    m_right_top_dof.push_back(2*m_right_top_nodes[i] + 1);
                }
                for (uint32_t i = 0; i < m_left_bottom_nodes.size(); ++i){
                    m_left_bottom_dof.push_back(2*m_left_bottom_nodes[i]);
                    m_left_bottom_dof.push_back(2*m_left_bottom_nodes[i] + 1);
                }
                for (uint32_t i = 0; i < m_corner_nodes.size(); ++i){
                    m_corner_dof.push_back(2*m_corner_nodes[i]);
                    m_corner_dof.push_back(2*m_corner_nodes[i] + 1);
                }
                for (uint32_t i = 0; i < m_interior_nodes.size(); ++i){
                    m_interior_dof.push_back(2*m_interior_nodes[i]);
                    m_interior_dof.push_back(2*m_interior_nodes[i] + 1);
                }
            }else{
                throw ("split_mesh_nodes_and_dof: error more than 2 dofxelem not implemented");
            }
        }
    }

    bool is_corner_node(const QPointF &coor,double xtoler,double ytoler,double xmin,double xmax,double ymin,double ymax)
    {
        bool bot_left_corner = isless(abs(xmin - coor.x()),xtoler) && isless(abs(ymin - coor.y()),ytoler);
        bool top_left_corner = isless(abs(xmin - coor.x()),xtoler) && isless(abs(ymax - coor.y()),ytoler);
        bool bot_right_corner = isless(abs(xmax - coor.x()),xtoler) && isless(abs(ymin - coor.y()),ytoler);
        bool top_right_corner = isless(abs(xmax - coor.x()),xtoler) && isless(abs(ymax - coor.y()),ytoler);
        return bot_left_corner || top_left_corner || bot_right_corner || top_right_corner;
    }
    bool is_right_edge_node(const QPointF &coor,double xtoler,double ytoler,double xmin,double xmax,double ymin,double ymax)
    {
        return isless(abs(xmax - coor.x()),xtoler) && !is_corner_node(coor,xtoler,ytoler,xmin,xmax,ymin,ymax);
    }
    bool is_top_edge_node(const QPointF &coor,double xtoler,double ytoler,double xmin,double xmax,double ymin,double ymax)
    {
        return isless(abs(ymax - coor.y()),ytoler) && !is_corner_node(coor,xtoler,ytoler,xmin,xmax,ymin,ymax);
    }
    bool is_bottom_edge_node(const QPointF &coor,double xtoler,double ytoler,double xmin,double xmax,double ymin,double ymax)
    {
        return isless(abs(ymin - coor.y()),ytoler) && !is_corner_node(coor,xtoler,ytoler,xmin,xmax,ymin,ymax);
    }
    bool is_left_edge_node(const QPointF &coor,double xtoler,double ytoler,double xmin,double xmax,double ymin,double ymax)
    {
        return isless(abs(xmin - coor.x()),xtoler) && !is_corner_node(coor,xtoler,ytoler,xmin,xmax,ymin,ymax);
    }
    /*
     * Lx: longitud en direccion x;
     * Ly: longitud en direccion y
     * Nx: numero de elementos rectangulares en direccion x
     * Ny: numero de elementos rectangulares en direccion y
     * diag: en caso de elementos triangulares define la diagonal, true = '/', false = '\'
    */
};

//extern template struct rectangular_mesh<ELEMENT_TYPE::QUAD4>;
//extern template struct rectangular_mesh<ELEMENT_TYPE::TRI3>;

#endif // MESH_H
