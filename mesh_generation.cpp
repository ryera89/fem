#include "mesh_generation.h"

std::vector<element<ELEMENT_TYPE::QUAD,4>> get_quad4_rect_mesh(double Lx, double Ly, uint32_t nelem_x, uint32_t nelem_y)
{
    double dx = Lx/nelem_x;
    double dy = Ly/nelem_y;
    uint32_t nelem = nelem_x*nelem_y;

    std::vector<element<ELEMENT_TYPE::QUAD,4>> elements(nelem);
    //elements.reserve(nelem);

    for (uint32_t i = 0; i < nelem_y; ++i){
        for (uint32_t j = 0 ; j < nelem_x; ++j){
            elements[i*nelem_x+j].nodes[0].id = i*nelem_x + j+i;
            elements[i*nelem_x+j].nodes[0].coord.setX(j*dx);
            elements[i*nelem_x+j].nodes[0].coord.setY(i*dy);

            elements[i*nelem_x+j].nodes[1].id = i*nelem_x + j+i+1;
            elements[i*nelem_x+j].nodes[1].coord.setX((j+1)*dx);
            elements[i*nelem_x+j].nodes[1].coord.setY(i*dy);

            elements[i*nelem_x+j].nodes[2].id = (i+1)*nelem_x + j+i+1;
            elements[i*nelem_x+j].nodes[2].coord.setX(j*dx);
            elements[i*nelem_x+j].nodes[2].coord.setY((i+1)*dy);

            elements[i*nelem_x+j].nodes[3].id = (i+1)*nelem_x + j+i+2;
            elements[i*nelem_x+j].nodes[3].coord.setX((j+1)*dx);
            elements[i*nelem_x+j].nodes[3].coord.setY((i+1)*dy);
        }
    }
    return elements;
}

template<ELEMENT_TYPE etype, uint8_t nnod>
std::vector<element<etype, nnod> > get_rectangular_mesh(double Lx, double Ly, uint32_t nelem_x, uint32_t nelem_y)
{
    if constexpr (etype == ELEMENT_TYPE::QUAD && nnod == 4){
        return get_quad4_rect_mesh(Lx,Ly,nelem_x,nelem_y);
    }
}


