#ifndef MESH_GENERATION_H
#define MESH_GENERATION_H

#include "element.h"
#include <vector>

std::vector<element<ELEMENT_TYPE::TRI,3>> get_tri3_rect_mesh(double Lx,double Ly,uint32_t nelem_x,uint32_t nelem_y,bool diag = true);
std::vector<element<ELEMENT_TYPE::QUAD,4>> get_quad4_rect_mesh(double Lx,double Ly,uint32_t nelem_x,uint32_t nelem_y);

template<ELEMENT_TYPE etype,uint8_t nnod>
std::vector<element<etype,nnod>> get_rectangular_mesh(double Lx,double Ly,uint32_t nelem_x,
                                                       uint32_t nelem_y);



#endif // MESH_GENERATION_H
