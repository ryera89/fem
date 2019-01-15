#include <QCoreApplication>

#include "mesh_generation.h"
#include "iostream"

int main(int argc, char *argv[])
{
QCoreApplication a(argc, argv);

uint32_t nelem_x = 4;
uint32_t nelem_y = 4;

//std::vector<element<ELEMENT_TYPE::QUAD,4>> elements = get_rectangular_mesh<ELEMENT_TYPE::QUAD,4>(1.0,1.0,nelem_x,nelem_y);
std::vector<element<ELEMENT_TYPE::QUAD,4>> elements_quad = get_quad4_rect_mesh(1.0,1.0,nelem_x,nelem_y);
std::vector<element<ELEMENT_TYPE::TRI,3>> elements_tri1 = get_tri3_rect_mesh(1.0,1.0,nelem_x,nelem_y);
std::vector<element<ELEMENT_TYPE::TRI,3>> elements_tri2 = get_tri3_rect_mesh(1.0,1.0,nelem_x,nelem_y,false);

std::cout << "QUAD4 MESH \n";
for (auto &elem:elements_quad){
    std::cout << elem.nodes[0].id << " (" << elem.nodes[0].coord.x() << "," << elem.nodes[0].coord.y() << ") " << "-- "
              << elem.nodes[1].id << " (" << elem.nodes[1].coord.x() << "," << elem.nodes[1].coord.y() << ") " << "-- "
              << elem.nodes[2].id << " (" << elem.nodes[2].coord.x() << "," << elem.nodes[2].coord.y() << ") " << "-- "
              << elem.nodes[3].id << " (" << elem.nodes[3].coord.x() << "," << elem.nodes[3].coord.y() << ") " << "\n";
}

std::cout << "TRI3_1 MESH \n";
for (auto &elem:elements_tri1){
    std::cout << elem.nodes[0].id << " (" << elem.nodes[0].coord.x() << "," << elem.nodes[0].coord.y() << ") " << "-- "
              << elem.nodes[1].id << " (" << elem.nodes[1].coord.x() << "," << elem.nodes[1].coord.y() << ") " << "-- "
              << elem.nodes[2].id << " (" << elem.nodes[2].coord.x() << "," << elem.nodes[2].coord.y() << ") " << "\n";
}

std::cout << "TRI3_2 MESH \n";
for (auto &elem:elements_tri2){
    std::cout << elem.nodes[0].id << " (" << elem.nodes[0].coord.x() << "," << elem.nodes[0].coord.y() << ") " << "-- "
              << elem.nodes[1].id << " (" << elem.nodes[1].coord.x() << "," << elem.nodes[1].coord.y() << ") " << "-- "
              << elem.nodes[2].id << " (" << elem.nodes[2].coord.x() << "," << elem.nodes[2].coord.y() << ") " << "\n";
}

return a.exec();
}
