#include <QCoreApplication>

#include "mesh_generation.h"
#include "iostream"
#include "mesh.h"

int main(int argc, char *argv[])
{
QCoreApplication a(argc, argv);

uint32_t nelem_x = 3;
uint32_t nelem_y = 3;

//std::vector<element<ELEMENT_TYPE::QUAD,4>> elements = get_rectangular_mesh<ELEMENT_TYPE::QUAD,4>(1.0,1.0,nelem_x,nelem_y);
std::vector<element<ELEMENT_TYPE::QUAD4,4>> elements_quad = get_quad4_rect_mesh(1.0,1.0,nelem_x,nelem_y);
std::vector<element<ELEMENT_TYPE::TRI3,3>> elements_tri1 = get_tri3_rect_mesh(1.0,1.0,nelem_x,nelem_y);
std::vector<element<ELEMENT_TYPE::TRI3,3>> elements_tri2 = get_tri3_rect_mesh(1.0,1.0,nelem_x,nelem_y,false);

rectangular_mesh<ELEMENT_TYPE::QUAD4> quad4_mesh(1.0,1.0,nelem_x,nelem_y);
rectangular_mesh<ELEMENT_TYPE::TRI3> tri3_mesh_1(1,1,nelem_x,nelem_y,true);
rectangular_mesh<ELEMENT_TYPE::TRI3> tri3_mesh_2(1,1,nelem_x,nelem_y,false);

std::cout << "QUAD4 MESH \n";
for (auto &elem:elements_quad){
    std::cout << elem.nodes[0].id << " (" << elem.nodes[0].coord.x() << "," << elem.nodes[0].coord.y() << ") " << "-- "
              << elem.nodes[1].id << " (" << elem.nodes[1].coord.x() << "," << elem.nodes[1].coord.y() << ") " << "-- "
              << elem.nodes[2].id << " (" << elem.nodes[2].coord.x() << "," << elem.nodes[2].coord.y() << ") " << "-- "
              << elem.nodes[3].id << " (" << elem.nodes[3].coord.x() << "," << elem.nodes[3].coord.y() << ") " << "\n";
}
std::cout << quad4_mesh.element_connect << "\n";

std::cout << "TRI3_1 MESH \n";
for (auto &elem:elements_tri1){
    std::cout << elem.nodes[0].id << " (" << elem.nodes[0].coord.x() << "," << elem.nodes[0].coord.y() << ") " << "-- "
              << elem.nodes[1].id << " (" << elem.nodes[1].coord.x() << "," << elem.nodes[1].coord.y() << ") " << "-- "
              << elem.nodes[2].id << " (" << elem.nodes[2].coord.x() << "," << elem.nodes[2].coord.y() << ") " << "\n";
}
std::cout << tri3_mesh_1.element_connect << "\n";

std::cout << "TRI3_2 MESH \n";
for (auto &elem:elements_tri2){
    std::cout << elem.nodes[0].id << " (" << elem.nodes[0].coord.x() << "," << elem.nodes[0].coord.y() << ") " << "-- "
              << elem.nodes[1].id << " (" << elem.nodes[1].coord.x() << "," << elem.nodes[1].coord.y() << ") " << "-- "
              << elem.nodes[2].id << " (" << elem.nodes[2].coord.x() << "," << elem.nodes[2].coord.y() << ") " << "\n";
}
std::cout << tri3_mesh_2.element_connect << "\n";

int i = 0;
for (auto &coor:quad4_mesh.nodes_coordinates) std::cout << i++ << "  " << coor.x() << "  " << coor.y() << "\n";

std::vector<MatDoub> X(5,MatDoub(4,2));
MatDoub dN = {{1,1,1,1},
              {2,2,2,2}};

std::iota(X[0].begin(),X[0].end(),0.0);
std::iota(X[1].begin(),X[1].end(),8);
std::iota(X[2].begin(),X[2].end(),16);
std::iota(X[3].begin(),X[3].end(),24);
std::iota(X[4].begin(),X[4].end(),32);

std::vector<MatDoub> v_JT = transpose_jacobian(dN,X);

std::cout << v_JT[0] << "\n";
std::cout << dN*X[0] << "\n";
std::cout << v_JT[1] << "\n";
std::cout << dN*X[1] << "\n";
std::cout << v_JT[2] << "\n";
std::cout << dN*X[2] << "\n";
std::cout << v_JT[3] << "\n";
std::cout << dN*X[3] << "\n";
std::cout << v_JT[4] << "\n";
std::cout << dN*X[4] << "\n";

std::complex<double> z(2,3);
double r = 10;
std::complex<double> z2 = r + z;

std::cout << z2 << "\n";

return a.exec();
}
