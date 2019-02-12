#include <QCoreApplication>

#include "mesh_generation.h"
#include "iostream"
#include "mesh.h"
#include "fononic_fem.h"

using namespace std;

int main(int argc, char *argv[])
{
QCoreApplication a(argc, argv);

uint32_t nelem_x = 3;
uint32_t nelem_y = 3;

//std::vector<element<ELEMENT_TYPE::QUAD,4>> elements = get_rectangular_mesh<ELEMENT_TYPE::QUAD,4>(1.0,1.0,nelem_x,nelem_y);
//std::vector<element<ELEMENT_TYPE::QUAD4,4>> elements_quad = get_quad4_rect_mesh(1.0,1.0,nelem_x,nelem_y);
//std::vector<element<ELEMENT_TYPE::TRI3,3>> elements_tri1 = get_tri3_rect_mesh(1.0,1.0,nelem_x,nelem_y);
//std::vector<element<ELEMENT_TYPE::TRI3,3>> elements_tri2 = get_tri3_rect_mesh(1.0,1.0,nelem_x,nelem_y,false);

rectangular_mesh<ELEMENT_TYPE::QUAD4> quad4_mesh(2,1.0,1.0,nelem_x,nelem_y);
quad4_mesh.split_mesh_nodes_and_dof();
//rectangular_mesh<ELEMENT_TYPE::TRI3> tri3_mesh_1(1,1,nelem_x,nelem_y,true);
//rectangular_mesh<ELEMENT_TYPE::TRI3> tri3_mesh_2(1,1,nelem_x,nelem_y,false);

//std::cout << "QUAD4 MESH \n";
//for (auto &elem:elements_quad){
//    std::cout << elem.nodes[0].id << " (" << elem.nodes[0].coord.x() << "," << elem.nodes[0].coord.y() << ") " << "-- "
//              << elem.nodes[1].id << " (" << elem.nodes[1].coord.x() << "," << elem.nodes[1].coord.y() << ") " << "-- "
//              << elem.nodes[2].id << " (" << elem.nodes[2].coord.x() << "," << elem.nodes[2].coord.y() << ") " << "-- "
//              << elem.nodes[3].id << " (" << elem.nodes[3].coord.x() << "," << elem.nodes[3].coord.y() << ") " << "\n";
//}
std::cout << quad4_mesh.m_element_connect << "\n";

//std::cout << "TRI3_1 MESH \n";
//for (auto &elem:elements_tri1){
//    std::cout << elem.nodes[0].id << " (" << elem.nodes[0].coord.x() << "," << elem.nodes[0].coord.y() << ") " << "-- "
//              << elem.nodes[1].id << " (" << elem.nodes[1].coord.x() << "," << elem.nodes[1].coord.y() << ") " << "-- "
//              << elem.nodes[2].id << " (" << elem.nodes[2].coord.x() << "," << elem.nodes[2].coord.y() << ") " << "\n";
//}
//std::cout << tri3_mesh_1.element_connect << "\n";

//std::cout << "TRI3_2 MESH \n";
//for (auto &elem:elements_tri2){
//    std::cout << elem.nodes[0].id << " (" << elem.nodes[0].coord.x() << "," << elem.nodes[0].coord.y() << ") " << "-- "
//              << elem.nodes[1].id << " (" << elem.nodes[1].coord.x() << "," << elem.nodes[1].coord.y() << ") " << "-- "
//              << elem.nodes[2].id << " (" << elem.nodes[2].coord.x() << "," << elem.nodes[2].coord.y() << ") " << "\n";
//}
//std::cout << tri3_mesh_2.element_connect << "\n";

int i = 0;
for (auto &coor:quad4_mesh.m_nodes_coordinates) std::cout << i++ << "  " << coor.x() << "  " << coor.y() << "\n";

cout << "left-bottom" << "\n";
for (auto &dof:quad4_mesh.m_left_bottom_dof) cout << dof << "\n";
cout << "right-top" << "\n";
for (auto &dof:quad4_mesh.m_right_top_dof) cout << dof << "\n";
cout << "corner" << "\n";
for (auto &dof:quad4_mesh.m_corner_dof) cout << dof << "\n";
cout << "interior" << "\n";
for (auto &dof:quad4_mesh.m_interior_dof) cout << dof << "\n";
//std::vector<MatDoub> X(5,MatDoub(4,2));
//MatDoub dN = {{1,1,1,1},
//              {2,2,2,2}};

//std::iota(X[0].begin(),X[0].end(),0.0);
//std::iota(X[1].begin(),X[1].end(),8);
//std::iota(X[2].begin(),X[2].end(),16);
//std::iota(X[3].begin(),X[3].end(),24);
//std::iota(X[4].begin(),X[4].end(),32);

//std::vector<MatDoub> v_JT = transpose_jacobian(dN,X);

//std::cout << v_JT[0] << "\n";
//std::cout << dN*X[0] << "\n";
//std::cout << v_JT[1] << "\n";
//std::cout << dN*X[1] << "\n";
//std::cout << v_JT[2] << "\n";
//std::cout << dN*X[2] << "\n";
//std::cout << v_JT[3] << "\n";
//std::cout << dN*X[3] << "\n";
//std::cout << v_JT[4] << "\n";
//std::cout << dN*X[4] << "\n";

//std::complex<double> z(2,3);
//double r = 10;
//std::complex<double> z2 = r + z;

//std::cout << z2 << "\n";

//std::cout << 7/2 << "\n";
double E1 = 1;
double E2 = 200;
double pois = 0.3;
//double lamda1 = E1*pois/((1+pois)*(1-2*pois));
//double mu1 = E1/(2*(1+pois));
double rho1 = 2;
//double lamda2 = E2*pois/((1+pois)*(1-2*pois));
//double mu2 = E2/(2*(1+pois));
double rho2 = 2;

isotropic_material mat1(E1,pois,rho1);
isotropic_material mat2(E2,pois,rho2);

rectangular_mesh<ELEMENT_TYPE::QUAD4> micro_cell_mesh(2,1.0,1.0,10,10);
micro_cell_mesh.split_mesh_nodes_and_dof();

//cout << micro_cell_mesh.m_element_connect << endl;

VecDoub k = {0,0};

VecDoub Xe(micro_cell_mesh.m_element_number);
Xe = 0.01;

gaussian_cuadrature gcuad(ELEMENT_TYPE::QUAD4);

auto start = std::chrono::high_resolution_clock::now();
Matrix<complexd,2,MATRIX_TYPE::CSR> K = fononic_elemental_stiffness_matrix(Xe,mat1,mat2,k,gcuad,micro_cell_mesh,quad4_master_element_shape_functions
                                                                           ,quad4_master_element_shape_functions_gradient);
auto end = std::chrono::high_resolution_clock::now();
auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end-start).count();
printf("ElapsedTime[sec] = %u \n",elapsed);
cout << "Assembly completed \n";
//K.printData();
cout << K.values().size() << "\n";
cout << K.rows() << "\n";

start = std::chrono::high_resolution_clock::now();
const Matrix<complexd,2,MATRIX_TYPE::CSR> Kred = fononic_reduced_system(K,micro_cell_mesh);
end = std::chrono::high_resolution_clock::now();
elapsed = std::chrono::duration_cast<std::chrono::seconds>(end-start).count();
printf("ElapsedTime[sec] = %u \n",elapsed);
cout << "reduction completed \n";
//Kred.printData();
cout << Kred.values().size() << "\n";
cout << Kred.rows() << "\n";

start = std::chrono::high_resolution_clock::now();
Matrix<complexd,2,MATRIX_TYPE::CSR> M = fononic_elemental_mass_matrix(Xe,mat1,mat2,k,gcuad,micro_cell_mesh,quad4_master_element_shape_functions
                                                                    ,quad4_master_element_shape_functions_gradient);
end = std::chrono::high_resolution_clock::now();
elapsed = std::chrono::duration_cast<std::chrono::seconds>(end-start).count();
printf("ElapsedTime[sec] = %u \n",elapsed);
cout << "Assembly completed \n";

start = std::chrono::high_resolution_clock::now();
const Matrix<complexd,2,MATRIX_TYPE::CSR> Mred = fononic_reduced_system(M,micro_cell_mesh);
end = std::chrono::high_resolution_clock::now();
elapsed = std::chrono::duration_cast<std::chrono::seconds>(end-start).count();
printf("ElapsedTime[sec] = %u \n",elapsed);
cout << "reduction completed \n";
cout << Mred.values().size() << "\n";
cout << Mred.rows() << "\n";

/*********************MKL Eigenvalues******************************/
MKL_INT fpm[128];
feastinit(fpm);
fpm[0] = 1;
//fpm[1] = 8;
//fpm[2] = 12;
//fpm[3] = 20;
//fpm[4] = 0;
//fpm[5] = 0;
//fpm[6] = 5;
//fpm[13] = 0;
fpm[26] = 1;
fpm[27] = 1;
//fpm[63] = 0;

std::vector<MKL_Complex16> A(Kred.nnz());
for (int i = 0; i < Kred.nnz(); ++i){
    A[i].real = Kred.values()[i].real();
    A[i].imag = Kred.values()[i].imag();
}
std::vector<MKL_Complex16> B(Mred.nnz());
for (int i = 0; i < Mred.nnz(); ++i){
    B[i].real = Mred.values()[i].real();
    B[i].imag = Mred.values()[i].imag();
}

const double emin = 0.0;
const double emax = 1500.0;
const char uplo = 'F';
const MKL_INT n = Kred.rows();

std::vector<MKL_INT> rowIndexA(n+1);
for (int i = 0; i < n; ++i) rowIndexA[i] = Kred.row_start()[i]+1;
rowIndexA[n] = Kred.nnz()+1;

std::vector<MKL_INT> ja(Kred.nnz());
for (int i = 0; i < Kred.nnz(); ++i) ja[i] = Kred.columns()[i]+1;


vector<MKL_INT> rowIndexB(n+1);
for (int i = 0; i < n; ++i) rowIndexB[i] = Mred.row_start()[i]+1;
rowIndexB[n] = Mred.nnz()+1;

vector<MKL_INT> jb(Mred.nnz());
for (int i = 0; i < Mred.nnz(); ++i) jb[i] = Mred.columns()[i]+1;

//cout << rowIndexA[0] << "  " << rowIndexA[1] << "\n";

MKL_INT loop;
MKL_INT m0 = 15;
MKL_INT m = 0;
MKL_INT info = 0;

double epsout;
vector<double> eigenvalues(m0);
vector<MKL_Complex16> eigenvectors(n*m0);
//Matrix<complexd,2> eigenvectors(n,m0);
vector<double> res(m0);
//double *eigenvalues = new double[m0];
//complexd *eigenvectors = new complexd[n*m0];
//double * res = new double[m0];
//std::vector<double> eigenvalues(m0);
//std::vector<complexd> eigenvectors(m0);
//std::vector<double> res(m0);

//zfeast_hcsrgv(&uplo,&n,Kred.valuesData(),rowIndexA.data(),ja.data(),Mred.valuesData(),rowIndexB.data(),jb.data(),
//              fpm,&epsout,&loop,&emin,&emax,&m0,eigenvalues.data(),eigenvectors.begin(),&m,res.data(),&info);

zfeast_hcsrgv(&uplo,&n,A.data(),rowIndexA.data(),ja.data(),B.data(),rowIndexB.data(),jb.data(),
              fpm,&epsout,&loop,&emin,&emax,&m0,eigenvalues.data(),eigenvectors.data(),&m,res.data(),&info);

cout << "Number of eigenvalues: " << m << "\n";

/******************************************************************/

//start = std::chrono::high_resolution_clock::now();
//Matrix<complexd,2,MATRIX_TYPE::CSR> SUM = K+K;
//end = std::chrono::high_resolution_clock::now();
//elapsed = std::chrono::duration_cast<std::chrono::seconds>(end-start).count();
//printf("ElapsedTime[sec] = %u \n",elapsed);
//cout << "sum completed \n";
//SUM.printData();

//M.printData();
//cout << M.values().size() << "\n";
//cout << M.rows() << "\n";



return a.exec();
}
