#ifndef FONONIC_FEM_H
#define FONONIC_FEM_H

#include "gaussian_cuadrature.h"
#include "fem.h"
#include "utilities.h"
#include <chrono>

typedef Matrix<complexd,2> MatComplexd;
typedef Matrix<complexd,2,MATRIX_TYPE::HER> HMatComplexd;
typedef Matrix<double,2,MATRIX_TYPE::SYMM> SMatDoub;

inline MatComplexd fononic_Bimag(const VecDoub &N,const VecDoub &k){ //imaginary part i(k x w)
    MatComplexd TMP(k.size(),N.size());
    for (size_t i = 0; i < TMP.rows(); ++i){
        for (size_t j = 0; j < TMP.cols(); ++j){
            TMP(i,j) = complexd(0,k(i)*N(j));
        }
    }
    return TMP;
}
inline Sparse_MatComplexd fononic_elemental_stiffness_matrix(const VecDoub &Xe,isotropic_material mat1,isotropic_material mat2,
                                                             const VecDoub &k,const gaussian_cuadrature &g_cuad,
                                                             const rectangular_mesh<ELEMENT_TYPE::QUAD4> &mesh,
                                                             VecDoub (*shapeFun)(const QPointF &p),
                                                             MatDoub (*shapeFunDer)(const QPointF &p)){
    /* Input
     * Xe: vector de variables de disenho
     * mat1 y mat2 materiales limites de disenho del material
     * k: vector de numero de onda
     * g_cuad: cuadratura gaussiane de integracion
     * mesh: malla
     * shapeFun: puntero a funciones de forma del elemento
     * shapeFunDer: puntero a derivada de funciones de forma
    */
    /* Output
     * dNnm + i*kn*Nm   n:(0,dim-1); m:(0-nodxelem-1) for every element
    */
    size_t nelem = Xe.size();
    size_t dim = k.size();
    size_t nodxelem = mesh.m_element_connect.cols(); //numero de nodos por elemento
    size_t dofxelem = mesh.m_dofxnode*nodxelem; // numero de dof por elemento
    size_t ipoints = g_cuad.gauss_points_number; //puntos de integracion
    //numeros de elementos y dimension 2D
    assert(mesh.m_element_connect.rows() == nelem && dim == 2 && nodxelem == g_cuad.nodxelem);
    std::vector<HMatComplexd> v_Kelem(nelem,HMatComplexd(dofxelem,complexd(0,0)));

    //HMatComplexd Kelem(dofxelem);
    //std::map<index_pair,complexd> table;
    //std::vector<uint32_t> vrowsdofsxelems(dofxelem);
    //uint32_t dofs = mesh.m_dofxnode*mesh.m_nodes_number;
    //std::vector<indexs_val<complexd>> v_index_val_table;
    //v_index_val_table.reserve(dofxelem*dofxelem*nelem);
    for (size_t i = 0; i < ipoints; ++i){
        QPointF p = g_cuad.gauss_points[i];
        double w = g_cuad.weights[i];
        VecDoub N = shapeFun(p); //shape functions
        MatComplexd imagB = fononic_Bimag(N,k); //imaginary part
        MatDoub dN = shapeFunDer(p); //shape functions derivatives
        assert(N.size() == nodxelem && dN.cols() == nodxelem && dN.rows() == 2);
        for (size_t j = 0; j < nelem; ++j){ //TODO: valorar paralelizar aca
            MatDoub MCOOR = element_vnodes_coordinates(mesh,j);
            MatDoub JT = dN*MCOOR ;
            auto [invJT,detJT] = inv_det_2x2_Matrix(JT);
            MatDoub R = invJT*dN;
            MatComplexd B = R + imagB;
            MatComplexd conjB = R - imagB;
            double wdetJ = w*detJT;
//            double C11 = Xe(j)*mat1.m_C11 + (1-Xe(j))*mat2.m_C11;
//            double C12 = Xe(j)*mat1.m_C12 + (1-Xe(j))*mat2.m_C12;
//            double C33 = Xe(j)*mat1.m_C33 + (1-Xe(j))*mat2.m_C33;
            double C11 = 266.5519;
            double C12 = 114.2365;
            double C33 = 76.1577;
            //******************************************************************
//            for (size_t m = 0; m < nodxelem; ++m){
//                //size_t m_indx = m/2;
//                uint32_t m0 = 2*mesh.m_element_connect(j,m);
//                uint32_t m1 = 2*mesh.m_element_connect(j,m) + 1;
//                for (size_t n = m; n < nodxelem; ++n){
//                    //size_t n_indx = n/2;
//                    uint32_t n0 = 2*mesh.m_element_connect(j,n);
//                    uint32_t n1 = 2*mesh.m_element_connect(j,n) + 1;
//                    complexd val00 = wdetJ*(C11*conjB(0,m)*B(0,n) + C33*conjB(1,m)*B(1,n)); //par-par m0-n0
//                    complexd val01 = wdetJ*(C12*conjB(0,m)*B(1,n) + C33*conjB(1,m)*B(0,n)); //par-impar m0-n1
//                    complexd val10 = wdetJ*(C12*conjB(1,m)*B(0,n) + C33*conjB(0,m)*B(1,n)); //impar-par m1-n0
//                    complexd val11 = wdetJ*(C11*conjB(1,m)*B(1,n) + C33*conjB(0,m)*B(0,n)); //impar-impar
//                    if (m!=n){
//                        v_index_val_table.emplace_back(indexs_val<complexd>(m0,n0,val00));
//                        v_index_val_table.emplace_back(indexs_val<complexd>(m0,n1,val01));
//                        v_index_val_table.emplace_back(indexs_val<complexd>(m1,n0,val10));
//                        v_index_val_table.emplace_back(indexs_val<complexd>(m1,n1,val11));
//                        v_index_val_table.emplace_back(indexs_val<complexd>(n0,m0,val00));
//                        v_index_val_table.emplace_back(indexs_val<complexd>(n0,m1,val01));
//                        v_index_val_table.emplace_back(indexs_val<complexd>(n1,m0,val10));
//                        v_index_val_table.emplace_back(indexs_val<complexd>(n1,m1,val11));
//                    }else{
//                        v_index_val_table.emplace_back(indexs_val<complexd>(m0,n0,val00));
//                        v_index_val_table.emplace_back(indexs_val<complexd>(m0,n1,val10));
//                        v_index_val_table.emplace_back(indexs_val<complexd>(m1,n0,val01));
//                        v_index_val_table.emplace_back(indexs_val<complexd>(m1,n1,val11));
//                    }
//                }
//            }
//                        index_pair p00(m0,n0);
//                        index_pair p01(m0,n1);
//                        index_pair p10(m1,n0);
//                        index_pair p11(m1,n1);
//                        index_pair sp00(n0,m0);
//                        index_pair sp01(n0,m1);
//                        index_pair sp10(n1,m0);
//                        index_pair sp11(n1,m1);
//                        table[p00] += val00;
//                        table[sp00] += val00;
//                        table[p01] += val01;
//                        table[sp01] += val01;
//                        table[p10] += val10;
//                        table[sp10] += val10;
//                        table[p11] += val11;
//                        table[sp11] += val11;

//                    }else{
//                        index_pair p00(m0,n0);
//                        index_pair p01(m0,n1);
//                        index_pair p10(m1,n0);
//                        index_pair p11(m1,n1);
//                        table[p00] += val00;
//                        table[p01] += val01;
//                        table[p10] += val10;
//                        table[p11] += val11;
//                    }

//                    if (n%2){ //n: impar
//                              //va un += si la idea de la tabla no funciona
//                        /*m: impar*/ if (m%2) v_Kelem[j](m,n) += wdetJ*(C11*conjB(1,m_indx)*B(1,n_indx) + C33*conjB(0,m_indx)*B(0,n_indx));
//                        /*m: par*/   else v_Kelem[j](m,n) += wdetJ*(C12*conjB(0,m_indx)*B(1,n_indx) + C33*conjB(1,m_indx)*B(0,n_indx));
//                    }else{ //n: par
//                        /*m: impar*/ if (m%2) v_Kelem[j](m,n) += wdetJ*(C12*conjB(1,m_indx)*B(0,n_indx) + C33*conjB(0,m_indx)*B(1,n_indx));
//                        /*m: par*/   else v_Kelem[j](m,n) += wdetJ*(C11*conjB(0,m_indx)*B(0,n_indx) + C33*conjB(1,m_indx)*B(1,n_indx));
//                    }
//                }
//            }
            //*********************************************************************
            for (size_t m = 0; m < dofxelem; ++m){
                size_t m_indx = m/2;
                for (size_t n = m; n < dofxelem; ++n){
                    size_t n_indx = n/2;
                    if (n%2){ //n: impar
                        //va un += si la idea de la tabla no funciona
                        /*m: impar*/ if (m%2) v_Kelem[j](m,n) += wdetJ*(C11*conjB(1,m_indx)*B(1,n_indx) + C33*conjB(0,m_indx)*B(0,n_indx));
                        /*m: par*/   else v_Kelem[j](m,n) += wdetJ*(C12*conjB(0,m_indx)*B(1,n_indx) + C33*conjB(1,m_indx)*B(0,n_indx));
                    }else{ //n: par
                        /*m: impar*/ if (m%2) v_Kelem[j](m,n) += wdetJ*(C12*conjB(1,m_indx)*B(0,n_indx) + C33*conjB(0,m_indx)*B(1,n_indx));
                        /*m: par*/   else v_Kelem[j](m,n) += wdetJ*(C11*conjB(0,m_indx)*B(0,n_indx) + C33*conjB(1,m_indx)*B(1,n_indx));
                    }
                }
            }
//            const HMatComplexd &Kconst = Kelem;

//            for (size_t ii = 0; ii < nodxelem; ++ii){
//                vrowsdofsxelems[2*ii] = 2*mesh.m_element_connect(j,ii);
//                vrowsdofsxelems[2*ii+1] = 2*mesh.m_element_connect(j,ii) + 1;
//            }
//            for (size_t ii = 0; ii < vrowsdofsxelems.size(); ++ii){
//                for (size_t jj = 0; jj < vrowsdofsxelems.size(); ++jj){
//                    index_pair ptmp(vrowsdofsxelems[ii],vrowsdofsxelems[jj]);
//                    complexd val = Kconst(ii,jj);
//                    table[ptmp] += val;
//                }
//            }

        }
    }
    return assembly(mesh,v_Kelem);
    //return Sparse(table,dofs,dofs);
    //return Sparse(v_index_val_table,dofs,dofs);
}
inline Matrix<double,2,MATRIX_TYPE::CSR> fononic_elemental_mass_matrix(const VecDoub &Xe,isotropic_material mat1,isotropic_material mat2,
                                                           const VecDoub &k,const gaussian_cuadrature &g_cuad,
                                                           const rectangular_mesh<ELEMENT_TYPE::QUAD4> &mesh,
                                                           VecDoub (*shapeFun)(const QPointF &p),
                                                           MatDoub (*shapeFunDer)(const QPointF &p)){
    //assert(rho.size() == v_ecoord.size());
    size_t nelem = Xe.size();
    size_t nodxelem = mesh.m_element_connect.cols(); //numero de nodos por elemento
    size_t dofxelem = mesh.m_dofxnode*nodxelem; // numero de dof por elemento
    size_t ipoints = g_cuad.gauss_points_number; //puntos de integracion
    uint32_t nrow = mesh.m_dofxnode*mesh.m_nodes_coordinates.size(); //numero de nodos x dofxnodo
    assert(mesh.m_element_connect.rows() == nelem && nodxelem == g_cuad.nodxelem);
    std::vector<SMatDoub> v_Melem(nelem,SMatDoub(dofxelem,0)); //TODO ver si esto se incializa a
     //SMatDoub Melem(dofxelem);
     std::vector<indexs_val<double>> v_index_val_table;
     for (uint32_t i = 0; i < ipoints; ++i){
        QPointF p = g_cuad.gauss_points[i];
        double w = g_cuad.weights[i];
        VecDoub N = shapeFun(p); //shape functions
        MatDoub dN = shapeFunDer(p); //shape functions derivatives
        assert(N.size() == nodxelem && dN.cols() == nodxelem && dN.rows() == 2);
        for (size_t j = 0; j < nelem; ++j){
            MatDoub MCOOR = element_vnodes_coordinates(mesh,j);
            MatDoub JT = dN*MCOOR;
            double detJT = det_2x2_Matrix(JT);
            double rho = Xe(j)*mat1.m_rho + (1-Xe(j))*mat2.m_rho;
            double tmp = w*rho*detJT;
            for (size_t m = 0; m < nodxelem; ++m){
                for (size_t n = m; n < nodxelem; ++n){
                    double tmpval = tmp*N(m)*N(n);
                    v_Melem[j](2*m,2*n) += tmpval;
                    v_Melem[j](2*m+1,2*n+1) += tmpval;
                }
            }
        }
    }
    return assembly(mesh,v_Melem);
}
template<typename T,ELEMENT_TYPE etype>
inline Matrix<T,2,MATRIX_TYPE::CSR> fononic_reduced_system(const Matrix<T,2,MATRIX_TYPE::CSR> &K,const rectangular_mesh<etype> &mesh){

    uint32_t dofxnode = mesh.m_dofxnode;

    std::vector<uint32_t> c0(dofxnode);
    std::vector<uint32_t> c1(dofxnode);
    std::vector<uint32_t> c2(dofxnode);
    std::vector<uint32_t> c3(dofxnode);

    switch (dofxnode) {
    case 1:
        assert(mesh.m_corner_dof.size() == 4);
        c0[0] = mesh.m_corner_dof[0];
        c1[0] = mesh.m_corner_dof[1];
        c2[0] = mesh.m_corner_dof[2];
        c3[0] = mesh.m_corner_dof[3];
        break;
    case 2:
        assert(mesh.m_corner_dof.size() == 8);
        c0[0] = mesh.m_corner_dof[0]; c0[1] = mesh.m_corner_dof[1];
        c1[0] = mesh.m_corner_dof[2]; c1[1] = mesh.m_corner_dof[3];
        c2[0] = mesh.m_corner_dof[4]; c2[1] = mesh.m_corner_dof[5];
        c3[0] = mesh.m_corner_dof[6]; c3[1] = mesh.m_corner_dof[7];
        break;
    case 3:
        assert(mesh.m_corner_dof.size() == 12);
        c0[0] = mesh.m_corner_dof[0]; c0[1] = mesh.m_corner_dof[1]; c0[2] = mesh.m_corner_dof[2];
        c1[0] = mesh.m_corner_dof[3]; c1[1] = mesh.m_corner_dof[4]; c1[2] = mesh.m_corner_dof[5];
        c2[0] = mesh.m_corner_dof[6]; c2[1] = mesh.m_corner_dof[7]; c2[2] = mesh.m_corner_dof[8];
        c3[0] = mesh.m_corner_dof[9]; c3[1] = mesh.m_corner_dof[10]; c3[2] = mesh.m_corner_dof[11];
        break;
    default:
        throw ("more than 3 dof per node has not been implemented");
        break;
    }

    auto start = std::chrono::high_resolution_clock::now();
    Matrix<T,2,MATRIX_TYPE::CSR> K00(K(mesh.m_interior_dof,mesh.m_interior_dof));
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end-start).count();
    printf("ElapsedTime[sec] = %u \n",elapsed);

    Matrix<T,2,MATRIX_TYPE::CSR> K01 = (K(mesh.m_interior_dof,mesh.m_left_bottom_dof)) +
                                       (K(mesh.m_interior_dof,mesh.m_right_top_dof));


    Matrix<T,2,MATRIX_TYPE::CSR> K02 = (K(mesh.m_interior_dof,c0)) + (K(mesh.m_interior_dof,c1)) +
                                       (K(mesh.m_interior_dof,c2)) + (K(mesh.m_interior_dof,c3));



    Matrix<T,2,MATRIX_TYPE::CSR> K10 = (K(mesh.m_left_bottom_dof,mesh.m_interior_dof)) +
                                       (K(mesh.m_right_top_dof,mesh.m_interior_dof));

    Matrix<T,2,MATRIX_TYPE::CSR> K11 = (K(mesh.m_left_bottom_dof,mesh.m_left_bottom_dof)) +
                                       (K(mesh.m_left_bottom_dof,mesh.m_right_top_dof));

    K11 = K11 + (K(mesh.m_right_top_dof,mesh.m_left_bottom_dof));
    K11 = K11 + (K(mesh.m_right_top_dof,mesh.m_right_top_dof));


    Matrix<T,2,MATRIX_TYPE::CSR> K12 = (K(mesh.m_left_bottom_dof,c0)) + (K(mesh.m_left_bottom_dof,c1));
    K12 = K12 + (K(mesh.m_left_bottom_dof,c2));
    K12 = K12 + (K(mesh.m_left_bottom_dof,c3));
    K12 = K12 + (K(mesh.m_right_top_dof,c0));
    K12 = K12 + (K(mesh.m_right_top_dof,c1));
    K12 = K12 + (K(mesh.m_right_top_dof,c2));
    K12 = K12 + (K(mesh.m_right_top_dof,c3));

    Matrix<T,2,MATRIX_TYPE::CSR> K20 = (K(c0,mesh.m_interior_dof)) + (K(c1,mesh.m_interior_dof));
    K20 = K20 + (K(c2,mesh.m_interior_dof));
    K20 = K20 + (K(c3,mesh.m_interior_dof));

    Matrix<T,2,MATRIX_TYPE::CSR> K21 = (K(c0,mesh.m_left_bottom_dof)) + (K(c1,mesh.m_left_bottom_dof));
    K21 = K21 + (K(c2,mesh.m_left_bottom_dof));
    K21 = K21 + (K(c3,mesh.m_left_bottom_dof));
    K21 = K21 + (K(c0,mesh.m_right_top_dof));
    K21 = K21 + (K(c1,mesh.m_right_top_dof));
    K21 = K21 + (K(c2,mesh.m_right_top_dof));
    K21 = K21 + (K(c3,mesh.m_right_top_dof));


    Matrix<T,2,MATRIX_TYPE::CSR> K22 = (K(c0,c0)) + (K(c0,c1)) + (K(c0,c2)) + (K(c0,c3)) +
                                       (K(c1,c0)) + (K(c1,c1)) + (K(c1,c2)) + (K(c1,c3)) +
                                       (K(c2,c0)) + (K(c2,c1)) + (K(c2,c2)) + (K(c2,c3)) +
                                       (K(c3,c0)) + (K(c3,c1)) + (K(c3,c2)) + (K(c3,c3));



    uint32_t dim1 = mesh.m_interior_dof.size();
    uint32_t dim2 = mesh.m_left_bottom_dof.size();
    uint32_t dim3 = c0.size();
    uint32_t dim = dim1+dim2+dim3;

    size_t n00 = K00.values().size();
    size_t n01 = K01.values().size();
    size_t n02 = K02.values().size();
    size_t n10 = K10.values().size();
    size_t n11 = K11.values().size();
    size_t n12 = K12.values().size();
    size_t n20 = K20.values().size();
    size_t n21 = K21.values().size();
    size_t n22 = K22.values().size();

    size_t nvals = n00+n01+n02+n10+n11+n12+n20+n21+n22;

    std::vector<T> vvals;
    vvals.reserve(nvals);
    std::vector<uint32_t> vcols;
    vcols.reserve(nvals);
    std::vector<uint32_t> row_start;
    row_start.reserve(dim);
    std::vector<uint32_t> row_end;
    row_end.reserve(dim);

    for (uint32_t i = 0; i < dim1; ++i){
        bool first_insertion = true;
        uint32_t beg = K00.row_start()[i];
        uint32_t end = K00.row_end()[i];
        uint32_t col;
        T value;
        for(;beg<end;++beg){
            col = K00.columns()[beg];
            value = K00.values()[beg];
            vvals.push_back(value);
            vcols.push_back(col);
            if (first_insertion){
               row_start.push_back(vvals.size()-1);
               first_insertion = false;
            }
        }
        beg = K01.row_start()[i];
        end = K01.row_end()[i];
        for(;beg<end;++beg){
            col = K01.columns()[beg] + dim1;
            value = K01.values()[beg];
            vvals.push_back(value);
            vcols.push_back(col);
            if (first_insertion){
               row_start.push_back(vvals.size()-1);
               first_insertion = false;
            }
        }
        beg = K02.row_start()[i];
        end = K02.row_end()[i];
        for(;beg<end;++beg){
            col = K02.columns()[beg] + dim1 + dim2;
            value = K02.values()[beg];
            vvals.push_back(value);
            vcols.push_back(col);
            if (first_insertion){
               row_start.push_back(vvals.size()-1);
               first_insertion = false;
            }
        }
        if (first_insertion){
            row_start.push_back(vvals.size());
            row_end.push_back(vvals.size());
        }else row_end.push_back(vvals.size());
    }
    for (uint32_t i = 0; i < dim2; ++i){
        bool first_insertion = true;
        uint32_t beg = K10.row_start()[i];
        uint32_t end = K10.row_end()[i];
        uint32_t col;
        T value;
        for(;beg<end;++beg){
            col = K10.columns()[beg];
            value = K10.values()[beg];
            vvals.push_back(value);
            vcols.push_back(col);
            if (first_insertion){
               row_start.push_back(vvals.size()-1);
               first_insertion = false;
            }
        }
        beg = K11.row_start()[i];
        end = K11.row_end()[i];
        for(;beg<end;++beg){
            col = K11.columns()[beg] + dim1;
            value = K11.values()[beg];
            vvals.push_back(value);
            vcols.push_back(col);
            if (first_insertion){
               row_start.push_back(vvals.size()-1);
               first_insertion = false;
            }
        }
        beg = K12.row_start()[i];
        end = K12.row_end()[i];
        for(;beg<end;++beg){
            col = K12.columns()[beg] + dim1 + dim2;
            value = K12.values()[beg];
            vvals.push_back(value);
            vcols.push_back(col);
            if (first_insertion){
               row_start.push_back(vvals.size()-1);
               first_insertion = false;
            }
        }
        if (first_insertion){
            row_start.push_back(vvals.size());
            row_end.push_back(vvals.size());
        }else row_end.push_back(vvals.size());
    }
    for (uint32_t i = 0; i < dim3; ++i){
        bool first_insertion = true;
        uint32_t beg = K20.row_start()[i];
        uint32_t end = K20.row_end()[i];
        uint32_t col;
        T value;
        for(;beg<end;++beg){
            col = K20.columns()[beg];
            value = K20.values()[beg];
            vvals.push_back(value);
            vcols.push_back(col);
            if (first_insertion){
               row_start.push_back(vvals.size()-1);
               first_insertion = false;
            }
        }
        beg = K21.row_start()[i];
        end = K21.row_end()[i];
        for(;beg<end;++beg){
            col = K21.columns()[beg] + dim1;
            value = K21.values()[beg];
            vvals.push_back(value);
            vcols.push_back(col);
            if (first_insertion){
               row_start.push_back(vvals.size()-1);
               first_insertion = false;
            }
        }
        beg = K22.row_start()[i];
        end = K22.row_end()[i];
        for(;beg<end;++beg){
            col = K22.columns()[beg] + dim1 + dim2;
            value = K22.values()[beg];
            vvals.push_back(value);
            vcols.push_back(col);
            if (first_insertion){
               row_start.push_back(vvals.size()-1);
               first_insertion = false;
            }
        }
        if (first_insertion){
            row_start.push_back(vvals.size());
            row_end.push_back(vvals.size());
        }else row_end.push_back(vvals.size());
    }
    return Matrix<T,2,MATRIX_TYPE::CSR>(dim,dim,row_start,row_end,vcols,vvals);
}

#endif // FONONIC_FEM_H
