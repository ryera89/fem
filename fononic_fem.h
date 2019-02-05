#ifndef FONONIC_FEM_H
#define FONONIC_FEM_H

#include "gaussian_cuadrature.h"
#include "fem.h"
#include "utilities.h"

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
inline std::vector<HMatComplexd> fononic_elemental_stiffness_matrix(const VecDoub &vC11,const VecDoub &vC12,const VecDoub &vC33,
                                                                 const VecDoub &k,const gaussian_cuadrature &g_cuad,
                                                                 const std::vector<MatDoub> &v_ecoord,
                                                                 VecDoub (*shapeFun)(QPointF p),MatDoub (*shapeFunDer)(QPointF p)){
    /* Input
     * v_dN: shape functions derivatives for every element.
     * N: shape functions.
     * k: wave number vector
    */
    /* Output
     * dNnm + i*kn*Nm   n:(0,dim-1); m:(0-nodxelem-1) for every element
    */
    size_t nelem = vC11.size();
    size_t dim = k.size();
    size_t nodxelem = g_cuad.nodxelem;
    size_t dofxelem = nodxelem*dim;
    size_t ipoints = g_cuad.gauss_points_number;
    assert(vC12.size() == nelem && vC33.size() == nelem && dim == 2);
    std::vector<HMatComplexd> v_Kelem(nelem,HMatComplexd(dofxelem)); //TODO ver si esto se incializa a (0,0)
    for (size_t i = 0; i < ipoints; ++i){
        QPointF p = g_cuad.gauss_points[i];
        double w = g_cuad.weights[i];
        VecDoub N = shapeFun(p); //shape functions
        MatComplexd imagB = fononic_Bimag(N,k); //imaginary part
        MatDoub dN = shapeFunDer(p); //shape functions derivatives
        assert(N.size() == nodxelem && dN.cols() == nodxelem && dN.rows() == 2);
        for (size_t j = 0; j < nelem; ++j){
            MatDoub JT = dN*v_ecoord[j];
            auto [invJT,detJT] = inv_det_2x2_Matrix(JT);
            MatComplexd B = invJT*dN + imagB;
            MatComplexd conjB = invJT*dN - imagB;
            double wdetJ = w*detJT;
            double C11 = vC11(j);
            double C12 = vC12(j);
            double C33 = vC33(j);
            for (size_t m = 0; m < dofxelem; ++m){
                size_t m_indx = m/2;
                for (size_t n = m; n < dofxelem; ++n){
                    size_t n_indx = n/2;
                    if (n%2){ //n: impar
           /*m: impar*/ if (m%2) v_Kelem[j](m,n) += wdetJ*(C11*conjB(m_indx,1)*B(n_indx,1) + C33*conjB(m_indx,0)*B(n_indx,0));
           /*m: par*/   else v_Kelem[j](m,n) += wdetJ*(C12*conjB(m_indx,0)*B(n_indx,1) + C33*conjB(m_indx,1)*B(n_indx,0));
                    }else{ //n: par
           /*m: impar*/ if (m%2) v_Kelem[j](m,n) += wdetJ*(C12*conjB(m_indx,1)*B(n_indx,0) + C33*conjB(m_indx,0)*B(n_indx,1));
           /*m: par*/   else v_Kelem[j](m,n) += wdetJ*(C11*conjB(m_indx,0)*B(n_indx,0) + C33*conjB(m_indx,1)*B(n_indx,1));
                    }
                }
            }
        }
    }
    return v_Kelem;
}
inline std::vector<SMatDoub> fononic_elemental_mass_matrix(const VecDoub &rho,const std::vector<MatDoub> &v_ecoord,
                                                          const gaussian_cuadrature &g_cuad,VecDoub (*shapeFun)(QPointF p),
                                                          MatDoub (*shapeFunDer)(QPointF p)){
    assert(rho.size() == v_ecoord.size());
    size_t nelem = rho.size();
    uint32_t nodxelem = g_cuad.nodxelem;
    uint32_t ipoints = g_cuad.gauss_points_number;
    std::vector<SMatDoub> v_Melem(nelem,SMatDoub(nodxelem)); //TODO ver si esto se incializa a
    for (uint32_t i = 0; i < ipoints; ++i){
        QPointF p = g_cuad.gauss_points[i];
        double w = g_cuad.weights[i];
        VecDoub N = shapeFun(p); //shape functions
        MatDoub dN = shapeFunDer(p); //shape functions derivatives
        assert(N.size() == nodxelem && dN.cols() == nodxelem && dN.rows() == 2);
        for (size_t j = 0; j < nelem; ++j){
            double erho = rho(j);
            MatDoub JT = dN*v_ecoord[j];
            double detJT = det_2x2_Matrix(JT);
            double tmp = w*erho*detJT;
            for (size_t m = 0; m < nodxelem; ++m){
                for (size_t n = m; n < nodxelem; ++n){
                    v_Melem[i](m,n) += tmp*N(m)*N(n);
                }
            }
        }
    }
    return v_Melem;
}
template<typename T,ELEMENT_TYPE etype>
inline Matrix<T,2,MATRIX_TYPE::CSR> fononic_reduced_system(const Matrix<T,2,MATRIX_TYPE::CSR> &K,const rectangular_mesh<etype> &mesh){

    uint32_t dofxnode = mesh.dofxnode;

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
        c3[0] = mesh.m_corner_dof[6]; c3[1] = mesh.m_corner_dof[8];
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

    Matrix<T,2,MATRIX_TYPE::CSR> K00(K(mesh.m_interior_dof,mesh.m_interior_dof));

    Matrix<T,2,MATRIX_TYPE::CSR> K01 = (K(mesh.m_interior_dof,mesh.m_left_bottom_dof)) +
                                       (K(mesh.m_interior_dof,mesh.m_right_top_dof));


    Matrix<T,2,MATRIX_TYPE::CSR> K02 = (K(mesh.m_interior_dof,c0)) + (K(mesh.m_interior_dof,c1)) +
                                       (K(mesh.m_interior_dof,c2)) + (K(mesh.m_interior_dof,c3));



    Matrix<T,2,MATRIX_TYPE::CSR> K10 = (K(mesh.m_left_bottom_dof,mesh.m_interior_dof)) +
                                       (K(mesh.m_right_top_dof,mesh.m_interior_dof));

    Matrix<T,2,MATRIX_TYPE::CSR> K11 = (K(mesh.m_left_bottom_dof,mesh.m_left_bottom_dof)) +
                                       (K(mesh.m_left_bottom_dof,mesh.m_right_top_dof)) +
                                       (K(mesh.m_right_top_dof,mesh.m_left_bottom_dof)) +
                                       (K(mesh.m_right_top_dof,mesh.m_right_top_dof));

    Matrix<T,2,MATRIX_TYPE::CSR> K12 = (K(mesh.m_left_bottom_dof,c0)) + (K(mesh.m_left_bottom_dof,c1)) +
                                       (K(mesh.m_left_bottom_dof,c2)) + (K(mesh.m_left_bottom_dof,c3)) +
                                       (K(mesh.m_right_top_dof,c0)) + (K(mesh.m_right_top_dof,c1)) +
                                       (K(mesh.m_right_top_dof,c2)) + (K(mesh.m_right_top_dof,c3));

    Matrix<T,2,MATRIX_TYPE::CSR> K20 = (K(c0,mesh.m_interior_dof)) + (K(c1,mesh.m_interior_dof)) +
                                       (K(c2,mesh.m_interior_dof)) + (K(c3,mesh.m_interior_dof));

    Matrix<T,2,MATRIX_TYPE::CSR> K21 = (K(c0,mesh.m_left_bottom_dof)) + (K(c1,mesh.m_left_bottom_dof)) +
                                       (K(c2,mesh.m_left_bottom_dof)) + (K(c3,mesh.m_left_bottom_dof)) +
                                       (K(c0,mesh.m_right_top_dof)) + (K(c1,mesh.m_right_top_dof)) +
                                       (K(c2,mesh.m_right_top_dof)) + (K(c3,mesh.m_right_top_dof));


    Matrix<T,2,MATRIX_TYPE::CSR> K22 = (K(c0,c0)) + (K(c0,c1)) + (K(c0,c2)) + (K(c0,c3)) +
                                       (K(c1,c0)) + (K(c1,c1)) + (K(c1,c2)) + (K(c1,c3)) +
                                       (K(c2,c0)) + (K(c2,c1)) + (K(c2,c2)) + (K(c2,c3)) +
                                       (K(c3,c0)) + (K(c3,c1)) + (K(c3,c2)) + (K(c3,c3));


    uint32_t dim1 = mesh.m_interior_dof.size();
    uint32_t dim2 = mesh.m_left_bottom_dof.size();
    uint32_t dim3 = c0.size();
    uint32_t dim = dim1+dim2+dim3;

    std::vector<T> vvals;
    std::vector<uint32_t> vcols;
    std::vector<uint32_t> row_start;
    std::vector<uint32_t> row_end;

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
    for (uint32_t i = 0; i < dim1; ++i){
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
