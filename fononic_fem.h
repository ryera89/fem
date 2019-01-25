#ifndef FONONIC_FEM_H
#define FONONIC_FEM_H

#include "gaussian_cuadrature.h"
#include "fem.h"
#include "utilities.h"


inline MatCompx fononic_Bimag(const VecDoub &N,const VecDoub &k){
    MatCompx TMP(k.size(),N.size());
    for (size_t i = 0; i < TMP.rows(); ++i){
        for (size_t j = 0; j < TMP.cols(); ++j){
            TMP(i,j) = complex(0,k(i)*N(j));
        }
    }
    return TMP;
}
inline std::vector<HMatCompx> fononic_elemental_stiffness_matrix(const VecDoub &vC11,const VecDoub &vC12,const VecDoub &vC33,
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
    std::vector<HMatCompx> v_Kelem(nelem,HMatCompx(dofxelem)); //TODO ver si esto se incializa a (0,0)
    for (size_t i = 0; i < ipoints; ++i){
        QPointF p = g_cuad.gauss_points[i];
        double w = g_cuad.weights[i];
        VecDoub N = shapeFun(p); //shape functions
        MatCompx imagB = fononic_Bimag(N,k); //imaginary part
        MatDoub dN = shapeFunDer(p); //shape functions derivatives
        assert(N.size() == nodxelem && dN.cols() == nodxelem && dN.rows() == 2);
        for (size_t j = 0; j < nelem; ++j){
            MatDoub JT = dN*v_ecoord[j];
            auto [invJT,detJT] = inv_det_2x2_Matrix(JT);
            MatCompx B = invJT*dN + imagB;
            MatCompx conjB = invJT*dN - imagB;
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

#endif // FONONIC_FEM_H
