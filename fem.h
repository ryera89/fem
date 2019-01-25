#ifndef FEM_H
#define FEM_H

#include "mesh.h"

//dN: matriz dN(dim,nodxelem): derivadas de funciones de forma elemento master
//X: vector de coordenadas por elemento X(nelem,COORD), COORD(nodxelem,dim): coordenadas de los nodos por elemento
inline std::vector<MatDoub> transpose_jacobian(const MatDoub &dN,const std::vector<MatDoub> &X){
    size_t nelem = X.size();
    size_t dim = dN.rows();
    assert(dim == 2 || dim == 3); //el problema es en 2 o 3 dimensiones
    std::vector<MatDoub> v_JT(nelem); //Medir si es mas optimo inicializar las matrices v_JT(nelem,MatDoub(dim,dim))
    std::transform(X.begin(),X.end(),v_JT.begin(),[&dN](const auto &x){assert(dN.rows() == x.cols());
                                                                       return dN*x;});
    return v_JT;
}
//dN(2,4); v_invJT(nelem); invJT(2,2)
inline std::vector<MatDoub> elemental_shape_functions_derivatives(const MatDoub &dN,const std::vector<MatDoub> &v_invJT){
    size_t nelem = v_invJT.size();
    std::vector<MatDoub> v_dN(nelem); //Medir si es mas optimo inicializar las matrices v_dN(nelem,MatDoub(dim,dim))
    std::transform(v_invJT.begin(),v_invJT.end(),v_dN.begin(),[&dN](const auto &invJT){assert(invJT.rows() == invJT.cols());
                                                                                       return invJT*dN;});
    return v_dN;
}

#endif // FEM_H
