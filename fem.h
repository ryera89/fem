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
typedef Matrix<double,2,MATRIX_TYPE::CSR> Sparse_MatDoub;
typedef Matrix<complexd,2,MATRIX_TYPE::CSR> Sparse_MatComplexd;

template<typename T>
struct indexs_val{
    uint32_t row;
    uint32_t col;
    T val;
    indexs_val() = default;
    indexs_val(uint32_t r,uint32_t c,T v):row(r),col(c),val(v){}
};
template<typename T>
inline bool operator == (const indexs_val<T> &iv1,const indexs_val<T> &iv2){return (iv1.row == iv2.row && iv1.col == iv2.col);}
template<typename T>
inline bool operator != (const indexs_val<T> &iv1,const indexs_val<T> &iv2){return !(iv1 == iv2);}
template<typename T>
inline bool operator < (const indexs_val<T> &iv1,const indexs_val<T> &iv2){
    return ((iv1.row < iv2.row) || (iv1.row == iv2.row && iv1.col < iv2.col));
}
template<typename T>
inline bool operator > (const indexs_val<T> &iv1,const indexs_val<T> &iv2){
    return ((iv1.row > iv2.row) || (iv1.row == iv2.row && iv1.col > iv2.col));
}
template<typename T>
inline bool operator <= (const indexs_val<T> &iv1,const indexs_val<T> &iv2){
    return (iv1 < iv2 || iv1 == iv2);
}
template<typename T>
inline bool operator >= (const indexs_val<T> &iv1,const indexs_val<T> &iv2){
    return (iv1 > iv2 || iv1 == iv2);
}

Sparse_MatComplexd Sparse(const std::vector<uint32_t> &vrows,const std::vector<uint32_t> &vcols,const std::vector<complexd> &vvals,
                 uint32_t nrow,uint32_t ncol){

    assert(vrows.size() == vcols.size() && vrows.size() == vvals.size());

    if (vrows.size() == 0) return Sparse_MatComplexd();

    std::vector<indexs_val<complexd>> vIndxValsTable(1,indexs_val<complexd>(vrows[0],vcols[0],vvals[0]));
    for (size_t ii = 1; ii < vrows.size(); ++ii){
        indexs_val<complexd> tmp(vrows[ii],vcols[ii],vvals[ii]);

        auto iter = std::find(vIndxValsTable.begin(),vIndxValsTable.end(),tmp);
        if (iter != vIndxValsTable.end()) iter->val+=tmp.val; //sumando aportes de distintos elementos al mismo nodo
        else vIndxValsTable.push_back(tmp); //sino esta el nodo se agrega
    }

    std::sort(vIndxValsTable.begin(),vIndxValsTable.end()); //se ordena la lista con privilegio de filas

    uint32_t nvals = vIndxValsTable.size();
    std::vector<complexd> values(nvals);
    std::vector<uint32_t> cols(nvals);
    std::vector<uint32_t> row_start(nrow);
    std::vector<uint32_t> row_end(nrow);
    uint32_t tmprow = std::numeric_limits<uint32_t>::max();
    uint32_t curr_row = 0;
    for (uint32_t i = 0; i < nvals; ++i){
        values[i] = vIndxValsTable[i].val;
        cols[i] = vIndxValsTable[i].col;
        curr_row = vIndxValsTable[i].row;
        if (curr_row != tmprow){ //primer elemento de la fila
            row_start[curr_row] = i;
            tmprow = curr_row;
            if (i != 0){ //fin de la fila anterior
                row_end[curr_row-1] = i;
            }
        }
    }
    row_end[nrow-1] = nvals;

    return Sparse_MatComplexd(nrow,ncol,row_start,row_end,cols,values);
}

Sparse_MatDoub Sparse(const std::vector<uint32_t> &vrows,const std::vector<uint32_t> &vcols,const std::vector<double> &vvals,
                 uint32_t nrow,uint32_t ncol){

    assert(vrows.size() == vcols.size() && vrows.size() == vvals.size());

    if (vrows.size() == 0) return Sparse_MatDoub();

    std::vector<indexs_val<double>> vIndxValsTable(1,indexs_val<double>(vrows[0],vcols[0],vvals[0]));
    for (size_t ii = 1; ii < vrows.size(); ++ii){
        indexs_val<double> tmp(vrows[ii],vcols[ii],vvals[ii]);

        auto iter = std::find(vIndxValsTable.begin(),vIndxValsTable.end(),tmp);
        if (iter != vIndxValsTable.end()) iter->val+=tmp.val; //sumando aportes de distintos elementos al mismo nodo
        else vIndxValsTable.push_back(tmp); //sino esta el nodo se agrega
    }

    std::sort(vIndxValsTable.begin(),vIndxValsTable.end()); //se ordena la lista con privilegio de filas

    uint32_t nvals = vIndxValsTable.size();
    std::vector<double> values(nvals);
    std::vector<uint32_t> cols(nvals);
    std::vector<uint32_t> row_start(nrow);
    std::vector<uint32_t> row_end(nrow);
    uint32_t tmprow = std::numeric_limits<uint32_t>::max();
    uint32_t curr_row = 0;
    for (uint32_t i = 0; i < nvals; ++i){
        values[i] = vIndxValsTable[i].val;
        cols[i] = vIndxValsTable[i].col;
        curr_row = vIndxValsTable[i].row;
        if (curr_row != tmprow){ //primer elemento de la fila
            row_start[curr_row] = i;
            tmprow = curr_row;
            if (i != 0){ //fin de la fila anterior
                row_end[curr_row-1] = i;
            }
        }
    }
    row_end[nrow-1] = nvals;
    return Sparse_MatDoub(nrow,ncol,row_start,row_end,cols,values);
}
template<ELEMENT_TYPE etype>
inline MatDoub element_vnodes_coordinates(const rectangular_mesh<etype> &mesh,size_t element_number){

    size_t nvnodes = mesh.element_connect.cols();
    MatDoub ecoord(nvnodes,2); //2D

    for (size_t i = 0; i < nvnodes; ++i){
        uint32_t node_id = mesh.element_connect(element_number,i);
        ecoord(i,0) = mesh.nodes_coordinates[node_id].x();
        ecoord(i,1) = mesh.nodes_coordinates[node_id].y();
    }

    return ecoord;
}

#endif // FEM_H
