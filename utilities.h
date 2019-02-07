#ifndef UTILITIES_H
#define UTILITIES_H

#include "ndimmatrix/matrix.h"
#include "mesh.h"
#include <map>

typedef Matrix<double,1> VecDoub;
typedef Matrix<double,2> MatDoub;

inline std::tuple<MatDoub,double> inv_det_2x2_Matrix(const MatDoub &A){
    assert(A.rows() == 2 && A.cols() == 2);
    const double a = A(0,0);
    const double b = A(0,1);
    const double c = A(1,0);
    const double d = A(1,1);
    MatDoub K(2,2);
    K(0,0) = d; K(0,1) = -b;
    K(1,0) = -c; K(1,1) = a;
    double det = (a*d-b*c);
    K/=det;
    return std::tuple<MatDoub,double>(K,det);
}
inline std::tuple<MatDoub,double> inv_det_3x3_Matrix(const MatDoub &A){
    assert(A.rows() == 3 && A.cols() == 3);
    double a = A(0,0);double b = A(0,1);double c = A(0,2);
    double d = A(1,0);double e = A(1,1);double f = A(1,2);
    double g = A(2,0);double h = A(2,1);double i = A(2,2);

    double aa = (e*i - f*h);
    double bb = -(d*i - f*g);
    double cc = (d*h - e*g);

    double dd = -(b*i - c*h);
    double ee = (a*i - c*g);
    double ff = -(a*h - b*g);

    double gg = (b*f - c*e);
    double hh = -(a*f - c*d);
    double ii = (a*e - b*d);

    double det = a*aa + b*bb + c*cc;

    MatDoub R(3,3);

    R(0,0) = aa/det; R(0,1) = dd/det; R(0,2) = gg/det;
    R(1,0) = bb/det; R(1,1) = ee/det; R(1,2) = hh/det;
    R(2,0) = cc/det; R(2,1) = ff/det; R(2,2) = ii/det;

    return std::tuple<MatDoub,double>(R,det);
}
inline double det_3x3_Matrix(const MatDoub &A){
    assert(A.rows() == 3 && A.cols() == 3);
    double a = A(0,0);double b = A(0,1);double c = A(0,2);
    double d = A(1,0);double e = A(1,1);double f = A(1,2);
    double g = A(2,0);double h = A(2,1);double i = A(2,2);

    double aa = (e*i - f*h);
    double bb = -(d*i - f*g);
    double cc = (d*h - e*g);

    return a*aa + b*bb + c*cc;
}
inline double det_2x2_Matrix(const MatDoub &A){
    assert(A.rows() == 2 && A.cols() == 2);
    return A(0,0)*A(1,1) - A(1,0)*A(0,1);
}

typedef Matrix<double,2,MATRIX_TYPE::CSR> Sparse_MatDoub;
typedef Matrix<complexd,2,MATRIX_TYPE::CSR> Sparse_MatComplexd;

struct index_pair{
    uint32_t row;
    uint32_t col;
    index_pair() = default;
    index_pair(uint32_t r,uint32_t c):row(r),col(c){}
};
inline bool operator == (const index_pair &iv1,const index_pair &iv2){return (iv1.row == iv2.row && iv1.col == iv2.col);}
inline bool operator != (const index_pair &iv1,const index_pair &iv2){return !(iv1 == iv2);}
inline bool operator < (const index_pair &iv1,const index_pair &iv2){
    return ((iv1.row < iv2.row) || (iv1.row == iv2.row && iv1.col < iv2.col));
}
inline bool operator > (const index_pair &iv1,const index_pair &iv2){
    return ((iv1.row > iv2.row) || (iv1.row == iv2.row && iv1.col > iv2.col));
}
inline bool operator <= (const index_pair &iv1,const index_pair &iv2){
    return (iv1 < iv2 || iv1 == iv2);
}
inline bool operator >= (const index_pair &iv1,const index_pair &iv2){
    return (iv1 > iv2 || iv1 == iv2);
}
template<typename T>
struct indexs_val{
    uint32_t row;
    uint32_t col;
    T val;
    indexs_val() = default;
    indexs_val(uint32_t r,uint32_t c,T v):row(r),col(c),val(v){}
    indexs_val(const indexs_val &other) = default;
    indexs_val(indexs_val &&other) = default;
    indexs_val& operator=(const indexs_val &other) = default;
    indexs_val& operator=(indexs_val &&other) = default;
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
//*********************************************************************************************************************
template<typename T>
inline Matrix<T,2,MATRIX_TYPE::CSR> Sparse(std::map<index_pair,T> &v_indx_table,uint32_t nrow,uint32_t ncol){
    //std::sort(v_indx_table.begin(),v_indx_table.end());
    size_t nvals = v_indx_table.size();
    //std::vector<T> values(nvals);
    //std::vector<uint32_t> cols(nvals);
    std::vector<uint32_t> row_start(nrow);
    std::vector<uint32_t> row_end(nrow);
    std::vector<T> values;
    values.reserve(nvals);
    std::vector<uint32_t> cols;
    cols.reserve(nvals);
    //std::vector<uint32_t> row_start(nrow);
    //std::vector<uint32_t> row_end(nrow);
    uint32_t tmprow = std::numeric_limits<uint32_t>::max();
    uint32_t curr_row = 0;
    size_t i = 0;
    for (auto &tmp:v_indx_table){
        auto [index,val] = tmp;
        if (val == T()) continue;
        //values[i] = val;
        values.push_back(val);
        //cols[i] = index.col;
        cols.push_back(index.col);
        curr_row = index.row;
        if (curr_row != tmprow){ //primer elemento de la fila
            row_start[curr_row] = i;
            tmprow = curr_row;
            if (i != 0){ //fin de la fila anterior
                row_end[curr_row-1] = i;
            }
        }
        ++i;
    }
    //row_end[nrow-1] = nvals;
    row_end[nrow-1] = values.size();

    return Matrix<T,2,MATRIX_TYPE::CSR>(nrow,ncol,row_start,row_end,cols,values);
}
template<typename T>
inline Matrix<T,2,MATRIX_TYPE::CSR> Sparse(std::vector<indexs_val<T>> &v_indx_table,uint32_t nrow,uint32_t ncol){
    std::sort(v_indx_table.begin(),v_indx_table.end());

    size_t nvals = v_indx_table.size();
    std::vector<T> values;
    values.reserve(nvals);
    std::vector<uint32_t> cols;
    cols.reserve(nvals);
    std::vector<uint32_t> row_start(nrow);
    std::vector<uint32_t> row_end(nrow);
    uint32_t tmprow = std::numeric_limits<uint32_t>::max();
    uint32_t curr_row = 0;

    indexs_val<T> tmp = v_indx_table[0];
    T val = tmp.val;
    //uint32_t row = 0;
    for (size_t i = 1; i < nvals; ++i){
        indexs_val<T> tmp1 = v_indx_table[i];
        if (tmp == tmp1){
            val += tmp1.val;
        }else{
            values.push_back(val);
            cols.push_back(tmp.col);
            curr_row = tmp.row;
            tmp = tmp1;
            val = tmp.val;
            if (curr_row != tmprow){
                if (tmprow != std::numeric_limits<uint32_t>::max()) row_end[curr_row-1] = values.size()-1;
                row_start[curr_row] = values.size()-1;
                tmprow = curr_row;
            }
        }
    }
    values.push_back(val);
    cols.push_back(tmp.col);
    curr_row = tmp.row;
//    if (curr_row != v_indx_table[nvals-1].row){

//    }
    row_end[nrow-1] = values.size();
//    for (uint32_t i = 0; i < nvals; ++i){
//        values[i] = v_indx_table[i].val;
//        cols[i] = v_indx_table[i].col;
//        curr_row = v_indx_table[i].row;
//        if (curr_row != tmprow){ //primer elemento de la fila
//            row_start[curr_row] = i;
//            tmprow = curr_row;
//            if (i != 0){ //fin de la fila anterior
//                row_end[curr_row-1] = i;
//            }
//        }
//    }
//    row_end[nrow-1] = nvals;

    return Matrix<T,2,MATRIX_TYPE::CSR>(nrow,ncol,row_start,row_end,cols,values);
}
template<typename T,MATRIX_TYPE matrix_type>
inline Matrix<T,2,MATRIX_TYPE::CSR> assembly(const rectangular_mesh<ELEMENT_TYPE::QUAD4> &mesh,const std::vector<Matrix<T,2,matrix_type>> &v_Kelem){

    assert(mesh.m_dofxnode == 2);
    uint32_t nelem = mesh.m_element_number;
    uint32_t nodxelem = mesh.m_element_connect.cols(); //numero nodos por elemento
    uint32_t dofxelem = mesh.m_dofxnode*nodxelem; //numero de dof por elemento
    //std::vector<indexs_val<T>> v_table; //(dofxelem*dofxelem*nelem);
    //v_table.reserve(dofxelem*dofxelem); //TODO implementar eso para casos simetrico y hermitico
    std::map<index_pair,T> table;

    std::vector<uint32_t> vrowsdofsxelems(dofxelem);
    //vrowsdofsxelems.reserve(dofxelem);

    for (size_t i = 0; i < nelem; ++i){
        for (size_t ii = 0; ii < nodxelem; ++ii){
            vrowsdofsxelems[2*ii] = 2*mesh.m_element_connect(i,ii);
            vrowsdofsxelems[2*ii+1] = 2*mesh.m_element_connect(i,ii) + 1;
        }
        for (size_t ii = 0; ii < vrowsdofsxelems.size(); ++ii){
            for (size_t jj = 0; jj < vrowsdofsxelems.size(); ++jj){
                index_pair ptmp(vrowsdofsxelems[ii],vrowsdofsxelems[jj]);
                T val = v_Kelem[i](ii,jj);
                table[ptmp] += val;

                //indexs_val<T> tmp(vrowsdofsxelems[ii],vrowsdofsxelems[jj],v_Kelem[i](ii,jj));

                //auto iter = std::find(v_table.begin(),v_table.end(),tmp);
                //if (iter != v_table.end()) (*iter).val+=tmp.val; //sumando aportes de distintos elementos al mismo nodo
                //else v_table.emplace_back(tmp);  //sino esta el nodo se agrega
            }
        }
    }
    uint32_t nrow = mesh.m_dofxnode*mesh.m_nodes_number; //numero de nodos x dofxnodo
    return Sparse(table,nrow,nrow);
}
//***********************************************************************************************************************
template<typename T,MATRIX_TYPE matrix_type>
std::vector<indexs_val<T>> index_val_table(const rectangular_mesh<ELEMENT_TYPE::QUAD4> &mesh,const Matrix<T,2,matrix_type> &Kelem,
                                           size_t elem_indx){
    assert(elem_indx < mesh.m_element_connect.rows());
    uint32_t nodxelem = mesh.m_element_connect.cols(); //numero nodos por elemento
    uint32_t dofxelem = mesh.m_dofxnode*nodxelem; //numero de dof por elemento
    assert(dofxelem == Kelem.rows());

    std::vector<indexs_val<T>> v_table;
    v_table.reserve(dofxelem*dofxelem); //TODO implementar eso para casos simetrico y hermitico

    std::vector<uint32_t> vrowsdofsxelems;
    vrowsdofsxelems.reserve(dofxelem);
    if (mesh.m_dofxnode == 2){
        for (size_t i = 0; i < nodxelem; ++i){
            vrowsdofsxelems.push_back(2*mesh.m_element_connect(elem_indx,i));
            vrowsdofsxelems.push_back(2*mesh.m_element_connect(elem_indx,i) + 1);
        }
    }else{
        throw ("index_val_table: error solo implementado para 2 dof por nodo");
    }
    for (size_t i = 0; i < vrowsdofsxelems.size(); ++i){
        for (size_t j = 0; j < vrowsdofsxelems.size(); ++j){
            indexs_val<T> tmp(vrowsdofsxelems[i],vrowsdofsxelems[j],Kelem(i,j));
            v_table.emplace_back(tmp);
        }
    }
    return v_table;
}
//inline Sparse_MatComplexd Sparse(const std::vector<uint32_t> &vrows,const std::vector<uint32_t> &vcols,const std::vector<complexd> &vvals,
//                          uint32_t nrow,uint32_t ncol){

//    assert(vrows.size() == vcols.size() && vrows.size() == vvals.size());

//    if (vrows.size() == 0) return Sparse_MatComplexd();

//    std::vector<indexs_val<complexd>> vIndxValsTable(1,indexs_val<complexd>(vrows[0],vcols[0],vvals[0]));
//    for (size_t ii = 1; ii < vrows.size(); ++ii){
//        indexs_val<complexd> tmp(vrows[ii],vcols[ii],vvals[ii]);

//        auto iter = std::find(vIndxValsTable.begin(),vIndxValsTable.end(),tmp);
//        if (iter != vIndxValsTable.end()) iter->val+=tmp.val; //sumando aportes de distintos elementos al mismo nodo
//        else vIndxValsTable.push_back(tmp); //sino esta el nodo se agrega
//    }

//    std::sort(vIndxValsTable.begin(),vIndxValsTable.end()); //se ordena la lista con privilegio de filas

//    uint32_t nvals = vIndxValsTable.size();
//    std::vector<complexd> values(nvals);
//    std::vector<uint32_t> cols(nvals);
//    std::vector<uint32_t> row_start(nrow);
//    std::vector<uint32_t> row_end(nrow);
//    uint32_t tmprow = std::numeric_limits<uint32_t>::max();
//    uint32_t curr_row = 0;
//    for (uint32_t i = 0; i < nvals; ++i){
//        values[i] = vIndxValsTable[i].val;
//        cols[i] = vIndxValsTable[i].col;
//        curr_row = vIndxValsTable[i].row;
//        if (curr_row != tmprow){ //primer elemento de la fila
//            row_start[curr_row] = i;
//            tmprow = curr_row;
//            if (i != 0){ //fin de la fila anterior
//                row_end[curr_row-1] = i;
//            }
//        }
//    }
//    row_end[nrow-1] = nvals;

//    return Sparse_MatComplexd(nrow,ncol,row_start,row_end,cols,values);
//}

//Sparse_MatDoub Sparse(const std::vector<uint32_t> &vrows,const std::vector<uint32_t> &vcols,const std::vector<double> &vvals,
//                      uint32_t nrow,uint32_t ncol){

//    assert(vrows.size() == vcols.size() && vrows.size() == vvals.size());

//    if (vrows.size() == 0) return Sparse_MatDoub();

//    std::vector<indexs_val<double>> vIndxValsTable(1,indexs_val<double>(vrows[0],vcols[0],vvals[0]));
//    for (size_t ii = 1; ii < vrows.size(); ++ii){
//        indexs_val<double> tmp(vrows[ii],vcols[ii],vvals[ii]);

//        auto iter = std::find(vIndxValsTable.begin(),vIndxValsTable.end(),tmp);
//        if (iter != vIndxValsTable.end()) iter->val+=tmp.val; //sumando aportes de distintos elementos al mismo nodo
//        else vIndxValsTable.push_back(tmp); //sino esta el nodo se agrega
//    }

//    std::sort(vIndxValsTable.begin(),vIndxValsTable.end()); //se ordena la lista con privilegio de filas

//    uint32_t nvals = vIndxValsTable.size();
//    std::vector<double> values(nvals);
//    std::vector<uint32_t> cols(nvals);
//    std::vector<uint32_t> row_start(nrow);
//    std::vector<uint32_t> row_end(nrow);
//    uint32_t tmprow = std::numeric_limits<uint32_t>::max();
//    uint32_t curr_row = 0;
//    for (uint32_t i = 0; i < nvals; ++i){
//        values[i] = vIndxValsTable[i].val;
//        cols[i] = vIndxValsTable[i].col;
//        curr_row = vIndxValsTable[i].row;
//        if (curr_row != tmprow){ //primer elemento de la fila
//            row_start[curr_row] = i;
//            tmprow = curr_row;
//            if (i != 0){ //fin de la fila anterior
//                row_end[curr_row-1] = i;
//            }
//        }
//    }
//    row_end[nrow-1] = nvals;
//    return Sparse_MatDoub(nrow,ncol,row_start,row_end,cols,values);
//}


#endif // UTILITIES_H
