#ifndef SPARSE_MATRIX_UTILS_H
#define SPARSE_MATRIX_UTILS_H

#include "matrix.h"

typedef Matrix<double,2,Matrix_Type::GEN,Matrix_Storage_Scheme::CSR3> SPMatDoub;
typedef Matrix<complex,2,Matrix_Type::GEN,Matrix_Storage_Scheme::CSR3> SPMatComplx;

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

SPMatComplx Sparse(const std::vector<uint32_t> &vrows,const std::vector<uint32_t> &vcols,const std::vector<complex> &vvals,
                 uint32_t nrow,uint32_t ncol){

    assert(vrows.size() == vcols.size() && vrows.size() == vvals.size());

    if (vrows.size() == 0) return SPMatComplx();

    std::vector<indexs_val<complex>> vIndxValsTable(1,indexs_val<complex>(vrows[0],vcols[0],vvals[0]));
    for (size_t ii = 1; ii < vrows.size(); ++ii){
        indexs_val<complex> tmp(vrows[ii],vcols[ii],vvals[ii]);

        auto iter = std::find(vIndxValsTable.begin(),vIndxValsTable.end(),tmp);
        if (iter != vIndxValsTable.end()) iter->val+=tmp.val;
        else vIndxValsTable.push_back(tmp);
    }

    std::sort(vIndxValsTable.begin(),vIndxValsTable.end());

    size_t nvals = vIndxValsTable.size();
    std::vector<complex> values(nvals);
    std::vector<uint32_t> cols(nvals);
    std::vector<uint32_t> rowIndex(1,0);
    uint32_t tmprow = 0;
    uint32_t curr_row = 0;
    for (size_t i = 0; i < nvals; ++i){
        values[i] = vIndxValsTable[i].val;
        cols[i] = vIndxValsTable[i].col;
        curr_row = vIndxValsTable[i].row;
        if (curr_row != tmprow){
            rowIndex.push_back(curr_row);
            tmprow = curr_row;
        }
    }
    rowIndex.push_back(nvals);
    return SPMatComplx(values,cols,rowIndex);
}

SPMatDoub Sparse(const std::vector<uint32_t> &vrows,const std::vector<uint32_t> &vcols,const std::vector<double> &vvals,
                 uint32_t nrow,uint32_t ncol){

    assert(vrows.size() == vcols.size() && vrows.size() == vvals.size());

    if (vrows.size() == 0) return SPMatDoub();

    std::vector<indexs_val<double>> vIndxValsTable(1,indexs_val<double>(vrows[0],vcols[0],vvals[0]));
    for (size_t ii = 1; ii < vrows.size(); ++ii){
        indexs_val<double> tmp(vrows[ii],vcols[ii],vvals[ii]);

        auto iter = std::find(vIndxValsTable.begin(),vIndxValsTable.end(),tmp);
        if (iter != vIndxValsTable.end()) iter->val+=tmp.val;
        else vIndxValsTable.push_back(tmp);
    }

    std::sort(vIndxValsTable.begin(),vIndxValsTable.end());

    size_t nvals = vIndxValsTable.size();
    std::vector<double> values(nvals);
    std::vector<uint32_t> cols(nvals);
    std::vector<uint32_t> rowIndex(1,0);
    uint32_t tmprow = 0;
    uint32_t curr_row = 0;
    for (size_t i = 0; i < nvals; ++i){
        values[i] = vIndxValsTable[i].val;
        cols[i] = vIndxValsTable[i].col;
        curr_row = vIndxValsTable[i].row;
        if (curr_row != tmprow){
            rowIndex.push_back(curr_row);
            tmprow = curr_row;
        }
    }
    rowIndex.push_back(nvals);
    return SPMatDoub(values,cols,rowIndex);
}


#endif // SPARSE_MATRIX_UTILS_H
