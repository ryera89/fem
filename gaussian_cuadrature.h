#ifndef GAUSSIAN_CUADRATURE_H
#define GAUSSIAN_CUADRATURE_H

#include "element.h"

struct gaussian_cuadrature{
    ELEMENT_TYPE element_type;
    uint32_t nodxelem;
    uint32_t gauss_points_number;
    std::vector<QPointF> gauss_points;
    std::vector<double> weights;

    gaussian_cuadrature(ELEMENT_TYPE etype):element_type(etype){
        switch (etype) {
        case ELEMENT_TYPE::QUAD4:{
            gauss_points_number = 4;
            nodxelem = 4;
            gauss_points = {QPointF(-1/std::sqrt(3),-1/std::sqrt(3)),
                                                 QPointF(1/std::sqrt(3),-1/std::sqrt(3)),
                                                 QPointF(1/std::sqrt(3),1/std::sqrt(3)),
                                                 QPointF(-1/std::sqrt(3),1/std::sqrt(3))};
            weights = {1.0,1.0,1.0,1.0};
            break;
        }
        default:
            throw ("gaussian_cuadrature: error invalid element type");
            //break;
        }
    }
    VecDoub integration_points_function_evaluation(double (*f)(const QPointF &x)){
        VecDoub vres(gauss_points_number);
        std::transform(gauss_points.begin(),gauss_points.end(),weights.begin(),
                       vres.begin(),[&f](auto &p,auto &w){return w*f(p);});
        return vres;
    }
    double integrate(double (*f)(const QPointF &x)) const{
        return std::inner_product(gauss_points.begin(),gauss_points.end(),weights.begin(),0.0,std::plus<double>(),
                                  [&f](const auto &p,const auto &w){return w*f(p);});
    }

    MatDoub integrate(MatDoub (*f)(const QPointF &x)) const{
        return std::inner_product(gauss_points.begin(),gauss_points.end(),weights.begin(),MatDoub(),std::plus<MatDoub>(),
                                  [&f](const auto &p,const auto &w){return w*f(p);});
    }
    template<typename F>
    MatDoub integrate(F fun) const{
        return std::inner_product(gauss_points.begin(),gauss_points.end(),weights.begin(),MatDoub(),std::plus<MatDoub>(),
                                  [&fun](const auto &p,const auto &w){return w*fun(p);});
    }
};


#endif // GAUSSIAN_CUADRATURE_H
