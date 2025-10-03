#pragma once
#ifndef CUBIC_SPLINE_H
#define CUBIC_SPLINE_H

#include <vector>
#include "Point.h"
#include "Spline.h"

namespace Com_Methods {

    class Cubic_Spline : public Spline {
    private:
        std::vector<Point> Points;
        std::vector<double> a, b, c, d;

    public:
        // Построение сплайна по точкам и значениям функции
        void Update_Spline(const std::vector<Point>& Points, 
                          const std::vector<double>& F_Value) override;
        
        //Вычисление значения сплайна и производных в точке P
        void Get_Value(const Point& P, double* Res) const override;
    };

}

#endif