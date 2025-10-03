#pragma once
#ifndef GRID_GENERATOR_H
#define GRID_GENERATOR_H

#include <vector>
#include "Point.h"

namespace Com_Methods {

    class Grid_Generator {
    public:
        // Регулярная равномерная сетка
        static std::vector<Point> uniform_grid(double a, double b, int segments);
        
        // Адаптивная сетка с коэффициентом разрядки
        static std::vector<Point> adaptive_grid(double a, double b, int segments, double r);
    };

}

#endif