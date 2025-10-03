#pragma once
#ifndef NUMERICAL_DIFFERENTIATOR_H
#define NUMERICAL_DIFFERENTIATOR_H

#include <vector>
#include <cmath>
#include <functional>

namespace Com_Methods {

    class Numerical_Differentiator {
    public:
    
        static double compute_derivative(const std::function<double(double)>& func, 
                                        double x, double epsilon);
        
        static double two_point_forward(const std::function<double(double)>& func, 
                                       double x, double h);
        
        static double three_point_central(const std::function<double(double)>& func, 
                                         double x, double h);
        
        static double five_point_central(const std::function<double(double)>& func, 
                                        double x, double h);
    };

}

#endif