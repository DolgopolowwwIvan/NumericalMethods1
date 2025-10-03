#include "Numerical_Differentiator.h"
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <cmath>

namespace Com_Methods {

    double Numerical_Differentiator::two_point_forward(const std::function<double(double)>& func, 
                                                     double x, double h) {
        return (func(x + h) - func(x)) / h;
    }

    double Numerical_Differentiator::three_point_central(const std::function<double(double)>& func, 
                                                       double x, double h) {
        return (func(x + h) - func(x - h)) / (2.0 * h);
    }

    double Numerical_Differentiator::five_point_central(const std::function<double(double)>& func, 
                                                      double x, double h) {
        return (-func(x + 2.0*h) + 8.0*func(x + h) - 8.0*func(x - h) + func(x - 2.0*h)) / (12.0 * h);
    }

    double Numerical_Differentiator::compute_derivative(const std::function<double(double)>& func, 
                                                      double x, double epsilon) {
        double h = 0.01;
        double derivative = two_point_forward(func, x, h);
        
        std::cout << "Метод: двухточечная разность" << std::endl;
        std::cout << "Шаг: h = " << h << std::endl;
        std::cout << "Вычисленное значение: " << derivative << std::endl;
        
        return derivative;
    }

}