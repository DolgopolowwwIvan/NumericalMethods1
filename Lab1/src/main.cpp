#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include "Grid_Generator.h"
#include "Cubic_Spline.h"
#include "Numerical_Differentiator.h"

using namespace Com_Methods;

void test_grid_generator()
{
    double a = 0.05, b = 0.30;
    int segments = 5;
    double ratio = 1.2;

    std::cout << "=== Тест генератора сеток ===" << std::endl;
    std::cout << "Отрезок: [" << a << ", " << b << "]" << std::endl;
    std::cout << "Сегментов: " << segments << std::endl;
    std::cout << "Коэффициент разрядки: " << ratio << std::endl
              << std::endl;

    // Тест равномерной сетки
    try
    {
        auto uniform = Grid_Generator::uniform_grid(a, b, segments);
        std::cout << "Равномерная сетка:" << std::endl;
        for (const auto &point : uniform)
        {
            std::cout << std::fixed << std::setprecision(6) << point.x() << " ";
        }
        std::cout << std::endl
                  << std::endl;

        // Тест адаптивной сетки
        auto adaptive = Grid_Generator::adaptive_grid(a, b, segments, ratio);
        std::cout << "Адаптивная сетка (r=" << ratio << "):" << std::endl;
        for (const auto &point : adaptive)
        {
            std::cout << std::fixed << std::setprecision(6) << point.x() << " ";
        }
        std::cout << std::endl;

        // Вывод шагов адаптивной сетки
        std::cout << "Шаги адаптивной сетки:" << std::endl;
        for (size_t i = 1; i < adaptive.size(); i++)
        {
            double step = adaptive[i].x() - adaptive[i - 1].x();
            std::cout << "h" << i << " = " << step << " (отношение: " << (i > 1 ? step / (adaptive[i - 1].x() - adaptive[i - 2].x()) : 1.0) << ")" << std::endl;
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << "Ошибка при генерации сетки: " << e.what() << std::endl;
    }
}

void test_spline_complete() {
    std::cout << "\n=== Исследование сплайна на вложенных сетках ===" << std::endl;
    
    auto func = [](double x) { return std::sin(x); };
    auto func_deriv = [](double x) { return std::cos(x); };
    auto func_deriv2 = [](double x) { return -std::sin(x); };
    
    double a = 0.05, b = 0.30;

    double base_h = (b - a) / 6.0;
    int base_segments = 6;
    
    std::cout << "Сетка h = " << std::fixed << std::setprecision(3) << base_h 
              << " (сегментов: " << base_segments << ")" << std::endl;
    
    std::vector<int> segments = {base_segments, 2 * base_segments, 4 * base_segments};
    std::vector<double> steps = {base_h, base_h/2, base_h/4};
    
    // Контрольные точки 
    std::vector<double> test_points = {
        0.06, 0.068, 0.087, 0.106, 0.125, 
        0.144, 0.163, 0.182, 0.201, 0.220
    };
    
    for (size_t grid_idx = 0; grid_idx < segments.size(); grid_idx++) {
        std::cout << "\nСетка h = " << std::fixed << std::setprecision(3) << steps[grid_idx] 
                  << " (сегментов: " << segments[grid_idx] << ")" << std::endl;
        
        auto nodes = Grid_Generator::uniform_grid(a, b, segments[grid_idx]);
        std::vector<double> values;
        for (const auto& node : nodes) {
            values.push_back(func(node.x()));
        }
        
        Cubic_Spline spline;
        spline.Update_Spline(nodes, values);
        
        // Заголовок таблицы
        std::cout << "x\t\tf(x)\t\tS(x)\t\tОшибка\t\tf'(x)\t\tS'(x)\t\tОшибка'\t\tf''(x)\t\tS''(x)\t\tОшибка''" << std::endl;
        std::cout << std::string(120, '-') << std::endl;
        
        double max_error[3] = {0, 0, 0};
        
        for (double x : test_points) {
            Point p(x, 0, 0);
            double spline_result[3];
            spline.Get_Value(p, spline_result);
            
            double analytic_f = func(x);
            double analytic_f1 = func_deriv(x);  
            double analytic_f2 = func_deriv2(x); 
            
            // Погрешности
            double errors[3] = {
                std::abs(analytic_f - spline_result[0]),
                std::abs(analytic_f1 - spline_result[1]), 
                std::abs(analytic_f2 - spline_result[2])
            };
            
            // Обновляем максимумы
            for (int i = 0; i < 3; i++) {
                if (errors[i] > max_error[i]) max_error[i] = errors[i];
            }
            
            // Вывод строки таблицы
            std::cout << std::fixed << std::setprecision(6) << x << "\t"
                      << std::setprecision(8) << analytic_f << "\t"
                      << spline_result[0] << "\t"
                      << std::scientific << std::setprecision(2) << errors[0] << "\t"
                      << std::fixed << std::setprecision(6) << analytic_f1 << "\t"
                      << spline_result[1] << "\t"
                      << std::scientific << std::setprecision(2) << errors[1] << "\t"
                      << std::fixed << std::setprecision(6) << analytic_f2 << "\t"
                      << spline_result[2] << "\t"
                      << std::scientific << std::setprecision(2) << errors[2] << std::endl;
        }
        
        // Вывод максимальных погрешностей
        std::cout << "\nМаксимальные погрешности:" << std::endl;
        std::cout << "Функции = " << std::scientific << std::setprecision(2) << max_error[0]
                  << ", производной 1-ой = " << max_error[1] 
                  << ", производной 2-ой = " << max_error[2] << std::endl;
    }
}

void test_numerical_differentiation() {
    std::cout << "\n=== Численное дифференцирование (пункт 4) ===" << std::endl;
    
    auto func = [](double x) { return std::sin(x); };
    auto func_deriv = [](double x) { return std::cos(x); };
    

    double center_point = (0.05 + 0.30) / 2.0;
    double exact_derivative = func_deriv(center_point);
    
    std::cout << "Центральная точка: x = " << std::fixed << std::setprecision(3) << center_point << std::endl;
    std::cout << "Точное значение: f'(" << center_point << ") = " << std::fixed << std::setprecision(6) << exact_derivative << std::endl;
    std::cout << "Требуемая точность: ε = 0.01" << std::endl << std::endl;
    
    // Заголовок таблицы
    std::cout << "h\t\tПравая разн.\tЛевая разн.\tЦентр. разн.\tОшибка прав.\tОшибка лев.\tОшибка центр." << std::endl;
    std::cout << std::string(100, '-') << std::endl;

    std::vector<double> test_h = {0.7, 0.5, 0.11, 0.01, 0.005, 0.0025, 0.00125};
    
    for (double h : test_h) {
        // Вычисляем разностные производные
        double right_diff = Numerical_Differentiator::two_point_forward(func, center_point, h);
        double left_diff = (func(center_point) - func(center_point - h)) / h;  // левая разность
        double central_diff = Numerical_Differentiator::three_point_central(func, center_point, h);
        
        // Вычисляем ошибки
        double error_right = std::abs(right_diff - exact_derivative);
        double error_left = std::abs(left_diff - exact_derivative);
        double error_central = std::abs(central_diff - exact_derivative);
        
        // Вывод строки таблицы
        std::cout << std::fixed << std::setprecision(3) << h << "\t\t"
                  << std::setprecision(6) << right_diff << "\t"
                  << left_diff << "\t"
                  << central_diff << "\t"
                  << std::scientific << std::setprecision(2) << error_right << "\t"
                  << error_left << "\t"
                  << error_central << std::endl;
        
    }
}

int main()
{
    test_grid_generator();            
    test_spline_complete();           
    test_numerical_differentiation(); 
}