#include "Grid_Generator.h"
#include <vector>
#include <cmath>
#include <stdexcept>

namespace Com_Methods
{

    std::vector<Point> Grid_Generator::uniform_grid(double a, double b, int segments)
    {
        std::vector<Point> grid_points;

        if (segments <= 0)
        {
            throw std::invalid_argument("Число сегментов, не может быть отрицательным!!!");
        }

        double step = (b - a) / segments;

        for (int i = 0; i <= segments; i++)
        {
            double x = a + i * step;
            grid_points.push_back(Point(x, 0, 0));
        }

        return grid_points;
    }

    std::vector<Point> Grid_Generator::adaptive_grid(double a, double b, int segments, double r)
    {
        std::vector<Point> grid_points;

        if (fabs(r - 1.0) < 1e-10)
        {
            double step = (b - a) / segments;
            for (int i = 0; i <= segments; i++)
            {
                grid_points.push_back(Point(a + i * step, 0, 0));
            }
            return grid_points;
        }

        double h1 = (b - a) * (1 - r) / (1 - pow(r, segments));
        grid_points.push_back(Point(a, 0, 0));

        double current = a;
        for (int i = 0; i < segments; i++)
        {
            current += h1 * pow(r, i);
            grid_points.push_back(Point(current, 0, 0));
        }

        return grid_points;
    }

}
