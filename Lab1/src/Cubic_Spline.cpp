#include "Cubic_Spline.h"
#include <vector>
#include <cmath>
#include <stdexcept>

namespace Com_Methods {

    void Cubic_Spline::Update_Spline(const std::vector<Point>& Nodes, 
                                                     const std::vector<double>& Node_Values) {
    
        if (Nodes.size() < 2) {
            throw std::invalid_argument("Для сплайна нужно минимум 2 точки");
        }
        if (Nodes.size() != Node_Values.size()) {
            throw std::invalid_argument("Число точек и значений должно совпадать");
        }

        Points = Nodes;
        int segment_count = Points.size() - 1;

        a.resize(segment_count);
        b.resize(segment_count);
        c.resize(segment_count);
        d.resize(segment_count);

        std::vector<double> h(segment_count);    
        std::vector<double> alpha(segment_count); 
        std::vector<double> l(segment_count + 1); 
        std::vector<double> mu(segment_count);    
        std::vector<double> z(segment_count + 1); 

        
        for (int i = 0; i < segment_count; i++) {
            h[i] = Points[i + 1].x() - Points[i].x();
        }

        for (int i = 1; i < segment_count; i++) {
            alpha[i] = 3.0 * ((Node_Values[i + 1] - Node_Values[i]) / h[i] - 
                             (Node_Values[i] - Node_Values[i - 1]) / h[i - 1]);
        }

        
        l[0] = 1.0;
        mu[0] = 0.0;
        z[0] = 0.0;

        for (int i = 1; i < segment_count; i++) {
            l[i] = 2.0 * (Points[i + 1].x() - Points[i - 1].x()) - h[i - 1] * mu[i - 1];
            mu[i] = h[i] / l[i];
            z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
        }

        l[segment_count] = 1.0;
        z[segment_count] = 0.0;

        for (int i = segment_count - 1; i >= 0; i--) {
            c[i] = z[i] - mu[i] * c[i + 1];
        }

        for (int i = 0; i < segment_count; i++) {
            a[i] = Node_Values[i];
            b[i] = (Node_Values[i + 1] - Node_Values[i]) / h[i] - 
                   h[i] * (c[i + 1] + 2.0 * c[i]) / 3.0;
            d[i] = (c[i + 1] - c[i]) / (3.0 * h[i]);
        }
    }

    void Cubic_Spline::Get_Value(const Point& Query_Point, double* Output) const {
        const double eps = 1e-10;
        int segment_count = Points.size() - 1;
        double x = Query_Point.x();

        for (int i = 0; i < segment_count; i++) {
            double x0 = Points[i].x();
            double x1 = Points[i + 1].x();

            if ((x >= x0 - eps && x <= x1 + eps)) {
                double dx = x - x0;
                Output[0] = a[i] + b[i] * dx + c[i] * dx * dx + d[i] * dx * dx * dx;
                Output[1] = b[i] + 2.0 * c[i] * dx + 3.0 * d[i] * dx * dx;
                Output[2] = 2.0 * c[i] + 6.0 * d[i] * dx;
                return;
            }
        }

        throw std::runtime_error("Точка находится вне области определения сплайна");
    }

}
