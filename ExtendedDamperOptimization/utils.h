#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <cmath>

inline double clamp(double x, double min, double max) {
    return std::max(min, std::min(max, x));
}
std::vector<double> computePSD(const std::vector<double>& signal, double dt, double f_min, double f_max);

#endif // UTILS_H
