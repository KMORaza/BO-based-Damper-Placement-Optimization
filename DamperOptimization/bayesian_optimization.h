#ifndef BAYESIAN_OPTIMIZATION_H
#define BAYESIAN_OPTIMIZATION_H

#include <vector>
#include <functional>
#include "utils.h"

struct ParameterBounds {
    double min;
    double max;
};
class BayesianOptimization {
public:
    BayesianOptimization(const std::vector<ParameterBounds>& bounds,
                        std::function<double(const std::vector<double>&)> objective);
    void initialize(int n_initial);
    void optimize(int n_iterations);
    std::pair<std::vector<double>, double> getBest() const;
private:
    std::vector<ParameterBounds> bounds_;
    std::function<double(const std::vector<double>&)> objective_;
    std::vector<std::vector<double>> samples_;
    std::vector<double> values_;
    double evaluateAcquisition(const std::vector<double>& x, double best_y) const;
    std::vector<double> optimizeAcquisition(double best_y) const;
};

#endif // BAYESIAN_OPTIMIZATION_H
