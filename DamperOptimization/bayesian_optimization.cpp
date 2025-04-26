#include "bayesian_optimization.h"
#include <Eigen/Dense>
#include <random>
#include <cmath>

namespace mock_gp {
struct GaussianProcess {
    std::vector<std::vector<double>> X_;
    std::vector<double> y_;
    double predict(const std::vector<double>& x, double& var) const {
        double mean = 0.0, var_est = 1.0;
        for (size_t i = 0; i < X_.size(); ++i) {
            double dist = 0.0;
            for (size_t j = 0; j < x.size(); ++j) {
                dist += (x[j] - X_[i][j]) * (x[j] - X_[i][j]);
            }
            mean += y_[i] * std::exp(-dist / 2.0);
        }
        return mean / (X_.size() + 1e-6);
    }
    void update(const std::vector<std::vector<double>>& X, const std::vector<double>& y) {
        X_ = X;
        y_ = y;
    }
};
}
BayesianOptimization::BayesianOptimization(const std::vector<ParameterBounds>& bounds,
                                         std::function<double(const std::vector<double>&)> objective)
    : bounds_(bounds), objective_(objective) {}
void BayesianOptimization::initialize(int n_initial) {
    std::random_device rd;
    std::mt19937 rng(rd());
    samples_.clear();
    values_.clear();
    //// Latin Hypercube Sampling
    for (int i = 0; i < n_initial; ++i) {
        std::vector<double> x(bounds_.size());
        for (size_t j = 0; j < bounds_.size(); ++j) {
            std::uniform_real_distribution<double> dist(bounds_[j].min, bounds_[j].max);
            x[j] = dist(rng);
        }
        samples_.push_back(x);
        values_.push_back(objective_(x));
    }
}
void BayesianOptimization::optimize(int n_iterations) {
    mock_gp::GaussianProcess gp;
    for (int iter = 0; iter < n_iterations; ++iter) {
        gp.update(samples_, values_);
        double best_y = *std::min_element(values_.begin(), values_.end());
        auto next_x = optimizeAcquisition(best_y);
        double y = objective_(next_x);
        samples_.push_back(next_x);
        values_.push_back(y);
    }
}
double BayesianOptimization::evaluateAcquisition(const std::vector<double>& x, double best_y) const {
    mock_gp::GaussianProcess gp;
    gp.update(samples_, values_);
    double mean, var;
    mean = gp.predict(x, var);
    double sigma = std::sqrt(var + 1e-6);
    double z = (best_y - mean) / sigma;
    return (best_y - mean) * 0.5 * (1.0 + std::erf(z / std::sqrt(2.0))) +
           sigma * std::exp(-z * z / 2.0) / std::sqrt(2.0 * M_PI);
}
std::vector<double> BayesianOptimization::optimizeAcquisition(double best_y) const {
    std::random_device rd;
    std::mt19937 rng(rd());
    std::vector<double> best_x(bounds_.size());
    double best_ei = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < 2000; ++i) { // Increased iterations
        std::vector<double> x(bounds_.size());
        for (size_t j = 0; j < bounds_.size(); ++j) {
            std::uniform_real_distribution<double> dist(bounds_[j].min, bounds_[j].max);
            x[j] = dist(rng);
        }
        double ei = evaluateAcquisition(x, best_y);
        if (ei > best_ei) {
            best_ei = ei;
            best_x = x;
        }
    }
    return best_x;
}
std::pair<std::vector<double>, double> BayesianOptimization::getBest() const {
    auto min_it = std::min_element(values_.begin(), values_.end());
    size_t idx = std::distance(values_.begin(), min_it);
    return {samples_[idx], *min_it};
}
