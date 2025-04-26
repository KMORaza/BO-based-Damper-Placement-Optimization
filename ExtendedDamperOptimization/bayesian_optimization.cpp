#include "bayesian_optimization.h"
#include <random>
#include <cmath>
#include <limits>

GaussianProcess::GaussianProcess(double length_scale, double sigma_f, double sigma_n)
    : length_scale_(length_scale), sigma_f_(sigma_f), sigma_n_(sigma_n) {}
double GaussianProcess::rbfKernel(const std::vector<double>& x1, const std::vector<double>& x2) const {
    double sum_sq = 0.0;
    for (size_t i = 0; i < x1.size(); ++i) {
        sum_sq += (x1[i] - x2[i]) * (x1[i] - x2[i]);
    }
    return sigma_f_ * sigma_f_ * std::exp(-sum_sq / (2.0 * length_scale_ * length_scale_));
}
void GaussianProcess::fit(const std::vector<std::vector<double>>& X, const std::vector<double>& y) {
    X_ = X;
    y_ = Eigen::VectorXd::Map(y.data(), y.size());
    //// Build covariance matrix K
    size_t n = X.size();
    K_.resize(n, n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            K_(i, j) = rbfKernel(X[i], X[j]);
            if (i == j) {
                K_(i, j) += sigma_n_ * sigma_n_;
            }
        }
    }
}
std::pair<double, double> GaussianProcess::predict(const std::vector<double>& x) const {
    if (X_.empty()) {
        return {0.0, sigma_f_ * sigma_f_};
    }
    //// k(x, X) vector
    size_t n = X_.size();
    Eigen::VectorXd k(n);
    for (size_t i = 0; i < n; ++i) {
        k(i) = rbfKernel(x, X_[i]);
    }
    //// Solve K^-1 * y using Cholesky decomposition
    Eigen::LLT<Eigen::MatrixXd> llt(K_);
    Eigen::VectorXd alpha = llt.solve(y_);
    //// Mean: k^T * K^-1 * y
    double mean = k.dot(alpha);
    //// Variance: sigma_f^2 - k^T * K^-1 * k
    Eigen::VectorXd v = llt.solve(k);
    double variance = sigma_f_ * sigma_f_ - k.dot(v);
    variance = std::max(variance, 1e-6);
    return {mean, variance};
}
BayesianOptimization::BayesianOptimization(const std::vector<ParameterBounds>& bounds,
                                         std::function<double(const std::vector<double>&)> objective)
    : bounds_(bounds), objective_(objective), gp_(1.0, 1.0, 0.1) {}
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
    for (int iter = 0; iter < n_iterations; ++iter) {
        gp_.fit(samples_, values_);
        double best_y = *std::min_element(values_.begin(), values_.end());
        auto next_x = optimizeAcquisition(best_y);
        double y = objective_(next_x);
        samples_.push_back(next_x);
        values_.push_back(y);
    }
}
double BayesianOptimization::evaluateAcquisition(const std::vector<double>& x, double best_y) const {
    std::pair<double, double> pred = gp_.predict(x);
    double mean = pred.first;
    double variance = pred.second;
    double sigma = std::sqrt(variance);
    if (sigma < 1e-6) {
        return 0.0;
    }
    double z = (best_y - mean) / sigma;
    return (best_y - mean) * 0.5 * (1.0 + std::erf(z / std::sqrt(2.0))) +
           sigma * std::exp(-z * z / 2.0) / std::sqrt(2.0 * M_PI);
}
std::vector<double> BayesianOptimization::optimizeAcquisition(double best_y) const {
    std::random_device rd;
    std::mt19937 rng(rd());
    std::vector<double> best_x(bounds_.size());
    double best_ei = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < 1000; ++i) {
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
