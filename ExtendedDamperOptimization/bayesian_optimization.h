#ifndef BAYESIAN_OPTIMIZATION_H
#define BAYESIAN_OPTIMIZATION_H
#include <vector>
#include <functional>
#include <Eigen/Dense>

struct ParameterBounds {
    double min;
    double max;
};
class GaussianProcess {
public:
    GaussianProcess(double length_scale = 1.0, double sigma_f = 1.0, double sigma_n = 0.1);
    void fit(const std::vector<std::vector<double>>& X, const std::vector<double>& y);
    std::pair<double, double> predict(const std::vector<double>& x) const;
private:
    double length_scale_; // Kernel length scale
    double sigma_f_;      // Signal variance
    double sigma_n_;      // Noise variance
    Eigen::MatrixXd K_;   // Covariance matrix
    Eigen::VectorXd y_;   // Training outputs
    std::vector<std::vector<double>> X_; // Training inputs
    double rbfKernel(const std::vector<double>& x1, const std::vector<double>& x2) const;
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
    GaussianProcess gp_;
    double evaluateAcquisition(const std::vector<double>& x, double best_y) const;
    std::vector<double> optimizeAcquisition(double best_y) const;
};

#endif // BAYESIAN_OPTIMIZATION_H
