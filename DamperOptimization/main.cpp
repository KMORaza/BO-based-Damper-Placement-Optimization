#include "bayesian_optimization.h"
#include "vehicle_dynamics.h"
#include <iostream>
#include <vector>

int main() {
    std::cout << "Starting DamperOptimization..." << std::endl;
    std::vector<ParameterBounds> bounds = {
        {500.0, 5000.0},   // C_c (Ns/m)
        {1000.0, 8000.0},  // C_r (Ns/m)
        {0.5, 2.0},        // Blow-off threshold (m/s)
        {5.0, 20.0},       // Gas preload pressure (bar)
        {0.5, 1.5},        // Motion ratio
        {0.0, 30.0},       // Inclination angle (degrees)
        {0.0, 2.0},        // Damping profile knee point (m/s)
        {0.05, 0.1},       // Damper stroke limit (m)
        {0.0, 1000.0},     // Spring preload (N)
        {50.0, 500.0},     // Tire damping (Ns/m)
        {1000.0, 10000.0}  // Anti-roll bar stiffness (N/rad)
    };
    std::cout << "Parameter bounds initialized." << std::endl;
    VehicleDynamics simulator;
    std::cout << "VehicleDynamics simulator initialized." << std::endl;
    BayesianOptimization optimizer(bounds, [&simulator](const std::vector<double>& params) {
        std::cout << "Evaluating objective for parameters..." << std::endl;
        double result = simulator.evaluateObjective(params);
        std::cout << "Objective value: " << result << std::endl;
        return result;
    });
    std::cout << "BayesianOptimization initialized." << std::endl;
    std::cout << "Starting initialization..." << std::endl;
    optimizer.initialize(15);
    std::cout << "Initialization complete. Starting optimization..." << std::endl;
    optimizer.optimize(50);
    std::cout << "Optimization complete." << std::endl;
    auto result = optimizer.getBest();
    auto best_params = result.first;
    auto best_value = result.second;
    std::cout << "Best Objective Value: " << best_value << "\nBest Parameters:\n";
    for (size_t i = 0; i < best_params.size(); ++i) {
        std::string param_names[] = {"C_c", "C_r", "Blow-off", "Gas Pressure", "Motion Ratio",
                                    "Inclination", "Knee Point", "Stroke Limit", "Spring Preload",
                                    "Tire Damping", "Anti-roll Stiffness"};
        std::cout << param_names[i] << ": " << best_params[i] << "\n";
    }
    std::cout << "Program finished successfully." << std::endl;
    return 0;
}
