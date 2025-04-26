#ifndef VEHICLE_DYNAMICS_H
#define VEHICLE_DYNAMICS_H
#include <vector>
#include <random>
#include "road_profile.h"

class VehicleDynamics {
public:
    VehicleDynamics();
    double evaluateObjective(const std::vector<double>& params) const;
private:
    //// Half-car
    double sprung_mass_ = 1200.0;   // kg (total vehicle mass)
    double unsprung_mass_ = 50.0;   // kg (per wheel)
    double tire_stiffness_ = 200000.0;  // N/m
    double cg_height_ = 0.5;        // m (center of gravity height)
    double gravity_ = 9.81;         // m/s^2
    double dt_ = 0.002;             // Time step (s)
    int sim_steps_ = 5000;          // Simulation steps
    RoadProfile road_profile_;      // Road profile generator
    //// State variables: [z_s, theta, phi, z_uf, z_ur, v_s, omega_y, omega_x, v_uf, v_ur]
    std::vector<double> initial_state_ = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double computeDampingForce(double v_rel, double C_c, double C_r, double blowoff, double knee_point,
                              double hysteresis, double temp_coeff, double damper_temp) const;
    void simulateHalfCar(const std::vector<double>& params,
                        std::vector<double>& accel, std::vector<double>& disp,
                        std::vector<double>& pitch, std::vector<double>& roll,
                        std::vector<double>& tire_force, double& max_stroke) const;
    double computeComfort(const std::vector<double>& accel) const;
    double computeVibration(const std::vector<double>& disp) const;
    double computeHandling(const std::vector<double>& pitch, const std::vector<double>& roll,
                          const std::vector<double>& tire_force) const;
    double computeConstraints(double max_stroke, double stroke_limit) const;
};

#endif // VEHICLE_DYNAMICS_H
