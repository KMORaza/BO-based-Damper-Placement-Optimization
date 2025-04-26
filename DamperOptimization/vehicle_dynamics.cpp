#include "vehicle_dynamics.h"
#include "utils.h"
#include <cmath>
#include <numeric>
#include <algorithm>

VehicleDynamics::VehicleDynamics() : road_profile_(sim_steps_, dt_) {}
double VehicleDynamics::evaluateObjective(const std::vector<double>& params) const {
    double C_c = params[0];           // Compression damping (Ns/m)
    double C_r = params[1];           // Rebound damping (Ns/m)
    double blowoff = params[2];       // Blow-off threshold (m/s)
    double gas_pressure = params[3];  // Gas preload (bar)
    double motion_ratio = params[4];  // Motion ratio
    double inclination = params[5];   // Inclination angle (degrees)
    double knee_point = params[6];    // Damping profile knee point (m/s)
    double stroke_limit = params[7];  // Stroke limit (m)
    double spring_preload = params[8]; // Spring preload (N)
    double tire_damping = params[9];  // Tire damping (Ns/m)
    double anti_roll_stiffness = params[10]; // Anti-roll bar stiffness (N/rad)
    double cos_theta = std::cos(inclination * M_PI / 180.0);
    double effective_C_c = C_c * motion_ratio * cos_theta;
    double effective_C_r = C_r * motion_ratio * cos_theta;
    std::vector<double> accel, disp, pitch, roll, tire_force;
    double max_stroke = 0.0;
    std::vector<double> adjusted_params = {effective_C_c, effective_C_r, blowoff, gas_pressure,
                                          motion_ratio, inclination, knee_point, stroke_limit,
                                          spring_preload, tire_damping, anti_roll_stiffness};
    simulateHalfCar(adjusted_params, accel, disp, pitch, roll, tire_force, max_stroke);
    double J_comfort = computeComfort(accel);
    double J_vibration = computeVibration(disp);
    double J_handling = computeHandling(pitch, roll, tire_force);
    double J_constraints = computeConstraints(max_stroke, stroke_limit);
    double w1 = 0.35, w2 = 0.25, w3 = 0.25, w4 = 0.15;
    return w1 * J_comfort + w2 * J_vibration + w3 * J_handling + w4 * J_constraints;
}
double VehicleDynamics::computeDampingForce(double v_rel, double C_c, double C_r, double blowoff, double knee_point) const {
    //// Nonlinear damping (digressive profile)
    double C_base = (v_rel >= 0) ? C_r : C_c;
    if (std::abs(v_rel) > blowoff) {
        C_base *= 0.5; // Blow-off reduction
    }
    if (std::abs(v_rel) > knee_point) {
        C_base *= knee_point / std::abs(v_rel); // Digressive taper
    }
    return -C_base * v_rel;
}
void VehicleDynamics::simulateHalfCar(const std::vector<double>& params,
                                     std::vector<double>& accel, std::vector<double>& disp,
                                     std::vector<double>& pitch, std::vector<double>& roll,
                                     std::vector<double>& tire_force, double& max_stroke) const {
    double C_c = params[0], C_r = params[1], blowoff = params[2], gas_pressure = params[3],
           motion_ratio = params[4], inclination = params[5], knee_point = params[6],
           spring_preload = params[8], tire_damping = params[9], anti_roll_stiffness = params[10];
    //// Gas pressure effect
    double adjusted_spring_stiffness = spring_stiffness_ * (1.0 + gas_pressure * 0.01);
    std::vector<double> state = initial_state_;
    accel.resize(sim_steps_);
    disp.resize(sim_steps_);
    pitch.resize(sim_steps_);
    roll.resize(sim_steps_);
    tire_force.resize(sim_steps_);
    max_stroke = 0.0;
    for (int t = 0; t < sim_steps_; ++t) {
        double z_s = state[0], theta = state[1], phi = state[2], z_uf = state[3], z_ur = state[4],
               v_s = state[5], omega_y = state[6], omega_x = state[7], v_uf = state[8], v_ur = state[9];
        //// Road inputs (front and rear)
        double z_rf = road_profile_.getDisplacement(t, 0.0);
        double z_rr = road_profile_.getDisplacement(t, wheelbase_ / 10.0);
        //// Suspension displacements
        double z_sf = z_s + wheelbase_ / 2 * theta - track_width_ / 2 * phi;
        double z_sr = z_s - wheelbase_ / 2 * theta + track_width_ / 2 * phi;
        //// Relative velocities
        double v_rel_f = v_s + wheelbase_ / 2 * omega_y - track_width_ / 2 * omega_x - v_uf;
        double v_rel_r = v_s - wheelbase_ / 2 * omega_y + track_width_ / 2 * omega_x - v_ur;
        //// Damping forces
        double F_damper_f = computeDampingForce(v_rel_f, C_c, C_r, blowoff, knee_point) * motion_ratio;
        double F_damper_r = computeDampingForce(v_rel_r, C_c, C_r, blowoff, knee_point) * motion_ratio;
        //// Spring forces
        double F_spring_f = -(adjusted_spring_stiffness * (z_sf - z_uf) + spring_preload);
        double F_spring_r = -(adjusted_spring_stiffness * (z_sr - z_ur) + spring_preload);
        //// Tire forces
        double F_tire_f = -tire_stiffness_ * (z_uf - z_rf) - tire_damping * (v_uf - road_profile_.getVelocity(t, 0.0));
        double F_tire_r = -tire_stiffness_ * (z_ur - z_rr) - tire_damping * (v_ur - road_profile_.getVelocity(t, wheelbase_ / 10.0));
        //// Anti-roll bar force (proportional to roll angle difference)
        double F_roll = -anti_roll_stiffness * phi;
        //// Equations of motion
        double I_y = sprung_mass_ * wheelbase_ * wheelbase_ / 12; // Pitch inertia
        double I_x = sprung_mass_ * track_width_ * track_width_ / 12; // Roll inertia
        double a_s = (F_spring_f + F_spring_r + F_damper_f + F_damper_r) / sprung_mass_;
        double alpha_y = ((F_spring_f + F_damper_f) * wheelbase_ / 2 - (F_spring_r + F_damper_r) * wheelbase_ / 2) / I_y;
        double alpha_x = ((F_spring_r - F_spring_f) * track_width_ / 2 + F_roll) / I_x;
        double a_uf = (-F_spring_f - F_damper_f + F_tire_f) / unsprung_mass_;
        double a_ur = (-F_spring_r - F_damper_r + F_tire_r) / unsprung_mass_;
        //// Euler integration
        state[5] += a_s * dt_;
        state[6] += alpha_y * dt_;
        state[7] += alpha_x * dt_;
        state[8] += a_uf * dt_;
        state[9] += a_ur * dt_;
        state[0] += state[5] * dt_;
        state[1] += state[6] * dt_;
        state[2] += state[7] * dt_;
        state[3] += state[8] * dt_;
        state[4] += state[9] * dt_;
        //// outputs
        accel[t] = a_s;
        disp[t] = z_s;
        pitch[t] = theta;
        roll[t] = phi;
        tire_force[t] = std::min(F_tire_f, F_tire_r); // Minimum contact force
        max_stroke = std::max(max_stroke, std::max(std::abs(z_sf - z_uf), std::abs(z_sr - z_ur)));
    }
}
double VehicleDynamics::computeComfort(const std::vector<double>& accel) const {
    //// RMS acceleration with frequency weighting (ISO 2631-1)
    double sum_sq = 0.0;
    for (double a : accel) {
        sum_sq += a * a * (std::abs(a) < 1.0 ? 1.0 : 1.4); // Weight higher accelerations
    }
    return std::sqrt(sum_sq / accel.size());
}
double VehicleDynamics::computeVibration(const std::vector<double>& disp) const {
    //// Power spectral density (approximate)
    std::vector<double> psd = computePSD(disp, dt_, 0.5, 20.0);
    double sum_psd = std::accumulate(psd.begin(), psd.end(), 0.0);
    return sum_psd / psd.size();
}
double VehicleDynamics::computeHandling(const std::vector<double>& pitch, const std::vector<double>& roll,
                                      const std::vector<double>& tire_force) const {
    //// Combined metric: max pitch, roll, and tire load variation
    double max_pitch = *std::max_element(pitch.begin(), pitch.end(), [](double a, double b) { return std::abs(a) < std::abs(b); });
    double max_roll = *std::max_element(roll.begin(), roll.end(), [](double a, double b) { return std::abs(a) < std::abs(b); });
    double min_tire_force = *std::min_element(tire_force.begin(), tire_force.end());
    return 0.4 * std::abs(max_pitch) + 0.4 * std::abs(max_roll) + 0.2 * (min_tire_force < 0 ? -min_tire_force : 0);
}
double VehicleDynamics::computeConstraints(double max_stroke, double stroke_limit) const {
    //// Penalty for stroke limit violation
    double penalty = 0.0;
    if (max_stroke > stroke_limit) {
        penalty += 1000.0 * (max_stroke - stroke_limit);
    }
    return penalty;
}
