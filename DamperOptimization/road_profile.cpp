#include "road_profile.h"
#include <cmath>

RoadProfile::RoadProfile(int steps, double dt)
    : profile_(steps), dt_(dt), rng_(std::random_device{}()), noise_dist_(0.0, 0.01) {
    generateProfile();
}

void RoadProfile::generateProfile() {
    // Simple ISO 8608-like roughness (random with frequency weighting)
    double freq = 0.1; // Base frequency (Hz)
    for (size_t i = 0; i < profile_.size(); ++i) {
        double t = i * dt_;
        profile_[i] = noise_dist_(rng_) * (1.0 / (1.0 + freq * t));
        if (i > 0) {
            profile_[i] += 0.95 * profile_[i - 1]; // Smooth transitions
        }
    }
}

double RoadProfile::getDisplacement(int t, double delay) const {
    int idx = static_cast<int>((t * dt_ + delay) / dt_);
    idx = std::max(0, std::min(static_cast<int>(profile_.size()) - 1, idx));
    return profile_[idx];
}

double RoadProfile::getVelocity(int t, double delay) const {
    int idx = static_cast<int>((t * dt_ + delay) / dt_);
    idx = std::max(0, std::min(static_cast<int>(profile_.size()) - 2, idx));
    return (profile_[idx + 1] - profile_[idx]) / dt_;
}
