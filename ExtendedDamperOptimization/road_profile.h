#ifndef ROAD_PROFILE_H
#define ROAD_PROFILE_H

#include <vector>
#include <random>

class RoadProfile {
public:
    RoadProfile(int steps, double dt);
    double getDisplacement(int t, double delay) const;
    double getVelocity(int t, double delay) const;

private:
    std::vector<double> profile_;
    double dt_;
    std::mt19937 rng_;
    std::normal_distribution<double> noise_dist_;
    void generateProfile();
};

#endif // ROAD_PROFILE_H
