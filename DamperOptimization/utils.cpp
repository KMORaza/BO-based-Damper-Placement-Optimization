#include "utils.h"
#include <vector>
#include <cmath>

std::vector<double> computePSD(const std::vector<double>& signal, double dt, double f_min, double f_max) {
    //// PSD estimation using autocorrelation and cosine transform approximation
    int N = signal.size();
    std::vector<double> autocorr(N, 0.0);
    double mean = 0.0;
    //// signal mean
    for (double x : signal) {
        mean += x;
    }
    mean /= N;

    //// autocorrelation
    int max_lag = N / 2;
    for (int lag = 0; lag < max_lag; ++lag) {
        for (int i = 0; i < N - lag; ++i) {
            autocorr[lag] += (signal[i] - mean) * (signal[i + lag] - mean);
        }
        autocorr[lag] /= (N - lag);
    }
    //// Approximate PSD using cosine transform over target frequency range
    std::vector<double> psd;
    double df = 1.0 / (N * dt); // Frequency resolution
    int n_freqs = static_cast<int>((f_max - f_min) / df) + 1;
    for (int k = 0; k < n_freqs; ++k) {
        double freq = f_min + k * df;
        if (freq > f_max) break;
        double power = 0.0;
        for (int lag = 0; lag < max_lag; ++lag) {
            double w = 2.0 * M_PI * freq * lag * dt;
            //// Blackman window to reduce noise
            double window = 0.42 - 0.5 * std::cos(2.0 * M_PI * lag / max_lag) + 0.08 * std::cos(4.0 * M_PI * lag / max_lag);
            power += autocorr[lag] * window * std::cos(w);
        }
        power *= 2.0 * dt;
        psd.push_back(std::max(power, 0.0));
    }
    return psd;
}
