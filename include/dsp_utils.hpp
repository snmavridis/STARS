///////////////////////////////////////////////////////////////////////////////
//FILE: dsp_utils.hpp
//
//Contains the dsp classes
//
//250412 Created by Savas N Mavridis
///////////////////////////////////////////////////////////////////////////////
#ifndef DSP_UTILS_H
#define DSP_UTILS_H

#include <deque> 
#include <vector>
#include <complex>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>
#include <iostream>
#include <iomanip> // For setprecision in dump functions
#include "kiss_fft.h"


using Complex = std::complex<double>;
using ComplexVec = std::vector<Complex>;
using ComplexMat = std::vector<ComplexVec>;
using DoubleVec = std::vector<double>;
using DoubleMat = std::vector<DoubleVec>;
using IntVec = std::vector<int>;
using IntMat = std::vector<IntVec>;

struct AoA_Estimate {
    int range_bin;
    int doppler_bin;
    double azimuth_est;      // radians
    double elevation_est;    // radians
    double azimuth_truth;    // radians
    double elevation_truth;  // radians
    double power;
};

ComplexVec fft(const ComplexVec& input);
ComplexVec kiss_fft_wrapper(const ComplexVec& input);
ComplexVec quantize_adc(
    const ComplexVec& input,
    int num_bits,
    double v_fullscale
);

void dump_range_profile(const ComplexVec& sweep, const std::string& filename);
ComplexVec ifft(const ComplexVec& input);
ComplexVec matched_filter_fft(const ComplexVec& signal, const ComplexVec& template_chirp);
void dump_rdm_binary(const std::vector<std::vector<int>>& binary_rdm, const std::string& filename);
void dump_rdm_power(const DoubleMat& rdm, const std::string& filename);
void dump_aoa_estimates(const std::vector<AoA_Estimate>& estimates, const std::string& filename);

class TemporalIntegrator {
private:
    std::deque<IntMat> history;
    int max_history = 7;   // N
    int threshold = 4;     // M

public:
    TemporalIntegrator(int N = 1, int M = 1);  // Constructor
    void add_detection_map(const IntMat& binary_map);
    IntMat integrate() const;
    void clear();
};

std::vector<AoA_Estimate> estimate_monopulse_aoa(
    const ComplexMat& ch1,
    const ComplexMat& ch2,
    const ComplexMat& ch3,
    const ComplexMat& ch4,
    const IntMat& detections,
    const DoubleMat& rdm_power,  // NEW
    double lambda,
    double power_threshold_db = -50.0, // dB default
    double los_az_rad=0.0,
    double los_az_el=0.0
);

bool kiss_fft_exec(kiss_fft_cfg cfg, kiss_fft_cpx* fin, kiss_fft_cpx* fout, size_t N, const ComplexVec& input);


struct DspBuffer4Ch {
    // --- Members ---
    ComplexMat ch1, ch2, ch3, ch4;
    int max_sweeps;       // Set via constructor
    int bins_per_sweep;   // Set via constructor

private: // Keep FFT resources private
    kiss_fft_cfg doppler_fft_cfg;
    kiss_fft_cpx* doppler_fft_in;
    kiss_fft_cpx* doppler_fft_out;
    bool doppler_fft_initialized;

public:
    // --- Declarations ---
    DspBuffer4Ch(int num_doppler_bins = 32, int num_range_bins = 256); // Constructor Declaration
    ~DspBuffer4Ch(); // Destructor Declaration
    
    void initialize_doppler_fft(); // Declaration
    void clear();                  // Declaration

    void add_sweep_4ch(const ComplexVec& s1, const ComplexVec& s2, const ComplexVec& s3, const ComplexVec& s4); // Declaration (or keep inline if simple)

    // Keep simple methods inline if preferred
    bool is_ready() const {
        return !ch1.empty() && ((int)ch1.size() == max_sweeps);
    }

    DoubleMat compute_rdm(const ComplexMat& channel); // Declaration
    ComplexMat get_complex_rdm(const ComplexMat& channel); // Declaration

    // Prevent copying (recommended for resource safety)
    DspBuffer4Ch(const DspBuffer4Ch&) = delete;
    DspBuffer4Ch& operator=(const DspBuffer4Ch&) = delete;
};


#endif