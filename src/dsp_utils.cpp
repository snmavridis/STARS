///////////////////////////////////////////////////////////////////////////////
//FILE: dsp_utils.cpp
//
//Contains DSP utility functions
//
//250411 Savas N Mavridis
///////////////////////////////////////////////////////////////////////////////

#include "dsp_utils.hpp"
#include <vector>
#include <complex>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>
#include <iostream>
#include "kiss_fft.h"

ComplexVec fft(const ComplexVec& input) {
    size_t N = input.size();
    if (N <= 1) return input;

    // Bit-reversal permutation
    ComplexVec output = input; // Start with a copy

    // Calculate log2N once if needed, or use std::log2 within loop s below
    // size_t log2N = static_cast<size_t>(std::round(std::log2(N))); // Ensure N is power of 2 elsewhere

    // Calculate bit-reversed indices and swap elements IN PLACE in the output vector
    for (size_t i = 0; i < N; ++i) {
        size_t j = 0; // Calculate bit-reversed index j for i
        size_t temp_i = i;
        // Efficient bit reversal (assuming N is power of 2)
        for (size_t bit = 0; (1UL << bit) < N; ++bit) { // Use unsigned long for shift safety
            if (temp_i & 1) { // Check the last bit of temp_i
                j |= (N >> 1) >> bit; // Set the corresponding top bit in j
            }
            temp_i >>= 1; // Shift temp_i right
        }
        // Alternative bit reversal (less efficient but maybe clearer):
        // size_t rev = 0;
        // size_t num_bits = static_cast<size_t>(std::round(std::log2(N)));
        // for (size_t k = 0; k < num_bits; ++k) {
        //     if ((i >> k) & 1) {
        //         rev |= 1 << (num_bits - 1 - k);
        //     }
        // }
        // size_t j = rev;


        if (j > i) {
            // Perform the swap only once (when j > i)
            std::swap(output[i], output[j]);
        }
    }

    // Cooley-Tukey FFT
    const double PI = std::acos(-1);
    for (size_t s = 1; s <= static_cast<size_t>(std::log2(N)); ++s) {
        size_t m = 1 << s;
        std::complex<double> wm = std::exp(std::complex<double>(0, -2.0 * PI / m));
        for (size_t k = 0; k < N; k += m) {
            std::complex<double> w = 1;
            for (size_t j = 0; j < m / 2; ++j) {
                std::complex<double> t = w * output[k + j + m / 2];
                std::complex<double> u = output[k + j];
                output[k + j] = u + t;
                output[k + j + m / 2] = u - t;
                w *= wm;
            }
        }
    }

    return output;
}

ComplexVec kiss_fft_wrapper(const ComplexVec& input) {
    size_t N = input.size();
    kiss_fft_cpx* in = new kiss_fft_cpx[N];
    kiss_fft_cpx* out = new kiss_fft_cpx[N];

    for (size_t i = 0; i < N; ++i) {
        in[i].r = static_cast<float>(input[i].real());
        in[i].i = static_cast<float>(input[i].imag());
    }

    kiss_fft_cfg cfg = kiss_fft_alloc(static_cast<int>(N), 0, nullptr, nullptr);
    kiss_fft(cfg, in, out);

    ComplexVec result(N);
    for (size_t i = 0; i < N; ++i) {
        result[i] = Complex(out[i].r, out[i].i);
    }

    free(cfg);
    delete[] in;
    delete[] out;
    return result;
}

ComplexVec quantize_adc(
    const ComplexVec& input,
    int num_bits,
    double v_fullscale)
{
    ComplexVec output;
    int levels = 1 << num_bits;
    double step = v_fullscale / levels;
    double half_scale = v_fullscale / 2.0;

    for (const auto& sample : input) {
        double real_clamped = std::min(std::max(std::real(sample), -half_scale), half_scale - step);
        double imag_clamped = std::min(std::max(std::imag(sample), -half_scale), half_scale - step);

        int real_code = static_cast<int>((real_clamped + half_scale) / step);
        int imag_code = static_cast<int>((imag_clamped + half_scale) / step);

        double real_q = (real_code * step) - half_scale;
        double imag_q = (imag_code * step) - half_scale;

        output.emplace_back(real_q, imag_q);
    }

    return output;
}


void dump_range_profile(const ComplexVec& sweep, const std::string& filename) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Failed to open dump file: " << filename << std::endl;
        return;
    }
    for (int i = 0; i < sweep.size(); ++i) {
        double power = std::norm(sweep[i]);  // magnitude^2
        out << i << " " << power << "\n";    // <-- space-separated
    }
    out.close();
}

ComplexVec ifft(const ComplexVec& input) {
    size_t N = input.size();
    ComplexVec conjugated(N);

    // Step 1: Conjugate input
    for (size_t i = 0; i < N; ++i)
        conjugated[i] = std::conj(input[i]);

    // Step 2: Forward FFT of conjugated input
    ComplexVec result = fft(conjugated);

    // Step 3: Conjugate again and scale
    for (auto& x : result)
        x = std::conj(x) / static_cast<double>(N);

    return result;
}

ComplexVec matched_filter_fft(const ComplexVec& signal, const ComplexVec& template_chirp) {
    int N = signal.size();
    int M = 2 * N;  // Zero-padding for linear convolution

    ComplexVec x(M, { 0, 0 });
    ComplexVec h(M, { 0, 0 });

    // Fill and zero-pad
    for (int i = 0; i < N; ++i) {
        x[i] = signal[i];
        h[i] = std::conj(template_chirp[N - 1 - i]);  // Conjugate time-reversed
    }

    ComplexVec X = fft(x);
    ComplexVec H = fft(h);

    // Multiply in frequency domain
    ComplexVec Y(M);
    for (int i = 0; i < M; ++i)
        Y[i] = X[i] * H[i];

    ComplexVec y = ifft(Y);

    // Return first N values (valid part of convolution)
    return ComplexVec(y.begin(), y.begin() + N);
}

void dump_rdm_binary(const std::vector<std::vector<int>>& binary_rdm, const std::string& filename) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Failed to open RDM binary dump file: " << filename << std::endl;
        return;
    }

    int doppler_bins = binary_rdm.size();
    int range_bins = binary_rdm[0].size();

    for (int i = 0; i < doppler_bins; ++i) {
        for (int j = 0; j < range_bins; ++j) {
            out << binary_rdm[i][j];
            if (j != range_bins - 1) out << " ";
        }
        out << "\n";
    }
    out.close();
}

void dump_rdm_power(const DoubleMat& rdm, const std::string& filename) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Failed to open RDM power dump file: " << filename << std::endl;
        return;
    }

    for (const auto& row : rdm) {
        for (size_t j = 0; j < row.size(); ++j) {
            out << row[j];
            if (j != row.size() - 1) out << " ";
        }
        out << "\n";
    }

    out.close();
}

void dump_aoa_estimates(const std::vector<AoA_Estimate>& estimates, const std::string& filename) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Failed to open AOA estimates dump file: " << filename << std::endl;
        return;
    }
    for (const auto& e : estimates) {
        out << e.range_bin << " " << e.doppler_bin << " "
            << e.azimuth_est << " " << e.azimuth_truth << " "
            << e.elevation_est << " " << e.elevation_truth << "\n";
    }

    out.close();
}

// Temporal integration
TemporalIntegrator::TemporalIntegrator(int N, int M)
    : max_history(N), threshold(M) {
}
void TemporalIntegrator::add_detection_map(const IntMat& binary_map) {
    if (history.size() >= max_history)
        history.pop_front();
    history.push_back(binary_map);
}

IntMat TemporalIntegrator::integrate() const {
    if (history.empty()) return {};

    int rows = history[0].size();
    int cols = history[0][0].size();
    IntMat result(rows, IntVec(cols, 0));

    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j) {
            int count = 0;
            for (const auto& map : history)
                count += map[i][j];
            if (count >= threshold)
                result[i][j] = 1;
        }

    return result;
}

void TemporalIntegrator::clear() {
    history.clear();
}
///////////////////////////////////////////////////////////////////////////////
//'Monopulse' module
//
// Significant simplifications were made to directly calculate phase while debugging
// Old implementation is below.
// Old implementation is closer to the real implementation. 
// Revisit and correct.
///////////////////////////////////////////////////////////////////////////////

std::vector<AoA_Estimate> estimate_monopulse_aoa(
    const ComplexMat& ch1,
    const ComplexMat& ch2,
    const ComplexMat& ch3,
    const ComplexMat& ch4,
    const IntMat& detections,
    const DoubleMat& rdm_power, // From ch1
    double lambda,
    double power_threshold_db, 
    double los_az_rad, 
    double los_el_rad)
{
    std::vector<AoA_Estimate> estimates;

    int range_bins = detections.size();
    int doppler_bins = detections[0].size();


    const double PI = acos(-1.0);
    std::ofstream debug_out("aoa_debug_dump.txt");
    debug_out << "# RangeBin DopplerBin EstAz(rad) EstEl(rad) TrueAz(rad) TrueEl(rad) RhoAz RhoEl SumMag\n";

    // Sanity check for dimensions
    if (ch1.size() != doppler_bins || ch1[0].size() != range_bins) {
        std::cerr << "*** ERROR: Mismatch between RDM dimensions and channel matrices ***" << std::endl;
        return estimates;
    }

    for (int r = 0; r < range_bins; ++r) {
        for (int d = 0; d < doppler_bins; ++d) {
            if (detections[r][d] != 1)
                continue;

            // --- Power Threshold Check ---
            double power = rdm_power[r][d];
            if (10.0 * log10(power + 1e-12) < power_threshold_db)
                continue; // Skip low-SNR detections

            // Access channel data as chX[doppler][range]
            Complex s1 = ch1[d][r]; // Top Right
            Complex s2 = ch2[d][r]; // Top Left
            Complex s3 = ch3[d][r]; // Bottom Right
            Complex s4 = ch4[d][r]; // Bottom Left
            // Calculate phase differences directly
            double phase1 = std::arg(s1);
            double phase2 = std::arg(s2);
            double phase3 = std::arg(s3);
            double phase4 = std::arg(s4);

            // Estimate delta_phi for Az (Right - Left average phase)
            double delta_phi_az_est = (phase1 + phase3) / 2.0 - (phase2 + phase4) / 2.0; // USE THIS            // Could also use (phase1+phase3)/2 - (phase2+phase4)/2
            // Previous was: double delta_phi_az_est = phase1 - phase2; // Approx TR-TL only - Less robust
            if (delta_phi_az_est > PI) delta_phi_az_est -= 2 * PI;
            if (delta_phi_az_est < -PI) delta_phi_az_est += 2 * PI;

            // Estimate delta_phi for El (Top - Bottom average phase)
            double delta_phi_el_est = (phase1 + phase2) / 2.0 - (phase3 + phase4) / 2.0; // USE THIS            // Could also use (phase1+phase2)/2 - (phase3+phase4)/2
            if (delta_phi_el_est > PI) delta_phi_el_est -= 2 * PI;
            if (delta_phi_el_est < -PI) delta_phi_el_est += 2 * PI;

            // Estimate angles using asin approximation (valid for d=lambda/2)
            // Clamp argument to asin to [-1, 1]
            double sin_arg_az = delta_phi_az_est / PI;
            double sin_arg_el = delta_phi_el_est / PI;
            sin_arg_az = std::max(-1.0, std::min(1.0, sin_arg_az));
            sin_arg_el = std::max(-1.0, std::min(1.0, sin_arg_el));

            double theta_az_est_rad = -asin(sin_arg_az);
            double theta_el_est_rad = asin(sin_arg_el);

            // Store estimate
            estimates.push_back({ r, d, theta_az_est_rad, theta_el_est_rad, los_az_rad, los_el_rad });

            // Write debug info using these new estimates
            debug_out << r << " " << d << " "
                << theta_az_est_rad << " " << theta_el_est_rad << " " // Using asin results
                << los_az_rad << " " << los_el_rad << " "
                << delta_phi_az_est << " " << delta_phi_el_est << " " // Raw phase diffs
                << std::abs(s1 + s2 + s3 + s4) << std::endl; // Sum magnitude
        }
    }
    debug_out.close();
    return estimates;
}

bool kiss_fft_exec(kiss_fft_cfg cfg, kiss_fft_cpx* fin, kiss_fft_cpx* fout, size_t N, const ComplexVec& input) {
    if (!cfg || !fin || !fout) {
        std::cerr << "*** Error: Invalid cfg or buffers passed to kiss_fft_exec ***" << std::endl;
        return false; // Indicate failure
    }

    // Copy input data to KissFFT input buffer
    // Add bounds check if input.size() might not match N?
    for (size_t i = 0; i < N; ++i) {
        // Check bounds if input vector size is not guaranteed externally
        if (i < input.size()) {
            fin[i].r = static_cast<float>(input[i].real());
            fin[i].i = static_cast<float>(input[i].imag());
        }
        else {
            fin[i].r = 0.0f; // Zero-pad if input is shorter than N
            fin[i].i = 0.0f;
        }
    }

    // Execute KissFFT
    kiss_fft(cfg, fin, fout);

    return true; // Indicate success
}
// --- DspBuffer4Ch Member Implementations ---

// Constructor Definition
DspBuffer4Ch::DspBuffer4Ch(int num_doppler_bins, int num_range_bins) :
    max_sweeps(num_doppler_bins),
    bins_per_sweep(num_range_bins),
    doppler_fft_cfg(nullptr),
    doppler_fft_in(nullptr),
    doppler_fft_out(nullptr),
    doppler_fft_initialized(false)
{
    std::cout << "DspBuffer4Ch constructed (FFT resources not allocated yet)." << std::endl;
    // IMPORTANT: Call initialize_doppler_fft() explicitly after creation
    // or ensure it's called before first use (e.g., in add_sweep_4ch).
}

// Destructor Definition
DspBuffer4Ch::~DspBuffer4Ch() {
    std::cout << "DspBuffer4Ch destructing..." << std::endl;
    if (doppler_fft_cfg) free(doppler_fft_cfg);
    delete[] doppler_fft_in;
    delete[] doppler_fft_out;
}

// initialize_doppler_fft Definition
void DspBuffer4Ch::initialize_doppler_fft() {
    if (doppler_fft_initialized) return;

    // Clean up just in case
    if (doppler_fft_cfg) free(doppler_fft_cfg);
    delete[] doppler_fft_in;
    delete[] doppler_fft_out;
    doppler_fft_cfg = nullptr; doppler_fft_in = nullptr; doppler_fft_out = nullptr;

    if (max_sweeps <= 0) {
        std::cerr << "*** Error: Invalid Doppler FFT size (max_sweeps) = " << max_sweeps << " ***" << std::endl;
        return;
    }

    std::cout << "Initializing Doppler FFT resources for size: " << max_sweeps << std::endl;
    try {
        doppler_fft_in = new kiss_fft_cpx[max_sweeps];
        doppler_fft_out = new kiss_fft_cpx[max_sweeps];
        doppler_fft_cfg = kiss_fft_alloc(max_sweeps, 0, nullptr, nullptr);

        if (!doppler_fft_cfg) { throw std::bad_alloc(); }
        doppler_fft_initialized = true;
    }
    catch (const std::bad_alloc& e) {
        std::cerr << "*** Error: Failed to allocate KissFFT resources for Doppler FFT: " << e.what() << " ***" << std::endl;
        delete[] doppler_fft_in; doppler_fft_in = nullptr;
        delete[] doppler_fft_out; doppler_fft_out = nullptr;
        if (doppler_fft_cfg) free(doppler_fft_cfg); doppler_fft_cfg = nullptr;
        doppler_fft_initialized = false;
    }
}
// clear Definition
void DspBuffer4Ch::clear() {
    ch1.clear();
    ch2.clear();
    ch3.clear();
    ch4.clear();
    // std::cout << "DspBuffer4Ch cleared." << std::endl;
}

// add_sweep_4ch Definition
void DspBuffer4Ch::add_sweep_4ch(const ComplexVec& s1, const ComplexVec& s2, const ComplexVec& s3, const ComplexVec& s4) {
    // Ensure FFT is initialized before adding data that might make it ready
    if (!doppler_fft_initialized) {
        initialize_doppler_fft();
        if (!doppler_fft_initialized) {
            std::cerr << "*** Error: Cannot add sweep, Doppler FFT not initialized. ***" << std::endl;
            return;
        }
    }

    if ((int)s1.size() != bins_per_sweep /*|| check others*/) { return; } // Basic check

    auto trim = [this](ComplexMat& mat) {
        if ((int)mat.size() >= max_sweeps)
            mat.erase(mat.begin());
        };
    trim(ch1); trim(ch2); trim(ch3); trim(ch4);

    ch1.push_back(s1);
    ch2.push_back(s2);
    ch3.push_back(s3);
    ch4.push_back(s4);
}

// compute_rdm Definition
DoubleMat DspBuffer4Ch::compute_rdm(const ComplexMat& channel) {
    if (!doppler_fft_initialized || !is_ready() || channel.empty() || channel[0].empty()) {
        std::cerr << "*** Error: Doppler FFT not ready or channel empty/invalid in compute_rdm ***" << std::endl;
        return DoubleMat();
    }
    // ... (rest of implementation using kiss_fft_exec as shown before) ...
    int current_bins_sweep = this->bins_per_sweep;
    int current_max_sweeps = this->max_sweeps;
    if ((int)channel.size() != current_max_sweeps || (int)channel[0].size() != current_bins_sweep) {
        std::cerr << "*** Error: Input channel dimensions mismatch buffer state in compute_rdm ***" << std::endl;
        return DoubleMat();
    }
    DoubleMat rdm(current_bins_sweep, std::vector<double>(current_max_sweeps, 0.0));
    ComplexVec slow_time_input(current_max_sweeps);
    const double PI = acos(-1.0);
    for (int bin = 0; bin < current_bins_sweep; ++bin) {
        for (int m = 0; m < current_max_sweeps; ++m) { slow_time_input[m] = channel[m][bin]; }
        ComplexVec windowed_slow_time(current_max_sweeps);
        for (int m = 0; m < current_max_sweeps; ++m) {
            double hann_mult = 0.5 * (1.0 - cos(2.0 * PI * m / (current_max_sweeps - 1)));
            windowed_slow_time[m] = slow_time_input[m] * hann_mult;
        }
        bool success = kiss_fft_exec(this->doppler_fft_cfg, this->doppler_fft_in, this->doppler_fft_out, current_max_sweeps, windowed_slow_time);
        if (success) {
            for (int d = 0; d < current_max_sweeps; ++d) {
                if (bin < (int)rdm.size() && d < (int)rdm[bin].size()) {
                    rdm[bin][d] = std::norm(Complex(this->doppler_fft_out[d].r, this->doppler_fft_out[d].i));
                }
            }
        }
        else { /* Handle error */ }
    }
    return rdm;
}

// get_complex_rdm Definition
ComplexMat DspBuffer4Ch::get_complex_rdm(const ComplexMat& channel) {
    if (!doppler_fft_initialized || !is_ready() || channel.empty() || channel[0].empty()) {
        std::cerr << "*** Error: Doppler FFT not ready or channel empty/invalid in get_complex_rdm ***" << std::endl;
        return ComplexMat();
    }
    // ... (rest of implementation using kiss_fft_exec as shown before) ...
    int current_bins_sweep = this->bins_per_sweep;
    int current_max_sweeps = this->max_sweeps;
    if ((int)channel.size() != current_max_sweeps || (int)channel[0].size() != current_bins_sweep) {
        std::cerr << "*** Error: Input channel dimensions mismatch buffer state in get_complex_rdm ***" << std::endl;
        return ComplexMat();
    }
    ComplexMat complex_rdm(current_bins_sweep, ComplexVec(current_max_sweeps));
    ComplexVec slow_time_input(current_max_sweeps);
    const double PI = acos(-1.0);
    for (int bin = 0; bin < current_bins_sweep; ++bin) {
        for (int m = 0; m < current_max_sweeps; ++m) { slow_time_input[m] = channel[m][bin]; }
        ComplexVec windowed_slow_time(current_max_sweeps);
        for (int m = 0; m < current_max_sweeps; ++m) {
            double hann_mult = 0.5 * (1.0 - cos(2.0 * PI * m / (current_max_sweeps - 1)));
            windowed_slow_time[m] = slow_time_input[m] * hann_mult;
        }
        bool success = kiss_fft_exec(this->doppler_fft_cfg, this->doppler_fft_in, this->doppler_fft_out, current_max_sweeps, windowed_slow_time);
        if (success) {
            for (int d = 0; d < current_max_sweeps; ++d) {
                if (bin < (int)complex_rdm.size() && d < (int)complex_rdm[bin].size()) {
                    complex_rdm[bin][d] = Complex(this->doppler_fft_out[d].r, this->doppler_fft_out[d].i);
                }
            }
        }
        else { /* Handle error */ }
    }
    return complex_rdm;
}

/* Near original implementation
std::vector<AoA_Estimate> estimate_monopulse_aoa(
    const ComplexMat& ch1,
    const ComplexMat& ch2,
    const ComplexMat& ch3,
    const ComplexMat& ch4,
    const IntMat& detections,
    const DoubleMat& rdm_power, // From ch1
    double lambda,
    double power_threshold_db,
    double los_az_rad,
    double los_el_rad)
{
    std::vector<AoA_Estimate> estimates;

    int range_bins = detections.size();
    int doppler_bins = detections[0].size();


    const double PI = acos(-1.0);
    std::ofstream debug_out("aoa_debug_dump.txt");
    debug_out << "# RangeBin DopplerBin EstAz(rad) EstEl(rad) TrueAz(rad) TrueEl(rad) RhoAz RhoEl SumMag\n";

    // Sanity check for dimensions
    if (ch1.size() != doppler_bins || ch1[0].size() != range_bins) {
        std::cerr << "*** ERROR: Mismatch between RDM dimensions and channel matrices ***" << std::endl;
        return estimates;
    }

    for (int r = 0; r < range_bins; ++r) {
        for (int d = 0; d < doppler_bins; ++d) {
            if (detections[r][d] != 1)
                continue;

            // --- Power Threshold Check ---
            double power = rdm_power[r][d];
            if (10.0 * log10(power + 1e-12) < power_threshold_db)
                continue; // Skip low-SNR detections

            // Access channel data as chX[doppler][range]
            Complex s1 = ch1[d][r]; // Top Right
            Complex s2 = ch2[d][r]; // Top Left
            Complex s3 = ch3[d][r]; // Bottom Right
            Complex s4 = ch4[d][r]; // Bottom Left

            // Monopulse sum and difference
            Complex sum = s1 + s2 + s3 + s4;
            if (std::abs(sum) < 1e-6)
                continue;

            // Changed az from s1 - s3 to s3 - s1 because of orientation
            Complex diff_az = (s1 + s3) - (s2 + s4);
            Complex diff_el = (s1 + s2) - (s3 + s4);

            double rho_az = std::real(diff_az / sum);
            double rho_el = std::real(diff_el / sum);
            // Added a negative to the theta az due to orientation

            //double theta_az = (4.0 / (4.0 *PI)) * atan(rho_az); // radians
            //double theta_el = (4.0 / (4.0 * PI)) * atan(rho_el); // radians

            // Using small angle appx and k_slope = 0.327
            double K = 1.6;
            double theta_hpbw = 30 * (PI / 180); // Radians
            double k_m_rad = K / theta_hpbw;
            double k_slope = 1 / k_m_rad;
            double theta_az = k_slope * atan(rho_az);
            double theta_el = k_slope * atan(rho_el);
            // === DEBUG: Log raw phase and magnitude ===
            double phase_s1 = std::arg(s1);
            double phase_s3 = std::arg(s3);
            //double delta_phase = phase_s1 - phase_s3;

            //Phase wrapping
            double delta_phase = phase_s1 - phase_s3;
            if (delta_phase > PI) delta_phase -= 2 * PI;
            if (delta_phase < -PI) delta_phase += 2 * PI;

            double mag_s1 = std::abs(s1);
            double mag_s3 = std::abs(s3);

            debug_out << r << " " << d << " "
                << theta_az << " " << theta_el << " "
                << los_az_rad << " " << los_el_rad << " "
                << delta_phase << " "
                << mag_s1 << " " << mag_s3 << std::endl;


            estimates.push_back({ r, d, theta_az, theta_el, los_az_rad, los_el_rad });
        }
    }
    debug_out.close();
    return estimates;
} */