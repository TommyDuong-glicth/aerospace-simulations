// lifting_line.cpp - C++ port of lifting_line.py (Prandtl/Glauert monoplane equation).
// No external libraries. Build:  g++ -O3 -std=c++17 lifting_line.cpp -o lifting_line
// Usage:  ./lifting_line verify          -> CL, e, CDi for rectangular wing AR=8, alpha=5 deg
//         ./lifting_line bench [reps]    -> 10,000-case sweep (100 alphas x 100 tapers), best time + checksum
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using Vec = std::vector<double>;
using Mat = std::vector<Vec>;
const double PI = 3.14159265358979323846;
const double A0 = 2.0 * PI;

struct Result { double CL, e, CDi; };

// Gaussian elimination with partial pivoting; M and b are taken by value (modified locally).
Vec gauss_solve(Mat M, Vec b) {
    const int n = static_cast<int>(b.size());
    for (int k = 0; k < n; ++k) {
        int piv = k;
        for (int r = k + 1; r < n; ++r)
            if (std::fabs(M[r][k]) > std::fabs(M[piv][k])) piv = r;
        std::swap(M[k], M[piv]);
        std::swap(b[k], b[piv]);
        for (int r = k + 1; r < n; ++r) {
            const double f = M[r][k] / M[k][k];
            for (int c = k; c < n; ++c) M[r][c] -= f * M[k][c];
            b[r] -= f * b[k];
        }
    }
    Vec x(n);
    for (int r = n - 1; r >= 0; --r) {            // back substitution
        double s = b[r];
        for (int c = r + 1; c < n; ++c) s -= M[r][c] * x[c];
        x[r] = s / M[r][r];
    }
    return x;
}

// Trapezoidal planform with taper ratio lambda (c/b). alpha measured from zero-lift.
Result solve_ll(double AR, double taper, double alpha_deg, double aL0_deg, int N) {
    Mat M(N, Vec(N));
    Vec rhs(N, (alpha_deg - aL0_deg) * PI / 180.0);
    const double c_root = 2.0 / (AR * (1.0 + taper));
    for (int i = 1; i <= N; ++i) {
        const double th = i * PI / (2.0 * N);
        const double c = c_root * (1.0 - (1.0 - taper) * std::cos(th));
        for (int j = 1; j <= N; ++j) {
            const double n = 2.0 * j - 1.0;
            M[i - 1][j - 1] = std::sin(n * th) * (4.0 / (A0 * c) + n / std::sin(th));
        }
    }
    const Vec A = gauss_solve(M, rhs);
    const double CL = PI * AR * A[0];
    double delta = 0.0;
    for (int j = 2; j <= N; ++j) {
        const double n = 2.0 * j - 1.0, r = A[j - 1] / A[0];
        delta += n * r * r;
    }
    const double e = 1.0 / (1.0 + delta);
    return {CL, e, CL * CL / (PI * e * AR)};
}

int main(int argc, char** argv) {
    const std::string mode = argc > 1 ? argv[1] : "verify";
    const int N = 40;
    if (mode == "verify") {
        const Result r = solve_ll(8.0, 1.0, 5.0, 0.0, N);
        std::cout << std::setprecision(12) << "CL " << r.CL << "\ne " << r.e << "\nCDi " << r.CDi << "\n";
        return 0;
    }
    const int reps = argc > 2 ? std::stoi(argv[2]) : 5;
    double best = 1e300, checksum = 0.0;
    for (int rep = 0; rep < reps; ++rep) {
        double sum = 0.0;
        const auto t0 = std::chrono::high_resolution_clock::now();
        for (int ia = 0; ia < 100; ++ia) {
            const double alpha = 0.1 + (10.0 - 0.1) * ia / 99.0;
            for (int it = 0; it < 100; ++it) {
                const double taper = 0.1 + (1.0 - 0.1) * it / 99.0;
                const Result r = solve_ll(8.0, taper, alpha, 0.0, N);
                sum += r.CL + r.e;               // keeps the optimizer from discarding the work
            }
        }
        const auto t1 = std::chrono::high_resolution_clock::now();
        best = std::min(best, std::chrono::duration<double>(t1 - t0).count());
        checksum = sum;
    }
    std::cout << std::setprecision(15) << "BENCH " << best << " " << checksum << "\n";
    return 0;
}
