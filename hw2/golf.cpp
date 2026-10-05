#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>

// We first define our physics parameters
const double G = 9.80;
const double MASS = 0.046;
const double RHO = 1.29;
const double AREA = 0.0014;
const double V0 = 70.0;
const double MAGNUS_CONST = 0.25;
const double pi = std::acos(-1.0);

enum ModelType { // These are the types of gold balls which we simulated. ModelType is restricted to these three options to prevent possible errors
    IDEAL = 1,
    SMOOTH = 2,
    DIMPLED = 3,
    DIMPLED_SPIN = 4
};

struct TrajectoryResult {
    double range;
    double flight_time;
};

// We compute drag coefficient C here.
double GetC(ModelType model, double v_mag) {
    if (model == IDEAL) return 0.0;
    if (model == SMOOTH || v_mag <= 14.0) return 0.5;
    return 7.0 / v_mag;
}

// We simulate the thing here
TrajectoryResult run_simulation(ModelType model, double theta_deg, double dt, std::ofstream* out_file = nullptr) {
    double rad = theta_deg * pi / 180.0;

    double x = 0.0;
    double y = 0.0;
    double vx = V0 * std::cos(rad);
    double vy = V0 * std::sin(rad);
    double t = 0.0;

    double drag_prefactor = (RHO * AREA) / MASS;

    while (true) {
        if (out_file && out_file->is_open()) {
            *out_file << std::fixed << std::setprecision(4) << t << ","
                      << std::setprecision(3) << x << "," << y << ","
                      << static_cast<int>(model) << "\n";
        }

        // Required order: x, y, v_mag, C, vx, vy
        double x_next = x + vx * dt;
        double y_next = y + vy * dt;

        if (y_next < 0.0 && t > 0.0) {
            // Ground intersection linear interpolation
            double fraction = -y / (y_next - y);
            double final_range = x + fraction * (x_next - x);
            double final_time = t + fraction * dt;
            if (out_file && out_file->is_open()) {
                *out_file << std::fixed << std::setprecision(4) << final_time << ","
                          << std::setprecision(3) << final_range << ",0.000,"
                          << static_cast<int>(model) << "\n";
            }
            return {final_range, final_time};
        }

        double v_mag = std::sqrt(vx * vx + vy * vy);
        double C = GetC(model, v_mag);

        // Compute accelerations by adding drag and spin terms
        double ax = -C * drag_prefactor * v_mag * vx;
        double ay = -G - C * drag_prefactor * v_mag * vy;

        if (model == DIMPLED_SPIN) {
            ax -= MAGNUS_CONST * vy;
            ay += MAGNUS_CONST * vx;
        }

        double vx_next = vx + ax * dt;
        double vy_next = vy + ay * dt;

        x = x_next;
        y = y_next;
        vx = vx_next;
        vy = vy_next;
        t += dt;
    }
}

// Helper to write CSV files directly to disk
void write_golf_csv(double theta_deg) {
    std::string fname = "golf_theta" + std::to_string(static_cast<int>(theta_deg)) + ".csv";
    std::ofstream out(fname);
    if (!out.is_open()) {
        std::cerr << "Error opening file: " << fname << "\n";
        return;
    }
    out << "t_s,x_m,y_m,model_id\n";
    const double dt = 0.005;
    for (int m = 1; m <= 4; ++m) {
        run_simulation(static_cast<ModelType>(m), theta_deg, dt, &out);
    }
    out.close();
    std::cout << "[I/O] Generated CSV: " << fname << "\n";
}

int main(int argc, char* argv[]) {
    double plot_theta = -1.0;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg.rfind("--plot=", 0) == 0) {
            plot_theta = std::stod(arg.substr(7));
        }
    }

    if (plot_theta > 0.0) {
        write_golf_csv(plot_theta);
        return 0;
    }

    std::cout << "=========================================================================\n";
    std::cout << "--- Problem 2 ---\n";

    // We compare now: Ideal/No-drag-or-spin trajectory vs analytic equations
    std::cout << "\n- Numerical Check: Ideal Trajectory vs Analytic Exact -\n";
    std::cout << std::setw(8) << "Theta"
              << std::setw(15) << "Simulated Range (m)"
              << std::setw(15) << "Analytic Range (m)"
              << std::setw(15) << "Simulated Time (s)"
              << std::setw(15) << "Analytic Time (s)\n";

    std::vector<double> angles = {45.0, 30.0, 15.0, 9.0};
    const double dt_standard = 0.01;

    for (double th : angles) {
        TrajectoryResult sim = run_simulation(IDEAL, th, dt_standard);
        double rad = th * pi / 180.0;
        double exact_range = (V0 * V0 * std::sin(2.0 * rad)) / G;
        double exact_time = (2.0 * V0 * std::sin(rad)) / G;

        std::cout << std::fixed << std::setprecision(1) << std::setw(7) << th << "°"
                  << std::setprecision(3) << std::setw(15) << sim.range
                  << std::setw(15) << exact_range
                  << std::setw(15) << sim.flight_time
                  << std::setw(15) << exact_time << "\n";
    }

    // Four models comparison across all four angles
    std::cout << "\n- Trajectory Ranges Depending on Model W/R to Angle (dt = 0.01 s) -\n";
    std::cout << std::setw(10) << "Angle"
              << std::setw(16) << "(1) Ideal (m)"
              << std::setw(16) << "(2) Smooth (m)"
              << std::setw(16) << "(3) Dimpled (m)"
              << std::setw(18) << "(4) Dimpled Ball w/ Spin\n";

    for (double th : angles) {
        TrajectoryResult r1 = run_simulation(IDEAL, th, dt_standard);
        TrajectoryResult r2 = run_simulation(SMOOTH, th, dt_standard);
        TrajectoryResult r3 = run_simulation(DIMPLED, th, dt_standard);
        TrajectoryResult r4 = run_simulation(DIMPLED_SPIN, th, dt_standard);

        std::cout << std::fixed << std::setprecision(1) << std::setw(9) << th << "°"
                  << std::setprecision(2) << std::setw(16) << r1.range
                  << std::setw(16) << r2.range
                  << std::setw(16) << r3.range
                  << std::setw(18) << r4.range << "\n";
    }

    // We now do the part b) here, Halving dt from 0.01 s to 0.005 s
    std::cout << "\n- Halfed Time-step (0.01 s -> 0.005 s) for Dimpled Ball w/ Spin -\n";
    std::cout << std::setw(10) << "Angle"
              << std::setw(18) << "Range when dt = 0.01s"
              << std::setw(18) << "Range when dt = 0.005s"
              << std::setw(16) << "Difference in Range\n";

    const double dt_half = 0.005;
    for (double th : angles) {
        TrajectoryResult r_std = run_simulation(DIMPLED_SPIN, th, dt_standard);
        TrajectoryResult r_half = run_simulation(DIMPLED_SPIN, th, dt_half);
        double delta = r_half.range - r_std.range;

        std::cout << std::fixed << std::setprecision(1) << std::setw(9) << th << "°"
                  << std::setprecision(3) << std::setw(18) << r_std.range
                  << std::setw(18) << r_half.range
                  << std::setw(16) << delta << " m\n";
    }
    std::cout << "=========================================================================\n";

    // Write all four trajectory CSVs for report figure generation
    for (double th : angles) {
        write_golf_csv(th);
    }

    return 0;
}