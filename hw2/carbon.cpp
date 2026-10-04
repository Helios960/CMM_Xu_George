#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>

int main(int argc, char* argv[]){
    // Physical Constants
    const double T_half = 5700.0;
    const double tau = T_half / std::log(2.0);
    const double mass_kg = 1.0e-12;
    const double mass_g = mass_kg * 1000;
    const double molar_mass = 14;
    const double N_A = 6.022e23;
    const double N0 = (mass_g / molar_mass) * N_A;
    const double t_max = 20000;

    double plot_dt = -1.0;
    for (int i = 1; i < argc; ++i){
        std::string arg = argv[i];
        if (arg.rfind("--plot=",0) == 0){ // For some reason I recall starts_with being a thing but I strangely cannot use it. Whatever :(. I'll use the older method.
            plot_dt = std::stod(arg.substr(7));
        }
    }
    
    if (plot_dt > 0.0) {
        std::cout << "t_years,activity_euler_Bq,activity_exact_Bq\n";
        const double s_p_y = 365.25 * 86400.0;
        double t = 0.0;
        double N = N0;
        while (t - t_max <= 1e-9) {
            double R_euler = (N / tau) / s_p_y;
            double R_exact = ((N0 * std::exp(-t / tau)) / tau) / s_p_y;
            std::cout << std::fixed << std::setprecision(2) << t << ","
                      << std::setprecision(5) << R_euler << "," << R_exact << "\n";
            N += (-N / tau) * plot_dt;
            t += plot_dt;
        }
        return 0;
    }

    std::cout << "========================================================\n";
    std::cout << "--- Problem 1: ---\n";
    std::cout << "Initial atoms N₀     : " << std::scientific << std::setprecision(5) << N0 << "\n";
    std::cout << "Decay constant tau   : " << std::fixed << std::setprecision(2) << tau << " years\n";
    std::cout << "Initial activity R₀  : " << (N0 / tau) / (365.25 * 86400.0) << " Bq\n\n";

    // We evaluate at all three time step sizes:
    std::vector<double> dt_list = {10.0, 100.0, 1000.0};
    for (double dt : dt_list) {
        double t = 0.0;
        double N = N0;
        while (t < t_max - 1e-9) {
            N += (-N / tau) * dt;
            t += dt;
        }
        double N_analytic = N0 * std::exp(-t / tau);
        double perc_err = (N - N_analytic) / N_analytic * 100.0;
        std::cout << "Step size dt = " << std::setw(4) << static_cast<int>(dt) << " yr" 
                  << " | N(20,000 yr) = " << std::scientific << std::setprecision(4) << N 
                  << " | Analytical Solution = " << N_analytic
                  << " | Percentage Error = " << std::fixed << std::setprecision(2) << perc_err << "%\n";
    }

    // Part (c): Error analysis at t = 2 * T_half = 11,400 years with dt = 1,000 years
    double t_target = 2.0 * T_half; // 11400 years
    double dt_c = 1000.0;
    double N_c = N0;
    double t_curr = 0.0;

    while (t_curr + dt_c - t_target <= 1e-9) {
        N_c += (-N_c / tau) * dt_c;
        t_curr += dt_c;
    }
    // We must also take the remainder of the 400 years per the instructions given by the professor 
    double remainder_dt = t_target - t_curr;
    if (remainder_dt > 1e-9) {
        N_c += (-N_c / tau) * remainder_dt;
        t_curr += remainder_dt;
    }

    double N_2 = N0 * 0.25;
    double actual_percent_dev = (N_c - N_2) / N_2 * 100.0;
    double order_error = -std::log(2.0) * (dt_c / tau) * 100.0;

    std::cout << "\n--- Part (c) Evaluation at 2 Half-Lives (t = 11,400 yr) ---\n";
    std::cout << "Euler result N(11,400 yr)   : " << std::scientific << N_c << "\n";
    std::cout << "Analytic result N₂        : " << N_2 << "\n";
    std::cout << "Actual percentage deviation : " << std::fixed << std::setprecision(3) << actual_percent_dev << " %\n";
    std::cout << "Expected second-order error : " << order_error << " %\n";
    std::cout << "========================================================\n";

    return 0;
}