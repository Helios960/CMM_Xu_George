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
    const double mass_g = mass_kg * 10000;
    const double molar_mass = 14;
    const double N_A = 6.022e23;
    const double N0 = (mass_g / molar_mass) * N_A;

    
    
}