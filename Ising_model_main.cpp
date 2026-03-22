#include <vector>
#include <cmath>
#include <iostream>
#include <numeric>
#include <fstream>
#include <random>
#include <string>
#include <iomanip>

using namespace std;

std::mt19937 gen(42);
std::uniform_real_distribution<double> dist_real(0.0, 1.0);

// Calculate Magnetization
double cal_mag(const vector<vector<int>>& mesh, int N) {
    double M = 0.0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            M += mesh[i][j];
        }
    }
    return fabs(M); 
}

// Calculate total energy
double cal_E(const vector<vector<int>>& mesh, int N, double J, bool FREE_BC) {
    double E = 0.0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            int spin = mesh[i][j];
            int right = 0, down = 0;

            if (FREE_BC) {
                if (j + 1 < N) right = mesh[i][j + 1];
                if (i + 1 < N) down  = mesh[i + 1][j];
            } else {
                right = mesh[i][(j + 1) % N];
                down  = mesh[(i + 1) % N][j];
            }

            E -= J * spin * (right + down);
        }
    }
    return E;
}

// Single Metropolis step: attempt to flip a random spin
void MC_steps(vector<vector<int>>& mesh, int N, double T, double J, bool FREE_BC) {
    std::uniform_int_distribution<int> dist_pos(0, N - 1);
    int i = dist_pos(gen);
    int j = dist_pos(gen);

    int spin = mesh[i][j];
    int up = 0, down = 0, left = 0, right = 0;

    if (FREE_BC) {
        if (i > 0)     up    = mesh[i - 1][j];
        if (i < N - 1) down  = mesh[i + 1][j];
        if (j > 0)     left  = mesh[i][j - 1];
        if (j < N - 1) right = mesh[i][j + 1];
    } else {
        up    = mesh[(i - 1 + N) % N][j];
        down  = mesh[(i + 1) % N][j];
        left  = mesh[i][(j - 1 + N) % N];
        right = mesh[i][(j + 1) % N];
    }

    double dE = 2.0 * J * spin * (up + down + left + right);
    if (dE <= 0.0 || dist_real(gen) < exp(-dE / T)) {
        mesh[i][j] = -spin;
    }
}

// Perform equilibration steps
void equilibration(vector<vector<int>>& mesh, int N, double T, double J, int steps, bool FREE_BC) {
    for (int i = 0; i < steps; i++) {
        for (int j = 0; j < N*N; j++) {
            MC_steps(mesh, N, T, J, FREE_BC);  // <-- add FREE_BC
        }
    }
}

vector<double> temp_grid_create(double T_min, double T_max, double T_step_coarse, double T_step_fine, double Tc, double Tc_range) {
    vector<double> temps;
    
    for (double T = T_min; T < Tc - Tc_range; T += T_step_coarse) {
        temps.push_back(T);
    }
    
    for (double T = Tc - Tc_range; T <= Tc + Tc_range; T += T_step_fine) {
        temps.push_back(T);
    }
    
    for (double T = Tc + Tc_range + T_step_coarse; T <= T_max; T += T_step_coarse) {
        if (T > Tc + Tc_range) {
            temps.push_back(T);
        }
    }
    
    return temps;
}

int main() {
    // Simulation parameters
    int N = 30;                  // System size
    double J = 1.0;             // Coupling constant ( +1.0 -> ferro; -1.0 -> anti-ferro)
    bool FREE_BC = false;      // true = free boundary, false = periodic
    double Tc = 2.269;        // Theoretical critical temperature for 2D Ising
    
    vector<double> temps = temp_grid_create(1.5, 3.5, 0.1, 0.02, Tc, 0.3);
    
    // Initialize simulation
    cout << "Starting 2D Ising model simulation" << endl;
    cout << "System size: " << N << "x" << N << endl;
    cout << "Temperature range: " << temps.front() << " to " << temps.back() << endl;
    
    int eq_steps = 5000;      // Equilibration sweeps
    int measure_steps = 20000; // Measurement sweeps
    
    // Create a single output file for all results
    string fname = "ising_results";
    if (J < 0)    fname += "_anti_ferro";
    else          fname += "_ferro";
    if (FREE_BC)  fname += "_FB";
    else          fname += "_PB";
    fname += ".txt";
    ofstream outfile(fname);
    outfile << "# Temperature Energy Magnetization SpecificHeat Susceptibility" << endl;
    
    // Run simulation
    for (double T : temps) {
        vector<vector<int>> mesh(N, vector<int>(N));
        
        // Initialize lattice with random spins
        std::uniform_int_distribution<int> dist_spin(0, 1);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                mesh[i][j] = (dist_spin(gen) == 0) ? -1 : 1;
            }
        }
        
        cout << "Temperature: " << T << endl;
        
        // Equilibration
        equilibration(mesh, N, T, J, eq_steps, FREE_BC);
        
        // Measurement
        vector<double> E_samples;
        vector<double> M_samples;
        
        for (int step = 0; step < measure_steps; step++) {
            for (int j = 0; j < N*N; j++) {
                MC_steps(mesh, N, T, J, FREE_BC);  // <-- add FREE_BC
            }
            if (step % 10 == 0) {
                double E = cal_E(mesh, N, J, FREE_BC);  // <-- add FREE_BC
                double M = cal_mag(mesh, N); 
                E_samples.push_back(E);
                M_samples.push_back(M);
            }
        }
        
        double E_avg = accumulate(E_samples.begin(), E_samples.end(), 0.0) / E_samples.size() / (N*N);
        double M_avg = accumulate(M_samples.begin(), M_samples.end(), 0.0) / M_samples.size() / (N*N);
        
        double E2_avg = 0.0;
        for (double E : E_samples) {
            E2_avg += (E/N/N) * (E/N/N);
        }
        E2_avg /= E_samples.size();
        
        double M2_avg = 0.0;
        for (double M : M_samples) {
            M2_avg += (M/N/N) * (M/N/N);
        }
        M2_avg /= M_samples.size();
        
        double Cv = N*N * (E2_avg - E_avg*E_avg) / (T*T);
        double chi = N*N * (M2_avg - M_avg*M_avg) / T;
        
        cout << "  E/N^2 = " << E_avg << ", M/N^2 = " << M_avg 
             << ", Cv = " << Cv << ", Chi = " << chi << endl;
             
        // Save results for this temperature
        outfile << fixed << setprecision(6)
                << T << " " << E_avg << " " << M_avg << " " << Cv << " " << chi << endl;
    }
    
    // Close the output file
    outfile.close();
    
    cout << "Results saved to " << fname << endl;

    return 0;
}