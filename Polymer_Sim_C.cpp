#include <iostream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <tuple>
#include <random>
#include <chrono>
#include <iomanip>


// ################# Simulation Parameter #################
double K_Bt = 1.0 ; 
double K_bond = 100.0 ; 
double b0 = 2.0 ;
double epsilon = 1.0 ; 
double sigma = 1.0 ; 
double dt = 0.01;
double m = 1.0;
double eta = 1.0;
double r = 1.0;
double alpha = 6.0 * 3.1456 * eta * r;
double a = (1.0 - alpha * dt / (2.0 * m)) / (1.0 + alpha * dt / (2.0 * m));
double b = 1.0 / (1.0 + alpha * dt / (2.0 * m));
double beta_std = std::sqrt(2.0 * alpha * K_Bt * dt);
double sigma_std = std::sqrt(K_Bt);

int n_steps = 500000;
int record_every = 200;


// ################# Random Number Initialization ############
std::mt19937 rng(std::random_device{}());
std::uniform_real_distribution<double> dist(0.0, 1.0);
std::normal_distribution<double> norm_dist(0.0, sigma_std);
std::normal_distribution<double> norm_dist_b(0.0, beta_std);
double x = dist(rng);


std::vector<std::vector<std::vector<double>>>
noise_3D(int n_frames, int n_part) {
    std::vector<std::vector<std::vector<double>>> noise(
        n_frames, std::vector<std::vector<double>>(
            n_part, std::vector<double>(3)));

    for (auto& frame : noise)
        for (auto& part : frame)
            for (auto& c : part)
                c = norm_dist_b(rng);

    return noise;
}

// ################# Position Initialization #################
int n_part = 150;
double bond_length = 2.0;
std::vector<double> Center = {0.0, 0.0, 0.0};

std::vector<std::vector<double>> initialize_polymer_chain(const int n_part, const double bond_length, const std::vector<double> Center){
    /*
    Initializes a linear polymer chain centered at a specified 3D coordinate.

    Parameters:
    - n_part: Number of particles in the polymer 
    - bond_length: Distance between adjacent monomers (default: 2.0)
    - Center: Array-like of shape (3,) specifying the center of the chain (default: [0.0, 0.0, 0.0])

    Returns:
    - positions: NumPy array of shape (n_part, 3) containing particle coordinates.
    */

    std::vector<std::vector<double>> positions(n_part, std::vector<double>(3));
    std::vector<double> offsets(n_part);

    // creating offsets
    for(int i = 0; i < n_part; i++){
        offsets[i] = (i - (n_part-1)/2.0)*(bond_length) ;
    }

    // creating positoin array in accord to offsets 
    for(int i = 0; i < n_part; i++){
        positions[i][0] = Center[0] + offsets[i];
        positions[i][1] = Center[1];
        positions[i][2] = Center[2];       
    }

    return positions;
}

// ################# Velocity Initialization #################
std::vector<std::vector<double>> Initialize_Velocity(const int n_part, const double K_Bt){
    std::vector<std::vector<double>> velocities(n_part, std::vector<double>(3));
    for (int i = 0; i < n_part; i++)
        for (int j = 0; j < 3; j++)
            velocities[i][j] = norm_dist(rng);
    
    return velocities;
}


// ################################# Forces Calculation ##############################
std::pair<std::vector<std::vector<double>>, double>
bond_force(const std::vector<std::vector<double>>& r,
           double k_bond, double b0) {

    std::vector<std::vector<double>> f(n_part, std::vector<double>(3, 0.0));
    double U_bond = 0.0;

    for (int i = 0; i < n_part - 1; i++) {
        double dr0 = r[i+1][0] - r[i][0];
        double dr1 = r[i+1][1] - r[i][1];
        double dr2 = r[i+1][2] - r[i][2];
        double dist = std::sqrt(dr0*dr0 + dr1*dr1 + dr2*dr2);
        U_bond += 0.5 * k_bond * (dist - b0) * (dist - b0);
        // avoiding collision 
        if (dist > 1e-10) {
            double fmag = k_bond * (dist - b0) / dist;
            f[i+1][0] -= fmag * dr0;
            f[i+1][1] -= fmag * dr1;
            f[i+1][2] -= fmag * dr2;
            f[i][0]   += fmag * dr0;
            f[i][1]   += fmag * dr1;
            f[i][2]   += fmag * dr2;
        }
    } 
    double U_bond_mean = U_bond/n_part; // will throw error here most prolly
    return {f,U_bond_mean};
}

std::vector<std::vector<double>> wca_force(const std::vector<std::vector<double>>& r, double epsilon, double sigma){
    std::vector<std::vector<double>> f(n_part, std::vector<double>(3, 0.0));
    double rcut = std::pow(2.0, 1.0/6.0) * sigma ;
    double rcut_sq = std::pow(rcut, 2) ;
    double sigma6 = std::pow(sigma, 6) ;
    // double sigma12 = std::pow(sigma, 12) ;
    for(int i = 0 ; i < n_part; i++){
        for(int j = i+1 ; j < n_part ; j++){
            double dx = r[j][0] - r[i][0];
            double dy = r[j][1] - r[i][1];
            double dz = r[j][2] - r[i][2];

            double dist_sq = dx*dx + dy*dy + dz*dz ;

            if (dist_sq < rcut_sq){
                double inv_r = 1.0/(std::sqrt(dist_sq));
                double inv_r6 = sigma6 * std::pow(inv_r, 6);
                double inv_r12 = inv_r6*inv_r6;

                double fmag = 24.0*epsilon*inv_r*(2.0*inv_r12-inv_r6);

                double fx = fmag*dx*inv_r;
                double fy = fmag*dy*inv_r;
                double fz = fmag*dz*inv_r;

                f[i][0] -= fx;
                f[i][1] -= fy;
                f[i][2] -= fz;
                f[j][0] += fx;
                f[j][1] += fy;
                f[j][2] += fz;
            }
        }
    }
    return f;
}


std::vector<std::vector<double>> total_force(const std::vector<std::vector<double>>& r, double b0, double epsilon, double sigma){
    std::vector<std::vector<double>> wca_forces_val = wca_force(r, epsilon, sigma) ;
    std::vector<std::vector<double>> bond_forces_val = bond_force(r, K_bond, b0).first ;
    std::vector<std::vector<double>> net_forces(
        wca_forces_val.size(),
        std::vector<double>(wca_forces_val[0].size(), 0.0));
    for (size_t i = 0; i < wca_forces_val.size(); i++)
        for (size_t j = 0; j < wca_forces_val[i].size(); j++)
            net_forces[i][j] = bond_forces_val[i][j] + wca_forces_val[i][j];

    return net_forces ;
}


// ##########################################################################
// ############ Grønbech-Jensen–Farago (GJF) Langevin integrator #############
// ##########################################################################
struct GjfResult {
    std::vector<std::vector<std::vector<double>>> r_hist;
    std::vector<std::vector<std::vector<double>>> u2gj_hist;
};

GjfResult run_2gj_integrator(
    const std::vector<std::vector<double>>& r0,
    const std::vector<std::vector<double>>& v0,
    int n_steps,
    const std::vector<std::vector<std::vector<double>>>& beta_array,
    double k_bond, double b0, double epsilon, double sigma,
    double dt, double m, double a, double b,
    int record_every)
{

    // copies of initial state
    std::vector<std::vector<double>> r = r0;
    std::vector<std::vector<double>> v = v0;
    int n_part = r0.size(); // redundant 

    int n_rec = n_steps / record_every;

    // allocate histories: (n_rec + 1, n_part, 3)
    std::vector<std::vector<std::vector<double>>> r_hist(
        n_rec + 1, std::vector<std::vector<double>>(n_part, std::vector<double>(3)));
    std::vector<std::vector<std::vector<double>>> u2gj_hist(
        n_rec + 1, std::vector<std::vector<double>>(n_part, std::vector<double>(3)));

    std::vector<std::vector<double>> f_n = total_force(r, b0, epsilon, sigma);

    // r_hist[0] = r.copy()
    for (int i = 0; i < n_part; i++)
        for (int j = 0; j < 3; j++)
            r_hist[0][i][j] = r[i][j];

    int rec_idx = 1;
    double sqrt_b = std::sqrt(b);

    for (int n = 0; n < n_steps; n++) {
        const auto& beta_np1 = beta_array[n];

        // Position update
        std::vector<std::vector<double>> r_new(n_part, std::vector<double>(3));
        for (int i = 0; i < n_part; i++) {
            for (int j = 0; j < 3; j++) {
                r_new[i][j] = r[i][j] + b * (
                    dt * v[i][j]
                    + (dt * dt / (2.0 * m)) * f_n[i][j]
                    + (dt / (2.0 * m)) * beta_np1[i][j]
                );
            }
        }

        std::vector<std::vector<double>> f_np1 =
            total_force(r_new, b0, epsilon, sigma);

        // Velocity update
        std::vector<std::vector<double>> v_new(n_part, std::vector<double>(3));
        for (int i = 0; i < n_part; i++) {
            for (int j = 0; j < 3; j++) {
                v_new[i][j] = a * v[i][j]
                    + (dt / (2.0 * m)) * (a * f_n[i][j] + f_np1[i][j])
                    + (b / m) * beta_np1[i][j];
            }
        }

        // Effective velocity u2gj = (r_new - r) / (sqrt_b * dt)
        std::vector<std::vector<double>> u2gj(n_part, std::vector<double>(3));
        for (int i = 0; i < n_part; i++) {
            for (int j = 0; j < 3; j++) {
                u2gj[i][j] = (r_new[i][j] - r[i][j]) / (sqrt_b * dt);
            }
        }

        r = std::move(r_new);
        v = std::move(v_new);
        f_n = std::move(f_np1);

        if ((n + 1) % record_every == 0) {
            for (int i = 0; i < n_part; i++) {
                for (int j = 0; j < 3; j++) {
                    r_hist[rec_idx][i][j] = r[i][j];
                    u2gj_hist[rec_idx][i][j] = u2gj[i][j];
                }
            }
            rec_idx++;
        }
    }

    return { std::move(r_hist), std::move(u2gj_hist) };
}

int main(){
    std::cout << "n_steps values is " << n_steps << std::endl ;
    std::cout << "n_part values is " << n_part << std::endl ;
    auto start_time = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<double>> r0 = initialize_polymer_chain(n_part, bond_length, Center);
    std::vector<std::vector<double>> v0 = Initialize_Velocity(n_part, K_Bt);
    std::vector<std::vector<double>> f0 = total_force(r0, b0, epsilon, sigma);
    std::vector<std::vector<std::vector<double>>> beta_array = noise_3D(n_steps, n_part);
    GjfResult out = run_2gj_integrator(r0, v0, n_steps, beta_array, K_bond, b0, epsilon, sigma, dt, m , a, b, record_every);

    // out is your GjfResult
    const auto& trajectories = out.r_hist;
    const auto& vel_frame    = out.u2gj_hist;

    int n_frames = trajectories.size();
    int n_part   = trajectories[0].size();

    int burn = n_frames / 5;

    // ---------- Kinetic energy ----------
    double KE_sum = 0.0;
    int    KE_count = 0;
    for (int t = burn; t < n_frames; t++) {
        double s = 0.0;
        for (int i = 0; i < n_part; i++)
            for (int j = 0; j < 3; j++)
                s += vel_frame[t][i][j] * vel_frame[t][i][j];
        KE_sum += 0.5 * m * s;
        KE_count++;
    }
    double mean_KE = (KE_sum / (KE_count)) / n_part;
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "<E_k> per DOF (2GJ velocity) = " << mean_KE / 3.0
            << "   (theory: " << 0.5 * K_Bt << ")\n";

    // ---------- End-to-end distance ----------
    double ee_sum = 0.0;
    double ee_count = 0;
    for (int t = burn; t < n_frames; t++) {
        double dx = trajectories[t][0][0] - trajectories[t][n_part-1][0];
        double dy = trajectories[t][0][1] - trajectories[t][n_part-1][1];
        double dz = trajectories[t][0][2] - trajectories[t][n_part-1][2];
        ee_sum += std::sqrt(dx*dx + dy*dy + dz*dz);
        ee_count++;
    }
    double end_end_dist_sq = ee_sum / ee_count;
    std::cout << "end_end_dist_sq " << end_end_dist_sq << "\n";

    // ---------- Radius of gyration ----------
    double ROG_sum = 0.0;
    int    ROG_count = 0;
    for (int t = burn; t < n_frames; t++) {
        double comx = 0.0, comy = 0.0, comz = 0.0;
        for (int i = 0; i < n_part; i++) {
            comx += trajectories[t][i][0];
            comy += trajectories[t][i][1];
            comz += trajectories[t][i][2];
        }
        comx /= n_part; comy /= n_part; comz /= n_part;

        double sum_sq = 0.0;
        for (int i = 0; i < n_part; i++) {
            double dx = trajectories[t][i][0] - comx;
            double dy = trajectories[t][i][1] - comy;
            double dz = trajectories[t][i][2] - comz;
            sum_sq += dx*dx + dy*dy + dz*dz;
        }
        ROG_sum += sum_sq / n_part;   // mean over particles
        ROG_count++;
    }
    double ROG_mean = ROG_sum / ROG_count;
    std::cout << "radius of gyration is " << ROG_mean << "\n";

    // ---------- Bond lengths ----------
    double bond_sum = 0.0;
    int    bond_count = 0;
    for (int t = burn; t < n_frames; t++) {
        for (int i = 1; i < n_part; i++) {
            double dx = trajectories[t][i][0] - trajectories[t][i-1][0];
            double dy = trajectories[t][i][1] - trajectories[t][i-1][1];
            double dz = trajectories[t][i][2] - trajectories[t][i-1][2];
            bond_sum += std::sqrt(dx*dx + dy*dy + dz*dz);
            bond_count++;
        }
    }
    double bond_mean = bond_sum / bond_count;
    std::cout << "bond_lengths " << bond_mean << "\n";
    std::cout << "no of particle is " << n_part << "\n";
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> execution_time = end_time - start_time;
    std::cout << "Process completed successfully." << std::endl;
    std::cout << "Execution time: " << std::fixed << std::setprecision(4) << execution_time.count() << std::endl;
}