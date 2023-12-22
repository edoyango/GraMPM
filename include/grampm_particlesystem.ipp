#ifndef GRAMPM_particlesystem_ipp
#define GRAMPM_particlesystem_ipp

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

constexpr double pi {3.14159265358979311599796346854};

namespace GraMPM {

    // vector-like api: at. Returns particle class
    template<typename F>
    particle<F> MPM_system<F>::p_at(const int &i) { 
        particle<F> p(p_x[i], p_y[i], p_z[i], p_vx[i], p_vy[i], p_vz[i], p_mass[i], p_rho[i], p_sigmaxx[i], p_sigmayy[i], p_sigmazz[i], 
            p_sigmaxy[i], p_sigmaxz[i], p_sigmayz[i], p_ax[i], p_ay[i], p_az[i], p_dxdt[i], p_dydt[i], p_dzdt[i], p_strainratexx[i], 
            p_strainrateyy[i], p_strainratezz[i], p_strainratexy[i], p_strainratexz[i], p_strainrateyz[i], p_spinratexy[i], 
            p_spinratexz[i], p_spinrateyz[i]);
        return p; 
    }

    // vector-like api: push_back. Takes particle class and appends its properties to particle_system member vectors
    template<typename F>
    void MPM_system<F>::p_push_back(const particle<F> &p) {
        p_x.push_back(p.x);
        p_y.push_back(p.y);
        p_z.push_back(p.z);
        p_vx.push_back(p.vx);
        p_vy.push_back(p.vy);
        p_vz.push_back(p.vz);
        p_ax.push_back(p.ax);
        p_ay.push_back(p.ay);
        p_az.push_back(p.az);
        p_dxdt.push_back(p.dxdt);
        p_dydt.push_back(p.dydt);
        p_dzdt.push_back(p.dzdt);
        p_mass.push_back(p.mass);
        p_rho.push_back(p.rho);
        p_sigmaxx.push_back(p.sigmaxx);
        p_sigmayy.push_back(p.sigmayy);
        p_sigmazz.push_back(p.sigmazz);
        p_sigmaxy.push_back(p.sigmaxy);
        p_sigmaxz.push_back(p.sigmaxz);
        p_sigmayz.push_back(p.sigmayz);
        p_strainratexx.push_back(p.strainratexx);
        p_strainrateyy.push_back(p.strainrateyy);
        p_strainratezz.push_back(p.strainratezz);
        p_strainratexy.push_back(p.strainratexy);
        p_strainratexz.push_back(p.strainratexz);
        p_strainrateyz.push_back(p.strainrateyz);
        p_spinratexy.push_back(p.spinratexy);
        p_spinratexz.push_back(p.spinratexz);
        p_spinrateyz.push_back(p.spinrateyz);
        p_grid_idx.push_back(
            ravel_grid_idx(
                calc_idxx(p.x),
                calc_idxy(p.y),
                calc_idxz(p.z)
            )
        );
        p_size++;
    }

    // vector-like api: clear. Makes size 0.
    template<typename F>
    void MPM_system<F>::p_clear() {
        p_x.clear();
        p_y.clear();
        p_z.clear();
        p_vx.clear();
        p_vy.clear();
        p_vz.clear();
        p_ax.clear();
        p_ay.clear();
        p_az.clear();
        p_dxdt.clear();
        p_dydt.clear();
        p_dzdt.clear();
        p_mass.clear();
        p_rho.clear();
        p_sigmaxx.clear();
        p_sigmayy.clear();
        p_sigmazz.clear();
        p_sigmaxy.clear();
        p_sigmaxz.clear();
        p_sigmayz.clear();
        p_strainratexx.clear();
        p_strainrateyy.clear();
        p_strainratezz.clear();
        p_strainratexy.clear();
        p_strainratexz.clear();
        p_strainrateyz.clear();
        p_spinratexy.clear();
        p_spinratexz.clear();
        p_spinrateyz.clear();
        p_grid_idx.clear();
        p_size = 0;
    }

    // vector-like api: empty. Checks whether all member vectors are empty. Also checks the m_size member
    template<typename F>
    bool MPM_system<F>::p_empty() {
        return p_x.empty() && p_y.empty() && p_z.empty() && p_vx.empty() && p_vy.empty() && p_vz.empty() && 
            p_ax.empty() && p_ay.empty() && p_az.empty() && p_dxdt.empty() && p_dydt.empty() && p_dzdt.empty() &&
            p_mass.empty() && p_rho.empty() && p_sigmaxx.empty() && p_sigmayy.empty() && p_sigmazz.empty() && 
            p_sigmaxy.empty() && p_sigmaxz.empty() && p_sigmayz.empty() && p_strainratexx.empty() && 
            p_strainrateyy.empty() && p_strainratezz.empty() && p_strainratexy.empty() && p_strainratexz.empty() && 
            p_strainrateyz.empty() && p_spinratexy.empty() && p_spinratexz.empty() && 
            p_spinrateyz.empty() &&p_grid_idx.empty() && p_size==0;
    }

    // vector-like api: resize. Shrinks/grows the member vectors. zeros the vectors.
    template<typename F>
    void MPM_system<F>::p_resize(const int n) {
        p_x.resize(n, 0.);
        p_y.resize(n, 0.);
        p_z.resize(n, 0.);
        p_vx.resize(n, 0.);
        p_vy.resize(n, 0.);
        p_vz.resize(n, 0.);
        p_ax.resize(n, 0.);
        p_ay.resize(n, 0.);
        p_az.resize(n, 0.);
        p_dxdt.resize(n, 0.);
        p_dydt.resize(n, 0.);
        p_dzdt.resize(n, 0.);
        p_mass.resize(n, 0.);
        p_rho.resize(n, 0.);
        p_sigmaxx.resize(n, 0.);
        p_sigmayy.resize(n, 0.);
        p_sigmazz.resize(n, 0.);
        p_sigmaxy.resize(n, 0.);
        p_sigmaxz.resize(n, 0.);
        p_sigmayz.resize(n, 0.);
        p_strainratexx.resize(n, 0.);
        p_strainrateyy.resize(n, 0.);
        p_strainratezz.resize(n, 0.);
        p_strainratexy.resize(n, 0.);
        p_strainratexz.resize(n, 0.);
        p_strainrateyz.resize(n, 0.);
        p_spinratexy.resize(n, 0.);
        p_spinratexz.resize(n, 0.);
        p_spinrateyz.resize(n, 0.);
        p_grid_idx.resize(n, 0);
        p_size = n;
    }

    // property getters
    template<typename F> 
    void MPM_system<F>::DP_params(F& phi, F& psi, F& coh) const { phi = m_phi; psi = m_psi; coh = m_coh; }
    template<typename F> 
    void MPM_system<F>::DP_params(F& phi, F& psi, F& coh, F& alpha_phi, F& alpha_psi, F& k_c) const {
         phi = m_phi; psi = m_psi; coh = m_coh; 
         alpha_phi = m_alphaphi; alpha_psi = m_alphapsi; k_c = m_kc;
    }

    template<typename F> 
    void MPM_system<F>::set_stress_update_function(std::function<void(MPM_system<F>&, const F&)> f) { 
        p_stress_update_function = f; 
    }
    
    template<typename F> void MPM_system<F>::set_DP_params(const F &phi, const F &psi, const F &coh) {
        m_phi = phi; m_psi = psi; m_coh = coh;
        // const F denom = std::sqrt(9.+12.*std::tan(phi)*std::tan(phi));
        // m_alphaphi = 3.*std::tan(phi)/denom;
        // m_alphapsi = 3.*std::tan(psi)/denom;
        // m_kc = 3.*coh/denom;
        const F denom = std::sqrt(3.)*(3.-std::sin(phi));
        m_alphaphi = 2.*std::sin(phi)/denom;
        m_alphapsi = 2.*std::sin(psi)/denom;
        m_kc = 6.*coh*std::cos(phi)/denom;
    }

    template<typename F> 
    void MPM_system<F>::p_update_stress(const F &dt) {
        p_stress_update_function(*this, dt);
    }

    template<typename F> 
    void MPM_system<F>::p_update_velocity(const F &dt) {
        #pragma omp for
        for (int i = 0; i < p_size; ++i) {
            p_vx[i] += dt*p_ax[i];
            p_vy[i] += dt*p_ay[i];
            p_vz[i] += dt*p_az[i];
        }
    }
    template<typename F> 
    void MPM_system<F>::p_update_position(const F &dt) {
        #pragma omp for
        for (int i = 0; i < p_size; ++i) {
            p_x[i] += dt*p_dxdt[i];
            p_y[i] += dt*p_dydt[i];
            p_z[i] += dt*p_dzdt[i];
        }
    }
    template<typename F> 
    void MPM_system<F>::p_update_density(const F &dt) {
        // update density using volumentric strain increment
        #pragma omp for
        for (int i = 0; i < p_size; ++i) {
            p_rho[i] /= 1. + dt*(p_strainratexx[i] + p_strainrateyy[i] + p_strainratezz[i]);
        }
    }

    template<typename F>
    std::array<int, 3> MPM_system<F>::p_unravelled_grid_idx(const int &i) const { 
        return unravel_grid_idx(p_grid_idx[i]); 
    }
}

#endif