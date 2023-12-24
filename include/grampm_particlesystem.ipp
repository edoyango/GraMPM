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
        particle<F> p(m_p_x[i], m_p_y[i], m_p_z[i], m_p_vx[i], m_p_vy[i], m_p_vz[i], m_p_mass[i], m_p_rho[i], 
            m_p_sigmaxx[i], m_p_sigmayy[i], m_p_sigmazz[i], m_p_sigmaxy[i], m_p_sigmaxz[i], m_p_sigmayz[i], m_p_ax[i], 
            m_p_ay[i], m_p_az[i], m_p_dxdt[i], m_p_dydt[i], m_p_dzdt[i], m_p_strainratexx[i], m_p_strainrateyy[i], 
            m_p_strainratezz[i], m_p_strainratexy[i], m_p_strainratexz[i], m_p_strainrateyz[i], m_p_spinratexy[i], 
            m_p_spinratexz[i], m_p_spinrateyz[i]);
        return p; 
    }

    // vector-like api: push_back. Takes particle class and appends its properties to particle_system member vectors
    template<typename F>
    void MPM_system<F>::p_push_back(const particle<F> &p) {
        m_p_x.push_back(p.x);
        m_p_y.push_back(p.y);
        m_p_z.push_back(p.z);
        m_p_vx.push_back(p.vx);
        m_p_vy.push_back(p.vy);
        m_p_vz.push_back(p.vz);
        m_p_ax.push_back(p.ax);
        m_p_ay.push_back(p.ay);
        m_p_az.push_back(p.az);
        m_p_dxdt.push_back(p.dxdt);
        m_p_dydt.push_back(p.dydt);
        m_p_dzdt.push_back(p.dzdt);
        m_p_mass.push_back(p.mass);
        m_p_rho.push_back(p.rho);
        m_p_sigmaxx.push_back(p.sigmaxx);
        m_p_sigmayy.push_back(p.sigmayy);
        m_p_sigmazz.push_back(p.sigmazz);
        m_p_sigmaxy.push_back(p.sigmaxy);
        m_p_sigmaxz.push_back(p.sigmaxz);
        m_p_sigmayz.push_back(p.sigmayz);
        m_p_strainratexx.push_back(p.strainratexx);
        m_p_strainrateyy.push_back(p.strainrateyy);
        m_p_strainratezz.push_back(p.strainratezz);
        m_p_strainratexy.push_back(p.strainratexy);
        m_p_strainratexz.push_back(p.strainratexz);
        m_p_strainrateyz.push_back(p.strainrateyz);
        m_p_spinratexy.push_back(p.spinratexy);
        m_p_spinratexz.push_back(p.spinratexz);
        m_p_spinrateyz.push_back(p.spinrateyz);
        m_p_grid_idx.push_back(
            ravel_grid_idx(
                calc_idxx(p.x),
                calc_idxy(p.y),
                calc_idxz(p.z)
            )
        );
        m_p_size++;
    }

    // vector-like api: clear. Makes size 0.
    template<typename F>
    void MPM_system<F>::p_clear() {
        m_p_x.clear();
        m_p_y.clear();
        m_p_z.clear();
        m_p_vx.clear();
        m_p_vy.clear();
        m_p_vz.clear();
        m_p_ax.clear();
        m_p_ay.clear();
        m_p_az.clear();
        m_p_dxdt.clear();
        m_p_dydt.clear();
        m_p_dzdt.clear();
        m_p_mass.clear();
        m_p_rho.clear();
        m_p_sigmaxx.clear();
        m_p_sigmayy.clear();
        m_p_sigmazz.clear();
        m_p_sigmaxy.clear();
        m_p_sigmaxz.clear();
        m_p_sigmayz.clear();
        m_p_strainratexx.clear();
        m_p_strainrateyy.clear();
        m_p_strainratezz.clear();
        m_p_strainratexy.clear();
        m_p_strainratexz.clear();
        m_p_strainrateyz.clear();
        m_p_spinratexy.clear();
        m_p_spinratexz.clear();
        m_p_spinrateyz.clear();
        m_p_grid_idx.clear();
        m_p_size = 0;
    }

    // vector-like api: empty. Checks whether all member vectors are empty. Also checks the m_size member
    template<typename F>
    bool MPM_system<F>::p_empty() {
        return m_p_x.empty() && m_p_y.empty() && m_p_z.empty() && m_p_vx.empty() && m_p_vy.empty() && m_p_vz.empty() && 
            m_p_ax.empty() && m_p_ay.empty() && m_p_az.empty() && m_p_dxdt.empty() && m_p_dydt.empty() && 
            m_p_dzdt.empty() && m_p_mass.empty() && m_p_rho.empty() && m_p_sigmaxx.empty() && m_p_sigmayy.empty() && 
            m_p_sigmazz.empty() && m_p_sigmaxy.empty() && m_p_sigmaxz.empty() && m_p_sigmayz.empty() && 
            m_p_strainratexx.empty() && m_p_strainrateyy.empty() && m_p_strainratezz.empty() && 
            m_p_strainratexy.empty() && m_p_strainratexz.empty() && m_p_strainrateyz.empty() && 
            m_p_spinratexy.empty() && m_p_spinratexz.empty() && m_p_spinrateyz.empty() &&m_p_grid_idx.empty() && 
            m_p_size==0;
    }

    // vector-like api: resize. Shrinks/grows the member vectors. zeros the vectors.
    template<typename F>
    void MPM_system<F>::p_resize(const int n) {
        m_p_x.resize(n, 0.);
        m_p_y.resize(n, 0.);
        m_p_z.resize(n, 0.);
        m_p_vx.resize(n, 0.);
        m_p_vy.resize(n, 0.);
        m_p_vz.resize(n, 0.);
        m_p_ax.resize(n, 0.);
        m_p_ay.resize(n, 0.);
        m_p_az.resize(n, 0.);
        m_p_dxdt.resize(n, 0.);
        m_p_dydt.resize(n, 0.);
        m_p_dzdt.resize(n, 0.);
        m_p_mass.resize(n, 0.);
        m_p_rho.resize(n, 0.);
        m_p_sigmaxx.resize(n, 0.);
        m_p_sigmayy.resize(n, 0.);
        m_p_sigmazz.resize(n, 0.);
        m_p_sigmaxy.resize(n, 0.);
        m_p_sigmaxz.resize(n, 0.);
        m_p_sigmayz.resize(n, 0.);
        m_p_strainratexx.resize(n, 0.);
        m_p_strainrateyy.resize(n, 0.);
        m_p_strainratezz.resize(n, 0.);
        m_p_strainratexy.resize(n, 0.);
        m_p_strainratexz.resize(n, 0.);
        m_p_strainrateyz.resize(n, 0.);
        m_p_spinratexy.resize(n, 0.);
        m_p_spinratexz.resize(n, 0.);
        m_p_spinrateyz.resize(n, 0.);
        m_p_grid_idx.resize(n, 0);
        m_p_size = n;
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
        for (size_t i = 0; i < m_p_size; ++i) {
            m_p_vx[i] += dt*m_p_ax[i];
            m_p_vy[i] += dt*m_p_ay[i];
            m_p_vz[i] += dt*m_p_az[i];
        }
    }
    template<typename F> 
    void MPM_system<F>::p_update_position(const F &dt) {
        #pragma omp for
        for (size_t i = 0; i < m_p_size; ++i) {
            m_p_x[i] += dt*m_p_dxdt[i];
            m_p_y[i] += dt*m_p_dydt[i];
            m_p_z[i] += dt*m_p_dzdt[i];
        }
    }
    template<typename F> 
    void MPM_system<F>::p_update_density(const F &dt) {
        // update density using volumentric strain increment
        #pragma omp for
        for (size_t i = 0; i < m_p_size; ++i) {
            m_p_rho[i] /= 1. + dt*(m_p_strainratexx[i] + m_p_strainrateyy[i] + m_p_strainratezz[i]);
        }
    }

    template<typename F>
    std::array<size_t, 3> MPM_system<F>::p_unravelled_grid_idx(const size_t &i) const { 
        return unravel_grid_idx(m_p_grid_idx[i]); 
    }
}

#endif