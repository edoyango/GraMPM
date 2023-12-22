#ifndef GRAMPM_particlesystem_ipp
#define GRAMPM_particlesystem_ipp

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

constexpr double pi {3.14159265358979311599796346854};

namespace GraMPM {

    // initialize a particle with size, but zerod properties
    template<typename F>
    particle_system<F>::particle_system(const long unsigned int size, std::array<F, 3> bf, grid<F> &ingrid, 
        kernel_base<F> &knl)
        : m_x(size, 0.)
        , m_y(size, 0.)
        , m_z(size, 0.)
        , m_vx(size, 0.)
        , m_vy(size, 0.)
        , m_vz(size, 0.)
        , m_ax(size, 0.)
        , m_ay(size, 0.)
        , m_az(size, 0.)
        , m_dxdt(size, 0.)
        , m_dydt(size, 0.)
        , m_dzdt(size, 0.)
        , m_mass(size, 0.)
        , m_rho(size, 0.)
        , m_sigmaxx(size, 0.)
        , m_sigmayy(size, 0.)
        , m_sigmazz(size, 0.)
        , m_sigmaxy(size, 0.)
        , m_sigmaxz(size, 0.)
        , m_sigmayz(size, 0.)
        , m_strainratexx(size, 0.)
        , m_strainrateyy(size, 0.)
        , m_strainratezz(size, 0.)
        , m_strainratexy(size, 0.)
        , m_strainratexz(size, 0.)
        , m_strainrateyz(size, 0.)
        , m_spinratexy(size, 0.)
        , m_spinratexz(size, 0.)
        , m_spinrateyz(size, 0.)
        , m_grid_idx(size, 0)
        , background_grid(ingrid)
        , m_capacity {size}
        , m_size {size}
        , m_body_force {bf}
        , m_knl {knl}
        , m_nneighbour_nodes_perp {static_cast<int>(8*std::ceil(knl.radius)*std::ceil(knl.radius)*std::ceil(knl.radius))}
    {
    }

    // initialize from std::vector of particle class
    template<typename F>
    particle_system<F>::particle_system(const std::vector<particle<F>> &pv, std::array<F, 3> bf, grid<F> &ingrid, kernel_base<F> &knl)
        : background_grid(ingrid)
        , m_capacity {pv.capacity()}
        , m_size {0}
        , m_body_force {bf}
        , m_knl {knl}
        , m_nneighbour_nodes_perp {static_cast<int>(8*std::ceil(knl.radius)*std::ceil(knl.radius)*std::ceil(knl.radius))}
    {
        for (int i = 0; i < pv.size(); ++i) push_back(pv[i]);
    }

    // initalize empty (no size)
    template<typename F>
    particle_system<F>::particle_system(grid<F> &ingrid, kernel_base<F> &knl)
        : background_grid(ingrid)
        , m_capacity {0}
        , m_size {0}
        , m_body_force {0., 0., 0.}
        , m_knl {knl}
        , m_nneighbour_nodes_perp {static_cast<int>(8*std::ceil(knl.radius)*std::ceil(knl.radius)*std::ceil(knl.radius))}
    {
    }

    // helper function to ravel idx
    template<typename F>
    int particle_system<F>::ravel_grid_idx(const int &idxx, const int &idxy, const int &idxz) const {
        return idxx*background_grid.m_ngridy*background_grid.m_ngridz + idxy*background_grid.m_ngridz + idxz;
    }
    
    // helper function to unravel index (return std::array)
    template<typename F>
    std::array<int, 3> particle_system<F>::unravel_grid_idx(const int &idx) const {
        std::array<int, 3> unravelled_idx;
        div_t tmp = std::div(idx, background_grid.m_ngridy*background_grid.m_ngridz);
        unravelled_idx[0] = tmp.quot;
        tmp = std::div(tmp.rem, background_grid.m_ngridz);
        unravelled_idx[1] = tmp.quot;
        unravelled_idx[2] = tmp.rem;
        return unravelled_idx;
    }

    // helper function to unravel index (modify args)
    template<typename F>
    void particle_system<F>::unravel_grid_idx(const int &idx, int &idxx, int &idxy, int &idxz) const {
        div_t tmp = std::div(idx, background_grid.m_ngridy*background_grid.m_ngridz);
        idxx = tmp.quot;
        tmp = std::div(tmp.rem, background_grid.m_ngridz);
        idxy = tmp.quot;
        idxz = tmp.rem;
    }

    // grid getter
    template<typename F>
    const grid<F>* particle_system<F>::grid_address() { return &background_grid; }

    // property getters
    template<typename F> void particle_system<F>::DP_params(F& phi, F& psi, F& coh) const { phi = m_phi; psi = m_psi; coh = m_coh; }
    template<typename F> void particle_system<F>::DP_params(F& phi, F& psi, F& coh, F& alpha_phi, F& alpha_psi, F& k_c) const {
         phi = m_phi; psi = m_psi; coh = m_coh; 
         alpha_phi = m_alphaphi; alpha_psi = m_alphapsi; k_c = m_kc;
    }
    template<typename F>
    std::array<int, 3> particle_system<F>::grid_idx(const int &i) const { return unravel_grid_idx(m_grid_idx[i]); }

    template<typename F>
    const int& particle_system<F>::ps_nn(const int i, const int j) const { 
        assert(j < m_nneighbour_nodes_perp); 
        return m_ps_nns[i*m_nneighbour_nodes_perp+j];
    }

    template<typename F>
    const F& particle_system<F>::ps_nn_dx(const int i, const int j) const {
        assert(j < m_nneighbour_nodes_perp);
        return m_ps_nns_dx[i*m_nneighbour_nodes_perp+j];
    }

    template<typename F>
    const F& particle_system<F>::ps_nn_dy(const int i, const int j) const {
        assert(j < m_nneighbour_nodes_perp);
        return m_ps_nns_dy[i*m_nneighbour_nodes_perp+j];
    }

    template<typename F>
    const F& particle_system<F>::ps_nn_dz(const int i, const int j) const {
        assert(j < m_nneighbour_nodes_perp);
        return m_ps_nns_dz[i*m_nneighbour_nodes_perp+j];
    }
    
    template<typename F>
    const F& particle_system<F>::ps_nn_w(const int i, const int j) const {
        assert(j < m_nneighbour_nodes_perp);
        return m_ps_nns_w[i*m_nneighbour_nodes_perp+j];
    }
    
    template<typename F>
    const F& particle_system<F>::ps_nn_dwdx(const int i, const int j) const {
        assert(j < m_nneighbour_nodes_perp);
        return m_ps_nns_dwdx[i*m_nneighbour_nodes_perp+j];
    }
    
    template<typename F>
    const F& particle_system<F>::ps_nn_dwdy(const int i, const int j) const {
        assert(j < m_nneighbour_nodes_perp);
        return m_ps_nns_dwdy[i*m_nneighbour_nodes_perp+j];
    }

    template<typename F>
    const F& particle_system<F>::ps_nn_dwdz(const int i, const int j) const {
        assert(j < m_nneighbour_nodes_perp);
        return m_ps_nns_dwdz[i*m_nneighbour_nodes_perp+j];
    }
    template<typename F> const long unsigned int& particle_system<F>::capacity() const { return m_capacity; }

    
    template<typename F> void particle_system<F>::set_stress_update_function(std::function<void(particle_system<F>&, const F&)> f) { m_stress_update_function = f; }
    template<typename F> void particle_system<F>::set_DP_params(const F &phi, const F &psi, const F &coh) {
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

    // vector-like api: at. Returns particle class
    template<typename F>
    particle<F> particle_system<F>::at(const int &i) { 
        particle<F> p(m_x[i], m_y[i], m_z[i], m_vx[i], m_vy[i], m_vz[i], m_mass[i], m_rho[i], m_sigmaxx[i], m_sigmayy[i], m_sigmazz[i], 
            m_sigmaxy[i], m_sigmaxz[i], m_sigmayz[i], m_ax[i], m_ay[i], m_az[i], m_dxdt[i], m_dydt[i], m_dzdt[i], m_strainratexx[i], 
            m_strainrateyy[i], m_strainratezz[i], m_strainratexy[i], m_strainratexz[i], m_strainrateyz[i], m_spinratexy[i], 
            m_spinratexz[i], m_spinrateyz[i]);
        return p; 
    }

    // vector-like api: push_back. Takes particle class and appends its properties to particle_system member vectors
    template<typename F>
    void particle_system<F>::push_back(const particle<F> &p) {
        m_x.push_back(p.x);
        m_y.push_back(p.y);
        m_z.push_back(p.z);
        m_vx.push_back(p.vx);
        m_vy.push_back(p.vy);
        m_vz.push_back(p.vz);
        m_ax.push_back(p.ax);
        m_ay.push_back(p.ay);
        m_az.push_back(p.az);
        m_dxdt.push_back(p.dxdt);
        m_dydt.push_back(p.dydt);
        m_dzdt.push_back(p.dzdt);
        m_mass.push_back(p.mass);
        m_rho.push_back(p.rho);
        m_sigmaxx.push_back(p.sigmaxx);
        m_sigmayy.push_back(p.sigmayy);
        m_sigmazz.push_back(p.sigmazz);
        m_sigmaxy.push_back(p.sigmaxy);
        m_sigmaxz.push_back(p.sigmaxz);
        m_sigmayz.push_back(p.sigmayz);
        m_strainratexx.push_back(p.strainratexx);
        m_strainrateyy.push_back(p.strainrateyy);
        m_strainratezz.push_back(p.strainratezz);
        m_strainratexy.push_back(p.strainratexy);
        m_strainratexz.push_back(p.strainratexz);
        m_strainrateyz.push_back(p.strainrateyz);
        m_spinratexy.push_back(p.spinratexy);
        m_spinratexz.push_back(p.spinratexz);
        m_spinrateyz.push_back(p.spinrateyz);
        m_grid_idx.push_back(
            ravel_grid_idx(
                background_grid.calc_idxx(p.x),
                background_grid.calc_idxy(p.y),
                background_grid.calc_idxz(p.z)
            )
        );
        m_size++;
    }

    // vector-like api: reserve. Reserves space for each member vector
    template<typename F>
    void particle_system<F>::reserve(const long unsigned int &n) {
        m_x.reserve(n);
        m_y.reserve(n);
        m_z.reserve(n);
        m_vx.reserve(n);
        m_vy.reserve(n);
        m_vz.reserve(n);
        m_ax.reserve(n);
        m_ay.reserve(n);
        m_az.reserve(n);
        m_dxdt.reserve(n);
        m_dydt.reserve(n);
        m_dzdt.reserve(n);
        m_mass.reserve(n);
        m_rho.reserve(n);
        m_sigmaxx.reserve(n);
        m_sigmayy.reserve(n);
        m_sigmazz.reserve(n);
        m_sigmaxy.reserve(n);
        m_sigmaxz.reserve(n);
        m_sigmayz.reserve(n);
        m_strainratexx.reserve(n);
        m_strainrateyy.reserve(n);
        m_strainratezz.reserve(n);
        m_strainratexy.reserve(n);
        m_strainratexz.reserve(n);
        m_strainrateyz.reserve(n);
        m_spinratexy.reserve(n);
        m_spinratexz.reserve(n);
        m_spinrateyz.reserve(n);
        m_grid_idx.reserve(n);
        m_capacity = std::max(m_capacity, n);
    }

    // vector-like api: clear. Makes size 0.
    template<typename F>
    void particle_system<F>::clear() {
        m_x.clear();
        m_y.clear();
        m_z.clear();
        m_vx.clear();
        m_vy.clear();
        m_vz.clear();
        m_ax.clear();
        m_ay.clear();
        m_az.clear();
        m_dxdt.clear();
        m_dydt.clear();
        m_dzdt.clear();
        m_mass.clear();
        m_rho.clear();
        m_sigmaxx.clear();
        m_sigmayy.clear();
        m_sigmazz.clear();
        m_sigmaxy.clear();
        m_sigmaxz.clear();
        m_sigmayz.clear();
        m_strainratexx.clear();
        m_strainrateyy.clear();
        m_strainratezz.clear();
        m_strainratexy.clear();
        m_strainratexz.clear();
        m_strainrateyz.clear();
        m_spinratexy.clear();
        m_spinratexz.clear();
        m_spinrateyz.clear();
        m_grid_idx.clear();
        m_ps_nns.clear();
        m_ps_nns_dx.clear();
        m_ps_nns_dy.clear();
        m_ps_nns_dz.clear();
        m_ps_nns_w.clear();
        m_ps_nns_dwdx.clear();
        m_ps_nns_dwdy.clear();
        m_ps_nns_dwdz.clear();
        m_size = 0;
    }

    // vector-like api: empty. Checks whether all member vectors are empty. Also checks the m_size member
    template<typename F>
    bool particle_system<F>::empty() {
        return m_x.empty() && m_y.empty() && m_z.empty() && m_vx.empty() && m_vy.empty() && m_vz.empty() && 
            m_ax.empty() && m_ay.empty() && m_az.empty() && m_dxdt.empty() && m_dydt.empty() && m_dzdt.empty() &&
            m_mass.empty() && m_rho.empty() && m_sigmaxx.empty() && m_sigmayy.empty() && m_sigmazz.empty() && 
            m_sigmaxy.empty() && m_sigmaxz.empty() && m_sigmayz.empty() && m_strainratexx.empty() && 
            m_strainrateyy.empty() && m_strainratezz.empty() && m_strainratexy.empty() && m_strainratexz.empty() && 
            m_strainrateyz.empty() && m_spinratexy.empty() && m_spinratexz.empty() && 
            m_spinrateyz.empty() &&m_grid_idx.empty() && m_ps_nns.empty() && 
            m_ps_nns_dx.empty() && m_ps_nns_dy.empty() && m_ps_nns_dz.empty() && 
            m_ps_nns_w.empty() && m_ps_nns_dwdx.empty() &&
            m_ps_nns_dwdy.empty() && m_ps_nns_dwdz.empty() && m_size==0;
    }

    // vector-like api: resize. Shrinks/grows the member vectors. zeros the vectors.
    template<typename F>
    void particle_system<F>::resize(const int n) {
        m_x.resize(n, 0.);
        m_y.resize(n, 0.);
        m_z.resize(n, 0.);
        m_vx.resize(n, 0.);
        m_vy.resize(n, 0.);
        m_vz.resize(n, 0.);
        m_ax.resize(n, 0.);
        m_ay.resize(n, 0.);
        m_az.resize(n, 0.);
        m_dxdt.resize(n, 0.);
        m_dydt.resize(n, 0.);
        m_dzdt.resize(n, 0.);
        m_mass.resize(n, 0.);
        m_rho.resize(n, 0.);
        m_sigmaxx.resize(n, 0.);
        m_sigmayy.resize(n, 0.);
        m_sigmazz.resize(n, 0.);
        m_sigmaxy.resize(n, 0.);
        m_sigmaxz.resize(n, 0.);
        m_sigmayz.resize(n, 0.);
        m_strainratexx.resize(n, 0.);
        m_strainrateyy.resize(n, 0.);
        m_strainratezz.resize(n, 0.);
        m_strainratexy.resize(n, 0.);
        m_strainratexz.resize(n, 0.);
        m_strainrateyz.resize(n, 0.);
        m_spinratexy.resize(n, 0.);
        m_spinratexz.resize(n, 0.);
        m_spinrateyz.resize(n, 0.);
        m_grid_idx.resize(n, 0);
        m_size = n;
    }

    // vector-like api: resize. Shrinks/grows the member vectors. Sets each vector value as per p properties
    template<typename F>
    void particle_system<F>::resize(const int n, const particle<F> p) {
        m_x.resize(n, p.x);
        m_y.resize(n, p.y);
        m_z.resize(n, p.z);
        m_vx.resize(n, p.vx);
        m_vy.resize(n, p.vy);
        m_vz.resize(n, p.vz);
        m_ax.resize(n, p.ax);
        m_ay.resize(n, p.ay);
        m_az.resize(n, p.az);
        m_dxdt.resize(n, p.ax);
        m_dydt.resize(n, p.ay);
        m_dzdt.resize(n, p.az);
        m_mass.resize(n, p.mass);
        m_rho.resize(n, p.rho);
        m_sigmaxx.resize(n, p.sigmaxx);
        m_sigmayy.resize(n, p.sigmayy);
        m_sigmazz.resize(n, p.sigmazz);
        m_sigmaxy.resize(n, p.sigmaxy);
        m_sigmaxz.resize(n, p.sigmaxz);
        m_sigmayz.resize(n, p.sigmayz);
        m_strainratexx.resize(n, p.strainratexx);
        m_strainrateyy.resize(n, p.strainrateyy);
        m_strainratezz.resize(n, p.strainratezz);
        m_strainratexy.resize(n, p.strainratexy);
        m_strainratexz.resize(n, p.strainratexz);
        m_strainrateyz.resize(n, p.strainrateyz);
        m_spinratexy.resize(n, p.spinratexy);
        m_spinratexz.resize(n, p.spinratexz);
        m_spinrateyz.resize(n, p.spinrateyz);
        m_grid_idx.resize(n, ravel_grid_idx(
            background_grid.calc_idxx(p.x),
            background_grid.calc_idxy(p.y),
            background_grid.calc_idxz(p.z)
        ));
        m_size = n;
    }

    template<typename F>
    void particle_system<F>::update_particle_to_cell_map(const int &start, const int &end) {
        #pragma omp for
        for (int i = start; i < end; ++i) {
            m_grid_idx[i] = ravel_grid_idx(
                background_grid.calc_idxx(m_x[i]),
                background_grid.calc_idxy(m_y[i]),
                background_grid.calc_idxz(m_z[i])
            );
        }
    }

    template<typename F>
    void particle_system<F>::update_particle_to_cell_map() {
        #pragma omp for
        for (int i = 0; i < m_size; ++i) {
            m_grid_idx[i] = ravel_grid_idx(
                background_grid.calc_idxx(m_x[i]),
                background_grid.calc_idxy(m_y[i]),
                background_grid.calc_idxz(m_z[i])
            );
        }
    }
    
    // NTS this could be faster
    template<typename F>
    void particle_system<F>::map_particles_to_grid() {

        // size output arrays (this should probably be done before omp parallel)
        #pragma omp single
        {
            m_ps_nns.resize(m_nneighbour_nodes_perp*m_size);
            m_ps_nns_dx.resize(m_nneighbour_nodes_perp*m_size);
            m_ps_nns_dy.resize(m_nneighbour_nodes_perp*m_size);
            m_ps_nns_dz.resize(m_nneighbour_nodes_perp*m_size);
            m_ps_nns_w.resize(m_nneighbour_nodes_perp*m_size);
            m_ps_nns_dwdx.resize(m_nneighbour_nodes_perp*m_size);
            m_ps_nns_dwdy.resize(m_nneighbour_nodes_perp*m_size);
            m_ps_nns_dwdz.resize(m_nneighbour_nodes_perp*m_size);
        }

        // update neighbour indices
        #pragma omp for
        for (int i = 0; i < m_size; ++i) {
            const int idx = m_grid_idx[i];
            int n = 0;
            for (int di=1-m_knl.radius; di <= m_knl.radius; ++di) {
                for (int dj = 1-m_knl.radius; dj <= m_knl.radius; ++dj) {
                    for (int dk = 1-m_knl.radius; dk <= m_knl.radius; ++dk) {
                        m_ps_nns[i*m_nneighbour_nodes_perp+n] = idx + ravel_grid_idx(di, dj, dk);
                        n++;
                    }
                }
            }
        }

        // update kernel and kernel gradient values
        #pragma omp for
        for (int i = 0; i < m_size; ++i) {
            int jstart = i*m_nneighbour_nodes_perp;
            for (int j = jstart; j < jstart+m_nneighbour_nodes_perp; ++j) {
                std::array<int, 3> idx;
                idx = unravel_grid_idx(m_ps_nns[j]);
                m_ps_nns_dx[j] = m_x[i] - (idx[0]*background_grid.m_dcell + background_grid.m_mingridx);
                m_ps_nns_dy[j] = m_y[i] - (idx[1]*background_grid.m_dcell + background_grid.m_mingridy);
                m_ps_nns_dz[j] = m_z[i] - (idx[2]*background_grid.m_dcell + background_grid.m_mingridz);
                m_knl.w_dwdx(
                    m_ps_nns_dx[j],
                    m_ps_nns_dy[j],
                    m_ps_nns_dz[j],
                    m_ps_nns_w[j],
                    m_ps_nns_dwdx[j],
                    m_ps_nns_dwdy[j],
                    m_ps_nns_dwdz[j]
                );
            }
        }
    }

    #pragma omp declare reduction(vec_plus : std::vector<double> : \
        std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
        initializer(omp_priv = std::vector<double>(omp_orig.size()))

    template<typename F> void particle_system<F>::map_mass_to_grid() { 

        // initialize grid data. Temporary arrays used because omp reduction needs class member or variable local to scope
        #pragma omp for
        for (int i = 0; i < background_grid.m_ncells; ++i) m_tmpgmass[i] = 0.;

        // map mass to grid
        #pragma omp for reduction(vec_plus:m_tmpgmass)
        for (int i = 0; i < m_size; ++i) {
            for (int j = 0; j < m_nneighbour_nodes_perp; ++j) {
                const int node_idx = ps_nn(i, j);
                m_tmpgmass[node_idx] += ps_nn_w(i, j)*m_mass[i];
            }
        }

        // copy data from temporary array to real location (needs to be improved)
        #pragma omp single
        background_grid.m_mass = m_tmpgmass;
    }

    template<typename F> void particle_system<F>::map_momentum_to_grid() { 
        
        // initialize grid data. Temporary arrays used because omp reduction needs class member or variable local to scope
        #pragma omp for
        for (int i = 0; i < background_grid.m_ncells; ++i) {
            m_tmpgmomentumx[i] = 0.;
            m_tmpgmomentumy[i] = 0.;
            m_tmpgmomentumz[i] = 0.;
        }

        // map momentum to grid
        #pragma omp for reduction(vec_plus:m_tmpgmomentumx,m_tmpgmomentumy,m_tmpgmomentumz)
        for (int i = 0; i < m_size; ++i) {
            for (int j = 0; j < m_nneighbour_nodes_perp; ++j) {
                const int node_idx = ps_nn(i, j);
                m_tmpgmomentumx[node_idx] += ps_nn_w(i, j)*m_mass[i]*m_vx[i];
                m_tmpgmomentumy[node_idx] += ps_nn_w(i, j)*m_mass[i]*m_vy[i];
                m_tmpgmomentumz[node_idx] += ps_nn_w(i, j)*m_mass[i]*m_vz[i];
            }
        }

        // copy data from temporary array to real location (needs to be improved)
        #pragma omp single
        {
            background_grid.m_momentumx = m_tmpgmomentumx;
            background_grid.m_momentumy = m_tmpgmomentumy;
            background_grid.m_momentumz = m_tmpgmomentumz;
        }
    }

    template<typename F> void particle_system<F>::map_force_to_grid() {

        // initialize grid data. Temporary arrays used because omp reduction needs class member or variable local to scope
        #pragma omp for
        for (int i = 0; i < background_grid.m_ncells; ++i) {
            m_tmpgforcex[i] = 0.;
            m_tmpgforcey[i] = 0.;
            m_tmpgforcez[i] = 0.;
        }

        // map force to grid (div sigma and body force together)
        #pragma omp for reduction (vec_plus:m_tmpgforcex,m_tmpgforcey,m_tmpgforcez)
        for (int i = 0; i < m_size; ++i) {
            for (int j = 0; j < m_nneighbour_nodes_perp; ++j) {
                const int node_idx = ps_nn(i, j);
                m_tmpgforcex[node_idx] += -m_mass[i]/m_rho[i]*(
                    m_sigmaxx[i]*ps_nn_dwdx(i, j) +
                    m_sigmaxy[i]*ps_nn_dwdy(i, j) +
                    m_sigmaxz[i]*ps_nn_dwdz(i, j)
                ) + m_body_force[0]*m_mass[i]*ps_nn_w(i, j);
                m_tmpgforcey[node_idx] += -m_mass[i]/m_rho[i]*(
                    m_sigmaxy[i]*ps_nn_dwdx(i, j) +
                    m_sigmayy[i]*ps_nn_dwdy(i, j) +
                    m_sigmayz[i]*ps_nn_dwdz(i, j)
                ) + m_body_force[1]*m_mass[i]*ps_nn_w(i, j);
                m_tmpgforcez[node_idx] += -m_mass[i]/m_rho[i]*(
                    m_sigmaxz[i]*ps_nn_dwdx(i, j) +
                    m_sigmayz[i]*ps_nn_dwdy(i, j) +
                    m_sigmazz[i]*ps_nn_dwdz(i, j)
                ) + m_body_force[2]*m_mass[i]*ps_nn_w(i, j);
            }
        }

        // copy data from temporary array to real location (needs to be improved)
        #pragma omp single
        {
            background_grid.m_forcex = m_tmpgforcex;
            background_grid.m_forcey = m_tmpgforcey;
            background_grid.m_forcez = m_tmpgforcez;
        }

    }

    template<typename F> void particle_system<F>::map_acceleration_to_particles() {

        // initialize particle data
        #pragma omp for
        for (int i = 0; i < m_size; ++i) {
            m_ax[i] = 0.;
            m_ay[i] = 0.;
            m_az[i] = 0.;
            m_dxdt[i] = 0.;
            m_dydt[i] = 0.;
            m_dzdt[i] = 0.;
        }

        // map acceleration and velocity to particles
        #pragma omp for
        for (int i = 0; i < m_size; ++i) {
            for (int j = 0; j < m_nneighbour_nodes_perp; ++j) {
                const int node_idx = ps_nn(i, j);
                m_ax[i] += ps_nn_w(i, j)*background_grid.m_forcex[node_idx]/background_grid.m_mass[node_idx];
                m_ay[i] += ps_nn_w(i, j)*background_grid.m_forcey[node_idx]/background_grid.m_mass[node_idx];
                m_az[i] += ps_nn_w(i, j)*background_grid.m_forcez[node_idx]/background_grid.m_mass[node_idx];
                m_dxdt[i] += ps_nn_w(i, j)*background_grid.m_momentumx[node_idx]/background_grid.m_mass[node_idx];
                m_dydt[i] += ps_nn_w(i, j)*background_grid.m_momentumy[node_idx]/background_grid.m_mass[node_idx];
                m_dzdt[i] += ps_nn_w(i, j)*background_grid.m_momentumz[node_idx]/background_grid.m_mass[node_idx];
            }
        }
    }

    template<typename F> void particle_system<F>::map_strainrate_to_particles() {

        // initialize particle data
        #pragma omp for
        for (int i = 0; i < m_size; ++i) {
            m_strainratexx[i] = 0.;
            m_strainrateyy[i] = 0.;
            m_strainratezz[i] = 0.;
            m_strainratexy[i] = 0.;
            m_strainratexz[i] = 0.;
            m_strainrateyz[i] = 0.;
            m_spinratexy[i] = 0.;
            m_spinratexz[i] = 0.;
            m_spinrateyz[i] = 0.;
        }

        // map strainrates to particles
        #pragma omp for
        for (int i = 0; i < m_size; ++i) {
            for (int j = 0; j < m_nneighbour_nodes_perp; ++j) {
                const int node_idx = ps_nn(i, j);
                m_strainratexx[i] += ps_nn_dwdx(i, j)*background_grid.m_momentumx[node_idx]/background_grid.m_mass[node_idx];
                m_strainrateyy[i] += ps_nn_dwdy(i, j)*background_grid.m_momentumy[node_idx]/background_grid.m_mass[node_idx];
                m_strainratezz[i] += ps_nn_dwdz(i, j)*background_grid.m_momentumz[node_idx]/background_grid.m_mass[node_idx];
                m_strainratexy[i] += 0.5*(ps_nn_dwdx(i, j)*background_grid.m_momentumy[node_idx]/background_grid.m_mass[node_idx] + 
                    ps_nn_dwdy(i, j)*background_grid.m_momentumx[node_idx]/background_grid.m_mass[node_idx]);
                m_strainratexz[i] += 0.5*(ps_nn_dwdx(i, j)*background_grid.m_momentumz[node_idx]/background_grid.m_mass[node_idx] + 
                    ps_nn_dwdz(i, j)*background_grid.m_momentumx[node_idx]/background_grid.m_mass[node_idx]);
                m_strainrateyz[i] += 0.5*(ps_nn_dwdy(i, j)*background_grid.m_momentumz[node_idx]/background_grid.m_mass[node_idx] + 
                    ps_nn_dwdz(i, j)*background_grid.m_momentumy[node_idx]/background_grid.m_mass[node_idx]);
                m_spinratexy[i] += 0.5*(ps_nn_dwdy(i, j)*background_grid.m_momentumx[node_idx]/background_grid.m_mass[node_idx] -
                    ps_nn_dwdx(i, j)*background_grid.m_momentumy[node_idx]/background_grid.m_mass[node_idx]);
                m_spinratexz[i] += 0.5*(ps_nn_dwdz(i, j)*background_grid.m_momentumx[node_idx]/background_grid.m_mass[node_idx] - 
                    ps_nn_dwdx(i, j)*background_grid.m_momentumz[node_idx]/background_grid.m_mass[node_idx]);
                m_spinrateyz[i] += 0.5*(ps_nn_dwdz(i, j)*background_grid.m_momentumy[node_idx]/background_grid.m_mass[node_idx] - 
                    ps_nn_dwdy(i, j)*background_grid.m_momentumz[node_idx]/background_grid.m_mass[node_idx]);
            }
        }
    }

    template<typename F> void particle_system<F>::update_stress(const F &dt) {
        m_stress_update_function(*this, dt);
    }

    template<typename F> void particle_system<F>::update_velocity(const F &dt) {
        #pragma omp for
        for (int i = 0; i < m_size; ++i) {
            m_vx[i] += dt*m_ax[i];
            m_vy[i] += dt*m_ay[i];
            m_vz[i] += dt*m_az[i];
        }
    }
    template<typename F> void particle_system<F>::update_position(const F &dt) {
        #pragma omp for
        for (int i = 0; i < m_size; ++i) {
            m_x[i] += dt*m_dxdt[i];
            m_y[i] += dt*m_dydt[i];
            m_z[i] += dt*m_dzdt[i];
        }
    }
    template<typename F> void particle_system<F>::update_density(const F &dt) {
        // update density using volumentric strain increment
        #pragma omp for
        for (int i = 0; i < m_size; ++i) {
            m_rho[i] /= 1. + dt*(m_strainratexx[i] + m_strainrateyy[i] + m_strainratezz[i]);
        }
    }

    template<typename F> void particle_system<F>::save_to_file(const std::string &prefix, const int &timestep) const {

        // convert timestep number to string (width of 7 chars, for up to 9,999,999,999 timesteps)
        std::string str_timestep = std::to_string(timestep);
        str_timestep = std::string(7-str_timestep.length(), '0') + str_timestep;

        std::string fname {prefix + str_timestep};

        std::ofstream outfile(fname);

        const int i_width = 7, f_width = 12, f_precision=10;

        outfile << std::setw(i_width) << "id" << ' '
                << std::setw(f_width) << "x" << ' '
                << std::setw(f_width) << "y" << ' '
                << std::setw(f_width) << "z" << ' '
                << std::setw(f_width) << "vx" << ' '
                << std::setw(f_width) << "vy" << ' '
                << std::setw(f_width) << "vz" << ' '
                << std::setw(f_width) << "mass" << ' '
                << std::setw(f_width) << "rho" << ' '
                << std::setw(f_width) << "sigmaxx" << ' '
                << std::setw(f_width) << "sigmayy" << ' '
                << std::setw(f_width) << "sigmazz" << ' '
                << std::setw(f_width) << "sigmaxy" << ' '
                << std::setw(f_width) << "sigmaxz" << ' '
                << std::setw(f_width) << "sigmayz" << ' '
                << std::setw(f_width) << "ax" << ' '
                << std::setw(f_width) << "ay" << ' '
                << std::setw(f_width) << "az" << ' '
                << std::setw(f_width) << "dxdt" << ' '
                << std::setw(f_width) << "dydt" << ' '
                << std::setw(f_width) << "dzdt" << ' '
                << std::setw(f_width) << "strainratexx" << ' '
                << std::setw(f_width) << "strainrateyy" << ' '
                << std::setw(f_width) << "strainratezz" << ' '
                << std::setw(f_width) << "strainratexy" << ' '
                << std::setw(f_width) << "strainratexz" << ' '
                << std::setw(f_width) << "strainrateyz" << ' '
                << std::setw(f_width) << "spinratexy" << ' '
                << std::setw(f_width) << "spinratexz" << ' '
                << std::setw(f_width) << "spinrateyz" << ' '
                << '\n';

        for (int i = 0; i < m_size; ++i) {
            outfile << std::setw(i_width) << i << ' ' << std::setprecision(f_precision)
                    << std::setw(f_width) << std::fixed << m_x[i] << ' '
                    << std::setw(f_width) << std::fixed << m_y[i] << ' '
                    << std::setw(f_width) << std::fixed << m_z[i] << ' '
                    << std::setw(f_width) << std::fixed << m_vx[i] << ' '
                    << std::setw(f_width) << std::fixed << m_vy[i] << ' '
                    << std::setw(f_width) << std::fixed << m_vz[i] << ' '
                    << std::setw(f_width) << std::fixed << m_mass[i] << ' '
                    << std::setw(f_width) << std::fixed << m_rho[i] << ' '
                    << std::setw(f_width) << std::fixed << m_sigmaxx[i] << ' '
                    << std::setw(f_width) << std::fixed << m_sigmayy[i] << ' '
                    << std::setw(f_width) << std::fixed << m_sigmazz[i] << ' '
                    << std::setw(f_width) << std::fixed << m_sigmaxy[i] << ' '
                    << std::setw(f_width) << std::fixed << m_sigmaxz[i] << ' '
                    << std::setw(f_width) << std::fixed << m_sigmayz[i] << ' '
                    << std::setw(f_width) << std::fixed << m_ax[i] << ' '
                    << std::setw(f_width) << std::fixed << m_ay[i] << ' '
                    << std::setw(f_width) << std::fixed << m_az[i] << ' '
                    << std::setw(f_width) << std::fixed << m_dxdt[i] << ' '
                    << std::setw(f_width) << std::fixed << m_dydt[i] << ' '
                    << std::setw(f_width) << std::fixed << m_dzdt[i] << ' '
                    << std::setw(f_width) << std::fixed << m_strainratexx[i] << ' '
                    << std::setw(f_width) << std::fixed << m_strainrateyy[i] << ' '
                    << std::setw(f_width) << std::fixed << m_strainratezz[i] << ' '
                    << std::setw(f_width) << std::fixed << m_strainratexy[i] << ' '
                    << std::setw(f_width) << std::fixed << m_strainratexz[i] << ' '
                    << std::setw(f_width) << std::fixed << m_strainrateyz[i] << ' '
                    << std::setw(f_width) << std::fixed << m_spinratexy[i] << ' '
                    << std::setw(f_width) << std::fixed << m_spinratexz[i] << ' '
                    << std::setw(f_width) << std::fixed << m_spinrateyz[i] << ' '
                    << '\n';
        }
    }

    // initalize empty (no size)
    template<typename F>
    particle_system<F>::particle_system(std::string fname, grid<F> &ingrid, kernel_base<F> &knl)
        : background_grid(ingrid)
        , m_capacity {0}
        , m_size {0}
        , m_body_force {0., 0., 0.}
        , m_knl {knl}
        , m_nneighbour_nodes_perp {static_cast<int>(8*std::ceil(knl.radius)*std::ceil(knl.radius)*std::ceil(knl.radius))}
    {
        std::ifstream file(fname);
        std::string line, header;

        // pull out header
        std::getline(file, header);

        while (std::getline(file, line)) {
            std::istringstream iss(line);
            GraMPM::particle<F> p;
            iss >> p.x >> p.y >> p.z >> p.vx >> p.vy >> p.vz >> p.mass >> p.rho >> p.sigmaxx >> p.sigmayy >> 
                p.sigmazz >> p.sigmaxy >> p.sigmaxz >> p.sigmayz >> p.ax >> p.ay >> p.az >> p.dxdt >> p.dydt >>
                p.dzdt >> p.strainratexx >> p.strainrateyy >> p.strainratezz >> p.strainratexy >> p.strainratexz >>
                p.strainrateyz >> p.spinratexy >> p.spinratexz >> p.spinrateyz;
                
            push_back(p);
            
        }
    }
}

#endif