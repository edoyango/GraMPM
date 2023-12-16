#ifndef GRAMPM_particlesystem_ipp
#define GRAMPM_particlesystem_ipp

#include <iostream>
#include <fstream>
#include <iomanip>
#include <immintrin.h>

static void m256d_reduce(__m256d va, double &b) {
    // adds all the elements of va to b
    __m128d vb = _mm_set1_pd(b);
    __m128d va_low = _mm256_castpd256_pd128(va);
    __m128d va_hi = _mm256_extractf128_pd(va, 1);
    va_low = _mm_add_pd(va_low, va_hi);
    va_low = _mm_hadd_pd(va_low, va_low);
    vb = _mm_add_sd(va_low, vb);
    _mm_store_sd(&b, vb);
}

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
        , m_momentumx(size, 0.)
        , m_momentumy(size, 0.)
        , m_momentumz(size, 0.)
        , m_forcex(size, 0.)
        , m_forcey(size, 0.)
        , m_forcez(size, 0.)
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
        , m_gax(ingrid.ncells())
        , m_gay(ingrid.ncells())
        , m_gaz(ingrid.ncells())
        , m_gdxdt(ingrid.ncells())
        , m_gdydt(ingrid.ncells())
        , m_gdzdt(ingrid.ncells())
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
        , m_gax(ingrid.ncells())
        , m_gay(ingrid.ncells())
        , m_gaz(ingrid.ncells())
        , m_gdxdt(ingrid.ncells())
        , m_gdydt(ingrid.ncells())
        , m_gdzdt(ingrid.ncells())
    {
    }

    // helper function to ravel idx
    template<typename F>
    int particle_system<F>::ravel_grid_idx(const int &idxx, const int &idxy, const int &idxz) const {
        return idxx*background_grid.ngridy()*background_grid.ngridz() + idxy*background_grid.ngridz() + idxz;
    }
    
    // helper function to unravel index (return std::array)
    template<typename F>
    std::array<int, 3> particle_system<F>::unravel_grid_idx(const int &idx) const {
        std::array<int, 3> unravelled_idx;
        div_t tmp = std::div(idx, background_grid.ngridy()*background_grid.ngridz());
        unravelled_idx[0] = tmp.quot;
        tmp = std::div(tmp.rem, background_grid.ngridz());
        unravelled_idx[1] = tmp.quot;
        unravelled_idx[2] = tmp.rem;
        return unravelled_idx;
    }

    // helper function to unravel index (modify args)
    template<typename F>
    void particle_system<F>::unravel_grid_idx(const int &idx, int &idxx, int &idxy, int &idxz) const {
        div_t tmp = std::div(idx, background_grid.ngridy()*background_grid.ngridz());
        idxx = tmp.quot;
        tmp = std::div(tmp.rem, background_grid.ngridz());
        idxy = tmp.quot;
        idxz = tmp.rem;
    }

    // grid getter
    template<typename F>
    const grid<F>* particle_system<F>::grid_address() { return &background_grid; }

    // property getters
    template<typename F>const F& particle_system<F>::x(const int &i) const { return m_x[i]; }
    template<typename F> std::vector<F>* particle_system<F>::x() { return &m_x; }
    template<typename F>const F& particle_system<F>::y(const int &i) const { return m_y[i]; }
    template<typename F> std::vector<F>* particle_system<F>::y() { return &m_y; }
    template<typename F>const F& particle_system<F>::z(const int &i) const { return m_z[i]; }
    template<typename F> std::vector<F>* particle_system<F>::z() { return &m_z; }
    template<typename F>const F& particle_system<F>::vx(const int &i) const { return m_vx[i]; }
    template<typename F> std::vector<F>* particle_system<F>::vx() { return &m_vx; }
    template<typename F>const F& particle_system<F>::vy(const int &i) const { return m_vy[i]; }
    template<typename F> std::vector<F>* particle_system<F>::vy() { return &m_vy; }
    template<typename F>const F& particle_system<F>::vz(const int &i) const { return m_vz[i]; }
    template<typename F> std::vector<F>* particle_system<F>::vz() { return &m_vz; }
    template<typename F>const F& particle_system<F>::ax(const int &i) const { return m_ax[i]; }
    template<typename F> std::vector<F>* particle_system<F>::ax() { return &m_ax; }
    template<typename F>const F& particle_system<F>::ay(const int &i) const { return m_ay[i]; }
    template<typename F> std::vector<F>* particle_system<F>::ay() { return &m_ay; }
    template<typename F>const F& particle_system<F>::az(const int &i) const { return m_az[i]; }
    template<typename F> std::vector<F>* particle_system<F>::az() { return &m_az; }
    template<typename F>const F& particle_system<F>::dxdt(const int &i) const { return m_dxdt[i]; }
    template<typename F> std::vector<F>* particle_system<F>::dxdt() { return &m_dxdt; }
    template<typename F>const F& particle_system<F>::dydt(const int &i) const { return m_dydt[i]; }
    template<typename F> std::vector<F>* particle_system<F>::dydt() { return &m_dydt; }
    template<typename F>const F& particle_system<F>::dzdt(const int &i) const { return m_dzdt[i]; }
    template<typename F> std::vector<F>* particle_system<F>::dzdt() { return &m_dzdt; }
    template<typename F>const F& particle_system<F>::mass(const int &i) const { return m_mass[i]; }
    template<typename F> std::vector<F>* particle_system<F>::mass() { return &m_mass; }
    template<typename F>const F& particle_system<F>::rho(const int &i) const { return m_rho[i]; }
    template<typename F> std::vector<F>* particle_system<F>::rho() { return &m_rho; }
    template<typename F>const F& particle_system<F>::sigmaxx(const int &i) const { return m_sigmaxx[i]; }
    template<typename F> std::vector<F>* particle_system<F>::sigmaxx() { return &m_sigmaxx; }
    template<typename F>const F& particle_system<F>::sigmayy(const int &i) const { return m_sigmayy[i]; }
    template<typename F> std::vector<F>* particle_system<F>::sigmayy() { return &m_sigmayy; }
    template<typename F>const F& particle_system<F>::sigmazz(const int &i) const { return m_sigmazz[i]; }
    template<typename F> std::vector<F>* particle_system<F>::sigmazz() { return &m_sigmazz; }
    template<typename F>const F& particle_system<F>::sigmaxy(const int &i) const { return m_sigmaxy[i]; }
    template<typename F> std::vector<F>* particle_system<F>::sigmaxy() { return &m_sigmaxy; }
    template<typename F>const F& particle_system<F>::sigmaxz(const int &i) const { return m_sigmaxz[i]; }
    template<typename F> std::vector<F>* particle_system<F>::sigmaxz() { return &m_sigmaxz; }
    template<typename F>const F& particle_system<F>::sigmayz(const int &i) const { return m_sigmayz[i]; }
    template<typename F> std::vector<F>* particle_system<F>::sigmayz() { return &m_sigmayz; }
    template<typename F>const F& particle_system<F>::strainratexx(const int &i) const { return m_strainratexx[i]; }
    template<typename F> std::vector<F>* particle_system<F>::strainratexx() { return &m_strainratexx; }
    template<typename F>const F& particle_system<F>::strainrateyy(const int &i) const { return m_strainrateyy[i]; }
    template<typename F> std::vector<F>* particle_system<F>::strainrateyy() { return &m_strainrateyy; }
    template<typename F>const F& particle_system<F>::strainratezz(const int &i) const { return m_strainratezz[i]; }
    template<typename F> std::vector<F>* particle_system<F>::strainratezz() { return &m_strainratezz; }
    template<typename F>const F& particle_system<F>::strainratexy(const int &i) const { return m_strainratexy[i]; }
    template<typename F> std::vector<F>* particle_system<F>::strainratexy() { return &m_strainratexy; }
    template<typename F>const F& particle_system<F>::strainratexz(const int &i) const { return m_strainratexz[i]; }
    template<typename F> std::vector<F>* particle_system<F>::strainratexz() { return &m_strainratexz; }
    template<typename F>const F& particle_system<F>::strainrateyz(const int &i) const { return m_strainrateyz[i]; }
    template<typename F> std::vector<F>* particle_system<F>::strainrateyz() { return &m_strainrateyz; }
    template<typename F>const F& particle_system<F>::spinratexy(const int &i) const { return m_spinratexy[i]; }
    template<typename F> std::vector<F>* particle_system<F>::spinratexy() { return &m_spinratexy; }
    template<typename F>const F& particle_system<F>::spinratexz(const int &i) const { return m_spinratexz[i]; }
    template<typename F> std::vector<F>* particle_system<F>::spinratexz() { return &m_spinratexz; }
    template<typename F>const F& particle_system<F>::spinrateyz(const int &i) const { return m_spinrateyz[i]; }
    template<typename F> std::vector<F>* particle_system<F>::spinrateyz() { return &m_spinrateyz; }
    template<typename F>const std::array<F, 3>& particle_system<F>::body_force() const { return m_body_force; }
    template<typename F>const F& particle_system<F>::body_force(const int &i) const { return m_body_force[i]; }
    template<typename F>const F& particle_system<F>::E() const { return m_E; }
    template<typename F>const F& particle_system<F>::v() const { return m_v; }
    template<typename F> void particle_system<F>::DP_params(F& phi, F& psi, F& coh) const { phi = m_phi; psi = m_psi; coh = m_coh; }
    template<typename F> void particle_system<F>::DP_params(F& phi, F& psi, F& coh, F& alpha_phi, F& alpha_psi, F& k_c) const {
         phi = m_phi; psi = m_psi; coh = m_coh; 
         alpha_phi = m_alphaphi; alpha_psi = m_alphapsi; k_c = m_kc;
    }
    template<typename F>const int& particle_system<F>::ravelled_grid_idx(const int &i) const { return m_grid_idx[i]; }
    template<typename F>
    std::array<int, 3> particle_system<F>::grid_idx(const int &i) const { return unravel_grid_idx(ravelled_grid_idx(i)); }

    template<typename F>
    const int& particle_system<F>::p2g_neighbour_node(const int i, const int j) const { 
        assert(j < m_nneighbour_nodes_perp); 
        return m_p2g_neighbour_nodes[i*m_nneighbour_nodes_perp+j];
    }

    template<typename F>
    const F& particle_system<F>::p2g_neighbour_node_dx(const int i, const int j) const {
        assert(j < m_nneighbour_nodes_perp);
        return m_p2g_neighbour_nodes_dx[i*m_nneighbour_nodes_perp+j];
    }

    template<typename F>
    const F& particle_system<F>::p2g_neighbour_node_dy(const int i, const int j) const {
        assert(j < m_nneighbour_nodes_perp);
        return m_p2g_neighbour_nodes_dy[i*m_nneighbour_nodes_perp+j];
    }

    template<typename F>
    const F& particle_system<F>::p2g_neighbour_node_dz(const int i, const int j) const {
        assert(j < m_nneighbour_nodes_perp);
        return m_p2g_neighbour_nodes_dz[i*m_nneighbour_nodes_perp+j];
    }
    
    template<typename F>
    const F& particle_system<F>::p2g_neighbour_node_w(const int i, const int j) const {
        assert(j < m_nneighbour_nodes_perp);
        return m_p2g_neighbour_nodes_w[i*m_nneighbour_nodes_perp+j];
    }
    
    template<typename F>
    const F& particle_system<F>::p2g_neighbour_node_dwdx(const int i, const int j) const {
        assert(j < m_nneighbour_nodes_perp);
        return m_p2g_neighbour_nodes_dwdx[i*m_nneighbour_nodes_perp+j];
    }
    
    template<typename F>
    const F& particle_system<F>::p2g_neighbour_node_dwdy(const int i, const int j) const {
        assert(j < m_nneighbour_nodes_perp);
        return m_p2g_neighbour_nodes_dwdy[i*m_nneighbour_nodes_perp+j];
    }

    template<typename F>
    const F& particle_system<F>::p2g_neighbour_node_dwdz(const int i, const int j) const {
        assert(j < m_nneighbour_nodes_perp);
        return m_p2g_neighbour_nodes_dwdz[i*m_nneighbour_nodes_perp+j];
    }
    template<typename F> const long unsigned int& particle_system<F>::capacity() const { return m_capacity; }
    template<typename F> const long unsigned int& particle_system<F>::size() const { return m_size; }

    
    template<typename F> void particle_system<F>::set_x(const int &i, const F &x) { m_x[i] = x;}
    template<typename F> void particle_system<F>::set_y(const int &i, const F &y) { m_y[i] = y;}
    template<typename F> void particle_system<F>::set_z(const int &i, const F &z) { m_z[i] = z;}
    template<typename F> void particle_system<F>::set_vx(const int &i, const F &vx) { m_vx[i] = vx;}
    template<typename F> void particle_system<F>::set_vy(const int &i, const F &vy) { m_vy[i] = vy;}
    template<typename F> void particle_system<F>::set_vz(const int &i, const F &vz) { m_vz[i] = vz;}
    template<typename F> void particle_system<F>::set_ax(const int &i, const F &ax) { m_ax[i] = ax;}
    template<typename F> void particle_system<F>::set_ay(const int &i, const F &ay) { m_ay[i] = ay;}
    template<typename F> void particle_system<F>::set_az(const int &i, const F &az) { m_az[i] = az;}
    template<typename F> void particle_system<F>::set_dxdt(const int &i, const F &dxdt) { m_dxdt[i] = dxdt;}
    template<typename F> void particle_system<F>::set_dydt(const int &i, const F &dydt) { m_dydt[i] = dydt;}
    template<typename F> void particle_system<F>::set_dzdt(const int &i, const F &dzdt) { m_dzdt[i] = dzdt;}
    template<typename F> void particle_system<F>::set_mass(const int &i, const F &m) { m_mass[i] = m; }
    template<typename F> void particle_system<F>::set_rho(const int &i, const F &rho) { m_rho[i] = rho; }
    template<typename F> void particle_system<F>::set_sigmaxx(const int &i, const F &sigmaxx) { m_sigmaxx[i] = sigmaxx; }
    template<typename F> void particle_system<F>::set_sigmayy(const int &i, const F &sigmayy) { m_sigmayy[i] = sigmayy; }
    template<typename F> void particle_system<F>::set_sigmazz(const int &i, const F &sigmazz) { m_sigmazz[i] = sigmazz; }
    template<typename F> void particle_system<F>::set_sigmaxy(const int &i, const F &sigmaxy) { m_sigmaxy[i] = sigmaxy; }
    template<typename F> void particle_system<F>::set_sigmaxz(const int &i, const F &sigmaxz) { m_sigmaxz[i] = sigmaxz; }
    template<typename F> void particle_system<F>::set_sigmayz(const int &i, const F &sigmayz) { m_sigmayz[i] = sigmayz; }
    template<typename F> void particle_system<F>::set_strainratexx(const int &i, const F &strainratexx) { m_strainratexx[i] = strainratexx; }
    template<typename F> void particle_system<F>::set_strainrateyy(const int &i, const F &strainrateyy) { m_strainrateyy[i] = strainrateyy; }
    template<typename F> void particle_system<F>::set_strainratezz(const int &i, const F &strainratezz) { m_strainratezz[i] = strainratezz; }
    template<typename F> void particle_system<F>::set_strainratexy(const int &i, const F &strainratexy) { m_strainratexy[i] = strainratexy; }
    template<typename F> void particle_system<F>::set_strainratexz(const int &i, const F &strainratexz) { m_strainratexz[i] = strainratexz; }
    template<typename F> void particle_system<F>::set_strainrateyz(const int &i, const F &strainrateyz) { m_strainrateyz[i] = strainrateyz; }
    template<typename F> void particle_system<F>::set_spinratexy(const int &i, const F &spinratexy) { m_spinratexy[i] = spinratexy; }
    template<typename F> void particle_system<F>::set_spinratexz(const int &i, const F &spinratexz) { m_spinratexz[i] = spinratexz; }
    template<typename F> void particle_system<F>::set_spinrateyz(const int &i, const F &spinrateyz) { m_spinrateyz[i] = spinrateyz; }
    template<typename F> void particle_system<F>::set_body_force(const std::array<F, 3> &bf) { m_body_force = bf; }
    template<typename F> void particle_system<F>::set_body_force(const F &bfx, const F &bfy, const F &bfz) { 
        m_body_force[0] = bfx;
        m_body_force[1] = bfy;
        m_body_force[2] = bfz;
    }
    template<typename F> void particle_system<F>::set_grid_index(const int &i, const int &idx) { m_grid_idx[i] = idx; }
    template<typename F> void particle_system<F>::set_stress_update_function(std::function<void(particle_system<F>&, const F&)> f) { m_stress_update_function = f; }
    template<typename F> void particle_system<F>::set_E(const F &E) { m_E = E; }
    template<typename F> void particle_system<F>::set_v(const F &v) { m_v = v; }
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
    template<typename F> void particle_system<F>::incrementNParticles() {m_size++;}

    // vector-like api: at. Returns particle class
    template<typename F>
    particle<F> particle_system<F>::at(const int &i) { 
        particle<F> p(x(i), y(i), z(i), vx(i), vy(i), vz(i), mass(i), rho(i), sigmaxx(i), sigmayy(i), sigmazz(i), 
            sigmaxy(i), sigmaxz(i), sigmayz(i), ax(i), ay(i), az(i), dxdt(i), dydt(i), dzdt(i), strainratexx(i), 
            strainrateyy(i), strainratezz(i), strainratexy(i), strainratexz(i), strainrateyz(i), spinratexy(i), 
            spinratexz(i), spinrateyz(i));
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
        m_momentumx.push_back(0.);
        m_momentumy.push_back(0.);
        m_momentumz.push_back(0.);
        m_forcex.push_back(0.);
        m_forcey.push_back(0.);
        m_forcez.push_back(0.);
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
        m_momentumx.reserve(n);
        m_momentumy.reserve(n);
        m_momentumz.reserve(n);
        m_forcex.reserve(n);
        m_forcey.reserve(n);
        m_forcez.reserve(n);
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
        m_momentumx.clear();
        m_momentumy.clear();
        m_momentumz.clear();
        m_forcex.clear();
        m_forcey.clear();
        m_forcez.clear();
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
        m_p2g_neighbour_nodes.clear();
        m_p2g_neighbour_nodes_dx.clear();
        m_p2g_neighbour_nodes_dy.clear();
        m_p2g_neighbour_nodes_dz.clear();
        m_p2g_neighbour_nodes_w.clear();
        m_p2g_neighbour_nodes_dwdx.clear();
        m_p2g_neighbour_nodes_dwdy.clear();
        m_p2g_neighbour_nodes_dwdz.clear();
        m_size = 0;
    }

    // vector-like api: empty. Checks whether all member vectors are empty. Also checks the m_size member
    template<typename F>
    bool particle_system<F>::empty() {
        return m_x.empty() && m_y.empty() && m_z.empty() && m_vx.empty() && m_vy.empty() && m_vz.empty() && 
            m_ax.empty() && m_ay.empty() && m_az.empty() && m_dxdt.empty() && m_dydt.empty() && m_dzdt.empty() &&
            m_momentumx.empty() && m_momentumy.empty() && m_momentumz.empty() && m_forcex.empty() && m_forcey.empty() &&
            m_forcez.empty() &&
            m_mass.empty() && m_rho.empty() && m_sigmaxx.empty() && m_sigmayy.empty() && m_sigmazz.empty() && 
            m_sigmaxy.empty() && m_sigmaxz.empty() && m_sigmayz.empty() && m_strainratexx.empty() && 
            m_strainrateyy.empty() && m_strainratezz.empty() && m_strainratexy.empty() && m_strainratexz.empty() && 
            m_strainrateyz.empty() && m_spinratexy.empty() && m_spinratexz.empty() && 
            m_spinrateyz.empty() &&m_grid_idx.empty() && m_p2g_neighbour_nodes.empty() && 
            m_p2g_neighbour_nodes_dx.empty() && m_p2g_neighbour_nodes_dy.empty() && m_p2g_neighbour_nodes_dz.empty() && 
            m_p2g_neighbour_nodes_w.empty() && m_p2g_neighbour_nodes_dwdx.empty() &&
            m_p2g_neighbour_nodes_dwdy.empty() && m_p2g_neighbour_nodes_dwdz.empty() && m_size==0;
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
        m_momentumx.reize(n, 0.);
        m_momentumy.reize(n, 0.);
        m_momentumz.reize(n, 0.);
        m_forcex.resize(n, 0.);
        m_forcey.resize(n, 0.);
        m_forcez.resize(n, 0.);
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
        m_momentumx.reize(n, 0.);
        m_momentumy.reize(n, 0.);
        m_momentumz.reize(n, 0.);
        m_forcex.resize(n, 0.);
        m_forcey.resize(n, 0.);
        m_forcez.resize(n, 0.);
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
    void particle_system<F>::resize_grid_scratchspace(const int n) {
        m_gax.resize(n);
        m_gay.resize(n);
        m_gaz.resize(n);
        m_gdxdt.resize(n);
        m_gdydt.resize(n);
        m_gdzdt.resize(n);
    }

    template<typename F>
    void particle_system<F>::update_particle_to_cell_map(const int &start, const int &end) {
        for (int i = start; i < end; ++i) {
            set_grid_index(i,
                ravel_grid_idx(
                    background_grid.calc_idxx(x(i)),
                    background_grid.calc_idxy(y(i)),
                    background_grid.calc_idxz(z(i))
                )
            );
        }
    }

    template<typename F>
    void particle_system<F>::update_particle_to_cell_map() {
        for (int i = 0; i < m_size; ++i) {
            set_grid_index( i,
                ravel_grid_idx(
                    background_grid.calc_idxx(x(i)),
                    background_grid.calc_idxy(y(i)),
                    background_grid.calc_idxz(z(i))
                )
            );
        }
    }
    
    // NTS this could be faster
    template<typename F>
    void particle_system<F>::map_particles_to_grid() {

        // size output arrays
        m_p2g_neighbour_nodes.resize(m_nneighbour_nodes_perp*m_size);
        m_p2g_neighbour_nodes_dx.resize(m_nneighbour_nodes_perp*m_size);
        m_p2g_neighbour_nodes_dy.resize(m_nneighbour_nodes_perp*m_size);
        m_p2g_neighbour_nodes_dz.resize(m_nneighbour_nodes_perp*m_size);
        m_p2g_neighbour_nodes_w.resize(m_nneighbour_nodes_perp*m_size);
        m_p2g_neighbour_nodes_dwdx.resize(m_nneighbour_nodes_perp*m_size);
        m_p2g_neighbour_nodes_dwdy.resize(m_nneighbour_nodes_perp*m_size);
        m_p2g_neighbour_nodes_dwdz.resize(m_nneighbour_nodes_perp*m_size);

        // update neighbour indices
        for (int i = 0; i < m_size; ++i) {
            const int idx = ravelled_grid_idx(i);
            int n = 0;
            for (int di=1-m_knl.radius; di <= m_knl.radius; ++di) {
                for (int dj = 1-m_knl.radius; dj <= m_knl.radius; ++dj) {
                    for (int dk = 1-m_knl.radius; dk <= m_knl.radius; ++dk) {
                        m_p2g_neighbour_nodes[i*m_nneighbour_nodes_perp+n] = idx + ravel_grid_idx(di, dj, dk);
                        n++;
                    }
                }
            }
        }

        // update kernel and kernel gradient values
        for (int i = 0; i < m_size; ++i) {
            for (int j = 0; j < m_nneighbour_nodes_perp; ++j) {
                std::array<int, 3> idx;
                idx = unravel_grid_idx(m_p2g_neighbour_nodes[i*m_nneighbour_nodes_perp+j]);
                m_p2g_neighbour_nodes_dx[i*m_nneighbour_nodes_perp+j] = x(i) - (idx[0]*background_grid.cell_size() + background_grid.mingridx());
                m_p2g_neighbour_nodes_dy[i*m_nneighbour_nodes_perp+j] = y(i) - (idx[1]*background_grid.cell_size() + background_grid.mingridy());
                m_p2g_neighbour_nodes_dz[i*m_nneighbour_nodes_perp+j] = z(i) - (idx[2]*background_grid.cell_size() + background_grid.mingridz());
                m_p2g_neighbour_nodes_w[i*m_nneighbour_nodes_perp+j] = m_knl.w(
                    m_p2g_neighbour_nodes_dx[i*m_nneighbour_nodes_perp+j],
                    m_p2g_neighbour_nodes_dy[i*m_nneighbour_nodes_perp+j],
                    m_p2g_neighbour_nodes_dz[i*m_nneighbour_nodes_perp+j]
                );
                m_knl.dwdx(
                    m_p2g_neighbour_nodes_dx[i*m_nneighbour_nodes_perp+j],
                    m_p2g_neighbour_nodes_dy[i*m_nneighbour_nodes_perp+j],
                    m_p2g_neighbour_nodes_dz[i*m_nneighbour_nodes_perp+j],
                    m_p2g_neighbour_nodes_dwdx[i*m_nneighbour_nodes_perp+j],
                    m_p2g_neighbour_nodes_dwdy[i*m_nneighbour_nodes_perp+j],
                    m_p2g_neighbour_nodes_dwdz[i*m_nneighbour_nodes_perp+j]
                );
            }
        }
    }

    template<typename F>
    void particle_system<F>::map2grid(const std::vector<F> &p_property, std::vector<F> *g_property) {

        // zero the grid
        std::fill(g_property->begin(), g_property->end(), 0.);

        for (int i = 0; i < m_size; ++i)
            for (int j = 0; j < m_nneighbour_nodes_perp; ++j) {
                const int node_idx = p2g_neighbour_node(i, j);
                (*g_property)[node_idx] += p2g_neighbour_node_w(i, j)*p_property[i];
            }
    }

    template<typename F>
    void particle_system<F>::map2particles(const std::vector<F> &g_property, std::vector<F> *p_property) {

        // zero the particles' array
        std::fill(p_property->begin(), p_property->end(), 0.);

        for (int i = 0; i < m_size; ++i) {
#ifdef __AVX2__
            __m256d p_property_v = _mm256_setzero_pd();
            for (int j = 0; j < m_nneighbour_nodes_perp; j+=4) {
                const int node_idx = p2g_neighbour_node(i, j);
                __m256d g_property_v = _mm256_loadu_pd(&g_property[node_idx]);
                __m256d p2g_w_v = _mm256_loadu_pd(&m_p2g_neighbour_nodes_w[i*m_nneighbour_nodes_perp+j]);
                p_property_v = _mm256_add_pd(_mm256_mul_pd(p2g_w_v, g_property_v), p_property_v);
            }
            m256d_reduce(p_property_v, (*p_property)[i]);
#else
            for (int j = 0; j < m_nneighbour_nodes_perp; ++j) {
                const int node_idx = p2g_neighbour_node(i, j);
                (*p_property)[i] += p2g_neighbour_node_w(i, j)*g_property[node_idx];
            }
#endif
        }
    }

    template<typename F> void particle_system<F>::map_mass_to_grid() { map2grid(m_mass, background_grid.mass()); }

    template<typename F> void particle_system<F>::map_momentum_to_grid() { 
        for (int i = 0; i < m_size; ++i) m_momentumx[i] = mass(i)*vx(i);
        map2grid(m_momentumx, background_grid.momentumx()); 
        for (int i = 0; i < m_size; ++i) m_momentumy[i] = mass(i)*vy(i);
        map2grid(m_momentumy, background_grid.momentumy()); 
        for (int i = 0; i < m_size; ++i) m_momentumz[i] = mass(i)*vz(i);
        map2grid(m_momentumz, background_grid.momentumz()); 
    }

    template<typename F> void particle_system<F>::map_force_to_grid() {

        // initialize grid force with body force
        for (int i = 0; i < m_size; ++i) m_forcex[i] = mass(i)*body_force(0);
        map2grid(m_forcex, background_grid.forcex());
        for (int i = 0; i < m_size; ++i) m_forcey[i] = mass(i)*body_force(1);
        map2grid(m_forcey, background_grid.forcey());
        for (int i = 0; i < m_size; ++i) m_forcez[i] = mass(i)*body_force(2);
        map2grid(m_forcez, background_grid.forcez());

        std::vector<F> *tmp_grid_forcex = background_grid.forcex(), 
            *tmp_grid_forcey = background_grid.forcey(),
            *tmp_grid_forcez = background_grid.forcez();

        // div sigma
        for (int i = 0; i < m_size; ++i)
            for (int j = 0; j < m_nneighbour_nodes_perp; ++j) {
                const int node_idx = p2g_neighbour_node(i, j);
                (*tmp_grid_forcex)[node_idx] -= m_mass[i]/m_rho[i]*(
                    sigmaxx(i)*p2g_neighbour_node_dwdx(i, j) +
                    sigmaxy(i)*p2g_neighbour_node_dwdy(i, j) +
                    sigmaxz(i)*p2g_neighbour_node_dwdz(i, j)
                );
                (*tmp_grid_forcey)[node_idx] -= m_mass[i]/m_rho[i]*(
                    sigmaxy(i)*p2g_neighbour_node_dwdx(i, j) +
                    sigmayy(i)*p2g_neighbour_node_dwdy(i, j) +
                    sigmayz(i)*p2g_neighbour_node_dwdz(i, j)
                );
                (*tmp_grid_forcez)[node_idx] -= m_mass[i]/m_rho[i]*(
                    sigmaxz(i)*p2g_neighbour_node_dwdx(i, j) +
                    sigmayz(i)*p2g_neighbour_node_dwdy(i, j) +
                    sigmazz(i)*p2g_neighbour_node_dwdz(i, j)
                );
            }

    }

    template<typename F> void particle_system<F>::map_acceleration_to_particles() {
        for (int i = 0; i < background_grid.ncells(); ++i) {
            m_gax[i] = background_grid.forcex(i)/background_grid.mass(i);
        }
        map2particles(m_gax, ax());
        for (int i = 0; i < background_grid.ncells(); ++i) {
            m_gay[i] = background_grid.forcey(i)/background_grid.mass(i);
        }
        map2particles(m_gay, ay());
        for (int i = 0; i < background_grid.ncells(); ++i) {
            m_gaz[i] = background_grid.forcez(i)/background_grid.mass(i);
        }
        map2particles(m_gaz, az());
        for (int i = 0; i < background_grid.ncells(); ++i) {
            m_gdxdt[i] = background_grid.momentumx(i)/background_grid.mass(i);
        }
        map2particles(m_gdxdt, dxdt());
        for (int i = 0; i < background_grid.ncells(); ++i) {
            m_gdydt[i] = background_grid.momentumy(i)/background_grid.mass(i);
        }
        map2particles(m_gdydt, dydt());
        for (int i = 0; i < background_grid.ncells(); ++i) {
            m_gdzdt[i] = background_grid.momentumz(i)/background_grid.mass(i);
        }
        map2particles(m_gdzdt, dzdt());
    }

    template<typename F> void particle_system<F>::map_strainrate_to_particles() {

        std::fill(m_strainratexx.begin(), m_strainratexx.end(), 0.);
        std::fill(m_strainrateyy.begin(), m_strainrateyy.end(), 0.);
        std::fill(m_strainratezz.begin(), m_strainratezz.end(), 0.);
        std::fill(m_strainratexy.begin(), m_strainratexy.end(), 0.);
        std::fill(m_strainratexz.begin(), m_strainratexz.end(), 0.);
        std::fill(m_strainrateyz.begin(), m_strainrateyz.end(), 0.);
        std::fill(m_spinratexy.begin(), m_spinratexy.end(), 0.);
        std::fill(m_spinratexz.begin(), m_spinratexz.end(), 0.);
        std::fill(m_spinrateyz.begin(), m_spinrateyz.end(), 0.);

        // calculating velocities at grid
        std::vector<F> tmp_gvx(background_grid.ncells()), tmp_gvy(background_grid.ncells()), tmp_gvz(background_grid.ncells());
        for (int i = 0; i < background_grid.ncells(); ++i) {
            tmp_gvx[i] = background_grid.momentumx(i)/background_grid.mass(i);
            tmp_gvy[i] = background_grid.momentumy(i)/background_grid.mass(i);
            tmp_gvz[i] = background_grid.momentumz(i)/background_grid.mass(i);
        }

        for (int i = 0; i < m_size; ++i) {
#ifdef __AVX2__
            __m256d strainratexx_cum_v = _mm256_setzero_pd();
            __m256d strainrateyy_cum_v = _mm256_setzero_pd();
            __m256d strainratezz_cum_v = _mm256_setzero_pd();
            __m256d strainratexy_cum_v = _mm256_setzero_pd();
            __m256d strainratexz_cum_v = _mm256_setzero_pd();
            __m256d strainrateyz_cum_v = _mm256_setzero_pd();
            __m256d spinratexy_cum_v = _mm256_setzero_pd();
            __m256d spinratexz_cum_v = _mm256_setzero_pd();
            __m256d spinrateyz_cum_v = _mm256_setzero_pd();
            for (int j = 0; j < m_nneighbour_nodes_perp; j += 4) {
                // __m128i node_idx_v = _mm_loadu_si128((__m128i*)&m_p2g_neighbour_nodes[i*m_nneighbour_nodes_perp+j]);
                // __m256d tmp_gvx_v = _mm256_i32gather_pd(tmp_gvx.data(), node_idx_v, 8);
                // __m256d tmp_gvy_v = _mm256_i32gather_pd(tmp_gvy.data(), node_idx_v, 8);
                // __m256d tmp_gvz_v = _mm256_i32gather_pd(tmp_gvz.data(), node_idx_v, 8);
                const int node_idx = m_p2g_neighbour_nodes[i*m_nneighbour_nodes_perp+j];
                __m256d tmp_gvx_v = _mm256_loadu_pd(&tmp_gvx[node_idx]);
                __m256d tmp_gvy_v = _mm256_loadu_pd(&tmp_gvy[node_idx]);
                __m256d tmp_gvz_v = _mm256_loadu_pd(&tmp_gvz[node_idx]);
                __m256d p2g_dwdx_v = _mm256_loadu_pd(&m_p2g_neighbour_nodes_dwdx[i*m_nneighbour_nodes_perp+j]);
                __m256d p2g_dwdy_v = _mm256_loadu_pd(&m_p2g_neighbour_nodes_dwdy[i*m_nneighbour_nodes_perp+j]);
                __m256d p2g_dwdz_v = _mm256_loadu_pd(&m_p2g_neighbour_nodes_dwdz[i*m_nneighbour_nodes_perp+j]);
                strainratexx_cum_v = _mm256_add_pd(_mm256_mul_pd(tmp_gvx_v, p2g_dwdx_v), strainratexx_cum_v);
                strainrateyy_cum_v = _mm256_add_pd(_mm256_mul_pd(tmp_gvy_v, p2g_dwdy_v), strainrateyy_cum_v);
                strainratezz_cum_v = _mm256_add_pd(_mm256_mul_pd(tmp_gvz_v, p2g_dwdz_v), strainratezz_cum_v);
                strainratexy_cum_v = _mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(0.5),_mm256_add_pd(_mm256_mul_pd(tmp_gvy_v, p2g_dwdx_v), _mm256_mul_pd(tmp_gvx_v, p2g_dwdy_v))), strainratexy_cum_v);
                strainratexz_cum_v = _mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(0.5),_mm256_add_pd(_mm256_mul_pd(tmp_gvz_v, p2g_dwdx_v), _mm256_mul_pd(tmp_gvx_v, p2g_dwdz_v))), strainratexz_cum_v);
                strainrateyz_cum_v = _mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(0.5),_mm256_add_pd(_mm256_mul_pd(tmp_gvz_v, p2g_dwdy_v), _mm256_mul_pd(tmp_gvy_v, p2g_dwdz_v))), strainrateyz_cum_v);
                spinratexy_cum_v = _mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(0.5),_mm256_sub_pd(_mm256_mul_pd(tmp_gvx_v, p2g_dwdy_v), _mm256_mul_pd(tmp_gvy_v, p2g_dwdx_v))), spinratexy_cum_v);
                spinratexz_cum_v = _mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(0.5),_mm256_sub_pd(_mm256_mul_pd(tmp_gvx_v, p2g_dwdz_v), _mm256_mul_pd(tmp_gvz_v, p2g_dwdx_v))), spinratexz_cum_v);
                spinrateyz_cum_v = _mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(0.5),_mm256_sub_pd(_mm256_mul_pd(tmp_gvy_v, p2g_dwdz_v), _mm256_mul_pd(tmp_gvz_v, p2g_dwdy_v))), spinrateyz_cum_v);
            }
            m256d_reduce(strainratexx_cum_v, m_strainratexx[i]);
            m256d_reduce(strainrateyy_cum_v, m_strainrateyy[i]);
            m256d_reduce(strainratezz_cum_v, m_strainratezz[i]);
            m256d_reduce(strainratexy_cum_v, m_strainratexy[i]);
            m256d_reduce(strainratexz_cum_v, m_strainratexz[i]);
            m256d_reduce(strainrateyz_cum_v, m_strainrateyz[i]);
            m256d_reduce(spinratexy_cum_v, m_spinratexy[i]);
            m256d_reduce(spinratexz_cum_v, m_spinratexz[i]);
            m256d_reduce(spinrateyz_cum_v, m_spinrateyz[i]);
        }
#else
            for (int j = 0; j < m_nneighbour_nodes_perp; ++j) {
                const int node_idx = p2g_neighbour_node(i, j);
                m_strainratexx[i] += p2g_neighbour_node_dwdx(i, j)*tmp_gvx[node_idx];
                m_strainrateyy[i] += p2g_neighbour_node_dwdy(i, j)*tmp_gvy[node_idx];
                m_strainratezz[i] += p2g_neighbour_node_dwdz(i, j)*tmp_gvz[node_idx];
                m_strainratexy[i] += 0.5*(p2g_neighbour_node_dwdx(i, j)*tmp_gvy[node_idx] + 
                    p2g_neighbour_node_dwdy(i, j)*tmp_gvx[node_idx]);
                m_strainratexz[i] += 0.5*(p2g_neighbour_node_dwdx(i, j)*tmp_gvz[node_idx] + 
                    p2g_neighbour_node_dwdz(i, j)*tmp_gvx[node_idx]);
                m_strainrateyz[i] += 0.5*(p2g_neighbour_node_dwdy(i, j)*tmp_gvz[node_idx] + 
                    p2g_neighbour_node_dwdz(i, j)*tmp_gvy[node_idx]);
                m_spinratexy[i] += 0.5*(p2g_neighbour_node_dwdy(i, j)*tmp_gvx[node_idx] -
                    p2g_neighbour_node_dwdx(i, j)*tmp_gvy[node_idx]);
                m_spinratexz[i] += 0.5*(p2g_neighbour_node_dwdz(i, j)*tmp_gvx[node_idx] - 
                    p2g_neighbour_node_dwdx(i, j)*tmp_gvz[node_idx]);
                m_spinrateyz[i] += 0.5*(p2g_neighbour_node_dwdz(i, j)*tmp_gvy[node_idx] - 
                    p2g_neighbour_node_dwdy(i, j)*tmp_gvz[node_idx]);
            }
        }
#endif
    }

    template<typename F> void particle_system<F>::update_stress(const F &dt) {
        m_stress_update_function(*this, dt);
    }

    template<typename F> void particle_system<F>::update_velocity(const F &dt) {
        for (int i = 0; i < m_size; ++i) {
            m_vx[i] += dt*m_ax[i];
            m_vy[i] += dt*m_ay[i];
            m_vz[i] += dt*m_az[i];
        }
    }
    template<typename F> void particle_system<F>::update_position(const F &dt) {
        for (int i = 0; i < m_size; ++i) {
            m_x[i] += dt*m_dxdt[i];
            m_y[i] += dt*m_dydt[i];
            m_z[i] += dt*m_dzdt[i];
        }
    }
    template<typename F> void particle_system<F>::update_density(const F &dt) {
        // update density using volumentric strain increment
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
                    << std::setw(f_width) << std::fixed << x(i) << ' '
                    << std::setw(f_width) << std::fixed << y(i) << ' '
                    << std::setw(f_width) << std::fixed << z(i) << ' '
                    << std::setw(f_width) << std::fixed << vx(i) << ' '
                    << std::setw(f_width) << std::fixed << vy(i) << ' '
                    << std::setw(f_width) << std::fixed << vz(i) << ' '
                    << std::setw(f_width) << std::fixed << mass(i) << ' '
                    << std::setw(f_width) << std::fixed << rho(i) << ' '
                    << std::setw(f_width) << std::fixed << sigmaxx(i) << ' '
                    << std::setw(f_width) << std::fixed << sigmayy(i) << ' '
                    << std::setw(f_width) << std::fixed << sigmazz(i) << ' '
                    << std::setw(f_width) << std::fixed << sigmaxy(i) << ' '
                    << std::setw(f_width) << std::fixed << sigmaxz(i) << ' '
                    << std::setw(f_width) << std::fixed << sigmayz(i) << ' '
                    << std::setw(f_width) << std::fixed << ax(i) << ' '
                    << std::setw(f_width) << std::fixed << ay(i) << ' '
                    << std::setw(f_width) << std::fixed << az(i) << ' '
                    << std::setw(f_width) << std::fixed << dxdt(i) << ' '
                    << std::setw(f_width) << std::fixed << dydt(i) << ' '
                    << std::setw(f_width) << std::fixed << dzdt(i) << ' '
                    << std::setw(f_width) << std::fixed << strainratexx(i) << ' '
                    << std::setw(f_width) << std::fixed << strainrateyy(i) << ' '
                    << std::setw(f_width) << std::fixed << strainratezz(i) << ' '
                    << std::setw(f_width) << std::fixed << strainratexy(i) << ' '
                    << std::setw(f_width) << std::fixed << strainratexz(i) << ' '
                    << std::setw(f_width) << std::fixed << strainrateyz(i) << ' '
                    << std::setw(f_width) << std::fixed << spinratexy(i) << ' '
                    << std::setw(f_width) << std::fixed << spinratexz(i) << ' '
                    << std::setw(f_width) << std::fixed << spinrateyz(i) << ' '
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