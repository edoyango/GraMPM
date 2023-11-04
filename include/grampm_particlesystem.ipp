#ifndef GRAMPM_particlesystem_ipp
#define GRAMPM_particlesystem_ipp

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
        , m_dvx(size, 0.)
        , m_dvy(size, 0.)
        , m_dvz(size, 0.)
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
    template<typename F>const F& particle_system<F>::dvx(const int &i) const { return m_dvx[i]; }
    template<typename F> std::vector<F>* particle_system<F>::dvx() { return &m_dvx; }
    template<typename F>const F& particle_system<F>::dvy(const int &i) const { return m_dvy[i]; }
    template<typename F> std::vector<F>* particle_system<F>::dvy() { return &m_dvy; }
    template<typename F>const F& particle_system<F>::dvz(const int &i) const { return m_dvz[i]; }
    template<typename F> std::vector<F>* particle_system<F>::dvz() { return &m_dvz; }
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
    template<typename F> void particle_system<F>::set_dvx(const int &i, const F &dvx) { m_dvx[i] = dvx;}
    template<typename F> void particle_system<F>::set_dvy(const int &i, const F &dvy) { m_dvy[i] = dvy;}
    template<typename F> void particle_system<F>::set_dvz(const int &i, const F &dvz) { m_dvz[i] = dvz;}
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
    template<typename F> void particle_system<F>::set_E(const F &E) { m_E = E; }
    template<typename F> void particle_system<F>::set_v(const F &v) { m_v = v; }
    template<typename F> void particle_system<F>::incrementNParticles() {m_size++;}

    // vector-like api: at. Returns particle class
    template<typename F>
    particle<F> particle_system<F>::at(const int &i) { 
        particle<F> p(x(i), y(i), z(i), vx(i), vy(i), vz(i), mass(i), rho(i), sigmaxx(i), sigmayy(i), sigmazz(i), 
            sigmaxy(i), sigmaxz(i), sigmayz(i), ax(i), ay(i), az(i), dvx(i), dvy(i), dvz(i), strainratexx(i), 
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
        m_dvx.push_back(p.dvx);
        m_dvy.push_back(p.dvy);
        m_dvz.push_back(p.dvz);
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
        m_dvx.reserve(n);
        m_dvy.reserve(n);
        m_dvz.reserve(n);
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
        m_dvx.clear();
        m_dvy.clear();
        m_dvz.clear();
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
            m_ax.empty() && m_ay.empty() && m_az.empty() && m_dvx.empty() && m_dvy.empty() && m_dvz.empty() &&
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
        m_dvx.resize(n, 0.);
        m_dvy.resize(n, 0.);
        m_dvz.resize(n, 0.);
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
        m_dvx.resize(n, p.ax);
        m_dvy.resize(n, p.ay);
        m_dvz.resize(n, p.az);
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

        for (int i = 0; i < m_size; ++i)
            for (int j = 0; j < m_nneighbour_nodes_perp; ++j) {
                const int node_idx = p2g_neighbour_node(i, j);
                (*p_property)[i] += p2g_neighbour_node_w(i, j)*g_property[node_idx];
            }
    }

    template<typename F> void particle_system<F>::map_mass_to_grid() { map2grid(m_mass, background_grid.mass()); }

    template<typename F> void particle_system<F>::map_momentum_to_grid() { 
        std::vector<F> tmp_momentum(m_size);
        for (int i = 0; i < m_size; ++i) tmp_momentum[i] = mass(i)*vx(i);
        map2grid(tmp_momentum, background_grid.momentumx()); 
        for (int i = 0; i < m_size; ++i) tmp_momentum[i] = mass(i)*vy(i);
        map2grid(tmp_momentum, background_grid.momentumy()); 
        for (int i = 0; i < m_size; ++i) tmp_momentum[i] = mass(i)*vz(i);
        map2grid(tmp_momentum, background_grid.momentumz()); 
    }

    template<typename F> void particle_system<F>::map_force_to_grid() {

        // initialize grid force with body force
        std::vector<F> tmp_force(m_size);
        for (int i = 0; i < m_size; ++i) tmp_force[i] = mass(i)*body_force(0);
        map2grid(tmp_force, background_grid.forcex());
        for (int i = 0; i < m_size; ++i) tmp_force[i] = mass(i)*body_force(1);
        map2grid(tmp_force, background_grid.forcey());
        for (int i = 0; i < m_size; ++i) tmp_force[i] = mass(i)*body_force(2);
        map2grid(tmp_force, background_grid.forcez());

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
        std::vector<F> tmp(background_grid.ncells());
        for (int i = 0; i < background_grid.ncells(); ++i) {
            tmp[i] = background_grid.forcex(i)/background_grid.mass(i);
        }
        map2particles(tmp, ax());
        for (int i = 0; i < background_grid.ncells(); ++i) {
            tmp[i] = background_grid.forcey(i)/background_grid.mass(i);
        }
        map2particles(tmp, ay());
        for (int i = 0; i < background_grid.ncells(); ++i) {
            tmp[i] = background_grid.forcez(i)/background_grid.mass(i);
        }
        map2particles(tmp, az());
        for (int i = 0; i < background_grid.ncells(); ++i) {
            tmp[i] = background_grid.momentumx(i)/background_grid.mass(i);
        }
        map2particles(tmp, dvx());
        for (int i = 0; i < background_grid.ncells(); ++i) {
            tmp[i] = background_grid.momentumy(i)/background_grid.mass(i);
        }
        map2particles(tmp, dvy());
        for (int i = 0; i < background_grid.ncells(); ++i) {
            tmp[i] = background_grid.momentumz(i)/background_grid.mass(i);
        }
        map2particles(tmp, dvz());
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

        for (int i = 0; i < m_size; ++i)
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

    template<typename F> void particle_system<F>::stress_update(const F &dt) {

        const F D0 {m_E/((1.+m_v)*(1.-2.*m_v))};

        std::vector<F> dsigmaxx(m_size), dsigmayy(m_size), dsigmazz(m_size), dsigmaxy(m_size), dsigmaxz(m_size),
            dsigmayz(m_size);

        // DE*dstrain
        for (int i = 0; i < m_size; ++i) {
            dsigmaxx[i] = dt*D0*((1.-m_v)*m_strainratexx[i] + m_v*m_strainrateyy[i] + m_v*m_strainratezz[i]);
            dsigmayy[i] = dt*D0*(m_v*m_strainratexx[i] + (1.-m_v)*m_strainrateyy[i] + m_v*m_strainratezz[i]);
            dsigmazz[i] = dt*D0*(m_v*m_strainratexx[i] + m_v*m_strainrateyy[i] + (1.-m_v)*m_strainratezz[i]);
            dsigmaxy[i] = dt*D0*m_strainratexy[i]*(1.-2.*m_v);
            dsigmaxz[i] = dt*D0*m_strainratexz[i]*(1.-2.*m_v);
            dsigmayz[i] = dt*D0*m_strainrateyz[i]*(1.-2.*m_v);
        }
        
        // jaumann stress rate
        for (int i = 0; i < m_size; ++i) {
            dsigmaxx[i] -= dt*2.*(m_spinratexy[i]*m_sigmaxy[i] + m_spinratexz[i]*m_sigmaxz[i]);
            dsigmayy[i] -= dt*2.*(-m_spinratexy[i]*m_sigmaxy[i] + m_spinrateyz[i]*m_sigmayz[i]);
            dsigmazz[i] += dt*2.*(m_spinratexz[i]*m_sigmaxz[i] + m_spinrateyz[i]*m_sigmayz[i]);
            dsigmaxy[i] += dt*(
                m_sigmaxx[i]*m_spinratexy[i] - m_sigmaxz[i]*m_spinrateyz[i] -
                m_spinratexy[i]*m_sigmayy[i] - m_spinratexz[i]*m_sigmayz[i]
            );
            dsigmaxz[i] += dt*(
                m_sigmaxx[i]*m_spinratexz[i] + m_sigmaxy[i]*m_spinrateyz[i] -
                m_spinratexy[i]*m_sigmayz[i] - m_spinratexz[i]*m_sigmazz[i]
            );
            dsigmayz[i] += dt*(
                m_sigmaxy[i]*m_spinratexz[i] + m_sigmayy[i]*m_spinrateyz[i] +
                m_spinratexy[i]*m_sigmaxz[i] - m_spinrateyz[i]*m_sigmazz[i]
            );
        }

        // update original stress states
        for (int i = 0; i < m_size; ++i) {
            m_sigmaxx[i] += dsigmaxx[i];
            m_sigmayy[i] += dsigmayy[i];
            m_sigmazz[i] += dsigmazz[i];
            m_sigmaxy[i] += dsigmaxy[i];
            m_sigmaxz[i] += dsigmaxz[i];
            m_sigmayz[i] += dsigmayz[i];
        }

    }
}

#endif