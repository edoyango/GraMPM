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
        , m_mass(size, 0.)
        , m_rho(size, 0.)
        , m_sigmaxx(size, 0.)
        , m_sigmayy(size, 0.)
        , m_sigmazz(size, 0.)
        , m_sigmaxy(size, 0.)
        , m_sigmaxz(size, 0.)
        , m_sigmayz(size, 0.)
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
    template<typename F>const F* particle_system<F>::x() const { return m_x.data(); }
    template<typename F>const F& particle_system<F>::y(const int &i) const { return m_y[i]; }
    template<typename F>const F* particle_system<F>::y() const { return m_y.data(); }
    template<typename F>const F& particle_system<F>::z(const int &i) const { return m_z[i]; }
    template<typename F>const F* particle_system<F>::z() const { return m_z.data(); }
    template<typename F>const F& particle_system<F>::vx(const int &i) const { return m_vx[i]; }
    template<typename F>const F* particle_system<F>::vx() const { return m_vx.data(); }
    template<typename F>const F& particle_system<F>::vy(const int &i) const { return m_vy[i]; }
    template<typename F>const F* particle_system<F>::vy() const { return m_vy.data(); }
    template<typename F>const F& particle_system<F>::vz(const int &i) const { return m_vz[i]; }
    template<typename F>const F* particle_system<F>::vz() const { return m_vz.data(); }
    template<typename F>const F& particle_system<F>::mass(const int &i) const { return m_mass[i]; }
    template<typename F>const F* particle_system<F>::mass() const { return &m_mass; }
    template<typename F>const F& particle_system<F>::rho(const int &i) const { return m_rho[i]; }
    template<typename F>const F* particle_system<F>::rho() const { return &m_rho; }
    template<typename F>const F& particle_system<F>::sigmaxx(const int &i) const { return m_sigmaxx[i]; }
    template<typename F>const F* particle_system<F>::sigmaxx() const { return &m_sigmaxx; }
    template<typename F>const F& particle_system<F>::sigmayy(const int &i) const { return m_sigmayy[i]; }
    template<typename F>const F* particle_system<F>::sigmayy() const { return &m_sigmayy; }
    template<typename F>const F& particle_system<F>::sigmazz(const int &i) const { return m_sigmazz[i]; }
    template<typename F>const F* particle_system<F>::sigmazz() const { return &m_sigmazz; }
    template<typename F>const F& particle_system<F>::sigmaxy(const int &i) const { return m_sigmaxy[i]; }
    template<typename F>const F* particle_system<F>::sigmaxy() const { return &m_sigmaxy; }
    template<typename F>const F& particle_system<F>::sigmaxz(const int &i) const { return m_sigmaxz[i]; }
    template<typename F>const F* particle_system<F>::sigmaxz() const { return &m_sigmaxz; }
    template<typename F>const F& particle_system<F>::sigmayz(const int &i) const { return m_sigmayz[i]; }
    template<typename F>const F* particle_system<F>::sigmayz() const { return &m_sigmayz; }
    template<typename F>const std::array<F, 3>& particle_system<F>::body_force() const { return m_body_force; }
    template<typename F>const F& particle_system<F>::body_force(const int &i) const { return m_body_force[i]; }
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
    template<typename F> void particle_system<F>::set_mass(const int &i, const F &m) { m_mass[i] = m; }
    template<typename F> void particle_system<F>::set_rho(const int &i, const F &rho) { m_rho[i] = rho; }
    template<typename F> void particle_system<F>::set_sigmaxx(const int &i, const F &sigmaxx) { m_sigmaxx[i] = sigmaxx; }
    template<typename F> void particle_system<F>::set_sigmayy(const int &i, const F &sigmayy) { m_sigmayy[i] = sigmayy; }
    template<typename F> void particle_system<F>::set_sigmazz(const int &i, const F &sigmazz) { m_sigmazz[i] = sigmazz; }
    template<typename F> void particle_system<F>::set_sigmaxy(const int &i, const F &sigmaxy) { m_sigmaxy[i] = sigmaxy; }
    template<typename F> void particle_system<F>::set_sigmaxz(const int &i, const F &sigmaxz) { m_sigmaxz[i] = sigmaxz; }
    template<typename F> void particle_system<F>::set_sigmayz(const int &i, const F &sigmayz) { m_sigmayz[i] = sigmayz; }
    template<typename F> void particle_system<F>::set_body_force(const std::array<F, 3> &bf) { m_body_force = bf; }
    template<typename F> void particle_system<F>::set_body_force(const F &bfx, const F &bfy, const F &bfz) { 
        m_body_force[0] = bfx;
        m_body_force[1] = bfy;
        m_body_force[2] = bfz;
    }
    template<typename F> void particle_system<F>::set_grid_index(const int &i, const int &idx) { m_grid_idx[i] = idx; }
    template<typename F> void particle_system<F>::incrementNParticles() {m_size++;}

    // vector-like api: at. Returns particle class
    template<typename F>
    particle<F> particle_system<F>::at(const int &i) { 
        particle<F> p(x(i), y(i), z(i), vx(i), vy(i), vz(i), mass(i), rho(i), sigmaxx(i), sigmayy(i), sigmazz(i),
            sigmaxy(i), sigmaxz(i), sigmayz(i));
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
        m_mass.push_back(p.mass);
        m_rho.push_back(p.rho);
        m_sigmaxx.push_back(p.sigmaxx);
        m_sigmayy.push_back(p.sigmayy);
        m_sigmazz.push_back(p.sigmazz);
        m_sigmaxy.push_back(p.sigmaxy);
        m_sigmaxz.push_back(p.sigmaxz);
        m_sigmayz.push_back(p.sigmayz);
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
        m_mass.reserve(n);
        m_rho.reserve(n);
        m_sigmaxx.reserve(n);
        m_sigmayy.reserve(n);
        m_sigmazz.reserve(n);
        m_sigmaxy.reserve(n);
        m_sigmaxz.reserve(n);
        m_sigmayz.reserve(n);
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
        m_mass.clear();
        m_rho.clear();
        m_sigmaxx.clear();
        m_sigmayy.clear();
        m_sigmazz.clear();
        m_sigmaxy.clear();
        m_sigmaxz.clear();
        m_sigmayz.clear();
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
            m_mass.empty() && m_grid_idx.empty() &&
            m_p2g_neighbour_nodes.empty() && m_p2g_neighbour_nodes_dx.empty() && m_p2g_neighbour_nodes_dy.empty() && 
            m_p2g_neighbour_nodes_dz.empty() && m_p2g_neighbour_nodes_w.empty() && m_p2g_neighbour_nodes_dwdx.empty() &&
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
        m_mass.resize(n, 0.);
        m_rho.resize(n, 0.);
        m_sigmaxx.resize(n, 0.);
        m_sigmayy.resize(n, 0.);
        m_sigmazz.resize(n, 0.);
        m_sigmaxy.resize(n, 0.);
        m_sigmaxz.resize(n, 0.);
        m_sigmayz.resize(n, 0.);
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
        m_mass.resize(n, p.mass);
        m_rho.resize(n, p.rho);
        m_sigmaxx.resize(n, p.rho);
        m_sigmayy.resize(n, p.rho);
        m_sigmazz.resize(n, p.rho);
        m_sigmaxy.resize(n, p.rho);
        m_sigmaxz.resize(n, p.rho);
        m_sigmayz.resize(n, p.rho);
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

        for (int i = 0; i < m_size; ++i) {
            for (int j = 0; j < m_nneighbour_nodes_perp; ++j) {
                const int node_idx = p2g_neighbour_node(i, j);
                (*g_property)[node_idx] += p2g_neighbour_node_w(i, j)*p_property[i];
            }
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
}

#endif