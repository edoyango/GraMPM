#ifndef GRAMPM_INTEGRATORS
#define GRAMPM_INTEGRATORS

#include <grampm.hpp>
#include <iostream>

namespace GraMPM {
    namespace integrators {
        template<typename F>
        void MUSL(GraMPM::particle_system<F> &p, const F &dt, const int &max_timestep, 
            const int &print_timestep_interval, const int &save_timestep_interval) {

            p.m_tmpgmass.resize(p.background_grid.ncells());
            p.m_tmpgmomentumx.resize(p.background_grid.ncells());
            p.m_tmpgmomentumy.resize(p.background_grid.ncells());
            p.m_tmpgmomentumz.resize(p.background_grid.ncells());
            p.m_tmpgforcex.resize(p.background_grid.ncells());
            p.m_tmpgforcey.resize(p.background_grid.ncells());
            p.m_tmpgforcez.resize(p.background_grid.ncells());

            #pragma omp parallel default(shared)
            {
            for (int itimestep = 1; itimestep < max_timestep+1; ++itimestep) {

                // print to terminal
                if (itimestep % print_timestep_interval == 0) {
                    #pragma omp single
                    std::cout << itimestep << ' ' << dt*itimestep << '\n';
                }

                // get which cell each particle is located within
                p.update_particle_to_cell_map();

                /* get all neighbour nodes of each particle, calculate dx between particles and adjacent nodes, and
                update kernel and kernel gradient values. */ 
                p.map_particles_to_grid();

                // map particles' mass to nodes
                // #pragma omp single
                p.map_mass_to_grid();

                // map particles' momentum to nodes
                // #pragma omp single
                p.map_momentum_to_grid();

                // apply user-defined momentum boundary conditions to grid
                // #pragma omp single
                p.background_grid.apply_momentum_boundary_conditions(itimestep, dt);

                // map particles' force to nodes
                // #pragma omp single
                p.map_force_to_grid();

                // apply user-defined force boundary conditions to grid
                // #pragma omp single
                p.background_grid.apply_force_boundary_conditions(itimestep, dt);

                // update nodal momentums
                // #pragma omp single
                p.background_grid.update_momentum(dt);

                // map nodal forces to particle accelerations
                // #pragma omp single
                p.map_acceleration_to_particles();

                // update particles' velocities with calculated accelerations
                // #pragma omp single
                p.update_velocity(dt);
                // update particles' position
                // #pragma omp single
                p.update_position(dt);

                // map particles' momentum to nodes, in preparation for updating stress
                // #pragma omp single
                p.map_momentum_to_grid();

                // apply user-defined momentum boundary conditions to grid
                // #pragma omp single
                p.background_grid.apply_momentum_boundary_conditions(itimestep, dt);

                // map nodal velocities to particle strain/spin rates
                // #pragma omp single
                p.map_strainrate_to_particles();

                // update particles' density
                // #pragma omp single
                p.update_density(dt);

                // update particles' stress
                // #pragma omp single
                p.update_stress(dt);

                if (itimestep % save_timestep_interval == 0) {
                    #pragma omp single
                    p.save_to_file("p_", itimestep);
                }

            }
            }
        }
    }
    
}
#endif