#ifndef GRAMPM_INTEGRATORS
#define GRAMPM_INTEGRATORS

#include <grampm.hpp>
#include <iostream>

namespace GraMPM {
    namespace integrators {
        template<typename F>
        void MUSL(GraMPM::MPM_system<F> &p, const F &dt, const size_t &max_timestep, 
            const size_t &print_timestep_interval, const size_t &save_timestep_interval) {

            #pragma omp parallel default(shared)
            {
            for (size_t itimestep = 1; itimestep < max_timestep+1; ++itimestep) {

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
                p.map_p2g_mass();

                // map particles' momentum to nodes
                p.map_p2g_momentum();

                // apply user-defined momentum boundary conditions to grid
                p.g_apply_momentum_boundary_conditions(itimestep, dt);

                // map particles' force to nodes
                p.map_p2g_force();

                // apply user-defined force boundary conditions to grid
                p.g_apply_force_boundary_conditions(itimestep, dt);

                // update nodal momentums
                p.g_update_momentum(dt);

                // map nodal forces to particle accelerations
                p.map_g2p_acceleration();

                // update particles' velocities with calculated accelerations
                p.p_update_velocity(dt);

                // update particles' position
                p.p_update_position(dt);

                // map particles' momentum to nodes, in preparation for updating stress
                p.map_p2g_momentum();

                // apply user-defined momentum boundary conditions to grid
                p.g_apply_momentum_boundary_conditions(itimestep, dt);

                // map nodal velocities to particle strain/spin rates
                p.map_g2p_strainrate();

                // update particles' density
                p.p_update_density(dt);

                // update particles' stress
                p.p_update_stress(dt);

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