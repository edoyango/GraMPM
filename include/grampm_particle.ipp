#ifndef GRAMPM_particle_ipp
#define GRAMPM_particle_ipp

namespace GraMPM {

    template<typename F>
    particle<F>::particle(const F &inx, const F &iny, const F &inz, const F &invx, const F &invy, const F &invz, 
        const F &inmass, const F &inrho)
        : x {inx}
        , y {iny}
        , z {inz}
        , vx {invx}
        , vy {invy}
        , vz {invz}
        , mass {inmass}
        , rho {inrho}
    {
    }
}

#endif