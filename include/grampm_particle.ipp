#ifndef GRAMPM_particle_ipp
#define GRAMPM_particle_ipp

namespace GraMPM {

    template<typename F>
    particle<F>::particle(const F inx, const F iny, const F inz, const F invx, const F invy, const F invz, 
        const F inmass, const F inrho, const F insigmaxx, const F insigmayy, const F insigmazz, const F insigmaxy,
        const F insigmaxz, const F insigmayz)
        : x {inx, iny, inz}
        , v {invx, invy, invz}
        , a {0., 0., 0.}
        , dxdt {0., 0., 0.}
        , sigma {insigmaxx, insigmayy, insigmazz, insigmaxy, insigmaxz, insigmayz}
        , strainrate {0., 0., 0., 0., 0., 0.}
        , spinrate {0., 0., 0.}
        , mass {inmass}
        , rho {inrho}
    {
    }

    template<typename F>
    particle<F>::particle(const F inx, const F iny, const F inz, const F invx, const F invy, const F invz, 
        const F inmass, const F inrho, const F insigmaxx, const F insigmayy, const F insigmazz, const F insigmaxy,
        const F insigmaxz, const F insigmayz, const F inax, const F inay, const F inaz, const F indxdt, 
        const F indydt, const F indzdt, const F instrainratexx, const F instrainrateyy, const F instrainratezz, 
        const F instrainratexy, const F instrainratexz, const F instrainrateyz, const F inspinratexy,
        const F inspinratexz, const F inspinrateyz)
        : x {inx, iny, inz}
        , v {invx, invy, invz}
        , a {inax, inay, inaz}
        , dxdt {indxdt, indydt, indzdt}
        , sigma {insigmaxx, insigmayy, insigmazz, insigmaxy, insigmaxz, insigmayz}
        , strainrate {instrainratexx, instrainrateyy, instrainratezz, instrainratexy, instrainratexz, instrainrateyz}
        , spinrate {inspinratexy, inspinratexz, inspinrateyz}
        , mass {inmass}
        , rho {inrho}
    {
    }

    template<typename F>
    particle<F>::particle(const std::array<F, 3> inx, const std::array<F, 3> inv, const F inmass, const F inrho, 
        const std::array<F, 6> insigma, const std::array<F, 3> ina, const std::array<F, 3> indxdt, 
        const std::array<F, 6> instrainrate, const std::array<F, 3> inspinrate)
        : x {inx}
        , v {inv}
        , a {ina}
        , dxdt {indxdt}
        , sigma {insigma}
        , strainrate {instrainrate}
        , spinrate {inspinrate}
        , mass {inmass}
        , rho {inrho}
    {
    }

    template<typename F> particle<F>::particle() {};
}

#endif