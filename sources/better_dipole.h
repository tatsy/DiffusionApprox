#ifndef _BETTER_DIPOLE_H_
#define _BETTER_DIPOLE_H_

#include <cstdio>

#include "common.h"

class BetterDipole : public DipoleBase {
public:
    // Public methods
    BetterDipole()
        : DipoleBase(0.0)
        , sigmap_s(0.0)
        , sigma_a(0.0) {
    }

    BetterDipole(double sigS, double sigA, double ior)
        : DipoleBase(ior)
        , sigmap_s(sigS)
        , sigma_a(sigA) {
        init();
    }

    ~BetterDipole() {
    }

    BetterDipole(const BetterDipole& dp) = default;
    BetterDipole& operator=(const BetterDipole& dp) = default;

    void init() {
        sigmap_t = sigmap_s + sigma_a;
        alphap = sigmap_s / sigmap_t;
        A = (1.0 + 3.0 * C2()) / (1.0 - 2.0 * C1());
        D = (2.0 * sigma_a + sigmap_s) / (3.0 * sigmap_t * sigmap_t);
        sigma_tr = std::sqrt(sigma_a / D);
        zr = 1.0 / sigmap_t;
        zb = 2.0 * A * D;
        zv = -zr - 2.0 * zb;
        Cphi = 0.25 * (1.0 - 2.0 * C1());
        CE   = 0.50 * (1.0 - 3.0 * C2());
    }

    double Rd(double r) const override {
        const double r2 = r * r;
        const double dr = std::sqrt(r2 + zr * zr);
        const double dv = std::sqrt(r2 + zv * zv);
        const double ar = zr * (sigma_tr * dr + 1.0) / (dr * dr);
        const double av = zv * (sigma_tr * dv + 1.0) / (dv * dv);
        const double br = std::exp(-sigma_tr * dr) / dr;
        const double bv = std::exp(-sigma_tr * dv) / dv;
        return (alphap * alphap / (4.0 * Pi)) *
               ((CE * ar + Cphi / D) * br -
                (CE * av + Cphi / D) * bv);
    }

private:
    double sigmap_s, sigma_a;
    double sigmap_t, alphap, A, D, sigma_tr;
    double zr, zb, zv, Cphi, CE;
};

#endif  // _BETTER_DIPOLE_H_
