#ifndef _DIPOLE_H_
#define _DIPOLE_H_

#include <cmath>

#include "common.h"

class Dipole : public DipoleBase {
public:
    Dipole()
        : DipoleBase(0.0)
        , sigmap_s(0.0)
        , sigma_a(0.0) {
    }

    Dipole(double sigS, double sigA, double ior)
        : DipoleBase(ior)
        , sigmap_s(sigS)
        , sigma_a(sigA) {
        init();
    }

    Dipole(const Dipole& dp) = default;
    Dipole& operator=(const Dipole& dp) = default;

    ~Dipole() {
    }

    double Rd(double r) const override {
        const double r2 = r * r;
        const double dr = std::sqrt(r2 + zr * zr);
        const double dv = std::sqrt(r2 + zv * zv);
        const double phi_r = zr * (dr * sigma_tr + 1.0) * std::exp(-sigma_tr * dr) / (dr * dr * dr);
        const double phi_v = zv * (dv * sigma_tr + 1.0) * std::exp(-sigma_tr * dv) / (dv * dv * dv);
        return (alphap / (4.0 * Pi)) * (phi_r - phi_v);
    }

private:
    // Private methods
    void init() {
        A = (1.0 + Fdr()) / (1.0 - Fdr());
        sigmap_t = sigmap_s + sigma_a;
        sigma_tr = std::sqrt(3.0 * sigma_a * sigmap_t);
        alphap = sigmap_s / sigmap_t;
        zr = 1.0 / sigmap_t;
        zv = -zr * (1.0 + (4.0 / 3.0) * A);
    }

    // Private parameters
    double sigmap_s, sigma_a;
    double sigmap_t, sigma_tr, alphap, A;
    double zr, zv;
};

#endif  // _DIPOLE_H_
