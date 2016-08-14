#ifndef _QUANTIZED_DIFFUSION_H_
#define _QUANTIZED_DIFFUSION_H_

#include <vector>

#include "common.h"

class QuantizedDiffusion : public DipoleBase {
public:
    QuantizedDiffusion()
        : DipoleBase(0.0)
        , sigmap_s(0.0)
        , sigma_a(0.0)
        , nPoles(0)
        , depth(0.0) {
    }

    QuantizedDiffusion(double sigS, double sigA, double ior,
                        int n = 0, double d = 1.0e5)
        : DipoleBase(ior)
        , sigmap_s(sigS)
        , sigma_a(sigA)
        , nPoles(n)
        , depth(d) {
        init();
    }

    ~QuantizedDiffusion() {}

    QuantizedDiffusion(const QuantizedDiffusion& qd) = default;
    QuantizedDiffusion& operator=(const QuantizedDiffusion& qd) = default;

    void init() {
        sigmap_t = sigma_a + sigmap_s;
        alphap = sigmap_s / sigmap_t;
        A = (1.0 + 3.0 * C2()) / (1.0 - 2.0 * C1());
        D = (2.0 * sigma_a + sigmap_s) / (3.0 * sigmap_t * sigmap_t);
        zb = 2.0 * A * D;
        Cphi = 0.25 * (1.0 - 2.0 * C1());
        CE   = 0.50 * (1.0 - 3.0 * C2());
        tau.resize(k + 1);
        tau[0] = t1;
        for (int i = 1; i <= k; i++) {
            tau[i] = tau[i - 1] * s;
        }
    }

    double Rd(double r) const override {
        double ret = 0.0;
        for (int i = 0; i < k; i++) {
            double vi = D * (tau[i] + tau[i + 1]);
            ret += w_R(vi) * w(i) * G2D(vi, r);
        }
        return alphap * ret;
    }

    double w(int i) const {
        return (std::exp(-tau[i] * sigma_a) - std::exp(-tau[i + 1] * sigma_a)) / sigma_a;
    }

    double G2D(double v, double r) const {
        return (1.0 / (2.0 * Pi * v)) * std::exp(- r * r / (2.0 * v));
    }

    double w_R(double vi) const {
        double ret = 0.0;
        for (int j = -nPoles; j <= nPoles; j++) {
            double m_r = 2.0 * j * (depth + 2.0 * zb);
            double m_v = 2.0 * j * (depth + 2.0 * zb) - 2.0 * zb;
            double w_phiR = w_phi(vi, 0.0, depth, m_r) + w_phi(vi, 0.0, depth, -m_v);
            double w_ER   = w_E(vi, 0.0, depth, m_r)   - w_E(vi, 0.0, depth, -m_v);
            ret += Cphi * w_phiR + CE * w_ER;
        }
        return ret;
    }

    double w_phi(double v, double z1, double z2, double m) const {
        double a1 = (m + sigmap_t * v + z1) / std::sqrt(2.0 * v);
        double a2 = (m + sigmap_t * v + z2) / std::sqrt(2.0 * v);
        return (alphap * sigmap_t * 0.5) *
                std::exp(m * sigmap_t + 0.5 * sigmap_t * sigmap_t * v) *
                (std::erf(a2) - std::erf(a1));
    }

    double w_E(double v, double z1, double z2, double m) const {
        double a1 = (m * m) / (2.0 * v) + (m + sigmap_t * v) * z1 / v + (z1 * z1) / (2.0 * v);
        double a2 = (m * m) / (2.0 * v) + (m + sigmap_t * v) * z2 / v + (z2 * z2) / (2.0 * v);
        return D * sigmap_t *
               (-w_phi(v, z1, z2, m) + alphap * (std::exp(-a1) - std::exp(-a2)) / std::sqrt(2.0 * Pi * v));
    }

private:
    double sigmap_s, sigma_a, depth;
    double sigmap_t, alphap, A, D, zr, zb;
    double Cphi, CE;
    int nPoles;

    static constexpr int k = 40;
    static constexpr double t1 = 1.0e-6;
    static constexpr double s = 1.618;
    std::vector<double> tau;
};

#endif  // _QUANTIZED_DIFFUSION_H_
