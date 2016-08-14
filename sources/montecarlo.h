#ifndef _MONTECARLO_H_
#define _MONTECARLO_H_

#include <cmath>
#include <ctime>
#include <vector>

#include "common.h"
#include "random.h"
#include "vec.h"
#include "ray.h"

class MonteCarlo {
public:
    MonteCarlo()
        : sigmap_s(0.0)
        , sigma_a(0.0)
        , eta(0.0) {
    }

    MonteCarlo(double sigS, double sigA, double ior)
        : sigmap_s(sigS)
        , sigma_a(sigA)
        , eta(ior) {
        sigmap_t = sigmap_s + sigma_a;
        alphap = sigmap_s / sigmap_t;
    }

    void compute(int samples, double dr, int divides) {
        hist.assign(divides, 0.0);

        Random rand((unsigned long)time(0));
        int accepted = 0;
        for (int s = 0; s < samples; s++) {
            double z  = rand.next01();
            double s2 = std::sqrt(1.0 - z * z);
            double phi = 2.0 * Pi * rand.next01();

            Vec dir(0.0, 0.0, -1.0);
            Ray ray(Vec(0.0, 0.0, 0.0), dir);
            double L = 1.0;
            for (int i = 0; i < 128; i++) {
                // Sample mean free path
                double t = -std::log((1.0 - rand.next01()) * eta) / sigmap_t;

                // Intersection test
                Vec next = ray.org + t * ray.dir;
                if (next.z >= 0.0) {
                    Vec po = Vec(0.0, 0.0, 0.0) - ray.org;
                    Vec n  = Vec(0.0, 0.0, -1.0);
                    double d = po.dot(n) / ray.dir.dot(n);
                    Vec p = ray.org + d * ray.dir;
                    if (std::abs(p.z) > 1.0e-8) {
                        std::cerr << "Error!!" << std::endl;
                        printf("%f %f %f\n", p.x, p.y, p.z);
                        std::exit(1);
                    }

                    double r = p.norm();
                    int index = (int)(r / dr);
                    if (index < divides) {
                        hist[index] += L;
                    }
                    accepted++;
                    break;
                }

                // Sample next direction
                double z = rand.next01() * 2.0 - 1.0;
                double s2 = std::sqrt(1.0 - z * z);
                double phi = 2.0 * Pi * rand.next01();

                Vec ndir = (std::cos(phi) * s2 * Vec(1.0, 0.0, 0.0) +
                            std::sin(phi) * s2 * Vec(0.0, 1.0, 0.0) +
                            z * Vec(0.0, 0.0, 1.0)).normalized();
                ray = Ray(next, ndir);
                L *= alphap;
            }

            if (s % 100 == 0) {
                printf("%6.2f %%\r", 100.0 * (s + 1) / samples);
            }
        }

        for (int i = 0; i < divides; i++) {
            double s = 2.0 * Pi * (i + 1) * dr * dr * accepted;
            hist[i] /= s;
        }
    }

    double operator[](int i) const {
        return hist[i];
    }

private:
    double sigmap_s, sigma_a, eta;
    double sigmap_t, alphap;
    std::vector<double> hist;
};

#endif  // _MONTECARLO_H_
