#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

#include "dipole.h"
#include "better_dipole.h"
#include "quantized_diffusion.h"
#include "montecarlo.h"

static const double sigmap_s = 1.00;
static const double sigma_a  = 0.1;
static const double eta      = 1.5;

static const int    divides  = 500;
static const double dr       = 0.03;

int main(int argc, char **argv) {
    // Compute Rd's and save the results to the file
    Dipole dp(sigmap_s, sigma_a, eta);
    BetterDipole btdp(sigmap_s, sigma_a, eta);
    QuantizedDiffusion qd(sigmap_s, sigma_a, eta);
    MonteCarlo mc(sigmap_s, sigma_a, eta);
    mc.compute(5000000, dr, divides);

    FILE *fp = fopen("output.dat", "w");
    for (int i = 0; i < divides; i++) {
        const double r = i * dr;
        const double dval = std::sqrt(2.0 * Pi * r * dp.Rd(r));
        const double bval = std::sqrt(2.0 * Pi * r * btdp.Rd(r));
        const double qval = std::sqrt(2.0 * Pi * r * qd.Rd(r));
        const double mval = std::sqrt(2.0 * Pi * r * mc[i]);
        fprintf(fp, "%12.8f %12.8f %12.8f %12.8f %12.8f\n", r, dval, bval, qval, mval);
    }
    fclose(fp);
}
