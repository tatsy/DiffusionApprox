#ifndef _COMMON_H_
#define _COMMON_H_

#include <cmath>

static const double Pi = 4.0 * std::atan(1.0);

class DipoleBase {
public:
    // Public methods
    explicit DipoleBase(double ior = 0.0) : eta(ior) {}
    virtual ~DipoleBase() {}

    virtual double Rd(double r) const = 0;

protected:
    // Protected methods
    double Fdr() const {
        const double eta2 = eta * eta;
        const double eta3 = eta2 * eta;
        if (eta >= 1.0) {
            return -1.4399 / eta2 + 0.7099 / eta + 0.6681 +
                    0.0636 * eta;
        } else {
            return -0.4399 + 0.7099 / eta - 0.3319 / eta2 +
                    0.0636 / eta3;
        }
    }

    double C1() const {
        const double eta2 = eta * eta;
        const double eta3 = eta2 * eta;
        const double eta4 = eta3 * eta;
        const double eta5 = eta4 * eta;
        if (eta < 1.0) {
            return (0.919317 - 3.4793 * eta + 6.75335 * eta2
                    -7.80989 * eta3 + 4.98554 * eta4 - 1.36881 * eta5) / 2.0;
        } else {
            return (-9.23372 + 22.2272 * eta - 20.9292 * eta2 + 10.2291 * eta3
                    - 2.54396 * eta4 + 0.254913 * eta5) / 2.0;
        }
    }

    double C2() const {
        const double eta2 = eta * eta;
        const double eta3 = eta2 * eta;
        const double eta4 = eta3 * eta;
        const double eta5 = eta4 * eta;
        if (eta < 1.0) {
            return (0.828421 - 2.62051 * eta + 3.36231 * eta2 - 1.95284 * eta3
                    + 0.236494 * eta4 + 0.145787 * eta5) / 3.0;
        } else {
            return (-1641.1 + 135.926 / eta3 - 656.175 / eta2 + 1376.53 / eta
                    + 1213.67 * eta - 568.556 * eta2 + 164.798 * eta3
                    - 27.0181 * eta4 + 1.91826 * eta5) / 3.0;
        }
    }

    // Protected parameters
    double eta;
};

#endif  // _COMMON_H_
