#ifndef _XCH_
#define _XCH_

#include <numbers> // For std::numbers::pi

inline double rs(double rho) {
    return std::pow(3 / (4 * std::numbers::pi * rho), 1.0 / 3);
}

class ExchangeCorrelation {
    int type;
    double A, C;

    static constexpr double alphax = 0.610887057710857; // (3/(2 * pi))^(2/3)
    static constexpr double Aw = 0.0311;
    static constexpr double Bw = -0.048;
    static constexpr double Cw = 0.002;
    static constexpr double D = -0.0116;
    static constexpr double gamma = -0.1423;
    static constexpr double beta1 = 1.0529;
    static constexpr double beta2 = 0.3334;
    static constexpr double Ap = 0.0621814;
    static constexpr double xp0 = -0.10498;
    static constexpr double bp = 3.72744;
    static constexpr double cp = 12.9352;
    static constexpr double Qp = 6.1519908;
    static constexpr double cp1 = 1.2117833;
    static constexpr double cp2 = 1.1435257;
    static constexpr double cp3 = -0.031167608;

public:
    explicit ExchangeCorrelation(int type_) : type(type_) {
        switch (type) {
        case 0:
            C = 0.0504;
            A = 30;
            break;
        case 1:
            C = 0.0666;
            A = 11.4;
            break;
        case 2:
            C = 0.045;
            A = 21;
            break;
        }
    }

    [[nodiscard]] double Vx(double rs) const {
        return -alphax / rs;
    }

    [[nodiscard]] double Ex(double rs) const {
        return -0.75 * alphax / rs;
    }

    [[nodiscard]] double ExVx(double rs) const {
        return 0.25 * alphax / rs;
    }

    [[nodiscard]] double Vc(double rs) const {
        if (type < 3) {
            double x = rs / A;
            return -0.5 * C * std::log(1 + 1 / x);
        } else if (type < 4) { // type=3 WVN
            double x = std::sqrt(rs);
            double xpx = x * x + bp * x + cp;
            double atnp = std::atan(Qp / (2 * x + bp));
            double ecp = 0.5 * Ap * (std::log(x * x / xpx) + cp1 * atnp - cp3 * (std::log((x - xp0) * (x - xp0) / xpx) + cp2 * atnp));
            return ecp - Ap / 6. * (cp * (x - xp0) - bp * x * xp0) / ((x - xp0) * xpx);
        } else {
            if (rs > 1) return gamma / (1 + beta1 * std::sqrt(rs) + beta2 * rs) * (1 + 7 / 6. * beta1 * std::sqrt(rs) + beta2 * rs) / (1 + beta1 * std::sqrt(rs) + beta2 * rs);
            else return Aw * std::log(rs) + Bw - Aw / 3. + 2 / 3. * Cw * rs * std::log(rs) + (2 * D - Cw) * rs / 3.;
        }
    }

    [[nodiscard]] double EcVc(double rs) const {
        if (type < 3) {
            double x = rs / A;
            double epsilon = -0.5 * C * ((1 + x * x * x) * std::log(1 + 1 / x) + 0.5 * x - x * x - 1 / 3.);
            return epsilon - Vc(rs);
        } else if (type < 4) { // type=3 VWN
            double x = std::sqrt(rs);
            return Ap / 6. * (cp * (x - xp0) - bp * x * xp0) / ((x - xp0) * (x * x + bp * x + cp));
        } else {
            if (rs > 1) return 2 * gamma / (1 + beta1 * std::sqrt(rs) + beta2 * rs) - Vc(rs);
            else return Aw * std::log(rs) + Bw + Cw * rs * std::log(rs) + D * rs - Vc(rs);
        }
    }

    [[nodiscard]] double Vxc(double rs) const {
        return Vx(rs) + Vc(rs);
    }
};

#endif //_XCH_
