#ifndef _XCH_
#define _XCH_
#include <cmath>
inline double rs(double rho)
{
    return pow(3 / (4 * M_PI * rho), 1 / 3.);
}

class ExchangeCorrelation {
    //******************************************************************************/
    //  Calculates Exchange&Correlation Energy and Potential                       */ 
    //  type=0 - due to U.von.Barth and L.Hedin, J.Phys.C5, 1629 (1972)            */
    //  type=1 - O.E.Gunnarsson and S.Lundqvist,  Phys.Rev.B                       */
    //  type=2 - V.L.Moruzzi, J.F.Janak, and A.R.Williams, Calculated              */
    //           Electronic Properties of Metals (New York, Pergamon Press, 1978)  */
    //  type=3 - S.H.Vosko, L.Wilk, and M.Nusair, Can.J.Phys.58, 1200 (1980)       */
    //  type=4 - Correlation of Perdew and Wang 1991                               */
    //******************************************************************************/
    int type;
    double A, C;
    //static const double alphax = 0.610887057710857;//(3/(2 Pi))^(2/3)
    static constexpr double alphax = 0.610887057710857;//(3/(2 Pi))^(2/3)
    //static const double Aw = 0.0311;
    static constexpr double Aw = 0.0311;
    //static const double Bw = -0.048;
    static constexpr double Bw = -0.048;
    //static const double Cw = 0.002;
    static constexpr double Cw = 0.002;
    //static const double D = -0.0116;
    static constexpr double D = -0.0116;
    //static const double gamma = -0.1423;
    static constexpr double gamma = -0.1423;
    //static const double beta1 = 1.0529;
    static constexpr double beta1 = 1.0529;
    //static const double beta2 = 0.3334;
    static constexpr double beta2 = 0.3334;
    //static const double Ap = 0.0621814;
    static constexpr double Ap = 0.0621814;
    //static const double xp0 = -0.10498;
    static constexpr double xp0 = -0.10498;
    //static const double bp = 3.72744;
    static constexpr double bp = 3.72744;
    //static const double cp = 12.9352;
    static constexpr double cp = 12.9352;
    //static const double Qp = 6.1519908;
    static constexpr double Qp = 6.1519908;
    //static const double cp1 = 1.2117833;
    static constexpr double cp1 = 1.2117833;
    //static const double cp2 = 1.1435257;
    static constexpr double cp2 = 1.1435257;
    //static const double cp3 = -0.031167608;
    static constexpr double cp3 = -0.031167608;
public:
    ExchangeCorrelation(int type_) : type(type_)
    {
        switch (type) {
        case 0: C = 0.0504; A = 30; break;
        case 1: C = 0.0666; A = 11.4; break;
        case 2: C = 0.045;  A = 21; break;
        }
    };
    double Vx(double rs) { return -alphax / rs; }
    double ExVx(double rs) { return 0.25 * alphax / rs; }
    double Ex(double rs) { return -0.75 * alphax / rs; }
    double Vc(double rs)
    {
        if (type < 3) {
            double x = rs / A;
            return -0.5 * C * log(1 + 1 / x);
        }
        else if (type < 4) {// type=3 WVN
            double x = sqrt(rs);
            double xpx = x * x + bp * x + cp;
            double atnp = atan(Qp / (2 * x + bp));
            double ecp = 0.5 * Ap * (log(x * x / xpx) + cp1 * atnp - cp3 * (log(sqr(x - xp0) / xpx) + cp2 * atnp));
            return ecp - Ap / 6. * (cp * (x - xp0) - bp * x * xp0) / ((x - xp0) * xpx);
        }
        else {
            if (rs > 1) return gamma / (1 + beta1 * sqrt(rs) + beta2 * rs) * (1 + 7 / 6. * beta1 * sqrt(rs) + beta2 * rs) / (1 + beta1 * sqrt(rs) + beta2 * rs);
            else return Aw * log(rs) + Bw - Aw / 3. + 2 / 3. * Cw * rs * log(rs) + (2 * D - Cw) * rs / 3.;
        }
    }
    double EcVc(double rs)
    {
        if (type < 3) {
            double x = rs / A;
            double epsilon = -0.5 * C * ((1 + x * x * x) * log(1 + 1 / x) + 0.5 * x - x * x - 1 / 3.);
            return epsilon - Vc(rs);
        }
        else if (type < 4) {// type=3 WVN
            double x = sqrt(rs);
            return Ap / 6. * (cp * (x - xp0) - bp * x * xp0) / ((x - xp0) * (x * x + bp * x + cp));
        }
        else {
            if (rs > 1) return 2 * gamma / (1 + beta1 * sqrt(rs) + beta2 * rs) - Vc(rs);
            else return Aw * log(rs) + Bw + Cw * rs * log(rs) + D * rs - Vc(rs);
        }
    }
    double Vxc(double rs) { return Vx(rs) + Vc(rs); }
};
#endif //_XCH_