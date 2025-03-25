/*#include <cstdio>

int main()
{
    printf("hello from LAPW!\n");
    return 0;
}*/
#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <list>
#include <algorithm>
#include "numerov.h"
#include "integrate.h"
#include "kmesh.h"
#include "bessel.h"
#include "zeroin.h"
//#include "erf.h"
#include <cmath>
#include "function.h"
#include "blas.h"
#include "exchcorr.h"
#include "core.h"
#include "util.h"

class FccLattice { // Class for storing reciprocal lattice
    double LatConst, Volume;
    dvector3 a0, a1, a2;    // Primitive vectors of fcc lattice
    dvector3 b0, b1, b2;    // Primitive vectors of reciprocal lattice
    dvector3 GammaPoint, LPoint, KPoint, XPoint, WPoint; // Special points in 1IRB
    vector<dvector3> Kmesh, kmesh;
    vector<double> wkp;     // weights for irreducible k-points
public:
    double Vol() const { return Volume; }
    int Ksize() const { return Kmesh.size(); }
    int ksize() const { return kmesh.size(); }
    double wk(int ik) const { return wkp[ik]; }
    const vector<double>& wk() const { return wkp; }
    const vector<dvector3>& Km() { return Kmesh; }
    const dvector3& K(int i) const { return Kmesh[i]; } // can not be changed, only read
    const dvector3& k(int i) const { return kmesh[i]; } // can not be changed, only read
    FccLattice(double LatConst_) : LatConst(LatConst_)
    {
        a0 = dvector3(0.5 * LatConst, 0.5 * LatConst, 0);
        a1 = dvector3(0.5 * LatConst, 0, 0.5 * LatConst);
        a2 = dvector3(0, 0.5 * LatConst, 0.5 * LatConst);
        Volume = fabs(Vproduct(a0, a1, a2));// Volume
        clog << "Volume is " << Volume << endl;
        b0 = (2 * M_PI / Volume) * cross(a1, a0);
        b1 = (2 * M_PI / Volume) * cross(a0, a2);
        b2 = (2 * M_PI / Volume) * cross(a2, a1);
        // Special points in Brillouin zone
        double brs = 2 * M_PI / LatConst;
        GammaPoint = dvector3(0, 0, 0);
        LPoint = dvector3(0.5 * brs, 0.5 * brs, 0.5 * brs);
        KPoint = dvector3(0.75 * brs, 0.75 * brs, 0);
        XPoint = dvector3(1 * brs, 0, 0);
        WPoint = dvector3(1 * brs, 0.5 * brs, 0);
    }
    void GenerateReciprocalVectors(int q, double CutOffK)
    {
        // Many reciprocal vectors are generated and later only the sortest are used
        list<dvector3> Kmesh0;
        for (int n = -q; n < q; n++) {
            for (int l = -q; l < q; l++) {
                for (int m = -q; m < q; m++) {
                    Kmesh0.push_back(n * b0 + l * b1 + m * b2);
                }
            }
        }
        Kmesh0.sort(cmp); // Sorting according to the length of vector. Shortest will be kept
        int Ksize = 0;
        for (list<dvector3>::const_iterator l = Kmesh0.begin(); l != Kmesh0.end(); l++, Ksize++) if (l->length() > CutOffK) break;
        Kmesh.resize(Ksize);
        int j = 0;
        for (list<dvector3>::const_iterator l = Kmesh0.begin(); l != Kmesh0.end() && j < Ksize; l++, j++) Kmesh[j] = *l;
        clog << "K-mesh size = " << Kmesh.size() << endl;
    }
    void ChoosePointsInFBZ(int nkp, int type = 0) {// Chooses the path in the 1BZ we will use
        if (type == 0) { // Choose mesh in the 1BZ to cover the whole space - for SC calculation
            list<dvector3> kp;
            Generate_All_K_Points(nkp, b0, b1, b2, kp);
            clog << "Number of all k-points = " << kp.size() << endl;
            ChooseIrreducible(kp, kmesh, wkp);
            clog << "Number of irreducible k points is " << kmesh.size() << endl;
        }
        else {        // Choose one particular path in the 1BZ - for plotting purposes
            int tm = 4 * static_cast<int>(nkp * nkp * nkp / 4.);
            int nkp = tm;
            clog << "nkp = " << nkp << endl;
            kmesh.resize(nkp);
            int N0 = kmesh.size() / 4;
            for (int i = 0; i < N0; i++) kmesh[i] = GammaPoint + (XPoint - GammaPoint) * i / (N0 - 1.);
            for (int i = 0; i < N0; i++) kmesh[N0 + i] = XPoint + (LPoint - XPoint) * i / (N0 - 1.);
            for (int i = 0; i < N0; i++) kmesh[N0 * 2 + i] = LPoint + (GammaPoint - LPoint) * i / (N0 - 1.);
            for (int i = 0; i < N0; i++) kmesh[N0 * 3 + i] = GammaPoint + (KPoint - GammaPoint) * i / (N0 - 1.);
        }

    }
    double RMuffinTin() const { return 0.5 * a0.length(); }// Spheres just touch
};

// This is the parametrization for the effective potential from experiment (for Cu)
double VeffP(double R)
{
    return 29 * exp(-2.3151241717834 * pow(R, 0.81266614122432) + 2.1984250222603e-2 * pow(R, 4.2246376280056))
        - 0.15595606773483 * R - 3.1350051440417e-3 * R * R + 5.1895222293006e-2 * pow(R, 3) - 2.8027608685637e-2 * pow(R, 4);
}

double extrapolate_f0(double f1, double f2, double x0, double x1, double x2)
{//extrapolation to zero
    return (f1 * (x2 - x0) - f2 * (x1 - x0)) / (x2 - x1);
}


class PartialWave {// Class for solving SCH equation
    int Z;
    vector<double> Rmesh;
    vector<double> rhs_MT, ur, urp, temp, inhom; // For solving SCH equation
    vector<double> dlogPsi, dlogPsip, Psi, Psip, PsipPsip;
    vector<vector<double> > Psi_l, Psip_l;
    vector<double> Veff0;    // Effective KS potential without centrifugal part
public:
    PartialWave(int N, double RMuffinTin, int Z_, int lMax) : Z(Z_), Rmesh(N), rhs_MT(N), ur(N), urp(N), temp(N), inhom(N),
        dlogPsi(lMax + 1), dlogPsip(lMax + 1), Psi(lMax + 1), Psip(lMax + 1), PsipPsip(lMax + 1),
        Psi_l(lMax + 1), Psip_l(lMax + 1), Veff0(N)
    {
        // Building equidistant radial mesh
        double dh = RMuffinTin / (N - 1.);
        for (int i = 0; i < N; i++) Rmesh[i] = i * dh;
    }
    int Rsize() const { return Rmesh.size(); }
    double R(int i) const { return Rmesh[i]; }
    const vector<double> R() const { return Rmesh; }
    double psi(int l) const { return Psi[l]; }
    double psip(int l) const { return Psip[l]; }
    double psi(int l, int r) const { return Psi_l[l][r]; }
    double psip(int l, int r) const { return Psip_l[l][r]; }
    double dlog(int l) const { return dlogPsi[l]; }
    double dlogp(int l) const { return dlogPsip[l]; }
    double PP(int l) const { return PsipPsip[l]; }
    double startSol(int Z, int l, double r) const// good choice for starting Numerov algorithm
    {
        return pow(r, l + 1) * (1 - Z * r / (l + 1));
    }
    // This function sets KS potential without centrifugal part
    void SetVeff0(const vector<double>& Veff)
    {
        for (int i = 0; i < Rmesh.size(); i++) Veff0[i] = Veff[i];
    }
    double SolveSCHEquation(const vector<double>& Enu)
    {
        for (int l = 0; l < Psi.size(); l++) {
            rhs_MT[0] = 0;
            for (int i = 1; i < Rmesh.size(); i++) {
                //	double Veff = -VeffP(Rmesh[i])/Rmesh[i] + 0.5*l*(l+1)/sqr(Rmesh[i]);
                double Veff = Veff0[i] + 0.5 * l * (l + 1) / sqr(Rmesh[i]);
                rhs_MT[i] = 2 * (Veff - Enu[l]);
            }
            double dh = Rmesh[1] - Rmesh[0];
            ur[0] = 0;
            ur[1] = startSol(Z, l, dh);
            // Solving radial SCH equation
            Numerov(rhs_MT, ur.size(), dh, ur);
            // Normalizing the result 
            for (int i = 0; i < Rmesh.size(); i++) temp[i] = sqr(ur[i]);
            double norm = 1. / sqrt(integrate4<double>(temp, dh, temp.size()));
            for (int i = 0; i < ur.size(); i++) ur[i] *= norm;

            Psi_l[l] = ur;// storing for future (density)

            for (int ir = 1; ir < Rmesh.size(); ir++) Psi_l[l][ir] /= Rmesh[ir]; // Psi=u/r!
            Psi_l[l][0] = extrapolate_f0(Psi_l[l][1], Psi_l[l][2], Rmesh[0], Rmesh[1], Rmesh[2]);
            //      Psi_l[l][0] = (Psi_l[l][1]*(Rmesh[2]-Rmesh[0])-Psi_l[l][2]*(Rmesh[1]-Rmesh[0]))/(Rmesh[2]-Rmesh[1]);//extrapolation to zero
            // Energy derivative of Psi
            for (int i = 0; i < Rmesh.size(); i++) inhom[i] = -2. * ur[i];
            urp[0] = 0;
            urp[1] = startSol(Z, l, dh);
            NumerovGen(rhs_MT, inhom, urp.size(), dh, urp);
            for (int i = 0; i < ur.size(); i++) temp[i] = ur[i] * urp[i];
            double alpha = integrate4<double>(temp, dh, temp.size());
            for (int i = 0; i < ur.size(); i++) urp[i] -= alpha * ur[i];

            Psip_l[l] = urp;      // Store it
            for (int ir = 1; ir < Rmesh.size(); ir++) Psip_l[l][ir] /= Rmesh[ir];// Psip=up/r!
            Psip_l[l][0] = extrapolate_f0(Psip_l[l][1], Psip_l[l][2], Rmesh[0], Rmesh[1], Rmesh[2]);
            //      Psip_l[l][0] = (Psip_l[l][1]*(Rmesh[2]-Rmesh[0])-Psip_l[l][2]*(Rmesh[1]-Rmesh[0]))/(Rmesh[2]-Rmesh[1]);// extrapolating

            for (int i = 0; i < urp.size(); i++) temp[i] = urp[i] * urp[i];
            PsipPsip[l] = integrate4<double>(temp, dh, temp.size()); // (Psip,Psip)

            // Here we estimate the derivatives at the Muffin-Tin boundary
            int N0 = ur.size() - 1;
            double RMuffinTin = Rmesh[N0];
            double v1 = rhs_MT[N0] * ur[N0];
            double v0 = rhs_MT[N0 - 1] * ur[N0 - 1];
            double w1 = rhs_MT[N0] * urp[N0] + inhom[N0];
            double w0 = rhs_MT[N0 - 1] * urp[N0 - 1] + inhom[N0 - 1];
            double dudr = (ur[N0] - ur[N0 - 1]) / dh + 0.125 * dh * (3 * v1 + v0);
            double dupdr = (urp[N0] - urp[N0 - 1]) / dh + 0.125 * dh * (3 * w1 + w0);
            dlogPsi[l] = RMuffinTin * dudr / ur[N0] - 1;
            dlogPsip[l] = RMuffinTin * dupdr / urp[N0] - 1;
            Psi[l] = ur[N0] / RMuffinTin;
            Psip[l] = urp[N0] / RMuffinTin;
        }
        return 0;
    }
};

// This is small class which is used as a functor for finding
// chemical potential. It can be used in connection with root finding routine
class FChemicalPotential {
    int Zval;
    const function2D<double>& epsk;
    const vector<double>& wk;
    double broad;
public:
    FChemicalPotential(int Zval_, const function2D<double>& epsk_, const vector<double>& wk_) :
        Zval(Zval_), epsk(epsk_), wk(wk_), broad(10 * pow(epsk.size_N(), 1 / 3.)) {}
    double operator()(double mu) const
    {
        double nn = 0;
        for (int ik = 0; ik < epsk.size_N(); ik++)
            for (int p = 0; p < epsk.size_Nd(); p++)
                //	if (epsk(ik,p)<=mu) nn += wk[ik];
                nn += wk[ik] * 0.5 * (1 + erf((mu - epsk(ik, p)) * broad));
        //    clog<<mu<<" "<<2*nn-Zval<<endl;
        return 2 * nn - Zval; // because of spin, factor of 2!
    }
};

void SolvePoisson(int Zq, const vector<double>& Rmesh, const function1D<double>& rho, vector<double>& Uhartree)
{// Given the input density rho, calculates the Hartree potential
  // The boundary conditions used are U(0)=0 and U(S)=Zq. The boundary condition at S is only a constant shift
  // of potential which is later readjusted by choosing MT zero. So, this boundary condition is not very relevant.
    static vector<double> RHS(Rmesh.size());
    for (int i = 0; i < Rmesh.size(); i++) RHS[i] = -4 * M_PI * Rmesh[i] * rho[i];
    Uhartree[0] = 0;  Uhartree[1] = (Rmesh[1] - Rmesh[0]);// Boundary condition for U_H=V_H/r
    NumerovInhom(RHS, RHS.size(), Rmesh[1] - Rmesh[0], Uhartree); // Solving the 2nd order differential equation
    // adding homogeneous solution to satisfay boundary conditions: U(0)=0, U(infinity)=Z
    int ilast = Uhartree.size() - 1;
    double U_last = Uhartree[ilast];
    double alpha = (Zq - U_last) / Rmesh[ilast];
    for (int i = 0; i < Rmesh.size(); i++) Uhartree[i] += alpha * Rmesh[i];
}

class LAPW {
    int Ksize, ksize, Rsize, lMax;
    double Vol, RMuffinTin;
    const vector<dvector3>& Km;
    function2D<double> omegal;             // constant takes care of continuity of the basis functions
    function2D<double>  C1;                // used for constracting C2
    vector<function2D<double> > C2l;       // MT constant C2
    function2D<double> C2_1, C2_2;         // C2 modification needed to calculate MT-charge
    vector<vector<vector<double> > > weigh0, weigh1, weigh2; // Weights for calculation of electronic charge
    vector<vector<double> > weighI;        // Weights for calculating carge in the interstitial region
    function2D<double> Olap_I;             // Overlap in the interstitials
    function2D<double> Olap, Ham;          // Total overlap and Hamiltonian
    function2D<double> temp0, temp1;       // and some temporary storage
public:
    LAPW(int Ksize_, int ksize_, int Rsize_, int lMax_, double Vol_, double RMuffinTin_, const vector<dvector3>& Km_) :
        Ksize(Ksize_), ksize(ksize_), Rsize(Rsize_), lMax(lMax_), Vol(Vol_), RMuffinTin(RMuffinTin_), Km(Km_),
        omegal(Ksize, lMax + 1), C1(Ksize, lMax + 1), C2l(lMax + 1), C2_1(Ksize, Ksize), C2_2(Ksize, Ksize),
        weigh0(ksize), weigh1(ksize), weigh2(ksize), weighI(ksize),
        Olap_I(Ksize, Ksize), Olap(Ksize, Ksize), Ham(Ksize, Ksize), temp0(Ksize, Ksize), temp1(Ksize, Ksize)

    {
        for (int l = 0; l <= lMax; l++) C2l[l].resize(Ksize, Ksize);
        for (int ik = 0; ik < ksize; ik++) {
            weighI[ik].resize(Ksize);
            weigh0[ik].resize(lMax + 1);
            weigh1[ik].resize(lMax + 1);
            weigh2[ik].resize(lMax + 1);
            for (int il = 0; il <= lMax; il++) {
                weigh0[ik][il].resize(Ksize);
                weigh1[ik][il].resize(Ksize);
                weigh2[ik][il].resize(Ksize);
            }
        }
    }
    void ComputeInterstitialOverlap()
    { // Overlap in the interstitials can be calculated outside the k-loop
        for (int i = 0; i < Ksize; i++) {
            Olap_I(i, i) = 1 - 4 * M_PI * sqr(RMuffinTin) * RMuffinTin / (3. * Vol);
            for (int j = i + 1; j < Ksize; j++) {
                double KKl = (Km[i] - Km[j]).length();
                Olap_I(i, j) = -4 * M_PI * sqr(RMuffinTin) * bessel_j(1, KKl * RMuffinTin) / (KKl * Vol);
                Olap_I(j, i) = Olap_I(i, j);
            }
        }
    }
    //void ComputeEigensystem(const dvector3& k, const PartialWave& wave, const vector<double>& Enu, double VKSi, function<double>& Energy)
    void ComputeEigensystem(const dvector3& k, const PartialWave& wave, const vector<double>& Enu, double VKSi, Function<double>& Energy)
    {// Hamiltonian and Overlap for the valence states is calculated
      // Bessel functions can be calculated only ones for each K-point
        for (int iK = 0; iK < Ksize; iK++)
            for (int il = 0; il <= lMax; il++) {
                double Dl, jl, jlDl;
                dlog_bessel_j(il, (k + Km[iK]).length() * RMuffinTin, Dl, jl, jlDl);
                omegal(iK, il) = -wave.psi(il) / wave.psip(il) * (Dl - wave.dlog(il)) / (Dl - wave.dlogp(il));
                //C1(iK,il) = sqrt(4*M_PI*(2*il+1)/Vol)*jl/(wave.psi(il)+omegal(iK,il)*wave.psip(il)); // This is less stable, but equivalent
                C1(iK, il) = sqrt(4 * M_PI * (2 * il + 1) / Vol) * (jlDl - jl * wave.dlogp(il)) / (wave.psi(il) * (-wave.dlogp(il) + wave.dlog(il)));
            }
        // Parts of the Hamiltonian matrix which do not depend on energy, are calculated
        // This part of the code needs most of the time. Should be very optimized
        for (int iK = 0; iK < Ksize; iK++) {
            for (int jK = 0; jK < Ksize; jK++) {
                dvector3 qi(k + Km[iK]);
                dvector3 qj(k + Km[jK]);

                double qi_len = qi.length();
                double qj_len = qj.length();
                double argv = (qi_len * qj_len == 0) ? 1. : qi * qj / (qi_len * qj_len);

                double olapMT = 0, hamMT = 0;
                for (int il = 0; il <= lMax; il++) {
                    double tC2l = C1(iK, il) * C1(jK, il) * Legendre(il, argv);
                    double toop = (1. + omegal(iK, il) * omegal(jK, il) * wave.PP(il));
                    olapMT += tC2l * toop;
                    hamMT += tC2l * (0.5 * (omegal(iK, il) + omegal(jK, il)) + toop * Enu[il]);
                    C2l[il](iK, jK) = tC2l;
                }
                Olap(iK, jK) = olapMT + Olap_I(iK, jK);
                Ham(iK, jK) = (0.25 * (qi * qi + qj * qj) + VKSi) * Olap_I(iK, jK) + hamMT;
            }
        }
        Eigensystem(Ksize, Energy, Olap, Ham); // When many K-points are used, this lapack call takes most of the time

    }
    void ComputeWeightsForDensity(int ik, const dvector3& k)
    { /// Calculation of valence density
        for (int il = 0; il <= lMax; il++) {
            for (int iK = 0; iK < Ksize; iK++)
                for (int jK = 0; jK < Ksize; jK++) {
                    C2_1(iK, jK) = C2l[il](iK, jK) * (omegal(iK, il) + omegal(jK, il));
                    C2_2(iK, jK) = C2l[il](iK, jK) * omegal(iK, il) * omegal(jK, il);
                }
            //temp0.Product("N", "N", Ham, C2l[il]);
            temp0.Product("N", "N", Ham, C2l[il], 1.0, 0.0);
            //temp1.Product("N", "T", temp0, Ham);
            temp1.Product("N", "T", temp0, Ham, 1.0, 0.0);
            for (int p = 0; p < Ksize; p++)	weigh0[ik][il][p] = temp1(p, p);

            //temp0.Product("N", "N", Ham, C2_1);
            temp0.Product("N", "N", Ham, C2_1, 1.0, 0.0);
            //temp1.Product("N", "T", temp0, Ham);
            temp1.Product("N", "T", temp0, Ham, 1.0, 0.0);
            for (int p = 0; p < Ksize; p++)	weigh1[ik][il][p] = temp1(p, p);

            //temp0.Product("N", "N", Ham, C2_2);
            temp0.Product("N", "N", Ham, C2_2, 1.0, 0.0);
            //temp1.Product("N", "T", temp0, Ham);
            temp1.Product("N", "T", temp0, Ham, 1.0, 0.0);
            for (int p = 0; p < Ksize; p++)	weigh2[ik][il][p] = temp1(p, p);
        }
        //temp0.Product("N", "N", Ham, Olap_I);
        temp0.Product("N", "N", Ham, Olap_I, 1.0, 0.0);
        //temp1.Product("N", "T", temp0, Ham);
        temp1.Product("N", "T", temp0, Ham, 1.0, 0.0);
        for (int p = 0; p < Ksize; p++)	weighI[ik][p] = temp1(p, p);
    }

    void ComputeMTDensity(function1D<double>& nMTRho, const function2D<double>& Energy, double mu, const vector<double>& wk, const PartialWave& wave)
    {
        nMTRho = 0;
        for (int il = 0; il <= lMax; il++) {
            double w0 = 0, w1 = 0, w2 = 0;
            for (int ik = 0; ik < ksize; ik++) {
                double sum0 = 0, sum1 = 0, sum2 = 0;
                for (int p = 0; p < Ksize; p++) {
                    if (Energy(ik, p) <= mu) {
                        sum0 += weigh0[ik][il][p];
                        sum1 += weigh1[ik][il][p];
                        sum2 += weigh2[ik][il][p];
                    }
                }
                w0 += sum0 * wk[ik];
                w1 += sum1 * wk[ik];
                w2 += sum2 * wk[ik];
            }

            //       cout.setf(ios::fixed,ios::floatfield);   // floatfield set to fixed
            //       cout.precision(10);
            //       cout<<setw(3)<<il<<" "<<setw(20)<<w0<<" "<<setw(20)<<w1<<" "<<setw(20)<<w2<<endl;
            for (int ir = 0; ir < Rsize; ir++)
                nMTRho[ir] += (w0 * sqr(wave.psi(il, ir)) + w1 * wave.psi(il, ir) * wave.psip(il, ir) + w2 * sqr(wave.psip(il, ir))) / (4 * M_PI);
        }
        for (int ir = 0; ir < Rsize; ir++) nMTRho[ir] *= 2; // due to spin
    }
    double ComputeInterstitialCharge(const function2D<double>& Energy, double mu, const vector<double>& wk)
    {
        double sIntRho = 0; // total interstitial charge
        for (int ik = 0; ik < ksize; ik++)
            for (int p = 0; p < Ksize; p++)
                if (Energy(ik, p) <= mu) sIntRho += weighI[ik][p] * wk[ik];
        sIntRho *= 2;// due to spin
        return sIntRho;
    }
    //void PrintBandStructure(int ik, const function<double>& Energy, ostream& out)
    void PrintBandStructure(int ik, const Function<double>& Energy, ostream& out)
    {
        out << setw(10) << ik / (ksize - 1.) << " ";
        for (int iK = 0; iK < Ksize; iK++) out << setw(12) << Energy[iK] << " ";
        out << endl;
    }
};

double IntegrateCharge(const vector<double>& R, const function1D<double>& rho)
{
    static function1D<double> temp(R.size());
    for (int i = 0; i < R.size(); i++) temp[i] = rho[i] * sqr(R[i]) * 4 * M_PI;
    return integrate4<double>(temp, R[1] - R[0], temp.size());
}

int main(int argc, char* argv[], char* env[])
{
    int Z = 29;                     // Number of electrons in the atom
    double LatConst = 6.8219117;  // Lattice constant
    int lMax = 5;                   // Maximum l considered in calculation
    int N = 1001;                 // Number of points in radial mesh
    int nkp = 8;                  // Number of k-points in 1BZ: (nkp x nkp x nkp)
    double CutOffK = 3.5;           // Largest length of reciprocal vectors K (only shorter vec. are taken into account)
    double dEz = 0.1;             // Step in serching for core states
    double mu_min = 0.0, mu_max = 1.5;// Interval where chemical potential is looking for
    double Mix = 0.2;             // Linear mixing parameter for charge
    double Vmix = 0.2;            // Linear mixing parameter for potential
    double precision = 1e-5;        // accuracy of energy
    int nitt = 200;               // maximum number of iterations
    bool read = false;              // weather to read input potential
    ///// Core states/////////////
    int lMaxCore = 2;
    vector<int> core(lMaxCore + 1);
    core[0] = 3; // 1s-3s is in core
    core[1] = 2; // 1p,2p is in core
    core[2] = 0; // no d in core

    int i = 0;
    while (++i < argc) {
        std::string str(argv[i]);
        if (str == "-Z" && i < argc - 1) Z = atoi(argv[++i]);
        if (str == "-dE" && i < argc - 1) dEz = atof(argv[++i]);
        if (str == "-lmax" && i < argc - 1) lMax = atoi(argv[++i]);
        if (str == "-N" && i < argc - 1) N = atoi(argv[++i]);
        if (str == "-nkp" && i < argc - 1) nkp = atoi(argv[++i]);
        if (str == "-CutOffK" && i < argc - 1) CutOffK = atof(argv[++i]);
        if (str == "-mumin" && i < argc - 1) mu_min = atof(argv[++i]);
        if (str == "-mumax" && i < argc - 1) mu_max = atof(argv[++i]);
        if (str == "-Mix" && i < argc - 1) Mix = atof(argv[++i]);
        if (str == "-Vmix" && i < argc - 1) Vmix = atof(argv[++i]);
        if (str == "-precision" && i < argc - 1) precision = atof(argv[++i]);
        if (str == "-nitt" && i < argc - 1) nitt = atoi(argv[++i]);
        if (str == "-read") read = true;
        if (str == "-h" || str == "--help") {
            std::clog << "**************** LAPW program for fcc lattice *********\n";
            std::clog << "**                                                  **\n";
            std::clog << "******************************************************\n";
            std::clog << "\n";
            std::clog << "lapw [-dE double] [] []\n";
            std::clog << "Options:   -Z          Number of electrons (" << Z << ")\n";
            std::clog << "           -dE         Step in searching for states (" << dEz << ")\n";
            std::clog << "           -lmax       Maximum l used in calculation (" << lMax << ")\n";
            std::clog << "           -N          Number of points in radial mesh (" << N << ")\n";
            std::clog << "           -nkp        Number of k-points from IRBZ used (nkp x nkp x nkp) (" << nkp << ")\n";
            std::clog << "           -CutOffK    Largest length of reciprocal vectors K (only shorter vec. are taken into account) (" << CutOffK << ")\n";
            std::clog << "           -mumin      Start looking for chemical potential (" << mu_min << ")\n";
            std::clog << "           -mumax      Stop looking for chemical potential (" << mu_max << ")\n";
            std::clog << "           -Mix        Linear mixing parameter for charge (" << Mix << ")\n";
            std::clog << "           -Vmix       Linear mixing parameter for potential (" << Vmix << ")\n";
            std::clog << "           -precision  Total energy accuracy required (" << precision << ")\n";
            std::clog << "           -nitt       Maximum number of iterations (" << nitt << ")\n";
            std::clog << "           -read       Weather to read from file Potential_input.dat (" << read << ")\n";
            std::clog << "*****************************************************\n";
            return 0;
        }
    }
    clog.precision(10);

    /////////////////////////////
    int Zcor = 0; // Core number of electrons
    for (int l = 0; l < core.size(); l++) Zcor += 2 * (2 * l + 1) * core[l];
    int Zval = Z - Zcor; // Valence number of electrons
    clog << "Z core = " << Zcor << " and Zval = " << Zval << endl;
    ///////////////////////////////
    // Generates and stores momentum points
    FccLattice fcc(LatConst);                  // Information about lattice
    double RMuffinTin = fcc.RMuffinTin();      // Muffin-Tin radius choosen such that spheres touch
    double VMT = 4 * M_PI * pow(RMuffinTin, 3) / 3.;  // Volume of MT
    double Vinter = fcc.Vol() - VMT;             // Volume of the interstitial region
    clog << "Muffin-Tin radius = " << RMuffinTin << endl;
    clog << "Volume of the MT sphere    = " << VMT << endl;
    clog << "Volume of the interstitial = " << Vinter << endl;
    //////////////////////////////
    // For solving SCH equation
    PartialWave wave(N, RMuffinTin, Z, lMax); // This is for partil waves in the valence band
    RadialWave wavec(wave.R());               // This is for core states - to solve SCH equation

    fcc.GenerateReciprocalVectors(4, CutOffK);// Reciprocal bravais lattice is builded, K points taken into account only for |K|<CutOff
    fcc.ChoosePointsInFBZ(nkp, 0);             // Chooses the path in the 1BZ or the k-points in the irreducible 1BZ
    ExchangeCorrelation XC(3);                // Exchange correlations class; VWN seems to be the best (look http://physics.nist.gov/PhysRefData/DFTdata/Tables/ptable.html)

    vector<double> Enu(lMax + 1);             // Linearization energies. Should be in the center of the occupied band
    Enu[0] = 0.11682;                       // Most of high energy partial waves should be centered around mu
    Enu[1] = 0.18794;                       // In general, it is a good idea to change them through iteration
    Enu[2] = 0.211145;
    for (int il = 3; il < Enu.size(); il++) Enu[il] = 0.3;
    double VKSi = 0;                       // Potential in the interstitials
    double mu = 0;                         // Chemical potential

    vector<double> Uhartree(wave.Rsize()), Vxc(wave.Rsize()), Veff(wave.Rsize()); // Hartree and exchange-correlation and KS potential
    function1D<double> TotRho(wave.Rsize()), nTotRho(wave.Rsize()), drho(wave.Rsize());// total density, input and output
    function2D<double> Energy(fcc.ksize(), fcc.Ksize());                               // E(k,p)- bands

    LAPW lapw(fcc.Ksize(), fcc.ksize(), wave.Rsize(), lMax, fcc.Vol(), RMuffinTin, fcc.Km()); //basic clas for LAPW calculation

    function1D<double>  MTRho(wave.Rsize()), coreRho(wave.Rsize()); // partial densities, core and valence-MT 

    // Starting guess for the Hartree and exchange correlation potential
    for (int i = 1; i < wave.Rsize(); i++) Veff[i] = -VeffP(wave.R(i)) / wave.R(i);
    Veff[0] = extrapolate_f0(Veff[1], Veff[2], wave.R(0), wave.R(1), wave.R(2));
    // Reads potential from file in read==true
    if (read) {
        ifstream inp("Potential_input.txt");
        double r, v, v0;
        int ii = 0;
        while (inp >> r && inp >> v && inp >> v0 && ii < wave.Rsize()) Veff[ii++] = -v / r;
        Veff[0] = extrapolate_f0(Veff[1], Veff[2], wave.R(0), wave.R(1), wave.R(2));
    }

    double zeroMT = Veff[wave.Rsize() - 1]; // adjusting MT zero
    for (int i = 0; i < wave.Rsize(); i++) Veff[i] -= zeroMT;



    double pEc = 0; // previous core energy - to calculate difference

    /////////////////// Starting calculation ////////////////////////////////
    lapw.ComputeInterstitialOverlap(); // Overlap in the interstitials can be calculated outside

    //////////////////// Main SC - iteration loop ///////////////////////////
    for (int itt = 0; itt < nitt; itt++) {
        clog << endl;
        clog << "****** Iteration number " << itt << " ************" << endl;
        /////////////////// Potential part /////////////////////////////
        if (itt > 0) {// We start with input potential rather than density. Calculation of potential skiped first time.
            SolvePoisson(Z, wave.R(), TotRho, Uhartree); // Poisson gives Hartree potential
            for (int i = 0; i < wave.Rsize(); i++) Vxc[i] = XC.Vxc(rs(TotRho[i])); // Adding exchange-correlation part
            // This is total KS effective potential
            for (int i = 1; i < wave.Rsize(); i++) Veff[i] = (1 - Vmix) * Veff[i] + Vmix * ((-Z + Uhartree[i]) / wave.R(i) + Vxc[i]);
            Veff[0] = extrapolate_f0(Veff[1], Veff[2], wave.R(0), wave.R(1), wave.R(2));
            // Zero of energy is choosen at each iteration such that potential vanishes on MT sphere
            double VMTzero = Veff[wave.Rsize() - 1]; // New MT zero
            for (int i = 0; i < wave.Rsize(); i++) Veff[i] -= VMTzero; // Shift potential to new MT zero
        }
        ///////////// Schroedinger equation for MT region //////////////////
        wave.SetVeff0(Veff);

        wave.SolveSCHEquation(Enu);// Energy Enu depends on l

        //////////////// Core part ////////////////////////////////////////
        wavec.SetVeff0(Veff);
        int Nc; double Ec;
        FindCoreStates(core, Z, dEz, wavec, Nc, Ec, coreRho);
        clog << "Core Z = " << Nc << " and Energy = " << Ec << endl;

        ///////////// Main LAPW loop over k points /////////////////////////
        for (int ik = 0; ik < fcc.ksize(); ik++) {
            dvector3 k = fcc.k(ik); // current k-point
            lapw.ComputeEigensystem(k, wave, Enu, VKSi, Energy[ik]);
            lapw.ComputeWeightsForDensity(ik, k);
        }
        /////////////////// New chemical potential ////////////////////////////
        FChemicalPotential chemp(Zval, Energy, fcc.wk());
        mu = zeroin(mu_min, mu_max, chemp, 1e-13);
        clog << "Chemical potential found at " << COLOR(BLUE, mu) << endl;
        //////////////////// New valence Density /////////////////////////////
        lapw.ComputeMTDensity(MTRho, Energy, mu, fcc.wk(), wave);

        double sIntRho = lapw.ComputeInterstitialCharge(Energy, mu, fcc.wk());

        double sMTRho = IntegrateCharge(wave.R(), MTRho);
        cout << "MTRho = " << sMTRho + sIntRho << endl;

        double scoreRho = IntegrateCharge(wave.R(), coreRho);
        ////////////////// New total charge /////////////////////////////////
        for (int i = 0; i < wave.Rsize(); i++) nTotRho[i] = MTRho[i] + coreRho[i];

        clog << "Weight in the MT sphere = " << sMTRho << " and in the interstitials = " << sIntRho << " and in core = " << scoreRho << endl;
        double renorm = Z / (sMTRho + sIntRho + scoreRho);
        clog << "Total charge found = " << scoreRho + sMTRho + sIntRho << ", should be " << Z << " -> renormalizing charge by " << renorm << endl;
        ///////////////// Renormalization of charge ////////////////////////
        for (int i = 0; i < wave.Rsize(); i++) nTotRho[i] *= renorm;
        //////////////////// Charge diference //////////////////////////////
        for (int i = 0; i < wave.Rsize(); i++) drho[i] = fabs(TotRho[i] - nTotRho[i]);
        double ChargeDifference = integrate4<double>(drho, wave.R(1) - wave.R(0), drho.size());
        if (itt == 0) Mix = 1; // Since we do not have charge in the first itteration, need to take new charge only
        //////////// Linear mixing. Could be improved with Broyden, Vanderbild or Johannson mixing /////////////
        for (int i = 0; i < wave.Rsize(); i++) TotRho[i] = TotRho[i] * (1 - Mix) + nTotRho[i] * Mix;
        ////////////// Convergence criteria ////////////////////////////////////////////////////
        clog << "Core energy difference = " << COLOR(YELLOW, fabs(Ec - pEc)) << ", Charge difference = " << COLOR(GREEN, ChargeDifference) << endl;
        if (fabs(Ec - pEc) < precision) break;
        pEc = Ec;
    }

    ////// Printing of band structure /////////////////////////
    ofstream bands("bands.txt");
    fcc.ChoosePointsInFBZ(nkp, 1);// Points for plotting band structure
    function1D<double> epsk(fcc.Ksize());
    if (nitt == 0) {// Schroedinger equation hasn't been solved yet
        wave.SetVeff0(Veff);
        wave.SolveSCHEquation(Enu);
    }
    for (int ik = 0; ik < fcc.ksize(); ik++) {
        dvector3 k = fcc.k(ik); // current k-point
        lapw.ComputeEigensystem(k, wave, Enu, VKSi, epsk);
        lapw.PrintBandStructure(ik, epsk, bands);
    }

    //////////////////// Some printing //////////////////////////////
    ofstream out("charge.txt");
    for (int ir = 0; ir < wave.Rsize(); ir++)
        out << wave.R(ir) << "  " << TotRho[ir] * 4 * M_PI * sqr(wave.R(ir)) << "  " << MTRho[ir] * 4 * M_PI * sqr(wave.R(ir)) << "  " << coreRho[ir] * 4 * M_PI * sqr(wave.R(ir)) << endl;

    ofstream outc("Charge.txt");
    for (int ir = 0; ir < wave.Rsize(); ir++)
        outc << wave.R(ir) << "  " << TotRho[ir] << "  " << MTRho[ir] << "  " << coreRho[ir] << endl;

    ofstream outu("potential.txt");
    for (int i = 0; i < wave.Rsize(); i++)
        outu << wave.R(i) << " " << Veff[i] << " " << Uhartree[i] << " " << Vxc[i] << endl;

    ofstream outp("Potential.txt");
    for (int i = 0; i < wave.Rsize(); i++)
        outp << wave.R(i) << " " << -Veff[i] * wave.R(i) << "  " << VeffP(wave.R(i)) << endl;

    return 0;
}
const double ExchangeCorrelation::alphax;
