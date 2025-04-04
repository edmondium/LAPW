#ifndef _CORE_
#define _CORE_

class RadialWave {
public:
    int N;
    vector<double> R;        // Radial mesh but reversed R[0]=Rmax and R.last=0
    vector<double> Solution; // Non-normalized solution of Schroedinger equation
    vector<double> rhs;      // Right-hand-site in solving Schroedinger equation
    vector<double> Veff;     // Effective KS potential
    vector<double> Veff0;    // Effective KS potential without centrifugal part
public:
    RadialWave(const vector<double>& Rmesh) : N(Rmesh.size()), Solution(N), R(N), rhs(N), Veff(N), Veff0(N)
    {
        for (int i = 0; i < N; i++) R[i] = Rmesh[N - 1 - i]; // The mesh is reversed
        Solution[0] = R[0] * exp(-R[0]); // Boundary (starting) points of integration by Numerov
        Solution[1] = R[1] * exp(-R[1]); // Boundary (starting) points of integration by Numerov
    }
    // This function-operator is used to find bound states. It is given to root-finding routine
    double operator()(double E)
    {
        double h = R[1] - R[0];
        for (int i = 0; i < N - 1; i++) rhs[i] = 2 * (Veff[i] - E); // The RHS of the SCHR-equation for choosen energy E
        rhs[R.size() - 1] = 0;                                // This is the zero frequency
        Numerov(rhs, Solution.size() - 1, h, Solution);     // Solving the radial SCH-equation
        int last = Solution.size() - 1;                     // The last point at zero frequency needs extrapolation
        Solution[last] = Solution[last - 1] * (2 + h * h * rhs[last - 1]) - Solution[last - 2];
        return Solution[last];                            // Value at zero frequency
    }
    // This function return the density of electrons (up+down per volume) of one (nl) state.
    void Density(vector<double>& rho)
    {
        rho.resize(Solution.size());
        int N = Solution.size();
        for (int i = 0; i < Solution.size(); i++) rho[i] = Solution[N - 1 - i] * Solution[N - 1 - i]; // The mesh outside this class is reversed!
        double norm = 1. / integrate4<double>(rho, R[0] - R[1], rho.size());                // Normalization constant 
        for (int i = 1; i < rho.size(); i++) rho[i] = rho[i] * norm / (4 * M_PI * sqr(R[N - i - 1]));   // rho_{nl}=u^2/(4*Pi*r^2)
        rho[0] = 2 * rho[1] - rho[2];                                                       // extrapolation to zero frequency
    }
    // This function sets KS potential without centrifugal part
    void SetVeff0(const vector<double>& Veff)
    {
        for (int i = 0; i < R.size(); i++) Veff0[i] = Veff[N - 1 - i];
    }
    void AddCentrifugal(int l)
    {
        for (int i = 0; i < R.size() - 1; i++) Veff[i] = Veff0[i] + 0.5 * l * (l + 1) / sqr(R[i]);
    }
    double V_KS0(int i) { return Veff0[N - 1 - i]; }
};

void FindCoreStates(const vector<int>& core, int Z, double dEz, RadialWave& wave, int& Nc, double& Ec, function1D<double>& coreRho)
{// Searches for bound state with given n,l. They are stored in vector<BState>

//   int l0=0;
//   wave.AddCentrifugal(l0);
//   double E = -310.0351139;
//   double dw = wave(E);

//   cout.setf(ios::fixed,ios::floatfield);   // floatfield set to fixed
//   cout.precision(10);
//   for (int i=0; i<wave.N; i++)
//     cout<<setw(20)<<wave.R[i]<<" "<<setw(20)<<wave.rhs[i]<<" "<<setw(20)<<wave.Solution[i]<<endl;

//   exit(1);



    static vector<double> drho(coreRho.size());
    for (int ir = 0; ir < coreRho.size(); ir++) coreRho[ir] = 0;
    Nc = 0;
    Ec = 0;
    for (int l = 0; l < core.size(); l++) {
        wave.AddCentrifugal(l);
        double x = -0.5 * Z * Z / sqr(l + 1) - 3.;// here is starts to look for zero
        double v0 = wave(x), v1 = v0;
        int n = 0, j = 0;
        while (n < core[l] && x < 10.) {
            x += dEz;                                     // Proceeding in small steps to bracket all zeros
            v1 = wave(x);                               // New value of radial function at origin
            if (v0 * v1 < 0) {                              // Change sign?
                double Energy = zeroin(x - dEz, x, wave, 1e-10); // Root-finder locates bound state very precisely
                //cout<<"Debugcore "<<l<<" "<<Energy<<endl;
                int dN = 2 * (2 * l + 1);  // degeneracy of each radial wave level
                wave.Density(drho);
                for (int ir = 0; ir < coreRho.size(); ir++) coreRho[ir] += drho[ir] * 2 * (2 * l + 1);
                Ec += 2 * (2 * l + 1) * Energy; // Sum of eigenvalues times degeneracy
                Nc += dN;
                clog << "Found core state for n = " << n + l << ", l = " << l << " at " << Energy << endl;
                n++;
                v0 = v1;
            }
        }
    }
}

#endif //_CORE_