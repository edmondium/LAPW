#ifndef _LBLAS_
#define _LBLAS_
using namespace std;
extern "C" {
    void dsygvd_(const int* ITYPE, const char* JOBZ, const char* UPLO, const int* N,
        double* A, const int* LDA, double* B, const int* LDB, double* W, double* WORK, const int* LWORK,
        int* IWORK, const int* LIWORK, int* INFO);
};
void dsygvd_(const int* ITYPE, const char* JOBZ, const char* UPLO, const int* N,
    double* A, const int* LDA, double* B, const int* LDB, double* W, double* WORK, const int* LWORK,
    int* IWORK, const int* LIWORK, int* INFO)
{

}
int Eigensystem(int N, function<double>& Energy, const function2D<double>& Olap, function2D<double>& F)
{// C++ wrapper function for general eigenvalue problem
    static function2D<double> tOlap(N, N);
    int lwork = 1 + 6 * N + 2 * N * N + 10;
    static function1D<double> work(lwork);
    int liwork = 3 + 5 * N + 10;
    static function1D<int> iwork(liwork);
    int itype = 1;
    int info = 0;
    tOlap = Olap;
    int lF = F.size_Nd();
    int lO = Olap.size_Nd();
    dsygvd_(&itype, "V", "U", &N, F.MemPt(), &lF, tOlap.MemPt(), &lO, Energy.MemPt(), work.MemPt(), &lwork, iwork.MemPt(), &liwork, &info);
    if (info) cerr << "Not sucessfull solving eigenvalue problem with dsygvd " << info << endl;
    return info;
}

#endif //_LBLAS_