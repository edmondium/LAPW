#ifndef FUNCTION_
#define FUNCTION_
#include <iostream>
#include <algorithm>
#include <string>
#include <memory.h>
#include "assert.h"
//////////////////////////// Hierarchy of classes: ///////////////////////////////////////
//                             
//                                     base clas  = function<>
//                   |                                                |
//               function1D<>                                     funProxy<>
//
//                                       function2D<funProxy>
//
template<class T> class Function;
//template<class T> class function;
template<class T> class function1D;
template<class T> class funProxy;
template<class T> class function2D;

//********************************************************************************//
// Base class for two derived classes: function1D<T> and function2D<T>.	  	  //
// It is also used as a proxy class for function2D. Function2D<T> consists of	  //
// arrays of function<T> rather than functions1D<T>.				  //
// Memory is allocated in a fortran-like fashion for better performance.	  //
// Linear interpolation is implemented with the operator() that takes one	  //
// argument (class intpar).							  //
//********************************************************************************//
template<class T>
//class function {
class Function {
protected:
    T* f;
    int N0, N;
public:
    T& operator[](int i) { Assert(i < N, "Out of range in function[]"); return f[i]; }
    const T& operator[](int i) const { Assert(i < N, "Out of range in function[]"); return f[i]; }
    const T& last() const { return f[N - 1]; }
    int size() const { return N; }
    int fullsize() const { return N0; }
    Function& operator+=(const Function& m);
    //function& operator+=(const function& m);
    //function& operator*=(const T& m);
    Function& operator*=(const T& m);
    T* MemPt() { return f; }
    const T* MemPt() const { return f; }
    Function& operator=(const T& c);
    //function& operator=(const T& c);
protected:
    //function() : f(NULL), N0(0), N(0) {};
    //explicit function(int N_) : N0(N_), N(N_) {};
    //~function() {};
    //function(const function&) {};
    Function() : f(NULL), N0(0), N(0) {};
    explicit Function(int N_) : N0(N_), N(N_) {};
    ~Function() {};
    Function(const Function&) {};
    template<class U> friend class function2D;
    //template <class U> friend U scalar_product(const function<U>& f1, const function<U>& f2);
};

//******************************************************************//
// One dimensional functions derived from function<T>. It has it's  //
// own constructors and destructors.				    //
//******************************************************************//
template <class T>
//class function1D : public function<T> {
class function1D : public Function<T> {
public:
    function1D() {};
    explicit function1D(int N_);
    ~function1D();
    function1D(const function1D& m);
    void resize(int N_);
    function1D& operator=(const function1D& m);
    function1D& operator=(const T& c) { Function<T>::operator=(c); return *this; }
    //function1D& operator=(const T& c) { function<T>::operator=(c); return *this; }
    //void Product(const function2D<T>& A, const function<T>& x, double alpha = 1., double beta = 0.);
};

template <class T>
//class funProxy : public function<T> {
class funProxy : public Function<T> {
public:
    void Initialize(int N_, T* f_);
    void ReInitialize(int N_, T* f_);
    void resize(int N_);
    funProxy& operator=(const function<T>& m);
    ~funProxy() {};
};

//**********************************************************************//
// Two dimentional function<T> derived from function<T>. It consists	//
// of an array of function<T> rather tham function1D<T>.		//
// Constructor calls operator new and aferwords placement new operator	//
// to allocate the whole memory in one single large peace. 		//
//**********************************************************************//
template<class T>
class function2D {
protected:
    void* memory;
    T* data;
    funProxy<T>* f;
    int N0, Nd0, N, Nd;
public:
    function2D() : memory(NULL), N0(0), Nd0(0), N(0), Nd(0) {};
    function2D(int N_, int Nd_);
    ~function2D();
    funProxy<T>& operator[](int i) { Assert(i < N, "Out of range in function2D[]"); return f[i]; }
    const funProxy<T>& operator[](int i) const { Assert(i < N, "Out of range in function2D[]"); return f[i]; }
    const T& operator()(int i, int j) const { Assert(i < N&& j < Nd, "Out of range in function2D(i,j)"); return f[i].f[j]; }
    T& operator()(int i, int j) { Assert(i < N&& j < Nd, "Out of range in function2D(i,j)"); return f[i].f[j]; }
    T* MemPt() { return data; }
    const T* MemPt() const { return data; }
    const int size_N() const { return N; }
    const int size_Nd() const { return Nd; }
    const int fullsize_N() const { return N0; }
    const int fullsize_Nd() const { return Nd0; }
    const int lda() const { return Nd0; }
    void resize(int N_, int Nd_);

    function2D& operator=(const function2D& m);
    function2D& operator+=(double x);
    function2D& operator+=(const function2D& m);
    function2D& operator-=(double x);
    function2D& operator-=(const function2D& m);
    function2D& operator=(const T& u);
    function2D& operator*=(const T& x);
    
    void Product(const std::string& transa, const std::string& transb, const function2D& A, const function2D& B, const T& alpha, const T& beta);
	void MProduct(const function2D& A, const function2D& B, const T& alpha, const T& beta);
};

// function ////////////////////////////////////////////////////////////////
template<class T>
//inline function<T>& function<T>::operator+=(const function& m)
inline Function<T>& Function<T>::operator+=(const Function& m)
{
    T_LOG(if (N != m.size()) cerr << "Functions not of equal length! Can't sum!" << std::endl;)
        for (int i = 0; i < N; i++) f[i] += m[i];
    return *this;
}

template<class T>
//inline function<T>& function<T>::operator*=(const T& m)
inline Function<T>& Function<T>::operator*=(const T& m)
{
    for (int i = 0; i < N; i++) f[i] *= m;
    return *this;
}

template <class T>
//inline function<T>& function<T>::operator=(const T& c)
inline Function<T>& Function<T>::operator=(const T& c)
{
    T_LOG(if (N <= 0) cerr << "Size of function is non positive! " << N << std::endl;)
        for (int i = 0; i < N; i++) f[i] = c;
    return *this;
}

// function1D ////////////////////////////////////////////////////////////
template<class T>
//inline function1D<T>::function1D(int N_) : function<T>(N_)
inline function1D<T>::function1D(int N_) : Function<T>(N_)
{
    this->f = new T[N_];
}

template<class T>
inline function1D<T>::~function1D()
{
    delete[] this->f;
    this->f = NULL;
}

template<class T>
inline void function1D<T>::resize(int n)
{
    if (n > this->N0) {
        if (this->f) delete[] this->f;
        this->f = new T[n];
        this->N0 = n;
    }
    this->N = n;
}

template<class T>
inline function1D<T>::function1D(const function1D& m)
{
    resize(m.N);
    std::copy(m.f, m.f + this->N, this->f);
}

template <class T>
inline function1D<T>& function1D<T>::operator=(const function1D<T>& m)
{
    resize(m.N);
    std::copy(m.f, m.f + this->N, this->f);
    return *this;
}

// funProxy ///////////////////////////////////////////////////////////////
template <class T>
inline void funProxy<T>::Initialize(int N_, T* f_)
{
    this->N = this->N0 = N_; this->f = f_;
}

template <class T>
inline void funProxy<T>::ReInitialize(int N_, T* f_)
{
    this->N = N_; this->f = f_;
}

template <class T>
inline void funProxy<T>::resize(int N_)
{
    if (N_ > this->N0) std::cerr << "Can't resize funProxy, to small funProxy!" << std::endl;
    else this->N = N_;
}

template <class T>
inline funProxy<T>& funProxy<T>::operator=(const function<T>& m)
{
    resize(m.size());
    std::copy(m.MemPt(), m.MemPt() + this->N, this->f);
    return *this;
}

#define HPoffset 8
// function2D ////////////////////////////////////////////////////////////
template<class T>
function2D<T>::function2D(int N_, int Nd_) : N0(N_), Nd0(Nd_), N(N_), Nd(Nd_)
{
    memory = operator new (sizeof(funProxy<T>) * N0 + sizeof(T) * Nd0 * N0 + HPoffset);

    Assert(memory != NULL, "Out of memory");

    f = new (memory) funProxy<T>[N0];

    int offset = sizeof(funProxy<T>) * N0 + HPoffset;
    data = reinterpret_cast<T*>(static_cast<char*>(memory) + offset);

    for (int i = 0; i < N0; i++) f[i].Initialize(Nd0, data + i * Nd0);
}

template<class T>
function2D<T>::~function2D()
{
    for (int i = 0; i < N0; i++) {
        f[i].~funProxy<T>();
    }
    operator delete(memory);
    memory = NULL;
}

template <class T>
inline function2D<T>& function2D<T>::operator=(const function2D& m)
{
    if (m.N <= N0 && m.Nd <= Nd0) {
        N = m.N; Nd = m.Nd;
        for (int i = 0; i < N; i++) memcpy(f[i].f, m.f[i].f, sizeof(T) * Nd);
    }
    else {
        int msize = sizeof(funProxy<T>) * m.N + sizeof(T) * m.Nd * m.N + HPoffset;
        operator delete(memory);
        memory = operator new (msize);
        Assert(memory != NULL, "Out of memory");
        memcpy(memory, m.memory, msize);
        N = N0 = m.N; Nd = Nd0 = m.Nd;
        f = new (memory) funProxy<T>[N];
        int offset = sizeof(funProxy<T>) * N + HPoffset;
        data = reinterpret_cast<T*>(static_cast<char*>(memory) + offset);
        for (int i = 0; i < N; i++) f[i].Initialize(Nd, data + i * Nd);
    }
    return *this;
}

template <class T>
inline void function2D<T>::resize(int N_, int Nd_)
{
    if (N_ > N0 || Nd_ > Nd0) {
        //    clog<<"Deleting function2D and resizing from "<<N0<<" "<<Nd0<<" to "<<N_<<" "<<Nd_<<std::endl;
        int msize = sizeof(funProxy<T>) * N_ + sizeof(T) * Nd_ * N_ + HPoffset;
        operator delete(memory);
        memory = operator new (msize);
        Assert(memory != NULL, "Out of memory");
        N = N0 = N_; Nd = Nd0 = Nd_;
        f = new (memory) funProxy<T>[N];
        int offset = sizeof(funProxy<T>) * N + HPoffset;
        data = reinterpret_cast<T*>(static_cast<char*>(memory) + offset);
        for (int i = 0; i < N; i++) f[i].Initialize(Nd, data + i * Nd);
    }
    else {
        N = N_; Nd = Nd_;
    }
}

template <class T>
inline function2D<T>& function2D<T>::operator+=(double x)
{
    if (N != Nd || !Nd || !N) {
        std::cerr << "Can't add number to non-square matrix!" << std::endl;
        return *this;
    }
    for (int i = 0; i < Nd; i++) f[i][i] += x;
    return *this;
}

template <class T>
inline function2D<T>& function2D<T>::operator+=(const function2D& m)
{
    if (N != m.N || Nd != m.Nd) {
        std::cerr << "Can't sum different matrices!" << std::endl;
        return *this;
    }
    for (int i = 0; i < N; i++)
        for (int j = 0; j < Nd; j++)
            f[i][j] += m[i][j];

    return *this;
}

template <class T>
inline function2D<T>& function2D<T>::operator-=(double x)
{
    if (N != Nd || !N || !Nd) {
        std::cerr << "Can't add number to non-square matrix!" << std::endl;
        return *this;
    }
    for (int i = 0; i < Nd; i++) f[i][i] -= x;
    return *this;
}

template <class T>
inline function2D<T>& function2D<T>::operator-=(const function2D& m)
{
    if (N != m.N || Nd != m.Nd) {
        std::cerr << "Can't sum different matrices!" << std::endl;
        return *this;
    }
    for (int i = 0; i < N; i++)
        for (int j = 0; j < Nd; j++)
            f[i][j] -= m[i][j];

    return *this;
}
template <class T>
inline function2D<T>& function2D<T>::operator=(const T& u)
{
    for (int i = 0; i < N; i++) for (int j = 0; j < Nd; j++) f[i].f[j] = u;
    return *this;
}
template <class T>
inline function2D<T>& function2D<T>::operator*=(const T& x)
{
    for (int i = 0; i < N; i++) for (int j = 0; j < Nd; j++) f[i][j] *= x;
    return *this;
}

template <class T>
std::ostream& operator<<(std::ostream& stream, const function<T>& f)
{
    int width = stream.width();
    for (int i = 0; i < f.size(); i++) stream << i << " " << std::setw(width) << f[i] << std::endl;
    return stream;
}

template <class T>
std::ostream& operator<<(std::ostream& stream, const function2D<T>& f)
{
    int width = stream.width();
    for (int i = 0; i < f.size_N(); i++) {
        for (int j = 0; j < f.size_Nd(); j++)
            stream << std::setw(width) << f[i][j] << " ";
        stream << std::endl;
    }
    return stream;
}

template <class T, class functor>
T accumulate(const function2D<T>& data, functor& f)
{
    T sum = 0;
    for (int i = 0; i < data.size_N(); i++)
        for (int j = 0; j < data.size_Nd(); j++)
            sum += f(data(i, j));
    return sum;
}

extern "C"
{
    void dgemm_(const char* transa, const char* transb, const int* m, const int* n, const int* k, const double* alpha, const double* A, const int* lda, const double* B, const int* ldb, const double* beta, double* C, const int* ldc);
    void dsymm_(const char* side, const char* uplo, const int* m, const int* n, const double* alpha, const double* A, const int* lda, const double* B, const int* ldb, const double* beta, double* C, const int* ldc);
}
void dgemm_(const char* transa, const char* transb, const int* m, const int* n, const int* k, const double* alpha, const double* A, const int* lda, const double* B, const int* ldb, const double* beta, double* C, const int* ldc)
{

}
inline void xgemm(const char* transa, const char* transb, const int m, const int n,
    const int k, const double alpha, const double* A,
    const int lda, const double* B, const int ldb, const double beta,
    double* C, const int ldc)
{
    dgemm_(transa, transb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
}

inline void xsymm(const char* side, const char* uplo, int m, int n, double alpha, const double* A, int lda, const double* B, int ldb, double beta, double* C, int ldc)
{
    dsymm_(side, uplo, &m, &n, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
}

// Test for non-quadratic matrices
template <class T>
inline void function2D<T>::Product(const std::string& transa, const std::string& transb, const function2D& A, const function2D& B, const T& alpha, const T& beta)
//inline void function2D<T>::Product(const std::string& transa, const std::string& transb, const function2D& A, const function2D& B, const T& alpha = 1, const T& beta = 0)
{
    if (transa != "N" && transa != "T") { cerr << "Did not recognize your task. Specify how to multiply matrices in dgemm!" << endl; return; }
    if (transa == "N" && transb == "N") {
        if (A.Nd != B.N || !B.Nd || !A.N || !A.Nd || !B.N || N0 < A.N || Nd0 < B.Nd)
            std::cerr << " Matrix sizes not correct" << std::endl;
        xgemm("N", "N", B.Nd, A.N, B.N, alpha, B.MemPt(), B.Nd0, A.MemPt(), A.Nd0, beta, MemPt(), Nd0);
        N = A.N; Nd = B.Nd;
    }
    else if (transa == "T" && transb == "N") {
        if (A.N != B.N || !B.Nd || !A.Nd || !A.N || !B.N || N0 < A.Nd || Nd0 < B.Nd)
            std::cerr << " Matrix sizes not correct" << std::endl;
        xgemm("N", "T", B.Nd, A.Nd, B.N, alpha, B.MemPt(), B.Nd0, A.MemPt(), A.Nd0, beta, MemPt(), Nd0);
        N = A.Nd; Nd = B.Nd;
    }
    else if (transa == "N" && transb == "T") {
        if (A.Nd != B.Nd || !B.N || !A.N || !A.Nd || !B.Nd || N0 < A.N || Nd0 < B.N)
            std::cerr << " Matrix sizes not correct" << std::endl;
        xgemm("T", "N", B.N, A.N, B.Nd, alpha, B.MemPt(), B.Nd0, A.MemPt(), A.Nd0, beta, MemPt(), Nd0);
        N = A.N; Nd = B.N;
    }
    else if (transa == "T" && transb == "T") {
        if (A.N != B.Nd || !B.N || !A.N || !A.Nd || !B.Nd || N0 < A.Nd || Nd0 < B.N)
            std::cerr << " Matrix sizes not correct" << std::endl;
        xgemm("T", "T", B.N, A.Nd, B.Nd, alpha, B.MemPt(), B.Nd0, A.MemPt(), A.Nd0, beta, MemPt(), Nd0);
        N = A.Nd; Nd = B.N;
    }
}

template <class T>
//inline void function2D<T>::MProduct(const function2D& A, const function2D& B, const T& alpha = 1, const T& beta = 0)
inline void function2D<T>::MProduct(const function2D& A, const function2D& B, const T& alpha, const T& beta)
{
    if (A.Nd != B.N || !B.Nd || !A.N || !A.Nd || !B.N || N0 < A.N || Nd0 < B.Nd)
        std::cerr << " Matrix sizes not correct" << std::endl;
    xgemm("N", "N", B.Nd, A.N, B.N, alpha, B.MemPt(), B.Nd0, A.MemPt(), A.Nd0, beta, MemPt(), Nd0);
    N = A.N; Nd = B.Nd;
}

template <class T>
inline T identity(const T& x) { return x; }

extern "C" void dgetrf_(int* n1, int* n2, double* a, int* lda, int* ipiv, int* info);

double Determinant(function2D<double>& A)
{
    if (A.size_Nd() != A.size_N()) { std::cerr << "Can't compute determinant of nonquadratic matrix!" << std::endl; return 0; }
    int info;
    int n = A.size_N();
    int lda = A.fullsize_Nd();
    function1D<int> ipiv(n);
    dgetrf_(&n, &n, A.MemPt(), &lda, ipiv.MemPt(), &info);
    if (info) { std::cerr << "LU factorization complains : " << info << std::endl; return 0; }
    double det = 1;
    for (int i = 0; i < n; i++) det *= ((ipiv[i] == i) ? 1 : -1) * A(i, i);
    return det;
}

#endif