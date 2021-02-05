#ifndef ED_H
#define ED_H

#include <complex>

class ED {
public:
    using Complex = std::complex<double>;
    ED(int Ns, double t, double U, double T, double mu);
    ED(const ED&) = delete;
    ED& operator=(const ED&) = delete;
    ~ED();

    void setT(double T);
    void setmu(double mu);

    double density() const;            // 1. density
    double Gtau(int k, double tau);    // 2. Green's function in imaginary time
    double DenDen(int k, double tau);  // 3. Density-Density correlation in imaginary time
    double SzSz(int k, double tau);    // 4. SzSZ correlation in imaginary time

private:
    const int Ns;            // Number of sites
    const int NN;            // NN = Ns * (Ns + 1)
    const int NNN;           // NNN = (Ns + 1) * NN
    const int DIM;           // DIM = 1 << (2 * Ns)
    const int onespindim;    // onespindim = 1 << Ns;
    const double t;          // Hopping Energy
    const double U;          // Onsite interaction
    static const double pi;  // pi = 3.1415926535897932385;
    static const Complex I;  // Imaginary unit
                             //
    double T;                // Temperature, default T = 1
    double mu;               // Chemical potential, default mu = U / 2

    int* ones;
    int* dim;
    int* index;
    int** base;

    double** H;         // Matrices of Hamiltonian
    double** egvalues;  // Eigen values of Hamiltonian
    double Z;           // Partition function

    double*** cd;
    int** cdmap;

    double*** Mup;
    double*** Mdn;
    int** Mmap;

    int onesbetween(int state, int k1, int k2) {
        return k1 < k2 ? ones[state >> (k1 + 1)] - ones[state >> k2] : ones[state >> (k2 + 1)] - ones[state >> k1];
    }

    void get_cd(int k);
    void get_Mup(int k);
    void get_Mdn(int k);
};

#endif  // ED_H