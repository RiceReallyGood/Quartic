#ifndef ED_H
#define ED_H

#include <complex>
#include <vector>
class ED
{
public:
    using Complex = std::complex<double>;
    ED(int Ns_, double t_, double U_, double T_, double delmu_);
    ED(const ED &) = delete;
    ED &operator=(const ED &) = delete;
    ~ED();
    Complex greensfunction(int k, int n) const;         //1. Green's function
    double nup(int k) const;                            //2. number of spinups
    double denscorr(int idx) const;                     //3. density-density correlation of equaltime
    double spincorr(int idx) const;                     //4. spin-spin correlation of equaltime
    Complex denscorr(int idx, int n) const;             //5. density-density correlation(with frequency)
    Complex spincorr(int idx, int n) const;             //6. spin-spin correlation(with frequency)
    double spectralfunction(int k, double omega) const; //7. spectral function
    std::vector<Complex> greensfunction(int k, std::vector<int> ntab) const;
    std::vector<double>  greensfunction(int k, std::vector<double> ttab) const;
    std::vector<Complex> spincorr(int idx, std::vector<int> ntab) const;
    std::vector<double> spectralfunction(int k, std::vector<double> omegatab) const;

private:
    //parameters
    const int Ns;           //number of sites
    const int NN;           //NN = Ns * (Ns + 1);
    const int NNN;          //NNN = (Ns + 1) * NN;
    const int Dim;          //Dim = 1 << (2 * Ns);
    const int onespindim;   //onespindim = 1 << Ns;
    const double t;         //hopping energy
    const double T;         //temperature
    const double U;         //onsite interaction
    const double mu;        //chemical potential mu = U/2 + delmu
    static const double pi; //pi = 3.1415926535897932385;
    static const Complex I; //imaginary unit
    std::vector<int> index;
    std::vector<int> dim;
    std::vector<std::vector<int>> base;
    std::vector<int> ones; //the number of '1' bits in the binary form of a integer
    int onesbetween(int state, int k1, int k2) const
    {
        if (k2 < k1)
            std::swap(k1, k2);
        return ones[state >> (k1 + 1)] - ones[state >> k2];
    }

    double **K;        //matrices of Hamiltonian
    double **egvalues; //eigen values of Hamiltonian
    double **w;        //exp(-egvalues/T)
    double Z;          //partition function
};

#endif