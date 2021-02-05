#include "ED.h"

#include <cmath>
#include <cstring>

#include "mkl.h"

const double ED::pi = 3.1415926535897932385;
const ED::Complex ED::I(0, 1);

ED::ED(int Ns_, double t_, double U_, double T_, double mu_)
    : Ns(Ns_),
      NN((Ns_ + 1) * Ns_),
      NNN((Ns_ + 1) * (Ns_ + 1) * Ns_),
      DIM(1 << (2 * Ns_)),
      onespindim(1 << Ns_),
      t(t_),
      U(U_),
      T(T_),
      mu(mu_),
      Z(0) {
    ones = new int[onespindim];
    ones[0] = 0;
    for (int n = 1; n < onespindim; n++) {
        ones[n] = ones[n & (n - 1)] + 1;
    }

    double TwoPiOverNs = 2. * pi / Ns;
    double* epsilon = new double[Ns];
    for (int k = 0; k < Ns; k++) {
        epsilon[k] = -2. * t * std::cos(TwoPiOverNs * k) - mu;
    }

    // sort states and initialize kinetic energy
    dim = new int[NNN];
    memset(dim, 0, NNN * sizeof(int));
    index = new int[DIM];
    int* block = new int[DIM];
    double* Ek = new double[DIM];
    for (int state = 0; state < DIM; state++) {
        int upstate = state % onespindim, downstate = state / onespindim;
        int nu = ones[upstate], nd = ones[downstate];
        int p = 0;
        double ek = 0;
        for (int k = 0; k < Ns; k++) {
            int nk = (upstate & 1) + (downstate & 1);
            p += nk * k;
            ek += nk * epsilon[k];
            upstate >>= 1, downstate >>= 1;
            if (upstate == 0 && downstate == 0) break;
        }

        Ek[state] = ek;
        p %= Ns;
        int blkid = NN * nu + Ns * nd + p;
        block[state] = blkid;
        index[state] = dim[blkid]++;
    }

    delete[] epsilon;

    base = new int*[NNN];
    for (int blkid = 0; blkid < NNN; blkid++) {
        if (dim[blkid] > 0) {
            base[blkid] = new int[dim[blkid]];
        }
    }
    for (int state = 0; state < DIM; state++) {
        base[block[state]][index[state]] = state;
    }

    delete[] block;

    H = new double*[NNN];
    egvalues = new double*[NNN];
    for (int blkid = 0; blkid < NNN; blkid++) {
        if (dim[blkid] > 0) {
            H[blkid] = new double[dim[blkid] * dim[blkid]];
            egvalues[blkid] = new double[dim[blkid]];
        }
    }

    // write H matrix and solve eigensystem
    double UOverNs = U / Ns;
    double ground_energy = 2 * Ns * U;
    for (int blkid = 0; blkid < NNN; blkid++) {
        if (dim[blkid] == 0) continue;
        memset(H[blkid], 0, dim[blkid] * dim[blkid] * sizeof(double));

#pragma omp parallel for
        for (int ridx = 0; ridx < dim[blkid]; ridx++) {
            int rstate = base[blkid][ridx];
            int rupstate = rstate % onespindim, rdownstate = rstate / onespindim;

            // 1. diagonal 2rd order
            H[blkid][ridx * dim[blkid] + ridx] = Ek[rstate];

            // 2. diagonal 4th order k1 == k2 && k3 == k4
            for (int k1 = 0; k1 < Ns; k1++) {
                if (((rupstate >> k1) & 1) == 0) continue;
                for (int k3 = 0; k3 < Ns; k3++) {
                    if (((rdownstate >> k3) & 1) == 0) continue;
                    H[blkid][ridx * dim[blkid] + ridx] += UOverNs;
                }
            }

            // 3. off diagonal
            for (int k1 = 0; k1 < Ns; k1++) {
                for (int k2 = 0; k2 < Ns; k2++) {
                    if (((rupstate >> k1) & 1) || !((rupstate >> k2) & 1)) continue;
                    int lupstate = (rupstate ^ (1 << k1)) ^ (1 << k2);
                    bool minusup = onesbetween(rupstate, k1, k2) & 1;

                    for (int k3 = 0; k3 < Ns; k3++) {
                        int k4 = (Ns + k1 + k3 - k2) % Ns;
                        if (((rdownstate >> k3) & 1) || !((rdownstate >> k4) & 1)) continue;

                        bool minusdown = onesbetween(rdownstate, k3, k4) & 1;
                        int ldownstate = (rdownstate ^ (1 << k3) ^ (1 << k4));
                        int lstate = ldownstate * onespindim + lupstate;
                        if (minusup == minusdown)
                            H[blkid][index[lstate] * dim[blkid] + ridx] += UOverNs;
                        else
                            H[blkid][index[lstate] * dim[blkid] + ridx] -= UOverNs;
                    }
                }
            }
        }

        LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', dim[blkid], H[blkid], dim[blkid], egvalues[blkid]);

        for (int idx = 0; idx < dim[blkid]; idx++) {
            ground_energy = std::min(ground_energy, egvalues[blkid][idx]);
        }
    }

    delete[] Ek;

    // shift eigen values by ground_energy in case Z to be too large
    for (int blkid = 0; blkid < NNN; blkid++) {
        for (int idx = 0; idx < dim[blkid]; idx++) {
            egvalues[blkid][idx] -= ground_energy;
            Z += std::exp(-egvalues[blkid][idx] / T);
        }
    }

    cd = new double**[Ns];
    cdmap = new int*[Ns];
    for (int k = 0; k < Ns; k++) {
        cdmap[k] = nullptr;
        cd[k] = nullptr;
    }

    Mup = new double**[Ns];
    Mdn = new double**[Ns];
    Mmap = new int*[Ns];
    for (int k = 0; k < Ns; k++) {
        Mmap[k] = nullptr;
        Mup[k] = nullptr;
        Mdn[k] = nullptr;
    }
}

ED::~ED() {
    for (int k = 0; k < Ns; k++) {
        if (Mmap[k] != nullptr) {
            for (int blkid = 0; blkid < NNN; blkid++) {
                if (Mup[k][blkid] != nullptr) {
                    delete[] Mup[k][blkid];
                }
                if (Mdn[k][blkid] != nullptr) {
                    delete[] Mdn[k][blkid];
                }
            }
            delete[] Mup[k];
            delete[] Mdn[k];
            delete[] Mmap[k];
        }
    }
    delete[] Mup;
    delete[] Mdn;
    delete[] Mmap;

    for (int k = 0; k < Ns; k++) {
        if (cdmap[k] != nullptr) {
            for (int blkid = 0; blkid < NNN; blkid++) {
                if (cd[k][blkid] != nullptr) {
                    delete[] cd[k][blkid];
                }
            }
            delete[] cdmap[k];
            delete[] cd[k];
        }
    }
    delete[] cdmap;
    delete[] cd;

    for (int blkid = 0; blkid < NNN; blkid++) {
        if (dim[blkid] > 0) {
            delete[] base[blkid];
            delete[] H[blkid];
            delete[] egvalues[blkid];
        }
    }

    delete[] base;
    delete[] H;
    delete[] egvalues;

    delete[] ones;
    delete[] dim;
    delete[] index;
}

void ED::setmu(double mu_) {
    Z = 0.;
    double ground_energy = 2 * Ns * U;
    for (int blkid = 0; blkid < NNN; blkid++) {
        double che_diff = (blkid / NN + blkid % NN / Ns) * (mu_ - mu);
        for (int idx = 0; idx < dim[blkid]; idx++) {
            egvalues[blkid][idx] -= che_diff;
            ground_energy = std::min(ground_energy, egvalues[blkid][idx]);
        }
    }
    mu = mu_;

    for (int blkid = 0; blkid < NNN; blkid++) {
        for (int idx = 0; idx < dim[blkid]; idx++) {
            egvalues[blkid][idx] -= ground_energy;
            Z += std::exp(-egvalues[blkid][idx] / T);
        }
    }
}

void ED::setT(double T_) {
    T = T_;
    Z = 0.;
    for (int blkid = 0; blkid < NNN; blkid++) {
        for (int idx = 0; idx < dim[blkid]; idx++) {
            Z += std::exp(-egvalues[blkid][idx] / T);
        }
    }
}

// 1. density
double ED::density() const {
    double den = 0.;
    for (int blkid = 0; blkid < NNN; blkid++) {
        double weight = 0.;
        for (int idx = 0; idx < dim[blkid]; idx++) {
            weight += std::exp(-egvalues[blkid][idx] / T);
        }
        den += weight * (blkid / NN + (blkid % NN) / Ns);
    }
    return den / (double(Ns) * Z);
}

void ED::get_cd(int k) {
    if (cd[k] != nullptr) return;
    cd[k] = new double*[NNN];
    if (cdmap[k] == nullptr) cdmap[k] = new int[NNN];

    for (int rblkid = 0; rblkid < NNN; rblkid++) {
        cd[k][rblkid] = nullptr;
        if (dim[rblkid] == 0) continue;
        int rp = rblkid % Ns;
        int lblkid = rblkid + NN - rp + (rp + k) % Ns;
        cdmap[k][rblkid] = lblkid;
        if (lblkid >= NNN || dim[lblkid] == 0) continue;

        double* temp = new double[dim[lblkid] * dim[rblkid]];
        memset(temp, 0, dim[lblkid] * dim[rblkid] * sizeof(double));

        for (int ridx = 0; ridx < dim[rblkid]; ridx++) {
            int rstate = base[rblkid][ridx];
            int rupstate = rstate % onespindim;

            if ((rupstate >> k) & 1) continue;
            int lidx = index[rstate ^ (1 << k)];
            memcpy(&temp[lidx * dim[rblkid]], &H[rblkid][ridx * dim[rblkid]], dim[rblkid] * sizeof(double));

            if ((ones[rupstate] - ones[rupstate >> k]) & 1) {
                cblas_dscal(dim[rblkid], -1., &temp[lidx * dim[rblkid]], 1);
            }
        }

        cd[k][rblkid] = new double[dim[lblkid] * dim[rblkid]];
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, dim[lblkid], dim[rblkid], dim[lblkid], 1., H[lblkid],
                    dim[lblkid], temp, dim[rblkid], 0, cd[k][rblkid], dim[rblkid]);
        delete[] temp;
    }
}

// 2. Green's function in imaginary time
double ED::Gtau(int k, double tau) {
    double res = 0.;
    get_cd(k);
    for (int rblkid = 0; rblkid < NNN; rblkid++) {
        if (cd[k][rblkid] == nullptr) continue;

        int lblkid = cdmap[k][rblkid];

        double temp = 0;
#pragma omp parallel for reduction(+ : temp)
        for (int lidx = 0; lidx < dim[lblkid]; lidx++) {
            for (int ridx = 0; ridx < dim[rblkid]; ridx++) {
                double ele = cd[k][rblkid][lidx * dim[rblkid] + ridx];
                temp +=
                    std::exp(-egvalues[rblkid][ridx] / T + tau * (egvalues[rblkid][ridx] - egvalues[lblkid][lidx])) *
                    ele * ele;
            }
        }
        res += temp;
    }
    return res / Z;
}

void ED::get_Mup(int k) {
    if (Mup[k] != nullptr) return;
    Mup[k] = new double*[NNN];
    if (Mmap[k] == nullptr) Mmap[k] = new int[NNN];

    // k == 0
    if (k == 0) {
        for (int blkid = 0; blkid < NNN; blkid++) {
            Mup[k][blkid] = nullptr;
            if (dim[blkid] == 0) continue;
            Mmap[k][blkid] = blkid;
            Mup[k][blkid] = new double[dim[blkid] * dim[blkid]];
            memset(Mup[k][blkid], 0, dim[blkid] * dim[blkid] * sizeof(double));
            int Nup = blkid / NN;
            for (int idx = 0; idx < dim[blkid]; idx++) {
                Mup[k][blkid][idx * dim[blkid] + idx] = Nup;
            }
        }
        return;
    }

    // k != 0
    for (int rblkid = 0; rblkid < NNN; rblkid++) {
        Mup[k][rblkid] = nullptr;
        if (dim[rblkid] == 0) continue;

        int rp = rblkid % Ns;
        int lblkid = rblkid - rp + (rp + k) % Ns;
        Mmap[k][rblkid] = lblkid;
        if (dim[lblkid] == 0) continue;

        Mup[k][rblkid] = new double[dim[lblkid] * dim[rblkid]];
        memset(Mup[k][rblkid], 0, dim[lblkid] * dim[rblkid] * sizeof(double));

        for (int ridx = 0; ridx < dim[rblkid]; ridx++) {
            int rstate = base[rblkid][ridx];
            int rupstate = rstate % onespindim;
            for (int k2 = 0; k2 < Ns; k2++) {
                int k1 = (k2 + k) % Ns;
                if (((rupstate >> k1) & 1) || !((rupstate >> k2) & 1)) continue;

                int lstate = rstate ^ (1 << k1) ^ (1 << k2);
                if (onesbetween(rupstate, k1, k2) & 1)
                    Mup[k][rblkid][index[lstate] * dim[rblkid] + ridx] -= 1.;
                else
                    Mup[k][rblkid][index[lstate] * dim[rblkid] + ridx] += 1.;
            }
        }

        double* temp = new double[dim[lblkid] * dim[rblkid]];
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, dim[lblkid], dim[rblkid], dim[lblkid], 1., H[lblkid],
                    dim[lblkid], Mup[k][rblkid], dim[rblkid], 0, temp, dim[rblkid]);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim[lblkid], dim[rblkid], dim[rblkid], 1., temp,
                    dim[rblkid], H[rblkid], dim[rblkid], 0, Mup[k][rblkid], dim[rblkid]);
        delete[] temp;
    }
}

void ED::get_Mdn(int k) {
    if (Mdn[k] != nullptr) return;
    Mdn[k] = new double*[NNN];
    if (Mmap[k] == nullptr) Mmap[k] = new int[NNN];

    // k == 0
    if (k == 0) {
        for (int blkid = 0; blkid < NNN; blkid++) {
            Mdn[k][blkid] = nullptr;
            if (dim[blkid] == 0) continue;
            Mmap[k][blkid] = blkid;
            Mdn[k][blkid] = new double[dim[blkid] * dim[blkid]];
            memset(Mdn[k][blkid], 0, dim[blkid] * dim[blkid] * sizeof(double));
            int Ndn = blkid % NN / Ns;
            for (int idx = 0; idx < dim[blkid]; idx++) {
                Mdn[k][blkid][idx * dim[blkid] + idx] = Ndn;
            }
        }
        return;
    }

    // k != 0
    for (int rblkid = 0; rblkid < NNN; rblkid++) {
        Mdn[k][rblkid] = nullptr;
        if (dim[rblkid] == 0) continue;

        int rp = rblkid % Ns;
        int lblkid = rblkid - rp + (rp + k) % Ns;
        Mmap[k][rblkid] = lblkid;
        if (dim[lblkid] == 0) continue;

        Mdn[k][rblkid] = new double[dim[lblkid] * dim[rblkid]];
        memset(Mdn[k][rblkid], 0, dim[lblkid] * dim[rblkid] * sizeof(double));

        for (int ridx = 0; ridx < dim[rblkid]; ridx++) {
            int rstate = base[rblkid][ridx];
            int rupstate = rstate % onespindim, rdnstate = rstate / onespindim;
            for (int k2 = 0; k2 < Ns; k2++) {
                int k1 = (k2 + k) % Ns;
                if (((rdnstate >> k1) & 1) || !((rdnstate >> k2) & 1)) continue;

                int lstate = (rdnstate ^ (1 << k1) ^ (1 << k2)) * onespindim + rupstate;
                if (onesbetween(rdnstate, k1, k2) & 1)
                    Mdn[k][rblkid][index[lstate] * dim[rblkid] + ridx] -= 1.;
                else
                    Mdn[k][rblkid][index[lstate] * dim[rblkid] + ridx] += 1.;
            }
        }

        double* temp = new double[dim[lblkid] * dim[rblkid]];
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, dim[lblkid], dim[rblkid], dim[lblkid], 1., H[lblkid],
                    dim[lblkid], Mdn[k][rblkid], dim[rblkid], 0, temp, dim[rblkid]);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim[lblkid], dim[rblkid], dim[rblkid], 1., temp,
                    dim[rblkid], H[rblkid], dim[rblkid], 0, Mdn[k][rblkid], dim[rblkid]);
        delete[] temp;
    }
}

// 3. Density-Density correlation in imaginary time
double ED::DenDen(int k, double tau) {
    double res = 0.;
    int invk = (Ns - k) % Ns;
    get_Mup(invk);
    get_Mup(k);
    get_Mdn(k);

    for (int rblkid = 0; rblkid < NNN; rblkid++) {
        if (dim[rblkid] == 0 || dim[Mmap[k][rblkid]] == 0) continue;
        int lblkid = Mmap[k][rblkid];

        double* ldiag = new double[dim[lblkid]];
        for (int lidx = 0; lidx < dim[lblkid]; lidx++) {
            ldiag[lidx] = std::exp(-tau * egvalues[lblkid][lidx]);
        }

        double* x = new double[dim[lblkid]];

        for (int ridx = 0; ridx < dim[rblkid]; ridx++) {
            vdMul(dim[lblkid], &Mup[invk][lblkid][ridx * dim[lblkid]], ldiag, x);
            double dot = cblas_ddot(dim[lblkid], x, 1, &Mup[k][rblkid][ridx], dim[rblkid]) +
                         cblas_ddot(dim[lblkid], x, 1, &Mdn[k][rblkid][ridx], dim[rblkid]);
            res += std::exp(-(1. / T - tau) * egvalues[rblkid][ridx]) * dot;
        }

        delete[] ldiag;
        delete[] x;
    }

    return 2. * res / (double(Ns) * Z);
}

// 4. SzSZ correlation in imaginary time
double ED::SzSz(int k, double tau) {
    double res = 0.;
    int invk = (Ns - k) % Ns;
    get_Mup(invk);
    get_Mup(k);
    get_Mdn(k);

    for (int rblkid = 0; rblkid < NNN; rblkid++) {
        if (dim[rblkid] == 0 || dim[Mmap[k][rblkid]] == 0) continue;
        int lblkid = Mmap[k][rblkid];

        double* ldiag = new double[dim[lblkid]];
        for (int lidx = 0; lidx < dim[lblkid]; lidx++) {
            ldiag[lidx] = std::exp(-tau * egvalues[lblkid][lidx]);
        }

        double* x = new double[dim[lblkid]];

        for (int ridx = 0; ridx < dim[rblkid]; ridx++) {
            vdMul(dim[lblkid], &Mup[invk][lblkid][ridx * dim[lblkid]], ldiag, x);
            double dot = cblas_ddot(dim[lblkid], x, 1, &Mup[k][rblkid][ridx], dim[rblkid]) -
                         cblas_ddot(dim[lblkid], x, 1, &Mdn[k][rblkid][ridx], dim[rblkid]);
            res += std::exp(-(1. / T - tau) * egvalues[rblkid][ridx]) * dot;
        }

        delete[] x;
    }

    return 2. * res / (double(Ns) * Z);
}