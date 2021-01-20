#include "ED.h"
#include <vector>
#include "mkl.h"
#include "mkl_spblas.h"
const double ED::pi = 3.1415926535897932385;
const ED::Complex ED::I(0, 1);

ED::ED(int Ns_, double t_, double U_, double T_, double delmu_ = 0)
    : Ns(Ns_),
      NN(Ns * (Ns + 1)),
      NNN((NN * (Ns + 1))),
      Dim(1 << (2 * Ns)),
      onespindim(1 << Ns),
      t(t_),
      U(U_),
      T(T_),
      mu(U / 2 + delmu_),
      index(Dim),
      dim(NNN),
      base(NNN, std::vector<int>()),
      ones(onespindim),
      Z(0)
{
    double TwoPiOverNs = 2. * pi / Ns;
    std::vector<double> epsilon(Ns);
    for (int k = 0; k < Ns; k++)
        epsilon[k] = -2. * t * cos(TwoPiOverNs * k) - mu;

    ones[0] = 0;
    for (int n = 1; n < ones.size(); n++)
        ones[n] = ones[n & (n - 1)] + 1;

    //sort states and initialize kinetic energy
    std::vector<double> Ek(Dim);
    for (int state = 0; state < Dim; state++)
    {
        int upstate = state % onespindim, downstate = state / onespindim;
        int nu = ones[upstate], nd = ones[downstate];
        int p = 0;
        double ek = 0;
        for (int k = 0; k < Ns; k++)
        {
            int nk = (upstate & 1) + (downstate & 1);
            p += nk * k;
            ek += nk * epsilon[k];
            upstate >>= 1, downstate >>= 1;
            if (upstate == 0 && downstate == 0)
                break;
        }

        Ek[state] = ek;
        p %= Ns;
        int block = NN * nu + Ns * nd + p;
        index[state] = base[block].size();
        base[block].push_back(state);
    }

    //record dimension of each block
    for (int block = 0; block < NNN; block++)
        dim[block] = base[block].size();

    //allocate memory
    K = new double *[NNN];
    egvalues = new double *[NNN];
    w = new double *[NNN];
    for (int block = 0; block < NNN; block++)
    {
        K[block] = new double[dim[block] * dim[block]];
        egvalues[block] = new double[dim[block]];
        w[block] = new double[dim[block]];
    }

    //write K matrix and solve eigensystem
    double UOverNs = U / Ns;
    double ground_energy = 2 * Ns * U;
    for (int block = 0; block < NNN; block++)
    {
        if (dim[block] == 0)
            continue;
        for (int i = 0; i < dim[block] * dim[block]; i++)
            K[block][i] = 0;
//int nu = block / NN, nd = (block % NN) / Ns, p = block % Ns;
#pragma omp parallel for
        for (int rindex = 0; rindex < dim[block]; rindex++)
        {
            int rstate = base[block][rindex];
            int rupstate = rstate % onespindim, rdownstate = rstate / onespindim;
            //1.diagonal 2rd order
            K[block][rindex * dim[block] + rindex] = Ek[rstate];
            //2. diagonal 4th order
            for (int k1 = 0; k1 < Ns; k1++)
            {
                if (!((rupstate >> k1) & 1))
                    continue;
                for (int k3 = 0; k3 < Ns; k3++)
                {
                    if (!((rdownstate >> k3) & 1))
                        continue;
                    K[block][rindex * dim[block] + rindex] += UOverNs;
                }
            }

            //3.off diagonal
            for (int k1 = 0; k1 < Ns; k1++)
            {
                for (int k2 = 0; k2 < Ns; k2++)
                {
                    if (((rupstate >> k1) & 1) || !((rupstate >> k2) & 1))
                        continue;
                    int lupstate = (rupstate ^ (1 << k1)) ^ (1 << k2);
                    bool minusup = onesbetween(rupstate, k1, k2) & 1;
                    for (int k3 = 0; k3 < Ns; k3++)
                    {
                        int k4 = (Ns + k1 + k3 - k2) % Ns;
                        if (((rdownstate >> k3) & 1) || !((rdownstate >> k4) & 1))
                            continue;
                        bool minusdown = onesbetween(rdownstate, k3, k4) & 1;
                        int ldownstate = (rdownstate ^ (1 << k3) ^ (1 << k4));
                        int lstate = ldownstate * onespindim + lupstate;
                        if (minusup == minusdown)
                            K[block][index[lstate] * dim[block] + rindex] += UOverNs;
                        else
                            K[block][index[lstate] * dim[block] + rindex] -= UOverNs;
                    }
                }
            }
        }

        LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', dim[block], K[block], dim[block], egvalues[block]);

        for (int index = 0; index < dim[block]; index++)
            if (egvalues[block][index] < ground_energy)
                ground_energy = egvalues[block][index];
    }

    //shift eigen values by ground_energy in case Z to be too large
    for (int block = 0; block < NNN; block++)
    {
        for (int index = 0; index < dim[block]; index++)
        {
            egvalues[block][index] -= ground_energy;
            w[block][index] = exp(-egvalues[block][index] / T);
            Z += w[block][index];
        }
    }
}

ED::~ED()
{
    for (int block = 0; block < NNN; block++)
    {
        delete[] K[block];
        delete[] egvalues[block];
    }

    delete[] K;
    delete[] egvalues;
}

//physical quantities

//1. Green's function
ED::Complex ED::greensfunction(int k, int n) const
{
    Complex iwn = I * pi * T * (2. * n + 1);
    Complex ret(0, 0);
    for (int rblock = 0; rblock < NNN; rblock++)
    {
        if (dim[rblock] == 0)
            continue;
        int rp = rblock % Ns;
        int lblock = rblock + NN + (rp + k) % Ns - rp;
        if (lblock >= NNN || dim[lblock] == 0)
            continue;
        //write cd matrix
        int nnz = 0;
        double *val = new double[dim[rblock]];
        int *row_index = new int[dim[rblock]];
        int *col_index = new int[dim[rblock]];
        for (int rindex = 0; rindex < dim[rblock]; rindex++)
        {
            int rstate = base[rblock][rindex];
            int rupstate = rstate % onespindim;
            if ((rstate >> k) & 1)
                continue;
            int lstate = rstate ^ (1 << k);
            if ((ones[rupstate] - ones[rupstate >> k]) & 1)
                val[nnz] = -1;
            else
                val[nnz] = 1;
            row_index[nnz] = index[lstate];
            col_index[nnz] = rindex;
            nnz++;
        }
        if (nnz == 0)
            continue;
        sparse_matrix_t cdold;
        mkl_sparse_d_create_coo(&cdold, SPARSE_INDEX_BASE_ZERO, dim[lblock], dim[rblock], nnz, row_index, col_index, val);
        struct matrix_descr descr = {SPARSE_MATRIX_TYPE_GENERAL, SPARSE_FILL_MODE_UPPER, SPARSE_DIAG_NON_UNIT};
        double *temp = new double[dim[lblock] * dim[rblock]];
        mkl_sparse_d_mm(SPARSE_OPERATION_NON_TRANSPOSE, 1., cdold, descr, SPARSE_LAYOUT_ROW_MAJOR, K[rblock],
                        dim[rblock], dim[rblock], 0., temp, dim[rblock]);
        delete[] val;
        delete[] row_index;
        delete[] col_index;
        mkl_sparse_destroy(cdold);

        double *cd = new double[dim[lblock] * dim[rblock]];
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, dim[lblock], dim[rblock], dim[lblock], 1., K[lblock],
                    dim[lblock], temp, dim[rblock], 0., cd, dim[rblock]);
        delete[] temp;
        double rtemp = 0, itemp = 0;
#pragma omp parallel for reduction(+ \
                                   : rtemp, itemp)
        for (int lindex = 0; lindex < dim[lblock]; lindex++)
        {
            for (int rindex = 0; rindex < dim[rblock]; rindex++)
            {
                cd[lindex * dim[rblock] + rindex] *= cd[lindex * dim[rblock] + rindex];
                Complex coeff = (w[lblock][lindex] + w[rblock][rindex]) / (iwn + egvalues[lblock][lindex] - egvalues[rblock][rindex]);
                rtemp -= coeff.real() * cd[lindex * dim[rblock] + rindex];
                itemp -= coeff.imag() * cd[lindex * dim[rblock] + rindex];
            }
        }
        ret += Complex(rtemp, itemp);
        delete[] cd;
    }
    return ret / Z;
}

//2. number of spinups
double ED::nup(int k) const
{
    double ret = 0;
    for (int block = 0; block < NNN; block++)
    {
        for (int lindex = 0; lindex < dim[block]; lindex++)
        {
            if ((base[block][lindex] >> k) & 1)
            {
                for (int rindex = 0; rindex < dim[block]; rindex++)
                {
                    ret += w[block][rindex] * K[block][lindex * dim[block] + rindex] * K[block][lindex * dim[block] + rindex];
                }
            }
        }
    }
    return ret / Z;
}

//3. density-density correlation
double ED::denscorr(int idx) const
{
    double TwoPiOverNs = 2. * pi / Ns;
    double ret = 0;
    for (int block = 0; block < NNN; block++)
    {
        if (dim[block] == 0)
            continue;
        double *Fold = new double[dim[block] * dim[block]];
        for (int i = 0; i < dim[block] * dim[block]; i++)
            Fold[i] = 0;
        //write matrix
        for (int rindex = 0; rindex < dim[block]; rindex++)
        {
            int rstate = base[block][rindex];
            int rupstate = rstate % onespindim, rdownstate = rstate / onespindim;
            //1.diagonal part k1==k2 && k3==k4
            for (int k1 = 0; k1 < Ns; k1++)
            {
                if (!((rupstate >> k1) & 1))
                    continue;
                for (int k3 = 0; k3 < Ns; k3++)
                    Fold[rindex * dim[block] + rindex] += ((rupstate >> k3) & 1) + ((rdownstate >> k3) & 1);
            }
            //2.off-diagonal part
            for (int k3 = 0; k3 < Ns; k3++)
            {
                for (int k4 = 0; k4 < Ns; k4++)
                {
                    //sigma = up
                    if (!((rupstate >> k3) & 1) && ((rupstate >> k4) & 1))
                    {
                        bool minus34 = onesbetween(rupstate, k3, k4) & 1;
                        int mupstate = rupstate ^ (1 << k3) ^ (1 << k4);
                        for (int k1 = 0; k1 < Ns; k1++)
                        {
                            int k2 = (Ns + k1 + k3 - k4) % Ns;
                            if (((mupstate >> k1) & 1) || !((mupstate >> k2) & 1))
                                continue;
                            bool minus12 = onesbetween(mupstate, k1, k2) & 1;
                            int lupstate = mupstate ^ (1 << k1) ^ (1 << k2);
                            int lstate = rdownstate * onespindim + lupstate;
                            if (minus34 == minus12)
                                Fold[index[lstate] * dim[block] + rindex] += cos(TwoPiOverNs * (k2 - k1) * idx);
                            else
                                Fold[index[lstate] * dim[block] + rindex] -= cos(TwoPiOverNs * (k2 - k1) * idx);
                        }
                    }
                    //sigma = down
                    if (!((rdownstate >> k3) & 1) && ((rdownstate >> k4) & 1))
                    {
                        bool minus34 = onesbetween(rdownstate, k3, k4) & 1;
                        int ldownstate = rdownstate ^ (1 << k3) ^ (1 << k4);
                        for (int k1 = 0; k1 < Ns; k1++)
                        {
                            int k2 = (Ns + k1 + k3 - k4) % Ns;
                            if (((rupstate >> k1) & 1) || !((rupstate >> k2) & 1))
                                continue;
                            bool minus12 = onesbetween(rupstate, k1, k2) & 1;
                            int lupstate = rupstate ^ (1 << k1) ^ (1 << k2);
                            int lstate = ldownstate * onespindim + lupstate;
                            if (minus34 == minus12)
                                Fold[index[lstate] * dim[block] + rindex] += cos(TwoPiOverNs * (k2 - k1) * idx);
                            else
                                Fold[index[lstate] * dim[block] + rindex] -= cos(TwoPiOverNs * (k2 - k1) * idx);
                        }
                    }
                }
            } //end of off-diagonal part
        }     //end of write matrix
        double *temp = new double[dim[block] * dim[block]];
        double *F = new double[dim[block] * dim[block]];
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim[block], dim[block], dim[block], 1., Fold,
                    dim[block], K[block], dim[block], 0., temp, dim[block]);
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, dim[block], dim[block], dim[block], 1., K[block],
                    dim[block], temp, dim[block], 0., F, dim[block]);
        delete[] Fold;
        delete[] temp;
        for (int i = 0; i < dim[block]; i++)
            ret += w[block][i] * F[i * dim[block] + i];
        delete[] F;
    }
    double num = 0;
    for (int k = 0; k < Ns; k++)
        num += nup(k);
    num = num * 2. / Ns;
    return 2 * ret / (Z * Ns * Ns) - num * num;
}

//4. spin-spin correlation
double ED::spincorr(int idx) const
{
    double TwoPiOverNs = 2. * pi / Ns;
    double ret = 0;
    for (int block = 0; block < NNN; block++)
    {
        if (dim[block] == 0)
            continue;
        double *Qold = new double[dim[block] * dim[block]];
        for (int i = 0; i < dim[block] * dim[block]; i++)
            Qold[i] = 0;
        //write matrix
        for (int rindex = 0; rindex < dim[block]; rindex++)
        {
            int rstate = base[block][rindex];
            int rupstate = rstate % onespindim, rdownstate = rstate / onespindim;
            //1.diagonal part k1==k2 && k3==k4
            for (int k1 = 0; k1 < Ns; k1++)
            {
                if (!((rupstate >> k1) & 1))
                    continue;
                for (int k3 = 0; k3 < Ns; k3++)
                    Qold[rindex * dim[block] + rindex] += ((rupstate >> k3) & 1) - ((rdownstate >> k3) & 1);
            }
            //2.off-diagonal part
            for (int k3 = 0; k3 < Ns; k3++)
            {
                for (int k4 = 0; k4 < Ns; k4++)
                {
                    //sigma = up
                    if (!((rupstate >> k3) & 1) && ((rupstate >> k4) & 1))
                    {
                        bool minus34 = onesbetween(rupstate, k3, k4) & 1;
                        int mupstate = rupstate ^ (1 << k3) ^ (1 << k4);
                        for (int k1 = 0; k1 < Ns; k1++)
                        {
                            int k2 = (Ns + k1 + k3 - k4) % Ns;
                            if (((mupstate >> k1) & 1) || !((mupstate >> k2) & 1))
                                continue;
                            bool minus12 = onesbetween(mupstate, k1, k2) & 1;
                            int lupstate = mupstate ^ (1 << k1) ^ (1 << k2);
                            int lstate = rdownstate * onespindim + lupstate;
                            if (minus34 == minus12)
                                Qold[index[lstate] * dim[block] + rindex] += cos(TwoPiOverNs * (k2 - k1) * idx);
                            else
                                Qold[index[lstate] * dim[block] + rindex] -= cos(TwoPiOverNs * (k2 - k1) * idx);
                        }
                    }
                    //sigma = down
                    if (!((rdownstate >> k3) & 1) && ((rdownstate >> k4) & 1))
                    {
                        bool minus34 = onesbetween(rdownstate, k3, k4) & 1;
                        int ldownstate = rdownstate ^ (1 << k3) ^ (1 << k4);
                        for (int k1 = 0; k1 < Ns; k1++)
                        {
                            int k2 = (Ns + k1 + k3 - k4) % Ns;
                            if (((rupstate >> k1) & 1) || !((rupstate >> k2) & 1))
                                continue;
                            bool minus12 = onesbetween(rupstate, k1, k2) & 1;
                            int lupstate = rupstate ^ (1 << k1) ^ (1 << k2);
                            int lstate = ldownstate * onespindim + lupstate;
                            if (minus34 == minus12)
                                Qold[index[lstate] * dim[block] + rindex] -= cos(TwoPiOverNs * (k2 - k1) * idx);
                            else
                                Qold[index[lstate] * dim[block] + rindex] += cos(TwoPiOverNs * (k2 - k1) * idx);
                        }
                    }
                }
            } //end of off-diagonal part
        }     //end of write matrix
        double *temp = new double[dim[block] * dim[block]];
        double *Q = new double[dim[block] * dim[block]];
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim[block], dim[block], dim[block], 1., Qold,
                    dim[block], K[block], dim[block], 0., temp, dim[block]);
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, dim[block], dim[block], dim[block], 1., K[block],
                    dim[block], temp, dim[block], 0., Q, dim[block]);
        delete[] Qold;
        delete[] temp;
        for (int i = 0; i < dim[block]; i++)
            ret += w[block][i] * Q[i * dim[block] + i];
        delete[] Q;
    }
    return 2 * ret / (Z * Ns * Ns);
}

//5. density-density correlation(with frequency)
ED::Complex ED::denscorr(int idx, int n) const
{
    Complex iwn = I * pi * T * (2. * n);
    Complex ret(0, 0);
    for (int q = 0; q < Ns; q++)
    {
        for (int lblock = 0; lblock < NNN; lblock++)
        {
            int lp = lblock % Ns;
            int rblock = lblock - lp + (lp + q) % Ns;
            if (dim[lblock] == 0 || dim[rblock] == 0)
                continue;

            //Matrix1 dim[rblock] * dim[lblock]
            double *M1 = new double[dim[rblock] * dim[lblock]];
            //Matrix2 dim[lblock] * dim[rblock]
            double *M2 = new double[dim[lblock] * dim[rblock]];
            for (int i = 0; i < dim[rblock] * dim[lblock]; i++)
                M1[i] = 0, M2[i] = 0;

            if (q != 0)
            {
                for (int lindex = 0; lindex < dim[lblock]; lindex++)
                {
                    int lstate = base[lblock][lindex];
                    for (int k2 = 0; k2 < Ns; k2++)
                    {
                        int k1 = (k2 + q) % Ns;
                        if (((lstate >> k1) & 1) || !((lstate >> k2) & 1))
                            continue;
                        int rstate = (lstate ^ (1 << k1)) ^ (1 << k2);
                        if (onesbetween(lstate % onespindim, k1, k2) & 1)
                            M1[index[rstate] * dim[lblock] + lindex] -= cos(2. * pi * q * idx / Ns);
                        else
                            M1[index[rstate] * dim[lblock] + lindex] += cos(2. * pi * q * idx / Ns);
                    }
                }

                for (int rindex = 0; rindex < dim[rblock]; rindex++)
                {
                    int rstate = base[rblock][rindex];
                    int rupstate = rstate % onespindim, rdownstate = rstate / onespindim;
                    for (int k4 = 0; k4 < Ns; k4++)
                    {
                        int k3 = (Ns + k4 - q) % Ns;

                        //sigma = up
                        if (!((rupstate >> k3) & 1) && ((rupstate >> k4) & 1))
                        {
                            int lupstate = (rupstate ^ (1 << k3)) ^ (1 << k4);
                            int lstate = rdownstate * onespindim + lupstate;
                            if (onesbetween(rupstate, k3, k4) & 1)
                                M2[index[lstate] * dim[rblock] + rindex] -= 1.;
                            else
                                M2[index[lstate] * dim[rblock] + rindex] += 1.;
                        }

                        //sigma = down
                        if (!((rdownstate >> k3) & 1) && ((rdownstate >> k4) & 1))
                        {
                            int ldownstate = (rdownstate ^ (1 << k3)) ^ (1 << k4);
                            int lstate = ldownstate * onespindim + rupstate;
                            if (onesbetween(rdownstate, k3, k4) & 1)
                                M2[index[lstate] * dim[rblock] + rindex] -= 1.;
                            else
                                M2[index[lstate] * dim[rblock] + rindex] += 1.;
                        }
                    }
                }
            }
            else //(q=0)diagonal
            {
                for (int index = 0; index < dim[lblock]; index++)
                {
                    int state = base[lblock][index];
                    int downstate = state / onespindim, upstate = state % onespindim;

                    for (int k = 0; k < Ns; k++)
                    {
                        M1[index * dim[lblock] + index] += (upstate >> k) & 1;
                        M2[index * dim[lblock] + index] += (upstate >> k) & 1;
                        M2[index * dim[lblock] + index] += (downstate >> k) & 1;
                    }
                }
            }

            double *temp = new double[dim[lblock] * dim[rblock]];
            for (int i = 0; i < dim[lblock] * dim[rblock]; i++)
                temp[i] = 0;

            cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, dim[rblock], dim[lblock], dim[rblock], 1., K[rblock],
                        dim[rblock], M1, dim[lblock], 0., temp, dim[lblock]);
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim[rblock], dim[lblock], dim[lblock], 1., temp,
                        dim[lblock], K[lblock], dim[lblock], 0., M1, dim[lblock]);

            for (int i = 0; i < dim[lblock] * dim[rblock]; i++)
                temp[i] = 0;
            cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, dim[lblock], dim[rblock], dim[lblock], 1., K[lblock],
                        dim[lblock], M2, dim[rblock], 0., temp, dim[rblock]);
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim[lblock], dim[rblock], dim[rblock], 1., temp,
                        dim[rblock], K[rblock], dim[rblock], 0., M2, dim[rblock]);

            for (int rindex = 0; rindex < dim[rblock]; rindex++)
            {
                for (int lindex = 0; lindex < dim[lblock]; lindex++)
                {
                    if (n == 0 && abs(egvalues[rblock][rindex] - egvalues[lblock][lindex]) < 1e-9)
                        ret += M1[rindex * dim[lblock] + lindex] * M2[lindex * dim[rblock] + rindex] *
                               (w[rblock][rindex]) / T;
                    else
                        ret += M1[rindex * dim[lblock] + lindex] * M2[lindex * dim[rblock] + rindex] *
                               (w[lblock][lindex] - w[rblock][rindex]) / (iwn + egvalues[rblock][rindex] - egvalues[lblock][lindex]);
                }
            }

            delete[] M1;
            delete[] M2;
            delete[] temp;
        }
    }
    if (n == 0)
    {
        double num = 0;
        for (int k = 0; k < Ns; k++)
            num += nup(k);
        num = num * 2. / Ns;
        return ret * 2. / (Z * Ns * Ns) - num * num / T;
    }
    else
    {
        return ret * 2. / (Z * Ns * Ns);
    }
}

//6. spin-spin correlation(with frequency)
ED::Complex ED::spincorr(int idx, int n) const
{
    Complex iwn = I * pi * T * (2. * n);
    Complex ret(0, 0);
    for (int q = 0; q < Ns; q++)
    {
        for (int lblock = 0; lblock < NNN; lblock++)
        {
            int lp = lblock % Ns;
            int rblock = lblock - lp + (lp + q) % Ns;
            if (dim[lblock] == 0 || dim[rblock] == 0)
                continue;

            //Matrix1 dim[rblock] * dim[lblock]
            double *M1 = new double[dim[rblock] * dim[lblock]];
            //Matrix2 dim[lblock] * dim[rblock]
            double *M2 = new double[dim[lblock] * dim[rblock]];
            for (int i = 0; i < dim[rblock] * dim[lblock]; i++)
                M1[i] = 0, M2[i] = 0;

            if (q != 0)
            {
                for (int lindex = 0; lindex < dim[lblock]; lindex++)
                {
                    int lstate = base[lblock][lindex];
                    for (int k2 = 0; k2 < Ns; k2++)
                    {
                        int k1 = (k2 + q) % Ns;
                        if (((lstate >> k1) & 1) || !((lstate >> k2) & 1))
                            continue;
                        int rstate = (lstate ^ (1 << k1)) ^ (1 << k2);
                        if (onesbetween(lstate % onespindim, k1, k2) & 1)
                            M1[index[rstate] * dim[lblock] + lindex] -= cos(2. * pi * q * idx / Ns);
                        else
                            M1[index[rstate] * dim[lblock] + lindex] += cos(2. * pi * q * idx / Ns);
                    }
                }

                for (int rindex = 0; rindex < dim[rblock]; rindex++)
                {
                    int rstate = base[rblock][rindex];
                    int rupstate = rstate % onespindim, rdownstate = rstate / onespindim;
                    for (int k4 = 0; k4 < Ns; k4++)
                    {
                        int k3 = (Ns + k4 - q) % Ns;

                        //sigma = up
                        if (!((rupstate >> k3) & 1) && ((rupstate >> k4) & 1))
                        {
                            int lupstate = (rupstate ^ (1 << k3)) ^ (1 << k4);
                            int lstate = rdownstate * onespindim + lupstate;
                            if (onesbetween(rupstate, k3, k4) & 1)
                                M2[index[lstate] * dim[rblock] + rindex] -= 1.;
                            else
                                M2[index[lstate] * dim[rblock] + rindex] += 1.;
                        }

                        //sigma = down
                        if (!((rdownstate >> k3) & 1) && ((rdownstate >> k4) & 1))
                        {
                            int ldownstate = (rdownstate ^ (1 << k3)) ^ (1 << k4);
                            int lstate = ldownstate * onespindim + rupstate;
                            if (onesbetween(rdownstate, k3, k4) & 1)
                                M2[index[lstate] * dim[rblock] + rindex] += 1.;
                            else
                                M2[index[lstate] * dim[rblock] + rindex] -= 1.;
                        }
                    }
                }
            }
            else //(q=0)diagonal
            {
                for (int index = 0; index < dim[lblock]; index++)
                {
                    int state = base[lblock][index];
                    int downstate = state / onespindim, upstate = state % onespindim;

                    for (int k = 0; k < Ns; k++)
                    {
                        M1[index * dim[lblock] + index] += (upstate >> k) & 1;
                        M2[index * dim[lblock] + index] += (upstate >> k) & 1;
                        M2[index * dim[lblock] + index] -= (downstate >> k) & 1;
                    }
                }
            }

            double *temp = new double[dim[lblock] * dim[rblock]];
            for (int i = 0; i < dim[lblock] * dim[rblock]; i++)
                temp[i] = 0;

            cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, dim[rblock], dim[lblock], dim[rblock], 1., K[rblock],
                        dim[rblock], M1, dim[lblock], 0., temp, dim[lblock]);
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim[rblock], dim[lblock], dim[lblock], 1., temp,
                        dim[lblock], K[lblock], dim[lblock], 0., M1, dim[lblock]);

            for (int i = 0; i < dim[lblock] * dim[rblock]; i++)
                temp[i] = 0;
            cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, dim[lblock], dim[rblock], dim[lblock], 1., K[lblock],
                        dim[lblock], M2, dim[rblock], 0., temp, dim[rblock]);
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim[lblock], dim[rblock], dim[rblock], 1., temp,
                        dim[rblock], K[rblock], dim[rblock], 0., M2, dim[rblock]);

            for (int rindex = 0; rindex < dim[rblock]; rindex++)
            {
                for (int lindex = 0; lindex < dim[lblock]; lindex++)
                {
                    if (n == 0 && abs(egvalues[rblock][rindex] - egvalues[lblock][lindex]) < 1e-9)
                        ret += M1[rindex * dim[lblock] + lindex] * M2[lindex * dim[rblock] + rindex] *
                               (w[rblock][rindex]) / T;
                    else
                        ret += M1[rindex * dim[lblock] + lindex] * M2[lindex * dim[rblock] + rindex] *
                               (w[lblock][lindex] - w[rblock][rindex]) / (iwn + egvalues[rblock][rindex] - egvalues[lblock][lindex]);
                }
            }

            delete[] M1;
            delete[] M2;
            delete[] temp;
        }
    }
    return ret * 2. / (Z * Ns * Ns);
}

//6. spectrum weight
double ED::spectralfunction(int k, double omega) const
{
    static double eta = 0.005;
    double ret = 0;
    for (int rblock = 0; rblock < NNN; rblock++)
    {
        if (dim[rblock] == 0)
            continue;
        int rp = rblock % Ns;
        int lblock = rblock + NN + (rp + k) % Ns - rp;
        if (lblock >= NNN || dim[lblock] == 0)
            continue;
        //write cd matrix
        int nnz = 0;
        double *val = new double[dim[rblock]];
        int *row_index = new int[dim[rblock]];
        int *col_index = new int[dim[rblock]];
        for (int rindex = 0; rindex < dim[rblock]; rindex++)
        {
            int rstate = base[rblock][rindex];
            int rupstate = rstate % onespindim;
            if ((rstate >> k) & 1)
                continue;
            int lstate = rstate ^ (1 << k);
            if ((ones[rupstate] - ones[rupstate >> k]) & 1)
                val[nnz] = -1;
            else
                val[nnz] = 1;
            row_index[nnz] = index[lstate];
            col_index[nnz] = rindex;
            nnz++;
        }
        if (nnz == 0)
            continue;
        sparse_matrix_t cdold;
        mkl_sparse_d_create_coo(&cdold, SPARSE_INDEX_BASE_ZERO, dim[lblock], dim[rblock], nnz, row_index, col_index, val);
        struct matrix_descr descr = {SPARSE_MATRIX_TYPE_GENERAL, SPARSE_FILL_MODE_UPPER, SPARSE_DIAG_NON_UNIT};
        double *temp = new double[dim[lblock] * dim[rblock]];
        mkl_sparse_d_mm(SPARSE_OPERATION_NON_TRANSPOSE, 1., cdold, descr, SPARSE_LAYOUT_ROW_MAJOR, K[rblock],
                        dim[rblock], dim[rblock], 0., temp, dim[rblock]);
        delete[] val;
        delete[] row_index;
        delete[] col_index;
        mkl_sparse_destroy(cdold);

        double *cd = new double[dim[lblock] * dim[rblock]];
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, dim[lblock], dim[rblock], dim[lblock], 1., K[lblock],
                    dim[lblock], temp, dim[rblock], 0., cd, dim[rblock]);
        delete[] temp;
        double sum = 0;
#pragma omp parallel for reduction(+ \
                                   : sum)
        for (int lindex = 0; lindex < dim[lblock]; lindex++)
        {
            for (int rindex = 0; rindex < dim[rblock]; rindex++)
            {
                sum += (w[lblock][lindex] + w[rblock][rindex]) * cd[lindex * dim[rblock] + rindex] * cd[lindex * dim[rblock] + rindex] *
                       2. * eta / ((-omega + egvalues[lblock][lindex] - egvalues[rblock][rindex]) * (-omega + egvalues[lblock][lindex] - egvalues[rblock][rindex]) + eta * eta);
            }
        }
        ret += sum;
        delete[] cd;
    }
    return ret / Z;
}

std::vector<ED::Complex> ED::greensfunction(int k, std::vector<int> ntab) const
{
    int sz = ntab.size();
    std::vector<Complex> iw(sz);
    std::vector<Complex> ret(sz, Complex(0, 0));
    for (int n = 0; n < sz; n++)
        iw[n] = I * pi * T * (2. * ntab[n] + 1.);
    for (int rblock = 0; rblock < NNN; rblock++)
    {
        if (dim[rblock] == 0)
            continue;
        int rp = rblock % Ns;
        int lblock = rblock + NN + (rp + k) % Ns - rp;
        if (lblock >= NNN || dim[lblock] == 0)
            continue;
        //write cd matrix
        int nnz = 0;
        double *val = new double[dim[rblock]];
        int *row_index = new int[dim[rblock]];
        int *col_index = new int[dim[rblock]];
        for (int rindex = 0; rindex < dim[rblock]; rindex++)
        {
            int rstate = base[rblock][rindex];
            int rupstate = rstate % onespindim;
            if ((rstate >> k) & 1)
                continue;
            int lstate = rstate ^ (1 << k);
            if ((ones[rupstate] - ones[rupstate >> k]) & 1)
                val[nnz] = -1;
            else
                val[nnz] = 1;
            row_index[nnz] = index[lstate];
            col_index[nnz] = rindex;
            nnz++;
        }
        if (nnz == 0)
            continue;
        sparse_matrix_t cdold;
        mkl_sparse_d_create_coo(&cdold, SPARSE_INDEX_BASE_ZERO, dim[lblock], dim[rblock], nnz, row_index, col_index, val);
        struct matrix_descr descr = {SPARSE_MATRIX_TYPE_GENERAL, SPARSE_FILL_MODE_UPPER, SPARSE_DIAG_NON_UNIT};
        double *temp = new double[dim[lblock] * dim[rblock]];
        mkl_sparse_d_mm(SPARSE_OPERATION_NON_TRANSPOSE, 1., cdold, descr, SPARSE_LAYOUT_ROW_MAJOR, K[rblock],
                        dim[rblock], dim[rblock], 0., temp, dim[rblock]);
        delete[] val;
        delete[] row_index;
        delete[] col_index;
        mkl_sparse_destroy(cdold);

        double *cd = new double[dim[lblock] * dim[rblock]];
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, dim[lblock], dim[rblock], dim[lblock], 1., K[lblock],
                    dim[lblock], temp, dim[rblock], 0., cd, dim[rblock]);
        delete[] temp;
        for (int lindex = 0; lindex < dim[lblock]; lindex++)
        {
            for (int rindex = 0; rindex < dim[rblock]; rindex++)
            {
                cd[lindex * dim[rblock] + rindex] *= cd[lindex * dim[rblock] + rindex];
                for (int n = 0; n < sz; n++)
                {
                    Complex coeff = (w[lblock][lindex] + w[rblock][rindex]) / (iw[n] + egvalues[lblock][lindex] - egvalues[rblock][rindex]);
                    ret[n] -= coeff * cd[lindex * dim[rblock] + rindex];
                }
            }
        }
        delete[] cd;
    }

    for (int n = 0; n < sz; n++)
        ret[n] /= Z;
    return ret;
}

std::vector<double> ED::greensfunction(int k, std::vector<double> ttab) const
{
    int sz = ttab.size();
    std::vector<double> ret(sz, 0);
    for (int rblock = 0; rblock < NNN; rblock++)
    {
        if (dim[rblock] == 0)
            continue;
        int rp = rblock % Ns;
        int lblock = rblock + NN + (rp + k) % Ns - rp;
        if (lblock >= NNN || dim[lblock] == 0)
            continue;
        //write cd matrix
        int nnz = 0;
        double *val = new double[dim[rblock]];
        int *row_index = new int[dim[rblock]];
        int *col_index = new int[dim[rblock]];
        for (int rindex = 0; rindex < dim[rblock]; rindex++)
        {
            int rstate = base[rblock][rindex];
            int rupstate = rstate % onespindim;
            if ((rstate >> k) & 1)
                continue;
            int lstate = rstate ^ (1 << k);
            if ((ones[rupstate] - ones[rupstate >> k]) & 1)
                val[nnz] = -1;
            else
                val[nnz] = 1;
            row_index[nnz] = index[lstate];
            col_index[nnz] = rindex;
            nnz++;
        }
        if (nnz == 0)
            continue;
        sparse_matrix_t cdold;
        mkl_sparse_d_create_coo(&cdold, SPARSE_INDEX_BASE_ZERO, dim[lblock], dim[rblock], nnz, row_index, col_index, val);
        struct matrix_descr descr = {SPARSE_MATRIX_TYPE_GENERAL, SPARSE_FILL_MODE_UPPER, SPARSE_DIAG_NON_UNIT};
        double *temp = new double[dim[lblock] * dim[rblock]];
        mkl_sparse_d_mm(SPARSE_OPERATION_NON_TRANSPOSE, 1., cdold, descr, SPARSE_LAYOUT_ROW_MAJOR, K[rblock],
                        dim[rblock], dim[rblock], 0., temp, dim[rblock]);
        delete[] val;
        delete[] row_index;
        delete[] col_index;
        mkl_sparse_destroy(cdold);

        double *cd = new double[dim[lblock] * dim[rblock]];
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, dim[lblock], dim[rblock], dim[lblock], 1., K[lblock],
                    dim[lblock], temp, dim[rblock], 0., cd, dim[rblock]);
        delete[] temp;
        for (int lindex = 0; lindex < dim[lblock]; lindex++)
        {
            for (int rindex = 0; rindex < dim[rblock]; rindex++)
            {
                cd[lindex * dim[rblock] + rindex] *= cd[lindex * dim[rblock] + rindex];
                for (int i = 0; i < sz; i++)
                {
                    double coeff = w[lblock][lindex] * std::exp(ttab[i] * (egvalues[lblock][lindex] - egvalues[rblock][rindex]));
                    ret[i] += coeff * cd[lindex * dim[rblock] + rindex];
                }
            }
        }
        delete[] cd;
    }

    for (int i = 0; i < sz; i++)
        ret[i] /= Z;
    return ret;
}

std::vector<ED::Complex> ED::spincorr(int idx, std::vector<int> ntab) const
{
    int sz = ntab.size();
    std::vector<Complex> iw(sz);
    std::vector<Complex> ret(sz, Complex(0, 0));
    for (int n = 0; n < sz; n++)
        iw[n] = I * pi * T * (2. * ntab[n]);
    for (int q = 0; q < Ns; q++)
    {
        for (int lblock = 0; lblock < NNN; lblock++)
        {
            int lp = lblock % Ns;
            int rblock = lblock - lp + (lp + q) % Ns;
            if (dim[lblock] == 0 || dim[rblock] == 0)
                continue;

            //Matrix1 dim[rblock] * dim[lblock]
            double *M1 = new double[dim[rblock] * dim[lblock]];
            //Matrix2 dim[lblock] * dim[rblock]
            double *M2 = new double[dim[lblock] * dim[rblock]];
            for (int i = 0; i < dim[rblock] * dim[lblock]; i++)
                M1[i] = 0, M2[i] = 0;

            if (q != 0)
            {
                for (int lindex = 0; lindex < dim[lblock]; lindex++)
                {
                    int lstate = base[lblock][lindex];
                    for (int k2 = 0; k2 < Ns; k2++)
                    {
                        int k1 = (k2 + q) % Ns;
                        if (((lstate >> k1) & 1) || !((lstate >> k2) & 1))
                            continue;
                        int rstate = (lstate ^ (1 << k1)) ^ (1 << k2);
                        if (onesbetween(lstate % onespindim, k1, k2) & 1)
                            M1[index[rstate] * dim[lblock] + lindex] -= cos(2. * pi * q * idx / Ns);
                        else
                            M1[index[rstate] * dim[lblock] + lindex] += cos(2. * pi * q * idx / Ns);
                    }
                }

                for (int rindex = 0; rindex < dim[rblock]; rindex++)
                {
                    int rstate = base[rblock][rindex];
                    int rupstate = rstate % onespindim, rdownstate = rstate / onespindim;
                    for (int k4 = 0; k4 < Ns; k4++)
                    {
                        int k3 = (Ns + k4 - q) % Ns;

                        //sigma = up
                        if (!((rupstate >> k3) & 1) && ((rupstate >> k4) & 1))
                        {
                            int lupstate = (rupstate ^ (1 << k3)) ^ (1 << k4);
                            int lstate = rdownstate * onespindim + lupstate;
                            if (onesbetween(rupstate, k3, k4) & 1)
                                M2[index[lstate] * dim[rblock] + rindex] -= 1.;
                            else
                                M2[index[lstate] * dim[rblock] + rindex] += 1.;
                        }

                        //sigma = down
                        if (!((rdownstate >> k3) & 1) && ((rdownstate >> k4) & 1))
                        {
                            int ldownstate = (rdownstate ^ (1 << k3)) ^ (1 << k4);
                            int lstate = ldownstate * onespindim + rupstate;
                            if (onesbetween(rdownstate, k3, k4) & 1)
                                M2[index[lstate] * dim[rblock] + rindex] += 1.;
                            else
                                M2[index[lstate] * dim[rblock] + rindex] -= 1.;
                        }
                    }
                }
            }
            else //(q=0)diagonal
            {
                for (int index = 0; index < dim[lblock]; index++)
                {
                    int state = base[lblock][index];
                    int downstate = state / onespindim, upstate = state % onespindim;

                    for (int k = 0; k < Ns; k++)
                    {
                        M1[index * dim[lblock] + index] += (upstate >> k) & 1;
                        M2[index * dim[lblock] + index] += (upstate >> k) & 1;
                        M2[index * dim[lblock] + index] -= (downstate >> k) & 1;
                    }
                }
            }

            double *temp = new double[dim[lblock] * dim[rblock]];
            for (int i = 0; i < dim[lblock] * dim[rblock]; i++)
                temp[i] = 0;

            cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, dim[rblock], dim[lblock], dim[rblock], 1., K[rblock],
                        dim[rblock], M1, dim[lblock], 0., temp, dim[lblock]);
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim[rblock], dim[lblock], dim[lblock], 1., temp,
                        dim[lblock], K[lblock], dim[lblock], 0., M1, dim[lblock]);

            for (int i = 0; i < dim[lblock] * dim[rblock]; i++)
                temp[i] = 0;
            cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, dim[lblock], dim[rblock], dim[lblock], 1., K[lblock],
                        dim[lblock], M2, dim[rblock], 0., temp, dim[rblock]);
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim[lblock], dim[rblock], dim[rblock], 1., temp,
                        dim[rblock], K[rblock], dim[rblock], 0., M2, dim[rblock]);

            for (int rindex = 0; rindex < dim[rblock]; rindex++)
            {
                for (int lindex = 0; lindex < dim[lblock]; lindex++)
                {
                    for (int n = 0; n < sz; n++)
                    {
                        if (ntab[n] == 0 && abs(egvalues[rblock][rindex] - egvalues[lblock][lindex]) < 1e-9)
                            ret[n] += M1[rindex * dim[lblock] + lindex] * M2[lindex * dim[rblock] + rindex] *
                                      (w[rblock][rindex]) / T;
                        else
                            ret[n] += M1[rindex * dim[lblock] + lindex] * M2[lindex * dim[rblock] + rindex] *
                                      (w[lblock][lindex] - w[rblock][rindex]) / (iw[n] + egvalues[rblock][rindex] - egvalues[lblock][lindex]);
                    }
                }
            }

            delete[] M1;
            delete[] M2;
            delete[] temp;
        }
    }
    for (int n = 0; n < sz; n++)
        ret[n] *= 2. / (Z * Ns * Ns);
    return ret;
}

std::vector<double> ED::spectralfunction(int k, std::vector<double> omegatab) const
{
    int sz = omegatab.size();
    static double eta = 0.005;
    std::vector<double> ret(sz, 0);
    for (int rblock = 0; rblock < NNN; rblock++)
    {
        if (dim[rblock] == 0)
            continue;
        int rp = rblock % Ns;
        int lblock = rblock + NN + (rp + k) % Ns - rp;
        if (lblock >= NNN || dim[lblock] == 0)
            continue;
        //write cd matrix
        int nnz = 0;
        double *val = new double[dim[rblock]];
        int *row_index = new int[dim[rblock]];
        int *col_index = new int[dim[rblock]];
        for (int rindex = 0; rindex < dim[rblock]; rindex++)
        {
            int rstate = base[rblock][rindex];
            int rupstate = rstate % onespindim;
            if ((rstate >> k) & 1)
                continue;
            int lstate = rstate ^ (1 << k);
            if ((ones[rupstate] - ones[rupstate >> k]) & 1)
                val[nnz] = -1;
            else
                val[nnz] = 1;
            row_index[nnz] = index[lstate];
            col_index[nnz] = rindex;
            nnz++;
        }
        if (nnz == 0)
            continue;
        sparse_matrix_t cdold;
        mkl_sparse_d_create_coo(&cdold, SPARSE_INDEX_BASE_ZERO, dim[lblock], dim[rblock], nnz, row_index, col_index, val);
        struct matrix_descr descr = {SPARSE_MATRIX_TYPE_GENERAL, SPARSE_FILL_MODE_UPPER, SPARSE_DIAG_NON_UNIT};
        double *temp = new double[dim[lblock] * dim[rblock]];
        mkl_sparse_d_mm(SPARSE_OPERATION_NON_TRANSPOSE, 1., cdold, descr, SPARSE_LAYOUT_ROW_MAJOR, K[rblock],
                        dim[rblock], dim[rblock], 0., temp, dim[rblock]);
        delete[] val;
        delete[] row_index;
        delete[] col_index;
        mkl_sparse_destroy(cdold);

        double *cd = new double[dim[lblock] * dim[rblock]];
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, dim[lblock], dim[rblock], dim[lblock], 1., K[lblock],
                    dim[lblock], temp, dim[rblock], 0., cd, dim[rblock]);
        delete[] temp;
        for (int lindex = 0; lindex < dim[lblock]; lindex++)
        {
            for (int rindex = 0; rindex < dim[rblock]; rindex++)
            {
                cd[lindex * dim[rblock] + rindex] *= cd[lindex * dim[rblock] + rindex];
#pragma omp parallel for
                for (int n = 0; n < sz; n++)
                {
                    double pole = -omegatab[n] + egvalues[lblock][lindex] - egvalues[rblock][rindex];
                    ret[n] += (w[lblock][lindex] + w[rblock][rindex]) * cd[lindex * dim[rblock] + rindex] *
                              2. * eta / (pole * pole + eta * eta);
                }
            }
        }
        delete[] cd;
    }
    for (int n = 0; n < sz; n++)
        ret[n] /= Z;
    return ret;
}