#include "QuarticApproximation2D.h"
#include "mkl.h"
#include <iostream>

const QuarticApproximation2D::Complex QuarticApproximation2D::I(0, 1);
const double QuarticApproximation2D::pi = 3.1415926535897932385;

QuarticApproximation2D::QuarticApproximation2D(int Nt_half, int Nx_, int Ny_, double T_, double t_, double U_, double delmu_)
    : Nt(Nt_half * 2), Nx(Nx_), Ny(Ny_), Ns(Nx * Ny), Nts(Nt * Ns), NN(Nts * Nts),
      T(T_), t(t_), U(U_), delmu(delmu_), f_5(0.5 * U * T / Ns)
{
    iw = new Complex[Nt * 3];
    h = new Complex[Nts * 3];
    g = new Complex[Nts];

    c = new Complex[NN];
    d1 = new Complex[NN];
    d2 = new Complex[NN];
    r = new Complex[NN];

    xid = new int[3 * Nts];
    yid = new int[3 * Nts];
    tid = new int[3 * Nts];

    for (int n = 0; n < 2 * Nt; n++)
        iw[n] = I * pi * T * (2. * (n - Nt_half) + 1);
    for (int n = 2 * Nt; n < 3 * Nt; n++)
        iw[n] = I * pi * T * (2. * (n - 3 * Nt - Nt_half) + 1);

    for (int alpha = 0; alpha < 3 * Nts; alpha++)
    {
        tid[alpha] = alpha / Ns;
        xid[alpha] = (alpha % Ns) / Ny;
        yid[alpha] = alpha % Ny;
    }

    for (int alpha = 0; alpha < Nts; alpha++)
        g[alpha] = -1. / (iw[tid[alpha]] - 2. * t * (cos(2. * pi * xid[alpha] / Nx) + cos(2. * pi * yid[alpha] / Ny)) - delmu);

    calculate();
}

void QuarticApproximation2D::calculate()
{
    double p;
    if (U / T < 8)
        p = 0.3;
    else if (U / T >= 8 && U / T < 15)
        p = 0.1;
    else
        p = 0.05;

    Complex *dd = new Complex[Nts];
    Complex *nc = new Complex[NN];
    Complex *nd1 = new Complex[NN];
    Complex *nd2 = new Complex[NN];
    Complex *nr = new Complex[NN];

    Complex *d11conv = new Complex[2 * NN];
    Complex *d12conv = new Complex[2 * NN];
    Complex *d21conv = new Complex[2 * NN];
    Complex *d22conv = new Complex[2 * NN];
    Complex *d11corr = new Complex[2 * NN];
    Complex *d12corr = new Complex[2 * NN];
    Complex *d21corr = new Complex[2 * NN];
    Complex *d22corr = new Complex[2 * NN];
    Complex *c1corr = new Complex[2 * NN];
    Complex *r1corr = new Complex[2 * NN];

    Complex *hco = new Complex[3 * Nt * (2 * Nx - 1) * (2 * Ny - 1)];
    Complex *hgconv = new Complex[2 * Nts];
    Complex *hgcorr = new Complex[2 * Nts];

    int hcospacedim = (2 * Nx - 1) * (2 * Ny - 1);

    int dims = 3;
    int hshape[] = {3 * Nt, 2 * Nx - 1, 2 * Ny - 1};
    int chainshape[] = {Nt, Nx, Ny};
    int resultshape[] = {2 * Nt, Nx, Ny};
    int hstride[] = {hcospacedim, 2 * Ny - 1, 1};
    int chain1stride[] = {Nts * Ns, Nts * Ny, Nts};
    int chain2stride[] = {Ns, Ny, 1};
    int resultstride[] = {Ns, Ny, 1};
    int convstart[] = {Nt, Nx - 1, Ny - 1};
    int corrstart[] = {-2 * Nt, -Nx + 1, -Ny + 1};

    int gshape[] = {Nt, Nx, Ny};
    int gstride[] = {Ns, Ny, 1};

    int Nit = 2000;
    for (int it = 0; it < Nit; it++)
    {
        std::cout << greensfunction(0, 0, 0) << std::endl;
        Complex nhalf = 0;
        for (int alpha = 0; alpha < 3 * Nts; alpha++)
            nhalf += g[alpha];
        nhalf *= T / Ns;

        for (int alpha = 0; alpha < 3 * Nts; alpha++)
            h[alpha] = -1. / (iw[tid[alpha]] - 2. * t * (cos(2. * pi * xid[alpha] / Nx) + cos(2. * pi * yid[alpha] / Ny)) - delmu /*- U / 2.*/ + U * nhalf);
        //std::cout << nhalf << std::endl;

        //write h to hco
        for (int ix = 0; ix < 2 * Nx - 1; ix++)
            for (int iy = 0; iy < 2 * Ny - 1; iy++)
                for (int n = 0; n < 3 * Nt; n++)
                    hco[n * hcospacedim + ix * (2 * Ny - 1) + iy] = h[(n + 2 * Nt) % (3 * Nt) * Ns + (ix + 1) % Nx * Ny + (iy + 1) % Ny];

        VSLConvTaskPtr convtask;
        vslzConvNewTaskX(&convtask, VSL_CONV_MODE_FFT, dims, hshape, chainshape, resultshape,
                         reinterpret_cast<MKL_Complex16 *>(hco), hstride);
        VSLCorrTaskPtr corrtask;
        vslzCorrNewTaskX(&corrtask, VSL_CONV_MODE_FFT, dims, hshape, chainshape, resultshape,
                         reinterpret_cast<MKL_Complex16 *>(hco), hstride);

        vslConvSetStart(convtask, convstart);
        vslCorrSetStart(corrtask, corrstart);

        vslzConvExecX(convtask, reinterpret_cast<MKL_Complex16 *>(g), gstride,
                      reinterpret_cast<MKL_Complex16 *>(hgconv), resultstride);
        vslzCorrExecX(corrtask, reinterpret_cast<MKL_Complex16 *>(g), gstride,
                      reinterpret_cast<MKL_Complex16 *>(hgcorr), resultstride);

#pragma omp parallel for
        for (int gamma = 0; gamma < Nts; gamma++)
        {
            vslzConvExecX(convtask, reinterpret_cast<MKL_Complex16 *>(&d1[gamma]), chain1stride,
                          reinterpret_cast<MKL_Complex16 *>(&d11conv[gamma * 2 * Nts]), resultstride);
            vslzConvExecX(convtask, reinterpret_cast<MKL_Complex16 *>(&d1[gamma * Nts]), chain2stride,
                          reinterpret_cast<MKL_Complex16 *>(&d12conv[gamma * 2 * Nts]), resultstride);
            vslzConvExecX(convtask, reinterpret_cast<MKL_Complex16 *>(&d2[gamma]), chain1stride,
                          reinterpret_cast<MKL_Complex16 *>(&d21conv[gamma * 2 * Nts]), resultstride);
            vslzConvExecX(convtask, reinterpret_cast<MKL_Complex16 *>(&d2[gamma * Nts]), chain2stride,
                          reinterpret_cast<MKL_Complex16 *>(&d22conv[gamma * 2 * Nts]), resultstride);
            vslzCorrExecX(corrtask, reinterpret_cast<MKL_Complex16 *>(&d1[gamma]), chain1stride,
                          reinterpret_cast<MKL_Complex16 *>(&d11corr[gamma * 2 * Nts]), resultstride);
            vslzCorrExecX(corrtask, reinterpret_cast<MKL_Complex16 *>(&d1[gamma * Nts]), chain2stride,
                          reinterpret_cast<MKL_Complex16 *>(&d12corr[gamma * 2 * Nts]), resultstride);
            vslzCorrExecX(corrtask, reinterpret_cast<MKL_Complex16 *>(&d2[gamma]), chain1stride,
                          reinterpret_cast<MKL_Complex16 *>(&d21corr[gamma * 2 * Nts]), resultstride);
            vslzCorrExecX(corrtask, reinterpret_cast<MKL_Complex16 *>(&d2[gamma * Nts]), chain2stride,
                          reinterpret_cast<MKL_Complex16 *>(&d22corr[gamma * 2 * Nts]), resultstride);
            vslzCorrExecX(corrtask, reinterpret_cast<MKL_Complex16 *>(&c[gamma]), chain1stride,
                          reinterpret_cast<MKL_Complex16 *>(&c1corr[gamma * 2 * Nts]), resultstride);
            vslzCorrExecX(corrtask, reinterpret_cast<MKL_Complex16 *>(&r[gamma]), chain1stride,
                          reinterpret_cast<MKL_Complex16 *>(&r1corr[gamma * 2 * Nts]), resultstride);
        }

        vslConvDeleteTask(&convtask);
        vslCorrDeleteTask(&corrtask);
#pragma omp parallel for
        for (int beta = 0; beta < Nts; beta++)
        {
            for (int gamma = 0; gamma < Nts; gamma++)
            {
                Complex gbetagamma = g[beta] * g[gamma];
                int betaplusgamma = (tid[beta] + tid[gamma]) * Ns + (xid[beta] + xid[gamma]) % Nx * Ny + (yid[beta] + yid[gamma]) % Ny;
                int betaminusgamma = (tid[beta] - tid[gamma] + Nt) * Ns + (xid[beta] - xid[gamma] + Nx) % Nx * Ny + (yid[beta] - yid[gamma] + Ny) % Ny;
                int gammaminusbeta = (tid[gamma] - tid[beta] + Nt) * Ns + (xid[gamma] - xid[beta] + Nx) % Nx * Ny + (yid[gamma] - yid[beta] + Ny) % Ny;
                //c
                nc[beta * Nts + gamma] = -4. * U * gbetagamma * hgconv[betaplusgamma];
                nc[beta * Nts + gamma] -= g[beta] * (d21conv[gamma * 2 * Nts + betaplusgamma] - d11conv[gamma * 2 * Nts + betaplusgamma]);
                nc[beta * Nts + gamma] -= g[gamma] * (d21conv[beta * 2 * Nts + betaplusgamma] - d11conv[beta * 2 * Nts + betaplusgamma]);
                nc[beta * Nts + gamma] = f_5 * nc[beta * Nts + gamma] / (1. + 2. * f_5 * hgconv[betaplusgamma]);

                //d1
                nd1[beta * Nts + gamma] = -g[beta] * (d11corr[gamma * 2 * Nts + gammaminusbeta] + d21corr[gamma * 2 * Nts + gammaminusbeta]);
                nd1[beta * Nts + gamma] -= g[gamma] * (d12corr[beta * 2 * Nts + betaminusgamma] + d22corr[beta * 2 * Nts + betaminusgamma]);
                nd1[beta * Nts + gamma] += (hgcorr[gammaminusbeta] + hgcorr[betaminusgamma]) * d2[beta * Nts + gamma];
                nd1[beta * Nts + gamma] = f_5 * nd1[beta * Nts + gamma] / (1. - f_5 * (hgcorr[gammaminusbeta] + hgcorr[betaminusgamma]));

                //d2
                nd2[beta * Nts + gamma] = -2. * U * gbetagamma * (hgcorr[betaminusgamma] + hgcorr[gammaminusbeta]);
                nd2[beta * Nts + gamma] += g[beta] * d11corr[gamma * 2 * Nts + gammaminusbeta];
                nd2[beta * Nts + gamma] += g[gamma] * d12corr[beta * 2 * Nts + betaminusgamma];
                nd2[beta * Nts + gamma] -= g[beta] * c1corr[gamma * 2 * Nts + betaminusgamma];
                nd2[beta * Nts + gamma] -= g[gamma] * r1corr[beta * 2 * Nts + gammaminusbeta];
                nd2[beta * Nts + gamma] = f_5 * nd2[beta * Nts + gamma] / (1. + f_5 * (hgcorr[gammaminusbeta] + hgcorr[betaminusgamma]));

                //r
                nr[beta * Nts + gamma] = -4. * U * gbetagamma * hgconv[betaplusgamma];
                nr[beta * Nts + gamma] -= g[beta] * (d22conv[gamma * 2 * Nts + betaplusgamma] - d12conv[gamma * 2 * Nts + betaplusgamma]);
                nr[beta * Nts + gamma] -= g[gamma] * (d22conv[beta * 2 * Nts + betaplusgamma] - d12conv[beta * 2 * Nts + betaplusgamma]);
                nr[beta * Nts + gamma] = f_5 * nr[beta * Nts + gamma] / (1. + 2. * f_5 * hgconv[betaplusgamma]);
            }
        }

        double chain_diff = 0.;
        for (int i = 0; i < NN; i++)
        {
            chain_diff = std::max(chain_diff, std::abs(nc[i] - c[i]));
            chain_diff = std::max(chain_diff, std::abs(nd1[i] - d1[i]));
            chain_diff = std::max(chain_diff, std::abs(nd2[i] - d2[i]));
            chain_diff = std::max(chain_diff, std::abs(nr[i] - r[i]));
        }

        for (int n = 0; n < NN; n++)
            c[n] = nc[n], d1[n] = nd1[n], d2[n] = nd2[n], r[n] = nr[n];

        for (int alpha = 0; alpha < Nts; alpha++)
        {
            dd[alpha] = 0;
            for (int nu = 0; nu < Nts; nu++)
                dd[alpha] += r[Nts * alpha + nu] + d1[Nts * alpha + nu] + 2. * d2[Nts * alpha + nu];
            dd[alpha] = (T / (6. * Ns)) * dd[alpha];
        }

        //get new g[n]
        double green_diff = 0.;
        for (int alpha = 0; alpha < Nts; alpha++)
        {
            //g[alpha] = h[alpha]*(1.+dd[alpha]);
            Complex newg_alpha = h[alpha] / (1. - h[alpha] * dd[alpha] / g[alpha]);
            //dd[alpha] = h[alpha] / (1. - h[alpha] * dd[alpha] / g[alpha]);
            newg_alpha = p * newg_alpha + (1. - p) * g[alpha];
            green_diff = std::max(green_diff, std::abs(newg_alpha - g[alpha]));
            g[alpha] = newg_alpha;
        }

        //std::cout << "g(0,0) = " << g[Nt / 2 * Ns].real() << "+I*" << g[Nt / 2 * Ns].imag() << std::endl;
        if (chain_diff < 1e-6 && green_diff < 1e-7)
            break;
    }

    delete[] hco;
    delete[] hgconv;
    delete[] hgcorr;
    delete[] dd;
    delete[] nc;
    delete[] nd1;
    delete[] nd2;
    delete[] nr;
    delete[] d11conv;
    delete[] d12conv;
    delete[] d21conv;
    delete[] d22conv;
    delete[] d11corr;
    delete[] d12corr;
    delete[] d21corr;
    delete[] d22corr;
    delete[] c1corr;
    delete[] r1corr;
}

void QuarticApproximation2D::setT(double T_)
{
    T = T_;
    f_5 = 0.5 * U * T / Ns;
    calculate();
}

void QuarticApproximation2D::sett(double t_)
{
    t = t_;
    calculate();
}

void QuarticApproximation2D::setdelmu(double delmu_)
{
    delmu = delmu_;
    calculate();
}

void QuarticApproximation2D::setparms(double U_, double delmu_)
{
    U = U_;
    delmu = delmu_;
    calculate();
}

QuarticApproximation2D::Complex QuarticApproximation2D::greensfunction(int kx, int ky, int n) const
{
    if (n > Nt / 2 - 1 || n < -Nt / 2)
        return 0;
    kx %= Nx, ky %= Ny;
    return g[(n + Nt / 2) * Ns + kx * Ny + ky];
}

QuarticApproximation2D::~QuarticApproximation2D()
{
    delete[] iw;
    delete[] h;
    delete[] g;
    delete[] c;
    delete[] d1;
    delete[] d2;
    delete[] r;
    delete[] tid;
    delete[] xid;
    delete[] yid;
}