#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "ED.h"
using namespace std;

int main() {
    int Ns = 8;
    double t = 1., U = 8.;
    ED mysystem(Ns, t, U, 1., U / 2.);

    double dmu = 0.025, eps = 1e-8;
    int nmu = int((U / 2 - 0. + eps) / dmu) + 1;
    vector<double> mutab(nmu);
    for (int i = 0; i < nmu; i++) {
        mutab[i] = i * dmu;
    }

    vector<double> betatab = {1, 2, 4, 6, 8, 16, 32};

    double dtau = 0.05;
    for (int bid = 0; bid < betatab.size(); bid++) {
        double beta = betatab[bid];
        int nt = (beta + eps) / dtau + 1;
        vector<double> ttab(nt);
        for (int i = 0; i < nt; i++) {
            ttab[i] = i * dtau;
        }

        ostringstream dirname("data/", ios::ate);
        dirname << "Ns" << Ns << "/beta" << beta << "/U" << U << "/";
        int status = system(("mkdir -p " + dirname.str()).c_str());
        if (status) {
            cerr << "Can not make diretory \"" << dirname.str() << "\"" << std::endl;
            exit(status);
        }

        ofstream denfile((dirname.str() + "DenVSmu.dat").c_str());
        if (!denfile) {
            cerr << "Can open file \"" << dirname.str() << "DenVSmu.dat"
                 << "\"" << std::endl;
            exit(1);
        }
        denfile.setf(ios::fixed);

        mysystem.setT(1. / beta);
        for (int muid = 0; muid < nmu; muid++) {
            mysystem.setmu(mutab[muid]);
            denfile << setprecision(3) << setw(5) << mutab[muid] << "\t" << setprecision(8) << setw(10)
                    << mysystem.density() << endl;

            if (muid % 20 == 0) {
                ostringstream subdir(dirname.str(), ios::ate);
                subdir << "mu" << mutab[muid] << "/";
                status = system(("mkdir -p " + subdir.str()).c_str());
                if (status) {
                    cerr << "Can not make diretory \"" << subdir.str() << "\"" << std::endl;
                    exit(status);
                }

                for (int k = 0; k < Ns; k++) {
                    ostringstream gfname(subdir.str(), ios::ate);
                    gfname << "GFtauk" << k << ".dat";
                    ofstream gfile(gfname.str().c_str());
                    if (!gfname) {
                        cerr << "Can open file \"" << gfname.str() << "\"" << std::endl;
                        exit(1);
                    }

                    for (int i = 0; i < ttab.size(); i++) {
                        gfile << resetiosflags(ios::floatfield) << setiosflags(ios::fixed) << setprecision(2) << setw(4)
                              << ttab[i] << "\t" << resetiosflags(ios::floatfield) << setiosflags(ios::scientific)
                              << setprecision(6) << setw(13) << mysystem.Gtau(k, ttab[i]) << std::endl;
                    }
                    gfile.close();
                }
            }
        }

        denfile.close();
    }
}