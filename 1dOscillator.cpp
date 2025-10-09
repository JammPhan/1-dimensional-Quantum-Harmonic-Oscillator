#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <fstream>

using Eigen::MatrixXd;
using namespace Eigen;
using namespace std;

int main() {

    string output = "wavefunctions.csv";
    
    const int N = 200;
    const double xrange = 8;
    const double hbar = 1.0;
    const double m = 1.0;
    const double k = 1.0;
    const double dx = (xrange*2)/(N-1);
    const double omega = 1.0;

    vector<double> x(N), V(N);

    // Potential Well
    for (int i = 0; i < N; i++) {
        x[i] = -xrange + i*dx;
        
        V[i] = 0.5 * x[i] * x[i] * m * omega * omega;

    }

    // Set up Hamiltonian

    MatrixXd H = MatrixXd::Zero(N,N);

    const double gamma = (hbar * hbar) / (2 * m * dx * dx);

    for (int i = 0; i < N; i++) {

        H(i,i) = V[i] + 2 * gamma;

        if (i > 0) {H(i, i-1) = -gamma;}
        if (i < N-1) { H(i, i+1) = -gamma;}
        
    }

    Eigen::SelfAdjointEigenSolver<MatrixXd> solver(H);

    if (solver.info() != Eigen::Success) {
        cout << "Failed to Compute Eigenvalues";
        return 1;
    }

    VectorXd energies = solver.eigenvalues();

    MatrixXd states = solver.eigenvectors();

    cout << setprecision(10);
    cout << "Lowest 6 energies:\n";
    for (int n = 0; n < 6; ++n) {
        // Analytical solution
        double E_analytic = hbar * omega * (n + 0.5);
        // relative error of numerical solution
        double rel_err = (energies[n] - E_analytic) / E_analytic;
        cout << "n=" << n << "  E_num=" << energies[n] << "  E_exact=" << E_analytic << "  rel_err=" << rel_err << endl;
    }

    int nwaves = 5; // number of wavestates to export

    ofstream fout(output);
    fout << "x";
    for (int n = 0; n < nwaves; ++n) fout << ",psi" << n;
    fout << "\n";


    for (int i = 0; i < N; i++) {
        fout << setprecision(10) << x[i];

        for (int n = 0; n < nwaves; n++) {
            VectorXd psi = states.col(n);
            double norm = sqrt(psi.array().square().sum() * dx);
            double psi_cont = states(i,n) / norm;
            fout << "," << psi_cont;
        }
        fout << "\n";
    }

    fout.close();
    std::cout << "Wrote " << output << " with first " << nwaves << " wavefunctions.\n";


    return 0;
}