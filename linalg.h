
#ifndef LINALG_H
#define LINALG_H

#include <Eigen/Dense>
#include <tuple>

typedef std::complex<double> cd;

extern "C" {
    void zggev_( // generalized eigen-solver for general matrices
                 // may be slower than generalized-self-adjoint!!
       char* JOBVL, char* JOBVR, int* N,
       cd* A ,    int* LDA, 
       cd* B,     int* LDB,
       cd* ALPHA, cd* BETA,
       cd* VL,    int* LDVL, 
       cd* VR,    int* LDVR,
       cd* WORK,  int* LWORK, double* RWORK, 
       int* INFO);
}

class Eigensolver {
    public:
        Eigensolver (double, Eigen::VectorXcd);
        double eigenvalue;
        Eigen::VectorXcd eigenvector;
};

Eigensolver::Eigensolver (double d, Eigen::VectorXcd v){
    eigenvalue = d;
    eigenvector = v;
}

void to_fortran_matrix(Eigen::MatrixXcd &A, cd* B)
{
    int N = A.rows();
    assert(A.rows() == A.cols());
    for(int i = 0; i < N; ++i)
        for(int j = 0; j < N; ++j)
            B[j*N + i] = A(i,j);
}

Eigensolver
lowest_eigenpair
(Eigen::MatrixXcd &H, Eigen::MatrixXcd &N, int MAX)
{

    cd f_H[MAX][MAX] = {0};
    cd f_N[MAX][MAX] = {0};

    to_fortran_matrix(H, &f_H[0][0]);
    to_fortran_matrix(N, &f_N[0][0]);

    /////////////////////// BEGIN SET DEFAULTS ///////////////////////

    int lda   = MAX;
    int ldb   = MAX;
    int ldc   = MAX;
    int ldvl  = MAX;
    int ldvr  = MAX;
    int lwork = MAX*3;
    int info  = 0;

    char nn = 'V';
    char yy = 'V';

    cd VL[MAX][MAX]; memset(VL, 0, sizeof(cd)*MAX*MAX);
    cd VR[MAX][MAX]; memset(VR, 0, sizeof(cd)*MAX*MAX);

    cd WORK[lwork];      memset (WORK, 0, sizeof(cd)*lwork);
    double RWORK[8*MAX]; memset (RWORK,0, sizeof(double)*8*MAX);

    cd alpha[MAX] = {0};
    cd beta [MAX] = {0};

    ///////////////////////  END  SET DEFAULTS ///////////////////////

    int n = H.rows();
    assert (H.rows() == N.rows());
    assert (H.cols() == N.cols());
    assert (H.rows() == H.cols());

    zggev_(&yy, &yy, &n,
           &f_H[0][0], &lda,
           &f_N[0][0], &ldb,
           &alpha[0], &beta[0],
           &VL[0][0], &ldvl,
           &VR[0][0], &ldvr,
           &WORK[0], &lwork, &RWORK[0],
           &info);

    Eigen::VectorXcd eigenvec(n);
    for(int i = 0; i < n; ++i)
        eigenvec(i) = VR[0][i];

    std::cout << "Here are the eigenvalues generated (alpha,beta): ";
    for(int i = 0; i < n; ++i)
        std::cout << "(" << alpha[i] << "," << beta[i] << "), ";
    std::cout << std::endl;

    return Eigensolver( std::real (alpha[0]/beta[0]), eigenvec );
}

#endif

