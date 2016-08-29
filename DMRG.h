
#ifndef MPS_DMRG_H
#define MPS_DMRG_H

#include "tensor.h"
#include "DMRG_contractions.h"

#include <Eigen/Dense>
#include <vector>
#include <array>
#include <complex>
#include <functional>
#include <utility>
#include "linalg.h"

tensor<cd,2> delta_tensor(int a, int b)
{
    typedef std::vector<int> vi;
    tensor<cd,2> delta(vi({{a,b}}));
    for(int i = 0; i < a; ++i) {
        for(int j = 0; j < b; ++j) {
            if (i == j)
                delta[i][j] = 1;
            else
                delta[i][j] = 0;
        }
    }
    return delta;
}

template<typename T>
Eigen::MatrixXcd
DMRG_matrixN(tensor<T,2> left, tensor<T,2> right, tensor<T,4> localHamiltonian)
{
    int a = localHamiltonian.shape()[1];
    int b = localHamiltonian.shape()[2];
    std::array<int,0> empty;
    tensor<T,4> T1 = contract(left,delta_tensor(a,b),empty,empty);
    tensor<T,6> T2 = contract(T1,  right,            empty,empty); // (l0 l1 d0 d1 r0 r1)
//std::array<int,3> rows = {{ T2.shape()[0], T2.shape()[2], T2.shape()[4] }}; // why are these the rows? shouldn't it be {0,2,4} & {135} instead of T2.shape()[..]?
//std::array<int,3> cols = {{ T2.shape()[1], T2.shape()[3], T2.shape()[5] }};
    std::array<int,3> rows = {{ 0,2,4 }};
    std::array<int,3> cols = {{ 1,3,5 }};

    return tensor_to_matrix(T2,rows,cols);
}

template<typename T>
Eigen::MatrixXcd
DMRG_matrixH(tensor<T,3> left, tensor<T,3> right, tensor<T,4> localHamiltonian)
{
    tensor<T,5> T1 = contract(left,localHamiltonian,ar::one,ar::zero); // (l0 l1 l2) (h0 h1 h2 h3)  ----->  (l0 l2 h1 h2 h3)
    tensor<T,6> T2 = contract(T1, right,            ar::four,ar::one);
// (l0 l2 h1 h2 h3) (r0 r1 r2) --->  (l0 l2 h1 h2 r0 r2) be very careful of
// the order of the middle two: this is analogous to (l0 l1 d1 d0 r0 r1)

//std::array<int,3> rows = {{ T2.shape()[0], T2.shape()[3], T2.shape()[4] }}; // why are these the rows? shouldn't it be {0,2,4} & {135} instead of T2.shape()[..]?
//std::array<int,3> cols = {{ T2.shape()[1], T2.shape()[2], T2.shape()[5] }};
    std::array<int,3> rows = {{ 0,3,4 }};
    std::array<int,3> cols = {{ 1,2,5 }};

    return tensor_to_matrix(T2,rows,cols);
}

tensor<std::complex<double>,3>
DMRG_eigenvector_to_mps(Eigen::VectorXcd vec, tensor<std::complex<double>,3> prevState)
{
    /* T = (a_{l-1}, sg_l, a_l)
     * v_T = M^{sg_l}_{a_{l-1}, a_l} == M_{a_{l-1} sg_l a_l} (3-tensor)
     * therefore: v[T] = M[ num_to_multi-index(T) ]
     */
    std::vector<int> shape(3);
    shape[0] = prevState.shape()[0];
    shape[1] = prevState.shape()[1];
    shape[2] = prevState.shape()[2];
    tensor<std::complex<double>,3> T(shape);

    std::vector<int> index;
    for(int i = 0; i < vec.size(); ++i){
        index = number_to_multi_index(i,shape);
        T(index) = vec[i];
    }
    return T;
}

template<typename T>
double DMRG_mps_site_update 
                          (int i, std::vector<tensor<T,3> > &mpsState,         std::vector<tensor<T,4> > &mpsHamiltonian,
                                  std::vector<tensor<T,3> > &tripleVectorLeft, std::vector<tensor<T,3> > &tripleVectorRight,
                                  std::vector<tensor<T,2> > &doubleVectorLeft, std::vector<tensor<T,2> > &doubleVectorRight)
{
    typedef std::complex<double> cd;

    int L   = mpsState.size();
    int vli = i;               // index for vectorLeft
    int vri = L-i-1;           // index for vectorRight

    Eigen::MatrixXcd H = DMRG_matrixH(tripleVectorLeft[vli], tripleVectorRight[vri], mpsHamiltonian[i]);
    Eigen::MatrixXcd N = DMRG_matrixN(doubleVectorLeft[vli], doubleVectorRight[vri], mpsHamiltonian[i]);

    std::cout << "DMRG_mps_site_update: The matrix H --" << std::endl
              << H << std::endl;
    std::cout << "DMRG_mps_site_update: The matrix N --" << std::endl
              << N << std::endl;


    int MAX = 100;                               // matrices should NOT be larger than MAX x MAX
    Eigensolver es = lowest_eigenpair(H,N,MAX);  // uses lapack directly
    double lowest_eigenvalue  = es.eigenvalue;
    //tensor<cd,3> eigenvector  = DMRG_eigenvector_to_mps(es.eigenvector, mpsState[i]);

/*
    std::cout << "DMRG_mps_site_update: eigenvector (matrix form): " << std::endl;
    std::cout << es.eigenvector << std::endl;

    std::cout << "DMRG_mps_site_update: mpsState[" << i << "] type -- ";
    output_tensortype(mpsState[i]);
    output_tensorfull(mpsState[i]);
*/

    auto eigenvector  = DMRG_eigenvector_to_mps(es.eigenvector, mpsState[i]);
//    std::cout << "DMRG_mps_site_update: eigenvector -- ";
//    output_tensortype(eigenvector);
//    output_tensorfull(eigenvector);

    // update states:
    mpsState[i] = eigenvector;
    DMRG_double_update_site(i, mpsState, doubleVectorLeft, doubleVectorRight);
    DMRG_triple_update_site(i, mpsState, mpsHamiltonian, tripleVectorLeft, tripleVectorRight);


    std::cout << "DMRG_mps_site_update: The lowest eigenvalue is: " << lowest_eigenvalue << std::endl;
    std::cout << "DMRG_mps_site_update: Here is the contraction <psi_new|H|psi_new>/<psi_new|psi_new> is: ";
    std::cout << DMRG_current_eigenvalue(mpsState, mpsHamiltonian) << std::endl;

    std::cout << "DMRG_mps_site_update: The corresponding eigenvector is:" << std::endl;
    output_tensorfull( contract(mpsState[0],contract(mpsState[1],contract(mpsState[2],mpsState[3],ar::two,ar::zero),ar::two,ar::zero),ar::two,ar::zero) ); 

    return lowest_eigenvalue;
}

template<typename T>
std::vector<double>
DMRG_sweep (std::vector<tensor<T,3> > &mpsState, std::vector<tensor<T,4> > &mpsHamiltonian)
{
    std::vector<tensor<T,2> > double_L = DMRG_double_left_recursive (mpsState);
    std::vector<tensor<T,2> > double_R = DMRG_double_right_recursive(mpsState);
    std::vector<tensor<T,3> > triple_L = DMRG_triple_left_recursive (mpsState,mpsHamiltonian);
    std::vector<tensor<T,3> > triple_R = DMRG_triple_right_recursive(mpsState,mpsHamiltonian);

    int L = mpsState.size();
    std::vector<double> eigenvalues_history;
    eigenvalues_history.reserve(2*L+1);

    for(int i = 0; i < L; ++i)
        eigenvalues_history.push_back( 
            DMRG_mps_site_update( i, mpsState, mpsHamiltonian, triple_L, triple_R, double_L, double_R));

    for(int i = L-2; i >=0; --i)
        eigenvalues_history.push_back( 
            DMRG_mps_site_update( i, mpsState, mpsHamiltonian, triple_L, triple_R, double_L, double_R));

    for(int i = 1; i < L; ++i)
        eigenvalues_history.push_back( 
            DMRG_mps_site_update( i, mpsState, mpsHamiltonian, triple_L, triple_R, double_L, double_R));

    for(int i = L-2; i >=0; --i)
        eigenvalues_history.push_back( 
            DMRG_mps_site_update( i, mpsState, mpsHamiltonian, triple_L, triple_R, double_L, double_R));

    return eigenvalues_history;
}


#endif //MPS_DMRG_H

