
#ifndef MPS_DMRG_H
#define MPS_DMRG_H

#include "tensor.h"
#include "DMRG_contractions.h"

#include <Eigen/Dense>
#include <vector>
#include <array>
#include <complex>

tensor<std::complex<double>,2> delta_tensor(int a, int b)
{
    typedef std::vector<int> vi;
    tensor<std::complex<double>,2> delta(vi({{a,b}}));
    for(int i = 0; i < a; ++i) {
        for(int j = 0; j < b; ++j) {
            if (i == j)
                delta(vi({{i,j}})) = 1;
            else
                delta(vi({{i,j}})) = 0;
        }
    }
}

template<typename T>
Eigen::MatrixXcd
DMRG_matrixN(tensor<T,2> left, tensor<T,2> right, tensor<T,4> localHamiltonian)
{
    int a = localHamiltonian.shape()[1];
    int b = localHamiltonian.shape()[2];
    std::array<int,0> empty;
    tensor<T,4> T1 = contract(left,delta_tensor(a,b),empty,empty);
    tensor<T,6> T2 = contract(T1,  right,            empty,empty);
    // (l0 l1 d0 d1 r0 r1)

    std::array<int,3> rows = {{ T2.shape()[0], T2.shape()[2], T2.shape()[4] }};
    std::array<int,3> cols = {{ T2.shape()[1], T2.shape()[3], T2.shape()[5] }};

    return tensor_to_matrix(T2,rows,cols);
}

template<typename T>
Eigen::MatrixXcd
DMRG_matrixH(tensor<T,3> left, tensor<T,3> right, tensor<T,4> localHamiltonian)
{
    typedef std::array<int,1> ar;
    tensor<T,5> T1 = contract(left,localHamiltonian,ar({{1}}),ar({{0}}));
    // (l0 l1 l2) (h0 h1 h2 h3)  ----->  (l0 l2 h1 h2 h3)

    tensor<T,6> T2 = contract(T1, right,            ar({{4}}),ar({{1}}));
    // (l0 l2 h1 h2 h3) (r0 r1 r2) --->  (l0 l2 h1 h2 r0 r2)
    //         this is NOT the same as the previous function!
    //         this is analogous to      (l0 l1 d1 d0 r0 r1)

    std::array<int,3> rows = {{ T2.shape()[0], T2.shape()[3], T2.shape()[4] }};
    std::array<int,3> cols = {{ T2.shape()[1], T2.shape()[2], T2.shape()[5] }};

    return tensor_to_matrix(T2,rows,cols);
}

tensor<std::complex<double>,3>
DMRG_eigenvector_to_mps(Eigen::VectorXcd vec, tensor<std::complex<double>,3> prevState)
{
    std::vector<int> tmp_shape(3);
    tmp_shape[0] = prevState.shape()[0];
    tmp_shape[1] = prevState.shape()[1];
    tmp_shape[2] = prevState.shape()[2];

    tensor<std::complex<double>,3> T(tmp_shape);
    std::vector<int> index;
    for(int i = 0; i < vec.size(); ++i){
        index = number_to_multi_index(i,tmp_shape);
        T(index) = vec[i];
    }
    return T;
}

/*
template<typename T>
void DMRG_mps_site_update (int i, std::vector<tensor<T,3> > &mpsState, std::vector<tensor<T,4> > &mpsHamiltonian,
                                  std::vector<tensor<T,3> > &tripleVectorLeft, std::vector<tensor<T,3> > &tripleVectorRight
                                  std::vector<tensor<T,3> > &doubleVectorLeft, std::vector<tensor<T,3> > &doubleVectorRight)
{
}
*/


#endif //MPS_DMRG_H

