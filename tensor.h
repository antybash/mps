
#ifndef MPS_TENSOR_H
#define MPS_TENSOR_H


#include <iostream>
#include <algorithm>
#include <boost/multi_array.hpp>
#include <array>
#include <iterator>
#include <vector>
#include <complex>
#include <Eigen/Dense>

#include "increment_indices.h"
#include "reindexing.h"

template<typename T, std::size_t N>
using tensor = boost::multi_array<T,N>;

template<std::size_t N, typename T>
void output_tensor(tensor<T,N> t) {
    copy(t.origin(),t.origin()+t.num_elements(), std::ostream_iterator<T>(std::cout, " "));
    std::cout << std::endl;
}

template<typename T, std::size_t N, std::size_t M, std::size_t dA, std::size_t dB>
    tensor<T,N+M-dA-dB> 
contract(tensor<T,N> A, tensor<T,M> B, std::array<int,dA> barA, std::array<int,dB> barB)
{

    assert(barA.size() == barB.size());
    for (int k = 0; k < barA.size(); ++k)
        assert(A.shape()[barA[k]] == B.shape()[barB[k]]);

    std::vector<int> vecA(A.shape(), A.shape()+N);       // vecA, vecB == shape of A, B
    std::vector<int> vecB(B.shape(), B.shape()+M);       //
    std::vector<int> vecC(N+M-dA-dB);                    // 
    std::vector<int> indA(vecA.size());                  // indA, indB == indices of A, B
    std::vector<int> indB(vecB.size());                  // these variables are looped over below
    std::vector<int> indC(N+M-dA-dB);                    // indC = [ Union(indA,indB) - Overlap ]

    combine(vecC,vecA,vecB,barA,barB);
    tensor<T,N+M-dA-dB> C(vecC);

    do {
        do {

            T sum = 0;

            do {

                sum += A(indA)*B(indB);
                
                increment_index_with_selection(indA, vecA, barA);
                increment_index_with_selection(indB, vecB, barB);

            }while(!check_all_zeros_with_selection(indA,barA));


            // set value of C
            combine(indC, indA, indB, barA, barB);
            C(indC) = sum;

            // increment the outside
            increment_index_with_barrier(indB,vecB,barB);

        } while(!check_all_zeros_with_barrier(indB,barB));

        increment_index_with_barrier(indA,vecA,barA);
    } while(!check_all_zeros_with_barrier(indA,barA));

    return C;
}

template<size_t R, size_t C, typename T>
    Eigen::MatrixXcd
tensor_to_matrix(tensor<T,R+C> A, std::array<int,R> row_indices, std::array<int,C> col_indices)
{
    //dimensions of Matrix
    int Row = 1;
    int Col = 1;
    for (auto i : row_indices) Row *= A.shape()[i];
    for (auto i : col_indices) Col *= A.shape()[i];

    Eigen::MatrixXcd M(Row,Col);
    std::vector<int> ind(R+C); // initially set to zero

    do {
        int a = multi_index_to_number(reindexing_subset(ind, row_indices), reindexing_subset(A.shape(), row_indices));
        int b = multi_index_to_number(reindexing_subset(ind, col_indices), reindexing_subset(A.shape(), col_indices));

        M(a,b) = A(ind);

        increment_index(ind, A.shape());
    } while(!check_all_zeros(ind));
    return M;
}


template<size_t R, size_t C>
    tensor<std::complex<double>,R+C>
matrix_to_tensor(Eigen::MatrixXcd M, std::array<int,R> row_indices, std::array<int,C> col_indices, std::vector<int> shape)
{

    tensor<std::complex<double>,R+C> A(shape);
    std::vector<int> ind(R+C); // initially set to zero

    do {
        int a = multi_index_to_number(reindexing_subset(ind, row_indices), reindexing_subset(A.shape(), row_indices));
        int b = multi_index_to_number(reindexing_subset(ind, col_indices), reindexing_subset(A.shape(), col_indices));

        A(ind) = M(a,b);

        increment_index(ind, A.shape());
    } while(!check_all_zeros(ind));
    return A;
}


#endif // MPS_TENSOR_H
