
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
#include <tuple>
#include <string>

#include "increment_indices.h"
#include "reindexing.h"

typedef std::complex<double> cd;
typedef std::vector<int> vi;

template<typename T, std::size_t N>
using tensor = boost::multi_array<T,N>;

template<std::size_t N, typename T>
void output_tensor(tensor<T,N> t)
{
    copy(t.origin(),t.origin()+t.num_elements(), std::ostream_iterator<T>(std::cout, " "));
    std::cout << std::endl;
}

template<typename T, std::size_t N, std::size_t M, std::size_t dA, std::size_t dB>
    tensor<T,N+M-dA-dB> 
contract(tensor<T,N> A, tensor<T,M> B, std::array<int,dA> barA, std::array<int,dB> barB, bool conjA=false, bool conjB=false)
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

                if (conjA && conjB)
                    sum += std::conj(A(indA)) * std::conj(B(indB));
                else if (!conjA && conjB)
                    sum += A(indA) * std::conj(B(indB));
                else if (conjA && !conjB)
                    sum += std::conj(A(indA)) * B(indB);
                else
                    sum += A(indA) * B(indB);
                
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
        // std::cout << "matrix_to_tensor: ";
        // oe(multi_index_to_number(ind,shape));
        // std::cout << "matrix_to_tensor (before increment_index)." << std::endl;
        increment_index(ind, A.shape());
        // std::cout << "matrix_to_tensor (after increment_index)." << std::endl;
    } while(!check_all_zeros(ind));
    return A;
}

template<typename T, size_t N>
std::vector<int> tensor_shape(tensor<T,N> t)
{
    std::vector<int> x(N);
    for(int i = 0; i < N; ++i)
        x[i] = t.shape()[i];
    return x;
}

int vector_trim_size(Eigen::VectorXcd S, double epsilon)
{
    for (int i = 0; i < S.size(); ++i)
        if(abs(S[i]) < epsilon)
            return (i-1);
    return S.size();
}

std::tuple<Eigen::MatrixXcd, Eigen::MatrixXcd, Eigen::MatrixXcd> custom_svd(Eigen::MatrixXcd M)
{
    Eigen::JacobiSVD<Eigen::MatrixXcd> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::VectorXcd vecS(svd.singularValues().size());
    for(int i = 0; i < svd.singularValues().size(); ++i)
        vecS[i] = svd.singularValues()[i];
    Eigen::MatrixXcd S = vecS.asDiagonal();
    int ur = svd.matrixU().rows();
    int uc = svd.matrixU().cols();
    int vr = svd.matrixV().rows();
    int vc = svd.matrixV().cols();

    if ( uc > vr )
        S.conservativeResize(uc, S.cols());
    else if ( uc < vr )
        S.conservativeResize(S.rows(), vr);

    return std::make_tuple(svd.matrixU(), S, svd.matrixV().conjugate().transpose());
}

std::tuple<int, Eigen::MatrixXcd, Eigen::MatrixXcd> svd_then_trim(Eigen::MatrixXcd M, double epsilon)
{
    Eigen::MatrixXcd U, S, V;
    std::tie(U,S,V) = custom_svd(M);
    V = S * V;
    int ts = vector_trim_size(S.diagonal(), epsilon);
    U.conservativeResize(U.rows(), ts);
    V.conservativeResize(ts, V.rows());
    return std::make_tuple(ts, U, V);
}


template<typename T,size_t N>
std::vector<tensor<T,3> > tensor_to_left_normalized_mps (tensor<T,N> A, double epsilon)
{
    std::array<int,2> fst = {{0,1}};
    std::array<int,1> snd = {{ 2 }};

    std::vector<tensor<T,3> > mpsState;

    // A = A_(i0... i_{n+1})
    // tmp = A (i0i1) (i2...i_{n+1})
    Eigen::MatrixXcd tmp       = tensor_to_matrix(A, fst, array_range<2,N-1>());
    std::vector<int> tmp_shape = tensor_shape(A);

    int trim;
    //std::cout << "Start left-normalized for-loop" << std::endl;
    for(int i = 0; i < N-3; ++i){
        Eigen::MatrixXcd U;
        Eigen::MatrixXcd V;
        std::tie (trim, U, V) = svd_then_trim(tmp, epsilon);

        //oe(i);
        //oe(U);
        mpsState.push_back( matrix_to_tensor(U, fst, snd, triple(tmp_shape[0],tmp_shape[1],trim)) );
        //oe("left-normalized: done adding U");

        // reshape(V)
        // A01234;     N = 5 = 1+3+1
        // A(01)(234) --> U(01)(s) V(s)(234)
        //            --> U[0s1]   V[s,2,3,4]

        tmp_shape.erase(tmp_shape.begin());
        tmp_shape[0] = trim;

        //oe("left-normalized: begin inplace swap");
        //oe(V);
        //out(tmp_shape);
        tmp          = inplace_index_swap_of_underlying_tensor(V, tmp_shape);
        //oe("left-normalized: end inplace swap");
    }

    // deal with last edge case
    //oe("left-normalized: last mps load");
    //oe(tmp);
    //out(tmp_shape);
    //oe(trim);
    //oe("left-normalized: last mps define");
    mpsState.push_back(  matrix_to_tensor(tmp, fst, snd, tmp_shape) );
    //oe("left-normalized: last mps done!");

    return mpsState;
}


#endif // MPS_TENSOR_H
