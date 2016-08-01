
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

        increment_index(ind, A.shape());
    } while(!check_all_zeros(ind));
    return A;
}

std::vector<int> vector_range(int a, int b)
    //returns a vector [a, a+1, ..., b-1, b] *inclusive*
{
    std::vector<int> x(b-a+1);
    for (int i = 0; i < b-a+1; ++i)
        x[i] = a+i;
    return x;
}

template<typename T, size_t N>
std::vector<int> tensor_shape(tensor<T,N> t)
{
    std::vector<int> x(N);
    for(int i = 0; i < N; ++i)
        x[i] = t.shape()[i];
    return x;
}

int middle_matrix_size_trim(Eigen::MatrixXcd diag_mat, double epsilon)
{
    Eigen::VectorXd vec = diag_mat.diagonal();
    for (int i = 0; i < vec.size(); ++i)
        if (abs(vec[i]) < epsilon) 
            return i-1;
    return vec.size();
}

void trim_matrix(Eigen::MatrixXcd &M, int index, int new_index_size)
{
    if (index == 0)
        M.conservativeResize(new_index_size, M.cols());
    else
        M.conservativeResize(M.rows(), new_index_size);
}

int vector_trim_size(Eigen::VectorXcd S, double epsilon)
{
    // if trim_size = -1,
    // then all abs(S[k]) > epsilon
    // else (S[k] < epsilon <=> k > trim_size)
    for (int i = 0; i < S.size(); ++i)
        if(abs(S[i]) < epsilon)
            return (i-1);
    return -1;
}

std::tuple<int, Eigen::MatrixXcd, Eigen::MatrixXcd> svd_then_trim(Eigen::MatrixXcd M, double epsilon)
{
    Eigen::JacobiSvd<MatrixXcd> svd(M);
    int ts = vector_trim_size(svd.singularValues(), epsilon);

    Eigen::MatrixXcd U = svd.matrixU();
    Eigen::MatrixXcd V = svd.singularValues().asDiagonal()*svd.matrixV();
    if(trim_size != -1){
        U.conservativeResize(svd.matrixU().rows(), trim_size);
        V.conservativeResize(trim_size, svd.matrixV().rows());
    }
    return std::make_tuple(ts, U, V);
}

template<typename T,size_t N>
std::vector<tensor<T,3> > tensor_to_left_normalized_mps (tensor<T,N> A, double epsilon)
{
    using std::begin;
    using std::advance;
    std::vector<tensor<T,3> > mpsState(N);
    // A = A_(i0... i_{n+1})
    // tmp = A (i0i1) (i2...i_{n+1})
    Eigen::MatrixXcd tmp       = tensor_to_matrix(A, vi({{0,1}}), vector_range(2,N-1));
    std::vector<int> tmp_shape = tensor_shape(A);

    int trim;
    for(int i = 0; i < N-2; ++i){
        Eigen::MatrixXcd U;
        Eigen::MatrixXcd V;
        std::tie (trim, U, V) = svd_then_trim(tmp, epsilon);

        vi first_three_elements(tmp_shape.begin(), tmp_shape.begin()+3);
        mpsState[i]   = matrix_to_tensor(U, vi({0,1}), vi({2}), vi({tmp_shape[0],tmp_shape[1],trim}));

        // reshape(V)
        // A01234;     N = 5 = 1+3+1
        // A(01)(234) --> U(01)(s) V(s)(234)
        //            --> U[0s1]   V[s,2,3,4]

        int beta = V.cols().size();
        int n    = tmp_shape.size();
        std::vector<int> tmp2_shape = number_to_multi_index(beta, vector_range(tmp_shape,2,n-1));
        tmp2_shape.insert(trim, tmp2_shape);
        std::vector<int> index (tmp2_shape.size());

        int rows_new = tmp2_shape[0] * tmp2_shape[1];
        int cols_new = 1;
        for(int i = 2; i < tmp2_shape.size(); ++i)
            cols_new *= tmp2_shape[i];

        Eigen::MatrixXcd tmp_new(rows_new, cols_new);

        do {
            // (s,a2),(a3,...,an)
            // this defines the new tmp_matrix;

            std::vector<int> vec_i1 = vector_range(index,0,0);
            std::vector<int> vec_j1 = vector_range(index,1,index.size()-1);

            std::vector<int> vec_a1 = vector_range(tmp2_shape,0,0);
            std::vector<int> vec_b1 = vector_range(tmp2_shape,1,tmp2_shape.size()-1);

            std::vector<int> vec_i2 = vector_range(index,0,1);
            std::vector<int> vec_j2 = vector_range(index,2,index.size()-1);

            std::vector<int> vec_a2 = vector_range(tmp2_shape,0,1);
            std::vector<int> vec_b2 = vector_range(tmp2_shape,2,tmp2_shape.size()-1);

            int i1 = multi_index_to_number(vec_i1, vec_a1);
            int j1 = multi_index_to_number(vec_j1, vec_b1);

            int i2 = multi_index_to_number(vec_i2, vec_a2);
            int j2 = multi_index_to_number(vec_j2, vec_b2);

            tmp_new[i2][j2] = V[i1][j1];

        } while (!check_all_zeros_with_selection(index));

        tmp = tmp_new;
    }
    // deal with last edge case
    mpsState[N-1] = matrix_to_tensor(tmp, vi({0,1}),vi({2}), vi({trim,tmp_shape[N-2],tmp_shape[N-1]}));
}


#endif // MPS_TENSOR_H
