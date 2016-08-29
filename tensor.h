
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
#include "utilities.h"

typedef std::complex<double> cd;
typedef std::vector<int> vi;

template<typename T, std::size_t N>
using tensor = boost::multi_array<T,N>;

template<std::size_t N, typename T>
void output_tensorfull(tensor<T,N> t)
{
    std::vector<int> ind(N);
    std::vector<int> shape(N);
    for(int i = 0; i < N; ++i)
        shape[i] = t.shape()[i];

    do{
        std::cout << "(";
        for(int i = 0; i < N; ++i)
            std::cout << ind[i] << ",";
        std::cout << ") = " << t(ind) << std::endl;

        increment_index(ind,shape);
    } while( !check_all_zeros(ind) );
}

template<std::size_t N, typename T>
void output_tensor(tensor<T,N> t)
{
    copy(t.origin(),t.origin()+t.num_elements(), std::ostream_iterator<T>(std::cout, " "));
    std::cout << std::endl;
}

template<typename T, size_t N>
void output_tensortype(tensor<T,N> t)
{
    std::cout << "(";
    for(int i = 0; i < N; ++i)
        std::cout << t.shape()[i] << ",";
    std::cout << ")" << std::endl;
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

/*
    std::cout << "Begin contraction:" << std::endl;
    std::cout << "A type -- ";       output_tensortype(A);
    std::cout << "A tensor --\n";    output_tensorfull(A);
    std::cout << std::endl;
    std::cout << "B type -- ";       output_tensortype(B);
    std::cout << "B tensor --\n";    output_tensorfull(B);
    std::cout << std::endl;
*/

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

/*
                std::cout << "\tA === (";
                for(int i = 0; i < N; ++i)
                    std::cout << indA[i] << ",";
                std::cout << ") = " << A(indA) << std::endl;
                std::cout << "\tB === (";
                for(int i = 0; i < M; ++i)
                    std::cout << indB[i] << ",";
                std::cout << ") = " << B(indB) << std::endl;
*/

                increment_index_with_selection(indA, vecA, barA);
                increment_index_with_selection(indB, vecB, barB);

                for (int k = 0; k < barA.size(); ++k)
                    assert(indA[barA[k]] == indB[barB[k]]);

            }while(!check_all_zeros_with_selection(indA,barA));
            combine(indC, indA, indB, barA, barB);        // set value of indC
            C(indC) = sum;                                // set value of C

/*
            std::cout << "\t\t C === (";
            for(int i = 0; i < N+M-dA-dB; ++i)
                std::cout << indC[i] << ",";
            std::cout << ") = " << C(indC) << std::endl;
*/

            increment_index_with_barrier(indB,vecB,barB); // increment the outside
        } while(!check_all_zeros_with_barrier(indB,barB));
        increment_index_with_barrier(indA,vecA,barA);
    } while(!check_all_zeros_with_barrier(indA,barA));

//    std::cout << "C type -- ";       output_tensortype(C);
//    std::cout << "C tensor --\n";    output_tensorfull(C);
//    std::cout << std::endl;
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
            return i;
    return S.size();
}

std::tuple<Eigen::MatrixXcd, Eigen::MatrixXcd, Eigen::MatrixXcd> custom_svd(Eigen::MatrixXcd M)
{
    Eigen::JacobiSVD<Eigen::MatrixXcd> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::VectorXcd vecS(svd.singularValues().size());
    for(int i = 0; i < svd.singularValues().size(); ++i)
        vecS[i] = svd.singularValues()[i];
    int ur = svd.matrixU().rows();
    int uc = svd.matrixU().cols();
    int vr = svd.matrixV().rows();
    int vc = svd.matrixV().cols();

    int r,c;
    if (uc > vr){
        r = uc;
        c = vecS.size();
    } else {
        r = vecS.size();
        c = vr;
    }

    //Eigen::MatrixXcd S(r,c) = Eigen::Matrix
    Eigen::MatrixXcd S(r,c);
    for(int i = 0; i < r; ++i)
        for(int j = 0; j < c; ++j)
            S(i,j) = cd(0,0);
    for(int i = 0; i < vecS.size(); ++i)
        S(i,i) = vecS[i];

    return std::make_tuple(svd.matrixU(), S, svd.matrixV().conjugate().transpose());
}

std::tuple<int, Eigen::MatrixXcd, Eigen::MatrixXcd> 
svd_then_trim (Eigen::MatrixXcd M, double epsilon)
{
    Eigen::MatrixXcd U, S, V;
    std::tie(U,S,V) = custom_svd(M);
    V = S * V;
    int ts = vector_trim_size(S.diagonal(), epsilon);
    U.conservativeResize(U.rows(), ts);                      //// CONSERVATIVE RESIZE !!!!! ////
    V.conservativeResize(ts, V.cols());
    return std::make_tuple(ts, U, V);
}


template<typename T,size_t N>
std::vector<tensor<T,3> > 
tensor_to_left_normalized_mps (tensor<T,N> A, double epsilon=1e-4)
{
    std::vector<tensor<T,3> > mpsState;
    // A = A_(i0... i_{n+1})
    // tmp = A (i0i1) (i2...i_{n+1})
    Eigen::MatrixXcd tmp       = tensor_to_matrix(A, ar::zeroone, array_range<2,N-1>());
    std::vector<int> tmp_shape = tensor_shape(A);

    int trim;
    for(int i = 0; i < N-3; ++i){ 
        Eigen::MatrixXcd U;
        Eigen::MatrixXcd V;
        std::tie (trim, U, V) = svd_then_trim(tmp, epsilon);
        mpsState.push_back( matrix_to_tensor(U, ar::zeroone, ar::two, triple(tmp_shape[0],tmp_shape[1],trim)) );
        // A(01)(234) --> U(01)(s) V(s)(234)
        //            --> U[01s]   V[s,2,3,4]
        tmp_shape.erase(tmp_shape.begin());
        tmp_shape[0] = trim;
        tmp          = inplace_index_swap_of_underlying_tensor(V, tmp_shape);
    }
    mpsState.push_back(  matrix_to_tensor(tmp, ar::zeroone, ar::two, tmp_shape) );
    return mpsState;
}


template<typename T, size_t N>
T simplify_constant_tensor (tensor<T,N> &t)
{
    std::vector<int> shape(N);
    for(int i = 0; i < N; ++i){
        assert(t.shape()[i] == 1);
        shape[i] = 0;
    }
    return t(shape);
}

void set_mpo(int r, int c, tensor<cd,4> &t, Eigen::MatrixXcd M)
{
    for (int sg_out = 0; sg_out < 2; ++sg_out)
        for (int sg_in = 0; sg_in < 2; ++sg_in)
            t[r][sg_in][sg_out][c] = M(sg_out,sg_in);     /// M_ab = <a|M|b>
}

namespace mpo {
    // pauli-matrices
    Eigen::Matrix2cd Z;
    Eigen::Matrix2cd I;
    Eigen::Matrix2cd zero;

    // spin-Hamiltonian
    tensor<cd,4> startH(vi({{1,2,2,3}}));
    tensor<cd,4> middleH(vi({{3,2,2,3}}));
    tensor<cd,4> endH(vi({{3,2,2,1}}));

    bool initialized_check_spin = false;

    void initialize(){
        //if (initialized_check_spin)
        //   break;
        Z(0,0) = 1; Z(0,1) =  0; Z(1,0) = 0; Z(1,1) = -1;
        I(0,0) = 1; I(0,1) = 0; I(1,0) = 0; I(1,1) = 1;
        zero(0,0) = 0; zero(0,1) = 0; zero(1,0) = 0; zero(1,1) = 0;

        set_mpo(0,0,startH, I);                     set_mpo(0,1,startH, Z);                     set_mpo(0,2,startH, zero);

        set_mpo(0,0, middleH, I);                   set_mpo(0,1, middleH, Z);                   set_mpo(0,2, middleH, zero);
        set_mpo(1,0, middleH, zero);               set_mpo(1,1, middleH, zero);               set_mpo(1,2, middleH, Z);
        set_mpo(2,0, middleH, zero);               set_mpo(2,1, middleH, zero);               set_mpo(2,2, middleH, I);

        set_mpo(0,0, endH,I);
        set_mpo(1,0, endH,Z);
        set_mpo(2,0, endH,I);
    }

};



#endif // MPS_TENSOR_H
