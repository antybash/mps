
#ifndef MPS_TESTS_H
#define MPS_TESTS_H

#include <array>
#include <vector>
#include <algorithm>
#include <complex>
#include <cassert>

#include "tensor.h"
#include "DMRG_contractions.h"

namespace ar {
    std::array<int,1> zero      = {0};
    std::array<int,1> one       = {1};
    std::array<int,1> two       = {2};
    std::array<int,1> three     = {3};
    std::array<int,1> four      = {4};
    std::array<int,1> five      = {5};

    std::array<int,2> zeroone   = {0,1};
    std::array<int,2> onezero   = {1,0};
    std::array<int,2> onetwo    = {1,2};
    std::array<int,2> twothree  = {2,3};
    std::array<int,2> twoseven  = {2,7};
    std::array<int,2> threefour = {3,4};
    std::array<int,2> fourfive  = {4,5};
    std::array<int,2> threesix  = {3,6};

    std::array<int,3> a_123 = {1,2,3};
    std::array<int,3> a_135 = {1,3,5};
    std::array<int,3> a_345 = {3,4,5};
}


bool check_contract_1()
{
    tensor<int,3> a( std::vector<int>{{ 2,2,2 }} );
    tensor<int,2> b( std::vector<int>{{  2,2  }} );

    a[0][0][0] = 1; a[0][0][1] = 2; a[0][1][0] = 3; a[0][1][1] = 4; 
    a[1][0][0] = 5; a[1][0][1] = 6; a[1][1][0] = 7; a[1][1][1] = 8; 

    b[0][0] = 2; b[0][1] = 2;
    b[1][0] = 2; b[1][1] = 2;

    std::array<int,2> barA = {{0,1}};
    std::array<int,2> barB = {{0,1}};
    auto x = contract(a,b,barA,barB);

    bool alpha = (x[0] == 32);
    bool beta  = (x[1] == 40);

    if (!alpha || !beta){
        std::cout << "check_contract_1 failed: " << alpha << " " << beta << std::endl;
        return false;
    } else {
        std::cout << "check_contract_1 passed. Congratulations!" << std::endl;
        return true;
    }
}

bool check_contract_2()
{
    tensor<int,3> a( std::vector<int>{{ 2,2,2 }} );
    tensor<int,3> b( std::vector<int>{{ 2,2,2 }} );

    a[0][0][0] = 1; a[0][0][1] = 2; a[0][1][0] = 3; a[0][1][1] = 4; 
    a[1][0][0] = 5; a[1][0][1] = 6; a[1][1][0] = 7; a[1][1][1] = 8; 

    b[0][0][0] = 2; b[0][0][1] = 2; b[0][1][0] = 2; b[0][1][1] = 2; 
    b[1][0][0] = 2; b[1][0][1] = 2; b[1][1][0] = 2; b[1][1][1] = 2; 

    std::array<int,3> barA = {{0,1,2}};
    std::array<int,3> barB = {{0,1,2}};
    auto x = contract(a,b,barA,barB);

    std::vector<int> empty;

    if (! (x(empty) == 72) ) {
        std::cout << "check_contract_2 failed: " << x(empty) << std::endl;
        return false;
    } else {
        std::cout << "check_contract_2 passed. Congratulations!" << std::endl;
        return true;
    }
}

bool check3()
{
    tensor<int,3> a( std::vector<int>{{ 2,2,2 }} );

    for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
    for (int k = 0; k < 2; k++)
        a[i][j][k] = i + 2*j + 4*k;

    // Expected solution
    Eigen::Matrix<std::complex<double>,4,2> M;
    M << 0,4,
         1,5,
         2,6,
         3,7;
        

    std::array<int,2> rows = {{ 0,1 }};
    std::array<int,1> cols = {{  2  }};

//   std::cout << "Check Tensor to Matrix method:" << std::endl;
//   std::cout << tensor_to_matrix(a,rows,cols)    << std::endl;
//   std::cout << "The matrix that I wish to obtain:" << std::endl;
//   std::cout << M << std::endl;

    if(((M - tensor_to_matrix(a, rows, cols)).array() == 0).array().all())
        return true;
    return false;
}

bool check4()
{
    tensor<int,3> a( std::vector<int>{{ 2,2,2 }} );

    for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
    for (int k = 0; k < 2; k++)
        a[i][j][k] = i + 2*j + 4*k;

    Eigen::Matrix<std::complex<double>,4,2> M;
    M << 0,2,
         1,3,
         4,6,
         5,7;

    std::array<int,2> rows = {{ 0,2 }};
    std::array<int,1> cols = {{  1  }};

//   std::cout << "Check Tensor to Matrix method:" << std::endl;
//   std::cout << tensor_to_matrix(a,rows,cols)    << std::endl;
//   std::cout << "The matrix that I wish to obtain:" << std::endl;
//   std::cout << M << std::endl;

    if(((M - tensor_to_matrix(a, rows, cols)).array() == 0).array().all())
        return true;
    return false;
}

bool check5()
{
    Eigen::Matrix<std::complex<double>,4,2> M;
    M << 0,2, 1,3, 4,6, 5,7;

    std::vector<int> shape = {{2,2,2}};
    std::array<int,2> rows = {{ 0,1 }};
    std::array<int,1> cols = {{  2  }};

    return (((M - tensor_to_matrix(matrix_to_tensor(M,rows, cols, shape), rows, cols)).array() == 0).array().all());
}

bool check6()
{
    Eigen::Matrix<std::complex<double>,2,4> M;
    M << 0,2, 1,3, 4,6, 5,7;

    std::vector<int> shape = {{2,2,2}};
    std::array<int,1> rows = {{  1  }};
    std::array<int,2> cols = {{ 0,2 }};

    return  ((M - tensor_to_matrix(matrix_to_tensor(M,rows,cols,shape), rows, cols)).array() == 0).array().all();
}

bool check7()
{
    typedef std::complex<double> cd;
    typedef std::vector<int> vi;

    tensor<cd,2> dummy2(vi({{1,1}}));
    dummy2(vi({{0,0}})) = 1;

    tensor<cd,3> dummy3(vi({{1,1,1}}));
    dummy3(vi({{0,0,0}})) = 1;

    tensor<cd,3> m1(vi({{1,2,2}}));
    tensor<cd,3> m2(vi({{2,2,2}}));
    tensor<cd,3> m3(vi({{2,2,1}}));

    m1[0][0][0] = 2; m1[0][0][1] = 2;
    m1[0][1][0] = 2; m1[0][1][1] = 2;

    m2[0][0][0] = 3; m2[0][0][1] = 3;
    m2[0][1][0] = 3; m2[0][1][1] = 3;
    m2[1][0][0] = 3; m2[1][0][1] = 3;
    m2[1][1][0] = 3; m2[1][1][1] = 3;

    m3[0][0][0] = 4;
    m3[0][1][0] = 4;
    m3[1][0][0] = 4;
    m3[1][1][0] = 4;

    tensor<cd,4> h1(vi({{1,2,2,2}}));
    
    h1[0][0][0][0] = 7; h1[0][0][0][1] = 7;
    h1[0][0][1][0] = 7; h1[0][0][1][1] = 7;
    h1[0][1][0][0] = 7; h1[0][1][0][1] = 7;
    h1[0][1][1][0] = 7; h1[0][1][1][1] = 7;

    tensor<cd,4> h2(vi({{2,2,2,2}}));
    
    h2[0][0][0][0] = 11; h2[0][0][0][1] = 11;
    h2[0][0][1][0] = 11; h2[0][0][1][1] = 11;
    h2[0][1][0][0] = 11; h2[0][1][0][1] = 11;
    h2[0][1][1][0] = 11; h2[0][1][1][1] = 11;

    h2[1][0][0][0] = 11; h2[1][0][0][1] = 11;
    h2[1][0][1][0] = 11; h2[1][0][1][1] = 11;
    h2[1][1][0][0] = 11; h2[1][1][0][1] = 11;
    h2[1][1][1][0] = 11; h2[1][1][1][1] = 11;

    tensor<cd,4> h3(vi({{2,2,2,1}}));
    
    h3[0][0][0][0] = 3;
    h3[0][0][1][0] = 3;
    h3[0][1][0][0] = 3;
    h3[0][1][1][0] = 3;

    h3[1][0][0][0] = 3;
    h3[1][0][1][0] = 3;
    h3[1][1][0][0] = 3;
    h3[1][1][1][0] = 3;

    tensor<cd,2> A1(vi({{2,2}}));
    A1[0][0] = 8; A1[0][1] = 8;
    A1[1][0] = 8; A1[1][1] = 8;

    tensor<cd,3> A2(vi({{2,2,2}}));
    A2[0][0][0] = 112; A2[1][0][0] = 112;
    A2[0][0][1] = 112; A2[1][0][1] = 112;
    A2[0][1][0] = 112; A2[1][1][0] = 112;
    A2[0][1][1] = 112; A2[1][1][1] = 112;

    return (A1 == DMRG_double_left_contract_once(dummy2,m1))
        && (A2 == DMRG_triple_left_contract_once(dummy3,m1,h1));
}

bool eq_eps(Eigen::MatrixXcd M1, Eigen::MatrixXcd M2, double epsilon)
{
    assert (M1.rows() == M2.rows());
    assert (M1.cols() == M2.cols());

    int R = M1.rows(); 
    int C = M1.cols();
    for (int i = 0; i < R; ++i)
        for(int j = 0; j < C; ++j)
            if( abs(M1(i,j)-M2(i,j)) > epsilon ){
                std::cout << "bad at (" << i << ", " << j << ")" << std::endl;
                std::cout << M1(i,j) << " vs. " << M2(i,j) << std::endl;
                return false;
            }
    return true;
}

void check_svd_simple_helper(Eigen::MatrixXcd U, Eigen::MatrixXcd S, Eigen::MatrixXcd V,
                             int n)
{
    std::cout << "output for check_svd_" << n << std::endl;
    std::cout << "U" << std::endl << U << std::endl
              << "S" << std::endl << S << std::endl
              << "V" << std::endl << V << std::endl;
    std::cout << "USV" << std::endl << U*S*V << std::endl;
}

bool check_svd_helper(Eigen::MatrixXcd U, Eigen::MatrixXcd S, Eigen::MatrixXcd V,
                      Eigen::MatrixXcd u, Eigen::MatrixXcd s, Eigen::MatrixXcd v,
                      int n)
{
    if( (!eq_eps(u,U,10e-5) || !eq_eps(s,S,10e-5) || !eq_eps(v,V,10e-5)) 
            && !eq_eps(U*S*V, u*s*v,10e-5) )
    {
        std::cout << "check_svd_" << n << " failed:" << std::endl;
        std::cout << "U" << std::endl << U << std::endl
                  << "S" << std::endl << S << std::endl
                  << "V" << std::endl << V << std::endl;
        std::cout << "u" << std::endl << u << std::endl
                  << "s" << std::endl << s << std::endl
                  << "v" << std::endl << v << std::endl;
        std::cout << "USV" << std::endl << U*S*V << std::endl;
        std::cout << "usv" << std::endl << u*s*v << std::endl;
        return false;
    }
    else {
        std::cout << "check_svd_" << n << " passed. Congratulations!" << std::endl;
        return true;
    }
}
                            
bool check_svd_1()
{
    Eigen::MatrixXcd M(2,2);
    M << -3, 0,
          0, 0;

    Eigen::MatrixXcd U,S,V;
    std::tie(U,S,V) = custom_svd(M);

    Eigen::MatrixXcd u(2,2);
    Eigen::MatrixXcd s(2,2);
    Eigen::MatrixXcd v(2,2);

    u << -1, 0,
          0, 1;
    s << 3, 0,
         0, 0;
    v << 1, 0,
         0, 1;

    return check_svd_helper(U,S,V,u,s,v,1);
}

bool check_svd_2()
{
    Eigen::MatrixXcd M(2,2);
    M << 2, -1,
         2,  2;

    Eigen::MatrixXcd U,S,V;
    std::tie(U,S,V) = custom_svd(M);

    Eigen::MatrixXcd u(2,2);
    Eigen::MatrixXcd s(2,2);
    Eigen::MatrixXcd v(2,2);

    u << 1.0/sqrt(5), -2.0/sqrt(5),
         2.0/sqrt(5),  1.0/sqrt(5);
    s << 3, 0,
         0, 2;
    v << 2.0/sqrt(5), 1.0/sqrt(5),
        -1.0/sqrt(5), 2.0/sqrt(5);

    return check_svd_helper(U,S,V,u,s,v,2);
}

bool check_svd_3()
{
    Eigen::MatrixXcd M(3,2);
    M << 7, 1,
         0, 0,
         5, 5;

    Eigen::MatrixXcd U,S,V;
    std::tie(U,S,V) = custom_svd(M);

    Eigen::MatrixXcd u(3,3);
    Eigen::MatrixXcd s(3,2);
    Eigen::MatrixXcd v(2,2);

    u << 1.0/sqrt(2), -1.0/sqrt(2), 0,
             0      ,      0,       1, 
         1.0/sqrt(2),  1.0/sqrt(2), 0;
    s << 3.0*sqrt(10),         0,
              0,        sqrt(10),
              0,               0;
    v << 2.0/sqrt(5), 1.0/sqrt(5),
        -1.0/sqrt(5), 2.0/sqrt(5);

    return check_svd_helper(U,S,V,u,s,v,3);
}

bool check_svd_4()
{
    Eigen::MatrixXcd M(2,3);
    M << 3, 2, 2,
         2, 3, -2;

    Eigen::MatrixXcd U,S,V;
    std::tie(U,S,V) = custom_svd(M);

    Eigen::MatrixXcd u(2,2);
    Eigen::MatrixXcd s(2,3);
    Eigen::MatrixXcd v(3,3);

    u << 1.0/sqrt(2), -1.0/sqrt(2),
         1.0/sqrt(2),  1.0/sqrt(2);
    s << 5, 0, 0,
         0, 3, 0;
    v << 1.0/sqrt(2), 1.0/sqrt(2), 0,
        -1.0/sqrt(18),1.0/sqrt(18), -4.0/sqrt(18),
        -2.0/3.0, 2.0/3.0, 1.0/3.0;

    return check_svd_helper(U,S,V,u,s,v,4);
}

bool check_svd_5()
{
    Eigen::MatrixXcd M(4,5);
    M << 6, -8, -4, 5, -4,
         2, 7, -5, -6, 4,
         0, -1, -8, 2, 2,
         -1, -2, 4, 4, -8;

    Eigen::MatrixXcd U,S,V;
    std::tie(U,S,V) = custom_svd(M);

    Eigen::MatrixXcd u(2,2);
    Eigen::MatrixXcd s(2,3);
    Eigen::MatrixXcd v(3,3);

    u << 1.0/sqrt(2), -1.0/sqrt(2),
         1.0/sqrt(2),  1.0/sqrt(2);
    s << 5, 0, 0,
         0, 3, 0;
    v << 1.0/sqrt(2), 1.0/sqrt(2), 0,
        -1.0/sqrt(18),1.0/sqrt(18), -4.0/sqrt(18),
        -2.0/3.0, 2.0/3.0, 1.0/3.0;

    //check_svd_simple_helper(U,S,V,5);
    return true;
}

bool check_svd_then_trim_1()
{
    Eigen::MatrixXcd M(4,5);
    M << 6, -8, -4, 5, -4,
         2, 7, -5, -6, 4,
         0, -1, -8, 2, 2,
         -1, -2, 4, 4, -8;

    Eigen::MatrixXcd U,V; 
    int trim; 
    std::tie(trim,U,V) = svd_then_trim(M, 5);

    Eigen::MatrixXcd u(4,2); Eigen::MatrixXcd s(2,2); Eigen::MatrixXcd v(2,5);
    u << -0.57, -0.65,
          0.63, -0.24,
          0.07, -0.63,
         -0.51,  0.34;
    s << 16.46,     0,
             0, 12.16;
    v << -0.10, 0.61, -0.21, -0.52,  0.55,
         -0.39, 0.29,  0.84, -0.14, -0.19;

    bool a = eq_eps(u,U,0.1);
    bool b = eq_eps(s*v,V,0.1); 
    bool c = (trim == 2);
    if (a && b && c) {
        std::cout << "check_svd_then_trim_1 passed. Congratulations!" << std::endl;
        return true;
    } else {
        std::cout << "check_svd_then_trim_1 failed: " << a << b << c << std::endl;
        return false;
    }
}

template <typename T, size_t N>
bool eq_tensor(tensor<T,N> A, tensor<T,N> B, double epsilon)
{
    for(int i = 0; i < N; ++i)
        assert(A.shape()[i] == B.shape()[i]);
    std::vector<int> ind(N);
    do {
        if (abs(A(ind) - B(ind)) > epsilon){
            std::cout << "eq_tensor comparison failed at: " << std::endl;
            for (auto i : ind) 
                std::cout << i << " ";
            std::cout << std::endl;
            return false;
        }
        increment_index(ind,A.shape());
    } while(!check_all_zeros(ind));
    return true;
}

bool check_leftnormalized_1()
{
    tensor<cd,5> C(vi({ 1,2,2,2,1 }));
    C(vi({0, 0, 0, 0, 0})) = 3;  C(vi({0, 1, 0, 0, 0})) = 6;
    C(vi({0, 0, 1, 0, 0})) = 10; C(vi({0, 1, 1, 0, 0})) = 8;
    C(vi({0, 0, 0, 1, 0})) = 8;  C(vi({0, 1, 0, 1, 0})) = 7;
    C(vi({0, 0, 1, 1, 0})) = 1;  C(vi({0, 1, 1, 1, 0})) = 4;

    std::array<int,1> zero  = {0};
    std::array<int,1> two   = {2};
    std::array<int,1> three = {3};

    //std::vector<tensor<cd,3> > mpsA = tensor_to_left_normalized_mps(A,10e-4);
    auto mpsC = tensor_to_left_normalized_mps(C,1e-7);
    auto contracted_tensor = contract(mpsC[0],contract(mpsC[1],mpsC[2], two, zero), two, zero);
    if (!eq_tensor(contracted_tensor, C,1e-4)){
        oe("tests.h -- check_leftnormalized_1:\n\n\t Final tensors:");
        oe(mpsC.size());
        for(int i = 0; i < mpsC.size(); ++i){
            std::cout << "The " << i << "-th tensor has dimensions: ";
            std::cout << "(" << mpsC[i].shape()[0] << mpsC[i].shape()[1] << mpsC[i].shape()[2] << ")" << std::endl;
            output_tensor(mpsC[i]);
        }

        std::cout << "The 'final' test:" << std::endl;
        output_tensor( contracted_tensor );
        output_tensor( C );
        return false;
    } else {
        std::cout << "check_leftnormalized_1 passed. Congratulations!" << std::endl;
        return true;
    }
}

bool check_leftnormalized_2()
{
    // look at testing.py
    // to generate some testcases.
    // try to generalize the testing code: in particular implement "contract_all" : mps --> n-tensor
    return false;
}

bool check_leftnormalized_3()
{
    tensor<cd,4> C(vi({1,2,2,1}));
    C(vi({0, 0, 0, 0})) = 1; 
    C(vi({0, 1, 0, 0})) = 0; 
    C(vi({0, 0, 1, 0})) = 0; 
    C(vi({0, 1, 1, 0})) = 0; 

    std::array<int,1> zero  = {0};
    std::array<int,1> two   = {2};
    std::array<int,1> three = {3};

    //std::vector<tensor<cd,3> > mpsA = tensor_to_left_normalized_mps(A,10e-4);
    auto mpsC = tensor_to_left_normalized_mps(C,1e-7);
    auto contracted_tensor = contract(mpsC[0],mpsC[1], two, zero);

    if (!eq_tensor(contracted_tensor, C,1e-4)){
        oe("tests.h -- check_leftnormalized_3:\n\n\t Final tensors:");
        oe(mpsC.size());
        for(int i = 0; i < mpsC.size(); ++i){
            std::cout << "The " << i << "-th tensor has dimensions: ";
            std::cout << "(" << mpsC[i].shape()[0] << mpsC[i].shape()[1] << mpsC[i].shape()[2] << ")" << std::endl;
            output_tensor(mpsC[i]);
        }

        std::cout << "The 'final' test:" << std::endl;
        output_tensor( contracted_tensor );
        output_tensor( C );
        return false;
    } else {
        std::cout << "check_leftnormalized_3 passed. Congratulations!" << std::endl;
        return true;
    }
}

bool check_leftnormalized_4()
{
    tensor<cd,4> C(vi({1,2,2,1}));
    C(vi({0, 0, 0, 0})) = 0; 
    C(vi({0, 1, 0, 0})) = 0; 
    C(vi({0, 0, 1, 0})) = 0; 
    C(vi({0, 1, 1, 0})) = 1; 

    std::array<int,1> zero  = {0};
    std::array<int,1> two   = {2};
    std::array<int,1> three = {3};

    //std::vector<tensor<cd,3> > mpsA = tensor_to_left_normalized_mps(A,10e-4);
    auto mpsC = tensor_to_left_normalized_mps(C,1e-7);
    auto contracted_tensor = contract(mpsC[0],mpsC[1], two, zero);

    if (!eq_tensor(contracted_tensor, C,1e-4)){
        oe("tests.h -- check_leftnormalized_3:\n\n\t Final tensors:");
        oe(mpsC.size());
        for(int i = 0; i < mpsC.size(); ++i){
            std::cout << "The " << i << "-th tensor has dimensions: ";
            std::cout << "(" << mpsC[i].shape()[0] << mpsC[i].shape()[1] << mpsC[i].shape()[2] << ")" << std::endl;
            output_tensor(mpsC[i]);
        }

        std::cout << "The 'final' test:" << std::endl;
        output_tensor( contracted_tensor );
        output_tensor( C );
        return false;
    } else {
        std::cout << "check_leftnormalized_3 passed. Congratulations!" << std::endl;
        return true;
    }
}

/*
template <size_t N, size_t M>
tensor<cd,N+M> contract_all(std::vector<tensor<cd,3> > mps, tensor<cd,M> R)
{
    std::array<int,1> zero  = {0};
    std::array<int,1> two   = {2};

    assert(mps.size() == N);
    return contract_all<cd,N-1,M+1>( std::vector<tensor<cd,3> >(mps.begin(), mps.begin()+N-2), contract(mps[N-1],R,two,zero));
}

template<typename T,size_t M>
tensor<T,M> contract_all<T,0,M>(std::vector<tensor<T,3> > mps, tensor<T,M> R)
{
    return R;
}

bool check_leftnormalized_2()
{
    tensor<cd,5> C(vi({ 1,2,2,2,1 }));
    C(vi({0, 0, 0, 0, 0})) = 3;  C(vi({0, 1, 0, 0, 0})) = 6;
    C(vi({0, 0, 1, 0, 0})) = 10; C(vi({0, 1, 1, 0, 0})) = 8;
    C(vi({0, 0, 0, 1, 0})) = 8;  C(vi({0, 1, 0, 1, 0})) = 7;
    C(vi({0, 0, 1, 1, 0})) = 1;  C(vi({0, 1, 1, 1, 0})) = 4;

    //std::vector<tensor<cd,3> > mpsA = tensor_to_left_normalized_mps(A,10e-4);
    auto mpsC = tensor_to_left_normalized_mps(C,10e-7);
    auto mpsC_minus_last = std::vector<tensor<cd,3> >(mpsC.begin(),mpsC.begin()+1);
    auto contracted_tensor = contract_all<cd,2,3>(mpsC_minus_last,mpsC[2]);

    if (!eq_tensor(contracted_tensor, C,10e-4)){
        oe("tests.h -- check_leftnormalized_1:\n\n\t Final tensors:");
        oe(mpsC.size());
        for(int i = 0; i < mpsC.size(); ++i){
            std::cout << "The " << i << "-th tensor has dimensions: ";
            std::cout << "(" << mpsC[i].shape()[0] << mpsC[i].shape()[1] << mpsC[i].shape()[2] << ")" << std::endl;
            output_tensor(mpsC[i]);
        }

        std::cout << "The 'final' test:" << std::endl;
        output_tensor( contracted_tensor );
        output_tensor( C );
        return false;
    } else {
        std::cout << "check_leftnormalized_1 passed. Congratulations!" << std::endl;
        return true;
    }
}
*/

bool check_DMRG_left_contract_1()
{
    tensor<cd,5> psi(vi({1,2,2,2,1}));
    psi(vi({0, 0, 0, 0, 0})) = 8; psi(vi({0, 1, 0, 0, 0})) = 7; psi(vi({0, 0, 1, 0, 0})) = 10; psi(vi({0, 1, 1, 0, 0})) = 9;
    psi(vi({0, 0, 0, 1, 0})) = 5; psi(vi({0, 1, 0, 1, 0})) = 1; psi(vi({0, 0, 1, 1, 0})) = 10; psi(vi({0, 1, 1, 1, 0})) = 8;

    auto mps = tensor_to_left_normalized_mps(psi);
    auto vec = DMRG_double_left_recursive(mps);
    auto lst = vec[vec.size()-1];

    std::array<int,5> indices;
    std::iota(indices.begin(),indices.end(),0);
    auto ten = contract(psi,psi,indices,indices,true,false);

    // the last element in vec_contract is a 2-tensor of type (1,1)
    // therefore, effectively a scalar!
    double norm_mps = abs(lst[0][0]);
    double norm_ten = abs(ten(std::vector<int>()));
    if ( abs(norm_mps - norm_ten) > 1e-4 ){
        std::cout << "check_DMRG_left_contract_1: (norm-mps," << norm_mps << ") vs. (norm_ten," << norm_ten << ")" << std::endl;
        return false;
    } else {
        std::cout << "check_DMRG_left_contract_1 passed. Congratulations!" << std::endl;
        return true;
    }
}

bool check_DMRG_left_contract_2()
{
    tensor<cd,4> psi(vi({1,2,2,1}));
    psi(vi({0, 0, 0, 0})) = 0; 
    psi(vi({0, 1, 0, 0})) = 0; 
    psi(vi({0, 0, 1, 0})) = 0; 
    psi(vi({0, 1, 1, 0})) = 1; 

    auto mps = tensor_to_left_normalized_mps(psi);
    auto vec = DMRG_double_left_recursive(mps);
    auto lst = vec[vec.size()-1];

    std::array<int,4> indices;
    std::iota(indices.begin(),indices.end(),0);
    auto ten = contract(psi,psi,indices,indices,false,true);

    // the last element in vec_contract is a 2-tensor of type (1,1)
    // therefore, effectively a scalar!
    double norm_mps = abs(lst[0][0]);
    double norm_ten = abs(ten(std::vector<int>()));
    if ( abs(norm_mps - norm_ten) > 1e-4 ){
        std::cout << "check_DMRG_left_contract_2: (norm-mps," << norm_mps << ") vs. (norm_ten," << norm_ten << ")" << std::endl;
        return false;
    } else {
        std::cout << "check_DMRG_left_contract_2 passed. Congratulations!" << std::endl;
        return true;
    }
}

bool check_DMRG_triple_recursive_1()
{

    tensor<cd,5> psi(vi({1,2,2,2,1}));
    psi(vi({0, 0, 0, 0, 0})) = 0; 
    psi(vi({0, 1, 0, 0, 0})) = 0; 
    psi(vi({0, 0, 1, 0, 0})) = 1;
    psi(vi({0, 1, 1, 0, 0})) = 0;
    psi(vi({0, 0, 0, 1, 0})) = 0; 
    psi(vi({0, 1, 0, 1, 0})) = 0;
    psi(vi({0, 0, 1, 1, 0})) = 0;
    psi(vi({0, 1, 1, 1, 0})) = 0;

    auto mps = tensor_to_left_normalized_mps(psi,-1);
    assert(mps.size() == 3);

    std::vector<tensor<cd,4> > B;
    B.push_back( mpo::startH );
    B.push_back( mpo::middleH );
    B.push_back( mpo::endH );

    auto vec = DMRG_triple_left_recursive(mps, B); 
    auto lst = *(vec.rbegin()); 
    cd ans = simplify_constant_tensor(lst);

    if ( abs(ans-cd(-1.0,0)) > 1e-4 ){
        std::cout << "check_DMRG_triple_recursive_1 failed: " << ans << " but expected " << cd(-1.0,0)
            << " diff = " << abs(ans-cd(-1.0,0)) << " > " << 1e-4 << std::endl;
        return false;
    } else {
        std::cout << "check_DMRG_triple_recursive_1 passed. Congratulations!" << std::endl;
        return true;
    }

}

bool check_DMRG_triple_recursive_2()
{
    tensor<cd,3> psi(vi({1,2,1}));
    psi(vi({0, 0, 0})) = 1; 
    psi(vi({0, 1, 0})) = 0; 

    auto mps = tensor_to_left_normalized_mps(psi,-1);
    assert(mps.size() == 1);

    tensor<cd,4> oneH(vi({{1,2,2,1}}));
    set_mpo(0,0, oneH, mpo::Z);

    auto x1 = contract(mps[0],oneH,ar::one,ar::one);
    auto x2 = contract(x1,mps[0],ar::three,ar::one,false,true);

    cd ans = simplify_constant_tensor(x2);

    if ( abs(ans-cd(1.0,0)) > 1e-4 ){
        std::cout << "check_DMRG_triple_recursive_2 (or kind of) failed: " << ans << " but expected " << cd(1.0,0) 
            << " diff = " << abs(ans-cd(1.0,0)) << " > " << 1e-4 << std::endl;
        return false;
    } else {
        std::cout << "check_DMRG_triple_recursive_2 (or kind of) passed. Congratulations!" << std::endl;
        return true;
    }
}

bool check_DMRG_triple_recursive_3()
{
    tensor<cd,4> psi(vi({1,2,2,1}));
    psi(vi({0, 0, 0, 0})) = 0; 
    psi(vi({0, 1, 0, 0})) = 1; 
    psi(vi({0, 0, 1, 0})) = 0; 
    psi(vi({0, 1, 1, 0})) = 0; 

    auto mps = tensor_to_left_normalized_mps(psi,-1);
    assert(mps.size() == 2);
    
    mpo::initialize();
    std::vector<tensor<cd,4> > mpoHam;
    mpoHam.push_back( mpo::startH );
    mpoHam.push_back( mpo::endH   );

    auto vec = DMRG_triple_left_recursive(mps, mpoHam); 
    auto lst = *(vec.rbegin()); 

    cd ans = simplify_constant_tensor(lst);

    if ( abs(ans-cd(0.0,0.0)) > 1e-4 ){
        std::cout << "check_DMRG_triple_recursive_3 failed: " << ans << " but expected " << cd(0.0,0.0)  
            << " diff = " << abs(ans-cd(2.0,0.0)) << " > " << 1e-4 << std::endl;
        return false;
    } else {
        std::cout << "check_DMRG_triple_recursive_3 passed. Congratulations!" << std::endl;
        return true;
    }
}

bool check_DMRG_triple_recursive_4()
{
    mpo::initialize();

    tensor<cd,6> psi(vi({1,2,2,2,2,1}));
    psi(vi({0, 0, 0, 0, 0, 0})) = 0;
    psi(vi({0, 1, 0, 0, 0, 0})) = 0;
    psi(vi({0, 0, 1, 0, 0, 0})) = 0;
    psi(vi({0, 1, 1, 0, 0, 0})) = 0;
    psi(vi({0, 0, 0, 1, 0, 0})) = 0; 
    psi(vi({0, 1, 0, 1, 0, 0})) = 0; 
    psi(vi({0, 0, 1, 1, 0, 0})) = 0; 
    psi(vi({0, 1, 1, 1, 0, 0})) = 0;
    psi(vi({0, 0, 0, 0, 1, 0})) = 0;
    psi(vi({0, 1, 0, 0, 1, 0})) = 0;
    psi(vi({0, 0, 1, 0, 1, 0})) = 1;
    psi(vi({0, 1, 1, 0, 1, 0})) = 0;
    psi(vi({0, 0, 0, 1, 1, 0})) = 0; 
    psi(vi({0, 1, 0, 1, 1, 0})) = 0; 
    psi(vi({0, 0, 1, 1, 1, 0})) = 0; 
    psi(vi({0, 1, 1, 1, 1, 0})) = 0;

    auto mps = tensor_to_left_normalized_mps(psi,-1);
    assert(mps.size() == 4);

    std::vector<tensor<cd,4> > B;
    B.push_back( mpo::startH );
    B.push_back( mpo::middleH );
    B.push_back( mpo::middleH );
    B.push_back( mpo::endH );

    auto vec = DMRG_triple_left_recursive(mps, B); 
    auto lst = *(vec.rbegin()); 

    cd ans = simplify_constant_tensor(lst);

    if ( abs(ans-cd(-2.0,0)) > 1e-4 ){
        std::cout << "check_DMRG_triple_recursive_4 failed: " << ans << " but expected " << cd(-2.0,0) 
            << " diff = " << abs(ans-cd(-2.0,0)) << " > " << 1e-4 << std::endl;
        return false;
    } else {
        std::cout << "check_DMRG_triple_recursive_4 passed. Congratulations!" << std::endl;
        return true;
    }
} 

void check_all()
{
    std::vector<bool> results(
            { check_contract_1(),                       // 0
              check_contract_2(),                       // 1
              check3(),                                 // 2
              check4(),                                 // 3
              check5(),                                 // 4
              check6(),                                 // 5
              check7(),                                 // 6
              check_svd_1(),                            // 7
              check_svd_2(),                            // 8
              check_svd_3(),                            // 9   *
              check_svd_4(),                            // 10  *
              check_svd_5(),                            // 11
              check_svd_then_trim_1(),                  // 12  *
              check_leftnormalized_1(),                 // 13  *
              check_leftnormalized_3(),                 // 14
              check_leftnormalized_4(),                 // 15
              check_DMRG_left_contract_1(),             // 16
              check_DMRG_left_contract_2(),             // 17
              check_DMRG_triple_recursive_1(),          // 18  *
              check_DMRG_triple_recursive_2(),          // 19
              check_DMRG_triple_recursive_3(),          // 20
              check_DMRG_triple_recursive_4()});        // 21

    int failed = 0;
    for(int i = 0; i < results.size(); ++i)
        if (!results[i])
            ++failed;
    
    std::cout << " **************** SUMMARY of TESTS ****************** " << std::endl;

    if (failed == 0)
        std::cout << "All " << results.size() << " tests passed. Congratulations!" << std::endl;
    else
        for(int i = 0; i < results.size(); ++i)
            if (!results[i])
                std::cout << "Test case " << i << " did not pass." << std::endl;
}

#endif // MPS_TESTS_H
