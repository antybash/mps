
#ifndef MPS_TESTS_H
#define MPS_TESTS_H

#include <array>
#include <vector>
#include <algorithm>
#include <complex>

#include "tensor.h"
#include "DMRG_contractions.h"

bool check1()
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

    if(x[0] == 32 && x[1] == 40) 
        return true;
    return false;
}

bool check2()
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
    if(x(empty) == 72) 
        return true;
    return false;
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
    dummy2(vi({{0,0}})) = cd(1,0);

    tensor<cd,3> dummy3(vi({{1,1,1}}));
    dummy3(vi({{0,0,0}})) = cd(1,0);

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
    if (M1.rows() != M2.rows() || M1.cols() != M2.cols())
        return false;
    int R = M1.rows(); 
    int C = M1.cols();
    for (int i = 0; i < R; ++i)
        for(int j = 0; j < C; ++j)
            if( abs(M1(i,j)-M2(i,j)) > epsilon ){
                std::cout << "bad at (" << i << ", " << j << ")" << std::endl;
                return false;
            }
    return true;
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

    if( !eq_eps(u,U,10e-5) || !eq_eps(s,S,10e-5) || !eq_eps(v,V,10e-5)
            && !eq_eps(U*S*V, u*s*v,10e-5) )
    {
        std::cout << "check_svd_1 failed:" << std::endl;
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
        std::cout << "check_svd_1 passed. Congratulations!" << std::endl;
        return true;
    }
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

    if( !eq_eps(u,U,10e-5) || !eq_eps(s,S,10e-5) || !eq_eps(v,V,10e-5) 
            && !eq_eps(U*S*V, u*s*v,10e-5) )
    {
        std::cout << "check_svd_2 failed:" << std::endl;
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
        std::cout << "check_svd_2 passed. Congratulations!" << std::endl;
        return true;
    }
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

    if( (!eq_eps(u,U,10e-5) || !eq_eps(s,S,10e-5) || !eq_eps(v,V,10e-5)) 
            && !eq_eps(U*S*V, u*s*v,10e-5) )
    {
        std::cout << "check_svd_3 failed:" << std::endl;
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
        std::cout << "check_svd_3 passed. Congratulations!" << std::endl;
        return true;
    }
}
#endif // MPS_TESTS_H
