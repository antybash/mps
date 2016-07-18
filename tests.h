
#ifndef MPS_TESTS_H
#define MPS_TESTS_H

#include <array>
#include <vector>
#include <algorithm>
#include <complex>

#include "tensor.h"

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

#endif // MPS_TESTS_H
