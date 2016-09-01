#include <iostream>
#include <vector>

#include "tensor.h"
#include "DMRG.h"
#include "utilities.h"

int main()
{
    mpo::initialize();

    tensor<cd,8> psi(vi({1,2,2,2,2,2,2,1}));
    //In [1]: import testing as T
    //In [2]: import numpy as np
    //In [3]: T.print_random_cpp_tensor(index_list=[1,2,2,2,2,2,2,1])
    psi(vi({0, 0, 0, 0, 0, 0, 0, 0})) = 1;
    psi(vi({0, 1, 0, 0, 0, 0, 0, 0})) = 6;
    psi(vi({0, 0, 1, 0, 0, 0, 0, 0})) = 2;
    psi(vi({0, 1, 1, 0, 0, 0, 0, 0})) = 5;
    psi(vi({0, 0, 0, 1, 0, 0, 0, 0})) = 10;
    psi(vi({0, 1, 0, 1, 0, 0, 0, 0})) = 6;
    psi(vi({0, 0, 1, 1, 0, 0, 0, 0})) = 5;
    psi(vi({0, 1, 1, 1, 0, 0, 0, 0})) = 9;
    psi(vi({0, 0, 0, 0, 1, 0, 0, 0})) = 5;
    psi(vi({0, 1, 0, 0, 1, 0, 0, 0})) = 7;
    psi(vi({0, 0, 1, 0, 1, 0, 0, 0})) = 2;
    psi(vi({0, 1, 1, 0, 1, 0, 0, 0})) = 6;
    psi(vi({0, 0, 0, 1, 1, 0, 0, 0})) = 3;
    psi(vi({0, 1, 0, 1, 1, 0, 0, 0})) = 1;
    psi(vi({0, 0, 1, 1, 1, 0, 0, 0})) = 7;
    psi(vi({0, 1, 1, 1, 1, 0, 0, 0})) = 2;
    psi(vi({0, 0, 0, 0, 0, 1, 0, 0})) = 7;
    psi(vi({0, 1, 0, 0, 0, 1, 0, 0})) = 3;
    psi(vi({0, 0, 1, 0, 0, 1, 0, 0})) = 5;
    psi(vi({0, 1, 1, 0, 0, 1, 0, 0})) = 9;
    psi(vi({0, 0, 0, 1, 0, 1, 0, 0})) = 2;
    psi(vi({0, 1, 0, 1, 0, 1, 0, 0})) = 1;
    psi(vi({0, 0, 1, 1, 0, 1, 0, 0})) = 2;
    psi(vi({0, 1, 1, 1, 0, 1, 0, 0})) = 5;
    psi(vi({0, 0, 0, 0, 1, 1, 0, 0})) = 7;
    psi(vi({0, 1, 0, 0, 1, 1, 0, 0})) = 5;
    psi(vi({0, 0, 1, 0, 1, 1, 0, 0})) = 9;
    psi(vi({0, 1, 1, 0, 1, 1, 0, 0})) = 5;
    psi(vi({0, 0, 0, 1, 1, 1, 0, 0})) = 8;
    psi(vi({0, 1, 0, 1, 1, 1, 0, 0})) = 5;
    psi(vi({0, 0, 1, 1, 1, 1, 0, 0})) = 2;
    psi(vi({0, 1, 1, 1, 1, 1, 0, 0})) = 5;
    psi(vi({0, 0, 0, 0, 0, 0, 1, 0})) = 10;
    psi(vi({0, 1, 0, 0, 0, 0, 1, 0})) = 2;
    psi(vi({0, 0, 1, 0, 0, 0, 1, 0})) = 1;
    psi(vi({0, 1, 1, 0, 0, 0, 1, 0})) = 9;
    psi(vi({0, 0, 0, 1, 0, 0, 1, 0})) = 5;
    psi(vi({0, 1, 0, 1, 0, 0, 1, 0})) = 2;
    psi(vi({0, 0, 1, 1, 0, 0, 1, 0})) = 10;
    psi(vi({0, 1, 1, 1, 0, 0, 1, 0})) = 2;
    psi(vi({0, 0, 0, 0, 1, 0, 1, 0})) = 8;
    psi(vi({0, 1, 0, 0, 1, 0, 1, 0})) = 8;
    psi(vi({0, 0, 1, 0, 1, 0, 1, 0})) = 6;
    psi(vi({0, 1, 1, 0, 1, 0, 1, 0})) = 3;
    psi(vi({0, 0, 0, 1, 1, 0, 1, 0})) = 9;
    psi(vi({0, 1, 0, 1, 1, 0, 1, 0})) = 6;
    psi(vi({0, 0, 1, 1, 1, 0, 1, 0})) = 10;
    psi(vi({0, 1, 1, 1, 1, 0, 1, 0})) = 5;
    psi(vi({0, 0, 0, 0, 0, 1, 1, 0})) = 7;
    psi(vi({0, 1, 0, 0, 0, 1, 1, 0})) = 3;
    psi(vi({0, 0, 1, 0, 0, 1, 1, 0})) = 10;
    psi(vi({0, 1, 1, 0, 0, 1, 1, 0})) = 3;
    psi(vi({0, 0, 0, 1, 0, 1, 1, 0})) = 10;
    psi(vi({0, 1, 0, 1, 0, 1, 1, 0})) = 4;
    psi(vi({0, 0, 1, 1, 0, 1, 1, 0})) = 1;
    psi(vi({0, 1, 1, 1, 0, 1, 1, 0})) = 1;
    psi(vi({0, 0, 0, 0, 1, 1, 1, 0})) = 9;
    psi(vi({0, 1, 0, 0, 1, 1, 1, 0})) = 1;
    psi(vi({0, 0, 1, 0, 1, 1, 1, 0})) = 8;
    psi(vi({0, 1, 1, 0, 1, 1, 1, 0})) = 5;
    psi(vi({0, 0, 0, 1, 1, 1, 1, 0})) = 7;
    psi(vi({0, 1, 0, 1, 1, 1, 1, 0})) = 10;
    psi(vi({0, 0, 1, 1, 1, 1, 1, 0})) = 6;
    psi(vi({0, 1, 1, 1, 1, 1, 1, 0})) = 8;
    
    //output_tensortype(psi);
    //output_tensorfull(psi);
    auto tmppp = contract(psi,psi,ar::a_123456,ar::a_123456);
    double norm_squared = std::real(simplify_constant_tensor(tmppp));
    for(int i = 0; i < 2; ++i)
    for(int j = 0; j < 2; ++j)
    for(int k = 0; k < 2; ++k)
    for(int l = 0; l < 2; ++l)
    for(int l1 = 0; l1 < 2; ++l1)
    for(int l2 = 0; l2 < 2; ++l2){
        psi(vi({0,i,j,k,l,l1,l2,0})) = psi(vi({0,i,j,k,l,l1,l2,0}))/sqrt(norm_squared);
    }
    //output_tensortype(psi);
    //output_tensorfull(psi);

    std::vector<tensor<cd,4> > mpoHam;
    mpoHam.push_back( mpo::startH );
    mpoHam.push_back( mpo::middleH );
    mpoHam.push_back( mpo::middleH );
    mpoHam.push_back( mpo::middleH );
    mpoHam.push_back( mpo::middleH );
    mpoHam.push_back( mpo::endH );

    DMRG<8> df(psi,mpoHam,(double)1e-4);
    df.left_sweep();
    df.right_sweep();
    df.left_sweep();
    df.output_eigenvalue_history();

    auto t = df.right_normalized_final_state();
    auto t1 = contract(t[0],
                contract(t[1],
                    contract(t[2],
                        contract(t[3],
                            contract(t[4],t[5],ar::two,ar::zero),
                        ar::two,ar::zero),
                    ar::two,ar::zero),
                ar::two,ar::zero),
               ar::two,ar::zero);
    output_tensortype(t1);
    output_tensorfull(t1);

    return 0;
}
/*
#include <iostream>

#include "tensor.h"
#include "DMRG.h"
#include "utilities.h"

int main()
{
    mpo::initialize();

    tensor<cd,6> psi(vi({1,2,2,2,2,1}));
    psi(vi({0, 0, 0, 0, 0, 0})) = 1;
    psi(vi({0, 1, 0, 0, 0, 0})) = 1;
    psi(vi({0, 0, 1, 0, 0, 0})) = 1;
    psi(vi({0, 1, 1, 0, 0, 0})) = 1;
    psi(vi({0, 0, 0, 1, 0, 0})) = 1;
    psi(vi({0, 1, 0, 1, 0, 0})) = 1;
    psi(vi({0, 0, 1, 1, 0, 0})) = 1;
    psi(vi({0, 1, 1, 1, 0, 0})) = 1;
    psi(vi({0, 0, 0, 0, 1, 0})) = 1;
    psi(vi({0, 1, 0, 0, 1, 0})) = 1;
    psi(vi({0, 0, 1, 0, 1, 0})) = 1;
    psi(vi({0, 1, 1, 0, 1, 0})) = 1;
    psi(vi({0, 0, 0, 1, 1, 0})) = 1;
    psi(vi({0, 1, 0, 1, 1, 0})) = 1;
    psi(vi({0, 0, 1, 1, 1, 0})) = 1;
    psi(vi({0, 1, 1, 1, 1, 0})) = 1;

    //output_tensortype(psi);
    //output_tensorfull(psi);
    auto tmppp = contract(psi,psi,ar::a_1234,ar::a_1234);
    double norm_squared = std::real(simplify_constant_tensor(tmppp));
    for(int i = 0; i < 2; ++i)
    for(int j = 0; j < 2; ++j)
    for(int k = 0; k < 2; ++k)
    for(int l = 0; l < 2; ++l){
        psi(vi({0,i,j,k,l,0})) = psi(vi({0,i,j,k,l,0}))/sqrt(norm_squared);
    }
    //output_tensortype(psi);
    //output_tensorfull(psi);

    std::vector<tensor<cd,4> > mpoHam;
    mpoHam.push_back( mpo::startH );
    mpoHam.push_back( mpo::middleH );
    mpoHam.push_back( mpo::middleH );
    mpoHam.push_back( mpo::endH );

    DMRG<6> df(psi,mpoHam,(double)1e-4);
    df.left_sweep();
    df.output_eigenvalue_history();

    auto t = df.right_normalized_final_state();
    auto t1 = contract(t[0],contract(t[1],contract(t[2],t[3],ar::two,ar::zero),ar::two,ar::zero),ar::two,ar::zero);
    output_tensortype(t1);
    output_tensorfull(t1);

    return 0;
}
*/


