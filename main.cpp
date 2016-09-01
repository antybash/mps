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

    /*

    auto mps = tensor_to_left_normalized_mps(psi,1e-4); // change to trimmed later
    assert(mps.size() == 4);

    for(auto k : mps){
        output_tensortype(k);
        output_tensorfull(k);
    }

    // check left-canonical
    int lch = check_left_canoncial(mps);
    int rch = check_right_canoncial(mps);

    if (lch == -1)
        std::cout << "Left canonical." << std::endl;
    else 
        std::cout << "Not left canonical at " << lch << std::endl;

    if (rch == -1)
        std::cout << "Right canonical." << std::endl;
    else 
        std::cout << "Not right canonical at " << rch << std::endl;
    */

//
//   auto eig = DMRG_sweep(mps,mpoHam);
//   for(auto i : eig)
//       std::cout << i << std::endl;
//
//   output_tensorfull( contract(mps[0],contract(mps[1],contract(mps[2],mps[3],ar::two,ar::zero),ar::two,ar::zero),ar::two,ar::zero) ); 

    return 0;
}

