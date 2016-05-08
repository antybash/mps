/*
 * Description: This is an implementation of tensors in a tensor network.
 * At the moment, there is only basic functionality. The end goal will be
 * to construct a working DMRG.
 *
 * Objects: 
 *      - Tensors are represented by   boost::multi_array
 *
 *     boost::multi_array<int,N>    ===    Tensor with N components
 *            T[a1][a2]..[aN]       ===         T_{a1,a2,...,aN}
 *
 * Methods:
 *      - Multiply two tensors
 *      - Contract a single tensor
 *      - Matrix <---> Tensor function 
 *
 * Compiling Notes:
 *      g++ -std=c++11 ma.cpp
 *      g++ -std=c++11 ma.cpp -larmadillo (if doing svd decomposition)
 * Dependencies:
 *      Armadillo ( http://arma.sourceforge.net/download.html              )
 *                ( Install procedure: untar, configure, make, make install)
 *
 * TODO:
 *      1) Implement an MPO
 *          (use multi_array with functions rather than constants?)
 *      2) Implement local variational procedure (for single site in DMRG)
 *      3) Write down the Sweeping procedure
 *      4) Done!?
 */

#include <iostream>
#include <vector>
#include <boost/multi_array.hpp>
#include <array>
#include <functional>
#include <utility>
#include <armadillo>

#include "tensor.h"

// Also *REMEMBER* to use     " g++ .... -larmadillo  "
//      for Singular Value Decomposition
/* 
 * TODO: 
 *       Sweeping algorithm for DMRG
 *       Need to be able to represent Hamiltonian in an MPO state
 *
 */

/*
 * 1. Inital right-normalized guess |psi>
 * 2. Calculate R-expressions iteratively $L-1 \to 1$.
 * 3. Right sweep (l=1 .. L-1)
 *      - solve EigProblem: M^{\sg_l}
 *      - left-normalize M^{\sg_l} (SVD)
 *      - multiply remaining matrices onto M^{\sg_{l+1}}
 *      - Build L iteratively
 * 4. Left sweep (l=L .. 2)
 */

int main(int argv, char** argc)
{
    // TESTS
    boost::multi_array<int,3> T2(boost::extents[2][2][3]);

    T2[0][0][0] = 1; T2[0][1][0] = 2; T2[1][0][0] = 3; T2[1][1][0] = 4;
    T2[0][0][1] = 1; T2[0][1][1] = 2; T2[1][0][1] = 3; T2[1][1][1] = 4; 
    T2[0][0][2] = 1; T2[0][1][2] = 2; T2[1][0][2] = 3; T2[1][1][2] = 4; 

    //output(contract(T2,0,1)); // It works! :D

    // Testing tensor_to_matrix function
    auto M = tensor_to_matrix(T2,std::vector<int>{0,1},std::vector<int>{2});
    copy(M.begin(),M.end(), std::ostream_iterator<int>(std::cout," "));
    std::cout << std::endl;

    return 0;
}
