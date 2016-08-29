
#ifndef MPS_DMRG_CONTRACTIONS_H
#define MPS_DMRG_CONTRACTIONS_H

#include "tensor.h"
#include <array>
#include <vector>
#include <complex>

// Contracting < psi | psi >
// from the left

template<typename T>
tensor<T,2>
DMRG_double_left_contract_once (tensor<T,2> &A, tensor<T,3> &P)
{
    std::array<int,1> x = {{ 0 }};
    tensor<T,3> T1 = contract(A,P,x,x,false,true);
    // (A0 A1) (P0* P1* P2*)  ---> (A1 P1* P2*)

    std::array<int,2> y = {{0,1}};
    tensor<T,2> T2 = contract(T1,P,y,y);
    // (A1 P1* P2*) (P0 P1 P2)---> (P2* P2)

    return T2;
}

template<typename T>
std::vector<tensor<T,2> >
DMRG_double_left_recursive (std::vector<tensor<T,3> > &mpsState)
{
    typedef std::vector<int> vi;
    typedef std::complex<double> cd;

    std::vector<tensor<T,2> > vectorLeft;

    vectorLeft.push_back( tensor<T,2>(vi({1,1})) );
    vectorLeft[0](vi({0,0})) = 1;

    for (int i = 1; i <= mpsState.size(); ++i)
        vectorLeft.push_back( DMRG_double_left_contract_once(vectorLeft[i-1],mpsState[i-1]) );
    
    return vectorLeft;
}

template<typename T>
tensor<T,2>
DMRG_double_right_contract_once (tensor<T,3> P, tensor<T,2> A)
{
    std::array<int,1> x = {{ 2 }};
    std::array<int,1> y = {{ 0 }};
    tensor<T,3> T1 = contract(P,A,x,y,false,false);
    // (P0 P1 P2) (A0 A1)    ---> T1 = (P0 P1 A1)

    std::array<int,2> z = {{1,2}};
    tensor<T,2> T2 = contract(P,T1,z,z,true,false);
    // (P0* P1* P2*) (P0 P1 A1) -> T2 = (P0* P0)

    return T2;
    //done
}

template<typename T>
std::vector<tensor<T,2> >
DMRG_double_right_recursive (std::vector<tensor<T,3> > mpsState)
{
    typedef std::vector<int> vi;
    typedef std::complex<double> cd;

    int L = mpsState.size();

    std::vector<tensor<T,2> > vectorRight;
    vectorRight.push_back( tensor<T,2>(vi ( {{1,1}} )) );
    vectorRight[0](vi({{0,0}})) = 1;

    for (int i = 1; i <= L; ++i)
        vectorRight.push_back( DMRG_double_right_contract_once(mpsState[L-i],vectorRight[i-1]) );

    return vectorRight;
}



template<typename T>
void
DMRG_double_update_site(int new_site, std::vector<tensor<T,3> > &mpsState, 
        std::vector<tensor<T,2> > &vectorLeft, std::vector<tensor<T,2> > &vectorRight)
{
    // Suppose mpsState[new_site] has just been updated.
    // This function will update the corresponding vectorLeft and
    // vectorRight by calling the DMRG_contract_once on each of them

    int L   = mpsState.size();
    int vli = new_site+1;
    int vri = L-new_site;

    vectorLeft [vli] = DMRG_double_left_contract_once  (vectorLeft[vli-1],  mpsState[new_site]); // check this with DMRG function
    vectorRight[vri] = DMRG_double_right_contract_once (mpsState[new_site], vectorRight[vri-1]);
    /*
     *     L = 6
     *    new_site   ==     0 1 2 3 4 5
     *   vectorLeft  ==   0 1 2 3 4 5 6
     *   vectorRight ==     6 5 4 3 2 1 0
     *
     *                    i      i      i = new_site
     *                   i+1    i+1
     *                   L-i    6-i
     */
}


/***************************************************************/


// Contracting < psi | H | psi >
// from the left

template<typename T>
tensor<T,3>
DMRG_triple_left_contract_once (tensor<T,3> A, tensor<T,3> P, tensor<T,4> H)
{

    /*
     *    |0       0 -------- 2
     *    |             1             2
     *    |1                    0 --------- 3
     *    |             1             1
     *    |2       0 -------- 2
     *
     */

    std::array<int,1> x1 = {{ 0 }};
    std::array<int,1> y1 = {{ 0 }};
    tensor<T,4> T1 = contract(A,P,x1,y1,false,true); // conjugate P
    
    std::array<int,2> x2 = {{0,2}};
    std::array<int,2> y2 = {{0,2}};
    tensor<T,4> T2 = contract(T1,H,x2,y2);

    std::array<int,2> x3 = {{0,2}};
    std::array<int,2> y3 = {{0,1}};
    tensor<T,3> T3 = contract(T2,P,x3,y3);

    return T3;
}

template<typename T>
std::vector<tensor<T,3> >
DMRG_triple_left_recursive (std::vector<tensor<T,3> > mpsState, std::vector<tensor<T,4> > mpsHamiltonian)
{
    typedef std::vector<int> vi;
    typedef std::complex<double> cd;

    std::vector<tensor<T,3> > vectorLeft;

    vectorLeft.push_back( tensor<T,3>(vi({{1,1,1}})) );
    vectorLeft[0](vi({{0,0,0}})) = 1;
    
    for (int i = 1; i <= mpsState.size(); ++i)
        vectorLeft.push_back( DMRG_triple_left_contract_once(vectorLeft[i-1],mpsState[i-1],mpsHamiltonian[i-1]) );
    
    return vectorLeft;
}

template<typename T>
tensor<T,3>
DMRG_triple_right_contract_once(tensor<T,3> P, tensor<T,4> H, tensor<T,3> A)
{
    /*
     *    |0       0 -------- 2                  0|
     *    |             1             2           |
     *    |1                    0 --------- 3    1|
     *    |             1             1           |
     *    |2       0 -------- 2                  2|
     *
     */

    std::array<int,1> x = {{ 2 }};
    tensor<T,4> T1 = contract(P,A,x,x);
    // T1 = (P0 P1 A0 A1)

    std::array<int,2> y = {{ 1,3 }};
    tensor<T,4> T2 = contract(H,T1,y,y);
    // (H0 H1 H2 H3) (P0 P1 A0 A1)  --> T2 = (H0 H2 P0 A0)

    std::array<int,2> z1 = {{ 1,2 }};
    std::array<int,2> z2 = {{ 1,3 }};
    tensor<T,3> T3 = contract(P,T2,z1,z2,true,false);
    // (P0* P1* P2*) (H0 H2 P0 A0)  --> T3 = (P0* H0 P0)

    return T3;
}

template<typename T>
std::vector<tensor<T,3> >
DMRG_triple_right_recursive (std::vector<tensor<T,3> > mpsState, std::vector<tensor<T,4> > mpsHamiltonian)
{
    typedef std::vector<int> vi;
    typedef std::complex<double> cd;

    int L = mpsState.size();

    std::vector<tensor<T,3> > vectorRight;
    vectorRight.push_back( tensor<T,3>(vi({{1,1,1}})) );
    vectorRight[0](vi({{0,0,0}})) = 1;

    for (int i = 1; i <= L; ++i)
        vectorRight.push_back( DMRG_triple_right_contract_once(mpsState[L-i], mpsHamiltonian[L-i],vectorRight[i-1]) );

    return vectorRight;
}

template<typename T>
void
DMRG_triple_update_site(int new_site, std::vector<tensor<T,3> > &mpsState, std::vector<tensor<T,4> > &mpsHamiltonian,
        std::vector<tensor<T,3> > &vectorLeft, std::vector<tensor<T,3> > &vectorRight)
{
    // Suppose mpsState[new_site] has just been updated.
    // This function will update the corresponding vectorLeft and
    // vectorRight by calling the DMRG_contract_once on each of them

    int L = mpsState.size();
    int vli = new_site+1;   // these are the vectorLeft & vectorRight indices corresponding to the mpsState[new_site]
    int vri = L-new_site;

//    std::cout << "DMRG_triple_update_site started --" << std::endl;
//
//    std::cout << "mpsState" << std::endl;       output_tensorfull(mpsState[new_site]);
//    std::cout << "mpsHamiltonian" << std::endl; output_tensorfull(mpsHamiltonian[new_site]);
//    std::cout << "vectorLeft" << std::endl;     output_tensorfull(vectorLeft[vli]);
//    std::cout << "vectorRight" << std::endl;    output_tensorfull(vectorRight[vri]);

      vectorLeft [vli] = DMRG_triple_left_contract_once (vectorLeft[vli-1],    mpsState[new_site],       mpsHamiltonian[new_site]); 
      vectorRight[vri] = DMRG_triple_right_contract_once(mpsState[new_site], mpsHamiltonian[new_site], vectorRight[vri-1]);

//    std::cout << "DMRG_triple_update_site finished!" << std::endl;

    /*
     *     L = 5
     *    new_site   ==     0 1 2 3 4 5
     *   vectorLeft  ==   0 1 2 3 4 5 6
     *   vectorRight ==     6 5 4 3 2 1 0
     *
     */
}

///////////////////////////////////////////////////////////////////////////////////////////////

template<typename T>
double
DMRG_current_eigenvalue (std::vector<tensor<T,3> > &mpsState, std::vector<tensor<T,4> > mpoHam)
{
    /*********************************
     *            <psi | H | psi>    *
     * Computes   ---------------    *
     *              <psi | psi>      *
     *********************************/
    auto numer = DMRG_triple_left_recursive(mpsState,mpoHam);
    auto denom = DMRG_double_left_recursive(mpsState);

    return std::real(simplify_constant_tensor(*numer.rbegin()) /
                     simplify_constant_tensor(*denom.rbegin())   );
}


#endif //MPS_DMRG_CONTRACTIONS_H
