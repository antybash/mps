
#ifndef DMRG_H
#define DMRG_H

#include <iostream>
#include <iterator>
#include <Eigen/Dense>
#include <vector>
#include <complex>
#include <stack>

#include "tensor.h"
#include "utilities.h"

typedef std::complex<double> cd;

template<size_t N>
class DMRG {
    private:
        double epsilon;

        std::stack<tensor<cd,3> > L1; // oL == 1L
        std::stack<tensor<cd,3> > R1; //     -----> left and right-normalized states

        std::stack<tensor<cd,2> > L2; // to compute the norm
        std::stack<tensor<cd,2> > R2;

        std::stack<tensor<cd,3> > L3; // tL == 3L  
        std::stack<tensor<cd,3> > R3; //     -----> L,R in LWR M = Psi_L Psi_R M

        std::vector<tensor<cd,4> > hamiltonian;

        std::vector<double> eigenvalue_history;

        std::tuple<Eigen::VectorXcd,double> eigen_solve(int);
        tensor<cd,3> eigenvector_to_mps(Eigen::VectorXcd, tensor<cd,3>);
        void left_push (Eigen::MatrixXcd, std::vector<int>, int);  // given a properly normalized eigentensor
        void right_push(Eigen::MatrixXcd, std::vector<int>, int);  // push the left/right states

        void right_sweep_once(); // doing one update loses information (and hence this method is private)
        void left_sweep_once();  // must update the entire sequence, and keep solving eigensolving!

    public:
        DMRG (tensor<cd,N>, std::vector<tensor<cd,4> >, double);

        void right_sweep(); 
        void left_sweep(); 

        void output_eigenvalue_history();
        std::vector<tensor<cd,3> > right_normalized_final_state();
};

template<size_t N>
void
DMRG<N>::output_eigenvalue_history()
{
    std::cout << "The eigenvalue history is: ";
    for(int i = 0; i < eigenvalue_history.size(); ++i)
        std::cout << eigenvalue_history[i] << ", ";
    std::cout << std::endl;
}

template<size_t N>
std::vector<tensor<cd,3> >
DMRG<N>::right_normalized_final_state()
{
    std::stack <tensor<cd,3> > right = R1;
    std::vector<tensor<cd,3> > x;
    for(int i = 0; i < N-2; ++i){
        x.push_back(R1.top());
        R1.pop();
    }
    return x;
}


template<size_t N>
std::tuple<Eigen::VectorXcd, double>
DMRG<N>::eigen_solve (int pos)
{
    tensor<cd,5> T1 = contract(L3.top(), hamiltonian[pos],   ar::one,  ar::zero); // (l0 l1 l2) (h0 h1 h2 h3)  ----->  (l0 l2 h1 h2 h3)
    tensor<cd,6> T2 = contract(T1,       R3.top(),           ar::four, ar::one);  // (l0 l2 h1 h2 h3) (r0 r1 r2) --->  (l0 l2 h1 h2 r0 r2)
    std::array<int,3> rows = {{ 0,3,4 }};
    std::array<int,3> cols = {{ 1,2,5 }};

    Eigen::MatrixXcd H = tensor_to_matrix(T2,rows,cols);  //TODO: check self-adjoint-ness
    std::cout << H << std::endl << std::endl;

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es;
    es.compute(H);

    return std::make_tuple(es.eigenvectors().col(0), es.eigenvalues()[0]);
}

template<size_t N>
tensor<cd,3>
DMRG<N>::eigenvector_to_mps(Eigen::VectorXcd vec, tensor<cd,3> prevState)
{
    std::vector<int> shape(3);
    shape[0] = prevState.shape()[0];
    shape[1] = prevState.shape()[1];
    shape[2] = prevState.shape()[2];
    tensor<cd,3> T(shape);

    std::vector<int> index;
    for(int i = 0; i < vec.size(); ++i){
        index = number_to_multi_index(i,shape);
        T(index) = vec[i];
    }
    return T;
}

template<size_t N>
void 
DMRG<N>::left_push (Eigen::MatrixXcd tmp, std::vector<int> tmp_shape, int ham_pos)
{
    auto H = hamiltonian[ham_pos];
    auto tU = matrix_to_tensor(tmp, ar::zeroone, ar::two, tmp_shape);
    auto x = L2.top();
    auto y = L3.top();

    auto x1 = contract(x,tU,ar::zero,ar::zero,false,true);
    auto x2 = contract(x1,tU,ar::zeroone,ar::zeroone,false,false);

    auto y1 = contract(y,tU,    ar::zero,   ar::zero,false,true);
    auto y2 = contract(y1,H,    ar::zerotwo,ar::zerotwo);
    auto y3 = contract(y2,tU,   ar::zerotwo,ar::zeroone);

    L1.push( tU );
    L2.push( x2 );
    L3.push( y3 );
}

template<size_t N>
void 
DMRG<N>::right_push(Eigen::MatrixXcd tmp, std::vector<int> tmp_shape, int ham_pos)
{
    auto H = hamiltonian[ham_pos];

    auto tU = matrix_to_tensor(tmp, ar::zero, ar::onetwo, tmp_shape); // TODO: check the indices
    auto x = R2.top();
    auto y = R3.top();

    auto x1 = contract(tU,x,ar::two,ar::zero,true,false);
    auto x2 = contract(x1,tU,ar::onetwo,ar::onetwo,false,false);

    auto y1 = contract(tU,y,    ar::two,      ar::two); 
    auto y2 = contract(H,y1,    ar::onethree, ar::onethree);
    auto y3 = contract(tU,y2,   ar::onetwo,   ar::onethree, true, false);

    R1.push( tU );
    R2.push( x2 );
    R3.push( y3 );
}


template<size_t N>
void DMRG<N>::right_sweep_once()
{
    int k = L1.size(); // this is the position that is of interest to us: mpoH[k] 

    auto prevState = R1.top();
    R1.pop();
    R2.pop();
    R3.pop();
    
    //TODO:check that it's an eigenvalue problem => Psi_L and Psi_R

    Eigen::VectorXcd vec;
    double eigenval;
    std::tie(vec,eigenval) = eigen_solve(k);
    eigenvalue_history.push_back(eigenval);

    tensor<cd,3> eigentensor = eigenvector_to_mps( vec, prevState );
    // CANNOT just add eigenvector to L1 and update L2, L3: we must left-normalize it first!

    auto M = tensor_to_matrix(eigentensor, ar::zeroone, ar::two);
    Eigen::MatrixXcd U,S,V;
    std::tie(U,S,V) = custom_svd(M);
    std::vector<int> prev_shape = { prevState.shape()[0], prevState.shape()[1], prevState.shape()[2] }; 
    left_push( U, prev_shape, k ); // This is the update of L1, L2, L3.
                                   // TODO: check whether the exist state
                                   // satisfies normalization conditions.
    assert( L1.size() == k+1 );
    assert( L2.size() == 1+L1.size() );
    assert( L3.size() == 1+L1.size() );
}

template<size_t N>
void DMRG<N>::right_sweep()
{
    for(int i = 0; i < N-2; ++i)
        right_sweep_once();
}

template<size_t N>
void 
DMRG<N>::left_sweep_once()
{
    int k = L1.size()-1; // this is the position that is of interest to us: mpoH[k] 
    auto prevState = L1.top();
    L1.pop();
    L2.pop();
    L3.pop();
    
    //TODO:check that it's an eigenvalue problem => Psi_L times Psi_R equals identity 
    Eigen::VectorXcd vec;
    double eigenval;
    std::tie(vec,eigenval) = eigen_solve(k);
    eigenvalue_history.push_back(eigenval);

    tensor<cd,3> eigentensor = eigenvector_to_mps( vec, prevState );

    // cannot just add eigenvector to R1 and update R2, R3: we must 
    // **RIGHT-NORMALIZE** it first!
    auto M = tensor_to_matrix(eigentensor, ar::zero, ar::onetwo);
    Eigen::MatrixXcd U,S,V;
    std::tie(U,S,V) = custom_svd(M);
    std::vector<int> prev_shape({(int)prevState.shape()[0], (int)prevState.shape()[1], (int)prevState.shape()[2] }); 
    right_push( V, prev_shape, k ); // This is the update of L1, L2, L3.
                                    // TODO: check whether the exist state
                                    // satisfies normalization conditions.
    assert( L1.size() + R1.size() == N-2 );
    assert( R1.size() == R2.size()-1 );
    assert( R1.size() == R3.size()-1 );
}


template<size_t N>
void DMRG<N>::left_sweep()
{
    for(int i = 0; i < N-2; ++i)
        left_sweep_once();
}

template<size_t N>
DMRG<N>::DMRG ( tensor<cd,N> psi, std::vector<tensor<cd,4> > H, double eps )
{
    epsilon = eps;
    hamiltonian = H;

    tensor<cd,2> f2(vi({1,1}));   f2(vi({0,0}))   = 1; L2.push(f2); R2.push(f2);  // initialize stacks
    tensor<cd,3> f3(vi({1,1,1})); f3(vi({0,0,0})) = 1; L3.push(f3); R3.push(f3);

    Eigen::MatrixXcd tmp       = tensor_to_matrix(psi, ar::zeroone, array_range<2,N-1>());
    std::vector<int> tmp_shape = tensor_shape(psi);

    int trim;
    for(int i = 0; i <= N-3; ++i){ 
        Eigen::MatrixXcd U;
        Eigen::MatrixXcd V;
        if (i < N-3) {
            std::tie (trim, U, V) = svd_then_trim(tmp, epsilon);
            auto tU = matrix_to_tensor(U, ar::zeroone, ar::two, triple(tmp_shape[0],tmp_shape[1],trim));
            auto x = L2.top();
            auto y = L3.top();

            auto x1 = contract(x,tU,ar::zero,ar::zero,false,true);
            auto x2 = contract(x1,tU,ar::zeroone,ar::zeroone,false,false);

            auto y1 = contract(y,tU,    ar::zero,   ar::zero,false,true); // conjugate P
            auto y2 = contract(y1,H[i], ar::zerotwo,ar::zerotwo);
            auto y3 = contract(y2,tU,   ar::zerotwo,ar::zeroone);

            L1.push( tU );
            L2.push( x2 );
            L3.push( y3 );

            tmp_shape.erase(tmp_shape.begin());          // 01234 -> 1234     A(01)(234) --> U(01)(s) V(s)(234)
            tmp_shape[0] = trim;                         // 1234  -> s234                --> U[01s]   V[s,2,3,4]
            tmp          = inplace_index_swap_of_underlying_tensor(V, tmp_shape); // (s,234) -> (s2,34)

        } else {
            auto tU = matrix_to_tensor(tmp, ar::zeroone, ar::two, tmp_shape);
            auto x = L2.top();
            auto y = L3.top();

            auto x1 = contract(x,tU,ar::zero,ar::zero,false,true);
            auto x2 = contract(x1,tU,ar::zeroone,ar::zeroone,false,false);

            auto y1 = contract(y,tU,    ar::zero,   ar::zero,false,true); // conjugate P
            auto y2 = contract(y1,H[i], ar::zerotwo,ar::zerotwo);
            auto y3 = contract(y2,tU,   ar::zerotwo,ar::zeroone);

            L1.push( tU );
            L2.push( x2 );
            L3.push( y3 );

        }
    }
}


#endif
