
#ifndef DMRG_H
#define DMRG_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <Eigen/Dense>
#include <vector>
#include <complex>
#include <stack>
#include <set>

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

        std::tuple<double, Eigen::VectorXcd, std::vector<int> > eigen_solve(int);
        tensor<cd,3> eigenvector_to_mps  (Eigen::VectorXcd, std::vector<int>);
        tensor<cd,3> left_normalize_mps  (tensor<cd,3>);
        tensor<cd,3> right_normalize_mps (tensor<cd,3>);
        tensor<cd,3> normalize_mps       (tensor<cd,3>, bool);

        void left_push (tensor<cd,3>, int);  // given a properly normalized eigentensor
        void right_push(tensor<cd,3>, int);  // push the left/right states

        void right_sweep_once(); // doing one update loses information (and hence this method is private)
        void left_sweep_once();  // must update the entire sequence, and keep solving eigensolving!

    public:
        //DMRG (tensor<cd,N>, std::vector<tensor<cd,4> >, double);                // (psi, mpo, cutoff)
        //DMRG (int, std::vector<tensor<cd,4> >, double);                           // (mpo, cutoff)

        DMRG (int, double); // assume heisenberg XYZ

        void right_sweep(); 
        void left_sweep(); 

        void sweep(int);

        void output_eigenvalue_history(std::ofstream &);
        void output_lowest_energy(std::ofstream & , bool );
        void output_lowest_energy_cout();


        std::vector<tensor<cd,3> > final_mps_state();

        void output_current_norm();

        //tensor<cd,N> contract_mps();
};

template<size_t N>
void
DMRG<N>::output_current_norm()
{
    if (R1.size() != 0) {
        std::cout << "The norm of the final_state is: "
                  << simplify_constant_tensor(R2.top())
                  << std::endl;
    } else {
        std::cout << "The norm of the final_state is: "
                  << simplify_constant_tensor(L2.top())
                  << std::endl;
    }
}

template<size_t N>
void
DMRG<N>::output_lowest_energy(std::ofstream & outfile, bool standard_output)
{
    outfile << "" << N-2 << " " << std::fixed << std::setprecision(7) << *(eigenvalue_history.rbegin()) << std::endl;
    if(standard_output)
        std::cout << "" << N-2 << " " << std::fixed << std::setprecision(7) << *(eigenvalue_history.rbegin()) << std::endl;
}

template<size_t N>
void
DMRG<N>::output_lowest_energy_cout()
{
    std::cout << "" << N-2 << " " << *(eigenvalue_history.rbegin()) << std::endl;
}

template<size_t N>
void
DMRG<N>::output_eigenvalue_history(std::ofstream & outfile)
{
    std::vector<double> t2(eigenvalue_history.begin(),eigenvalue_history.end());
    std::unique(t2.begin(), t2.end());

    outfile << "The eigenvalue history is: " << std::endl;
    for(int i = 0; i < t2.size(); ++i)
        outfile << std::fixed << std::setprecision(10) << t2[i] << std::endl;
    outfile << std::endl;
}

template<size_t N>
std::vector<tensor<cd,3> >
DMRG<N>::final_mps_state()
{
    if (R1.size() != 0) {
        std::stack <tensor<cd,3> > right = R1;
        std::vector<tensor<cd,3> > x;
        for(int i = 0; i < N-2; ++i){
            x.push_back(right.top());
            right.pop();
        }
        return x;
    } else {
        std::stack <tensor<cd,3> > left = L1;
        std::vector<tensor<cd,3> > x;
        for(int i = 0; i < N-2; ++i){
            x.push_back(left.top());
            left.pop();
        }
        return x;
    }
}

template<size_t N>
std::tuple<double, Eigen::VectorXcd, std::vector<int> >
DMRG<N>::eigen_solve (int pos)
{
    tensor<cd,5> T1 = contract(L3.top(), hamiltonian[pos],   ar::one,  ar::zero); // (l0 l1 l2) (h0 h1 h2 h3)  ----->  (l0 l2 h1 h2 h3)
    tensor<cd,6> T2 = contract(T1,       R3.top(),           ar::four, ar::one);  // (l0 l2 h1 h2 h3) (r0 r1 r2) --->  (l0 l2 h1 h2 r0 r2)

    std::array<int,3> rows = {{ 0,3,4 }};
    std::array<int,3> cols = {{ 1,2,5 }};

    Eigen::MatrixXcd H = tensor_to_matrix(T2,rows,cols);  //TODO: check self-adjoint-ness

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es;
    es.compute(H);

    return std::make_tuple(
            es.eigenvalues()[0],
            es.eigenvectors().col(0), 
            vi({ (int)T2.shape()[1], (int)T2.shape()[3], (int)T2.shape()[5]  }) );
}

////////////////////////////////////////////////////////////////////////
// This section is to take a solution to the eigenvalue equation      //
// and produce the correct mps that is will replace the old one.      //
////////////////////////////////////////////////////////////////////////

template<size_t N>
tensor<cd,3>
DMRG<N>::eigenvector_to_mps(Eigen::VectorXcd vec, std::vector<int> shape)
{
    tensor<cd,3> T(shape);

    std::vector<int> index;
    for(int i = 0; i < vec.size(); ++i){
        index = number_to_multi_index(i,shape);
        T(index) = vec[i];
    }
    return T;
}

template<size_t N>
tensor<cd,3>
DMRG<N>::left_normalize_mps(tensor<cd,3> mps)
{
    Eigen::MatrixXcd U,S,V;
    auto M = tensor_to_matrix(mps, ar::zeroone, ar::two);   // mps_(sσ)(s') = M_αβ
    std::tie(U,S,V) = custom_svd(M);

    int ts = vector_trim_size(S.diagonal(), epsilon);
    U.conservativeResize(U.rows(), ts);                      //// CONSERVATIVE RESIZE !!!!! ////

    std::vector<int> shape = { (int) mps.shape()[0], (int) mps.shape()[1], ts };
    return matrix_to_tensor(U, ar::zeroone, ar::two, shape);
}

template<size_t N>
tensor<cd,3>
DMRG<N>::right_normalize_mps(tensor<cd,3> mps)
{
    Eigen::MatrixXcd U,S,V;
    auto M = tensor_to_matrix(mps, ar::zero, ar::onetwo);   // mps_(s)(σs') = M_αβ
    std::tie(U,S,V) = custom_svd(M);

    int ts = vector_trim_size(S.diagonal(), epsilon);
    V.conservativeResize(ts, V.cols());                      //// CONSERVATIVE RESIZE !!!!! ////

    std::vector<int> shape = { ts, (int) mps.shape()[1], (int) mps.shape()[2] };
    return matrix_to_tensor(V, ar::zero, ar::onetwo, shape);
}

template<size_t N>
tensor<cd,3> 
DMRG<N>::normalize_mps(tensor<cd,3> mps, bool left_or_right)
{
    if (left_or_right == 0) // left
        return left_normalize_mps (mps);
    else                    // right
        return right_normalize_mps(mps);
}



template<size_t N>
void 
DMRG<N>::left_push (tensor<cd,3> mps, int ham_pos) // for right-sweep
{
    auto H = hamiltonian[ham_pos];
    auto x = L2.top();
    auto y = L3.top();

    auto x1 = contract(x, mps,ar::zero,ar::zero,false,true);
    auto x2 = contract(x1,mps,ar::zeroone,ar::zeroone,false,false);

    auto y1 = contract(y,mps,    ar::zero,   ar::zero,false,true);
    auto y2 = contract(y1,H,     ar::zerotwo,ar::zerotwo);
    auto y3 = contract(y2,mps,   ar::zerotwo,ar::zeroone);

    L1.push( mps);
    L2.push( x2 );
    L3.push( y3 );
}

template<size_t N>
void 
DMRG<N>::right_push(tensor<cd,3> mps, int ham_pos) // for left-sweep
{
    auto H = hamiltonian[ham_pos];

    auto x = R2.top();
    auto y = R3.top();

    auto x1 = contract(mps,x,ar::two,ar::zero,true,false);
    auto x2 = contract(x1,mps,ar::onetwo,ar::onetwo,false,false);

    auto y1 = contract(mps,y,    ar::two,      ar::two); 
    auto y2 = contract(H,y1,     ar::onethree, ar::onethree);
    auto y3 = contract(mps,y2,   ar::onetwo,   ar::onethree, true, false);

    R1.push( mps);
    R2.push( x2 );
    R3.push( y3 );
}


template<size_t N>
void DMRG<N>::right_sweep_once()
{
    int k = L1.size(); // this is the position that is of interest to us: mpoH[k] 
    R1.pop();
    R2.pop();
    R3.pop();

    //TODO:check that it's an eigenvalue problem => Psi_L and Psi_R
    Eigen::VectorXcd  vec;
    double            eigenval;
    std::vector<int>  new_shape;

    std::tie(eigenval,vec,new_shape) = eigen_solve(k);
    eigenvalue_history.push_back(eigenval);

    tensor<cd,3> eigentensor_not_normalized = eigenvector_to_mps( vec, new_shape );
    tensor<cd,3> eigentensor_normalized     = left_normalize_mps(eigentensor_not_normalized);

    left_push( eigentensor_normalized, k ); // This is the update of L1, L2, L3.  TODO: check whether the exit state satisfies normalization conditions.

    assert( L1.size() == k+1 ); assert( L2.size() == 1+L1.size() ); assert( L3.size() == 1+L1.size() );
}

template<size_t N>
void DMRG<N>::right_sweep()
{
    for(int i = 0; i < N-2; ++i){
        right_sweep_once();
    }
}

template<size_t N>
void DMRG<N>::left_sweep_once()
{
    int k = L1.size()-1; // this is the position that is of interest to us: mpoH[k] 
    L1.pop();
    L2.pop();
    L3.pop();

    //TODO:check that it's an eigenvalue problem => Psi_L times Psi_R equals identity 
    Eigen::VectorXcd vec;
    double eigenval;
    std::vector<int> new_shape;
    std::tie(eigenval,vec,new_shape) = eigen_solve(k);
    eigenvalue_history.push_back(eigenval);

    tensor<cd,3> eigentensor_not_normalized = eigenvector_to_mps ( vec, new_shape );
    tensor<cd,3> eigentensor_normalized     = right_normalize_mps(eigentensor_not_normalized);

    right_push( eigentensor_normalized, k ); // This is the update of R1, R2, R3. TODO: check whether the exit state satisfies normalization conditions.

    assert( L1.size() + R1.size() == N-2 ); assert( R1.size() == R2.size()-1 ); assert( R1.size() == R3.size()-1 );
}

template<size_t N>
void DMRG<N>::left_sweep()
{
    for(int i = 0; i < N-2; ++i){
        left_sweep_once();
    }
}

template<size_t N>
void DMRG<N>::sweep(int num_sweeps)
{
    for(int i = 0; i < num_sweeps; ++i){
        left_sweep();
        right_sweep();
    }
}

template<size_t N>
DMRG<N>::DMRG ( int D, double eps )
{

    epsilon = eps;
    
    ///////////// DECLARE MPO /////////////////

    std::vector<tensor<cd,4> > H;
    H.push_back( mpo_heis::startH );
    for(int i = 0; i < N-4; ++i)
        H.push_back( mpo_heis::middleH );
    H.push_back( mpo_heis::endH );

    hamiltonian = H;

    ///////////// DECLARE MPS /////////////////

    tensor<cd,2> f2(vi({1,1}));   f2(vi({0,0}))   = 1; L2.push(f2); R2.push(f2);  // initialize stacks
    tensor<cd,3> f3(vi({1,1,1})); f3(vi({0,0,0})) = 1; L3.push(f3); R3.push(f3);

    //Eigen::MatrixXcd tmp       = tensor_to_matrix(psi, ar::zeroone, array_range<2,N-1>());
    //std::vector<int> tmp_shape = tensor_shape(psi);

    std::vector<int> tmp_shape;
    int trim;
    Eigen::MatrixXcd U;
    Eigen::MatrixXcd V;

    for(int i = 0; i <= N-3; ++i){ 
        if (i < N-3) {
            if (i == 0)
                tmp_shape = triple(1,2,D);
            else
                tmp_shape = triple(trim,2,D);

            auto M = tensor_to_matrix(random_mps_site(tmp_shape), ar::zeroone, ar::two);
            std::tie (trim, U, V) = svd_then_trim(M, epsilon);
            
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

            //tmp_shape.erase(tmp_shape.begin());          // 01234 -> 1234     A(01)(234) --> U(01)(s) V(s)(234)
            //tmp_shape[0] = trim;                         // 1234  -> s234                --> U[01s]   V[s,2,3,4]
            //tmp          = inplace_index_swap_of_underlying_tensor(V, tmp_shape); // (s,234) -> (s2,34)

        } else {
            tmp_shape = triple(trim,2,1);
            //std::tie (trim, U, V) = svd_then_trim(random_mps_site(tmp_shape), epsilon);
            auto tU = random_mps_site(tmp_shape);
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

/*
template<size_t N>
tensor<cd,N-2> contract_mps_recursive<N,N-3>( tensor<cd,N-3> old, std::vector<tensor<cd,3> > rest )
{
    assert(rest.size() == 1);
    return contract(old, rest[0], ar::two, ar::zero);
}

template<size_t N, size_t M>
tensor<cd,M+1> contract_mps_recursive( tensor<cd,M> old, std::vector<tensor<cd,3> > rest )
{
    assert(N-M+1 == rest.size());
//   if (M+1 == N-2){ // redundant
//       assert(rest.size() == 1);
//       return contract(old, rest[0], ar::two, ar::zero);
//   }
//   else{ 
        tensor<cd,3> tmp = rest[N-M];
        rest.erase(rest.begin()+N-M);
        return contract_mps_recursive<N,M+1>( contract(tmp, old, ar::two, ar::zero ), rest ); 
//    }
}


template<size_t N>
tensor<cd,N> 
DMRG<N>::contract_mps()
{
    std::vector<tensor<cd,3> > rest = final_mps_state(); // size = N-2
    tensor<cd,3> first = rest[N-3];
    rest.erase(rest.begin()+N-3);
    return contract_mps_recursive<N,3>( first, rest );
}
*/



#endif
