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
 *      - Multiply two tensors:
 *              C_{a... b...} = A_{a...} B_{b...}
 *
 *      - Contract a single tensor:
 *              F'_{a,b...}   = F_{mu,mu,a,b,...}
 *
 *      - Represent a tensor in a matrix form (this is where SVD comes in!)
 *
 *              M_{ab} = F_{ (a1,a2,...) (b1,b2,...) }
 *
 *        the a_i's don't have to come before the b_j's, as long as there
 *        is a bipartition of the indices of F, this method handles the
 *        general case.
 *
 * Compiling Notes:
 *      g++ -std=c++11 ma.cpp
 *      g++ -std=c++11 ma.cpp -larmadillo (if doing svd decomposition)
 *
 *
 * TODO:
 *      0) Matrix -> Tensor method
 *          (Needed after the SVD procedure)
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

template<std::size_t N>  // N-dimensional tensor
void output(boost::multi_array<int,N> a)
{
    auto i = a.origin();
    for (int k = 1; k < a.num_elements()+1; ++k)
        std::cout << *(i++) << " "; 
    std::cout << std::endl;
}

std::vector<int> int_to_index(int k, std::vector<int> shape)
{
    // "i" enumerates the index abc
    // the following does: "i" -> abc
    int prod = std::accumulate(shape.begin(),shape.end(), 1, std::multiplies<int>());
    std::vector<int> vec_index(shape.size());
    for(int j = 0; j < shape.size(); ++j){
        vec_index[j] = (k% shape[j]);
        k -= vec_index[j]; 
        k /= shape[j];
    }
    return vec_index;
}

template<class Type, std::size_t N, std::size_t M> 
boost::multi_array<Type, N+M> 
combine_tensors(boost::multi_array<Type,N> t1, boost::multi_array<Type,M> t2)
{
    // declares the components of tensor
    std::vector<Type> index; index.reserve(N+M);
    for (auto p = t1.shape(); p != t1.shape()+N; ++p)
        index.push_back(*p);
    for (auto p = t2.shape(); p != t2.shape()+M; ++p)
        index.push_back(*p);

    boost::multi_array<Type, N+M> new_tensor(index);

    // eg. Z_abc = X_a*Y_bc
    for (int i = 0; i < new_tensor.num_elements(); ++i){
        // "i" enumerates the index abc
        // the following does: "i" -> abc
        int k = i;
        std::vector<int> i_index(N+M);
        for(int j = 0; j < N+M; ++j){
            i_index[j] = (k% (new_tensor.shape())[j]);
            k -= i_index[j];
            k /= (new_tensor.shape())[j];
        }
        new_tensor(i_index) = 
            t1( std::vector<int>(i_index.begin(),i_index.begin()+N) )
          * t2( std::vector<int>(i_index.begin()+N,i_index.end()) );
    }
    return new_tensor;
}


template<std::size_t N>
boost::multi_array<int,N-2> 
contract(boost::multi_array<int,N> t1, int x, int y){
    /*  eg. F_abca 
     *  This method contracts one pair of indices within a single tensor.
     */

    if ((t1.shape())[x] != (t1.shape())[y]) throw;

    std::vector<int> new_index;
    new_index.reserve(N-2);
    new_index.insert(new_index.end(), t1.shape(), t1.shape()+x);
    new_index.insert(new_index.end(), t1.shape()+x+1,   t1.shape()+x+y);
    new_index.insert(new_index.end(), t1.shape()+x+y+1, t1.shape()+t1.num_dimensions());

    boost::multi_array<int,N-2> ans(new_index);
    //
    // now we use new_index and old_index to generate the new tensor
    std::vector<int> old_index(N);
    for(int i = 0; i < ans.num_elements(); ++i){
        new_index = int_to_index(i, std::vector<int>(ans.shape(),ans.shape()+ans.num_dimensions()));
        ans(new_index) = 0;
        for (int j = 0; j < (t1.shape())[x]; ++j){
            old_index[x] = j; old_index[y] = j;
            copy(new_index.begin(),     new_index.begin()+x,   old_index.begin());
            copy(new_index.begin()+x,   new_index.begin()+y-1, old_index.begin()+x+1);
            copy(new_index.begin()+y-1, new_index.end(),       old_index.begin()+y+1);

            ans(new_index) += t1(old_index);
        }
    }
    return ans;
}

template<std::size_t N>
arma::mat tensor_to_matrix(boost::multi_array<int,N> t, 
        std::vector<int> first_index, 
        std::vector<int> second_index)
{
    int a = std::accumulate(first_index.begin(), first_index.end(), 1,
            [&t](int x, int y){ return x * (t.shape())[y]; });
    int b = std::accumulate(second_index.begin(), second_index.end(), 1,
            [&t](int x, int y){ return x * (t.shape())[y]; });

    arma::mat M(a,b);

    std::vector<int> a_shape(first_index.size());
    std::vector<int> b_shape(second_index.size());
    for(int i = 0; i < first_index.size(); ++i)    // a_shape definition = pull back of t_shape by first_index
        a_shape[i] = (t.shape())[first_index[i]];
    for(int j = 0; j < second_index.size(); ++j)
        b_shape[j] = (t.shape())[first_index[j]];

    std::vector<int> a_index(first_index.size());
    std::vector<int> b_index(second_index.size());
    std::vector<int> t_index(t.num_dimensions());
    for (int i = 0; i < a; ++i){
        a_index = int_to_index(i,a_shape);         // a_index < a_shape
        for (int k = 0; k < first_index.size(); ++k)
            t_index[first_index[k]] = a_index[k];  // t_index == pushforward of (a_index) by (first_index)
        for(int j = 0; j < b; ++j){
            b_index = int_to_index(j,b_shape);     // b_index < b_shape
            for (int k = 0; k < second_index.size(); ++k)
                t_index[second_index[k]] = b_index[k];//t_index = pushforward 

            M(i,j) = t(t_index);
        }
    }
    return M;
}

// SVD
// arma::svd( (mat) U, (vec) s, (mat) V, (mat) X)
// I assume this stores the singular value decomposition in
// the undetermined entries U,s,V
// Also *******REMEMBER******* to use                      g++ .... -larmadillo
// in case you actually use this function! 
/* Next: 
 *       Sweeping algorithm for DMRG
 *       Need to be able to represent Hamiltonian in an MPO state
 */

int main(int argv, char** argc)
{
    // TESTS
    boost::multi_array<int,1> T1(boost::extents[2]);
    boost::multi_array<int,3> T2(boost::extents[2][2][3]);
    boost::multi_array<int,4> T4(boost::extents[2][2][4][3]);

    T1[0] = 1;
    T1[1] = 2;
    T2[0][0][0] = 1;
    T2[0][1][0] = 2;
    T2[1][0][0] = 3;
    T2[1][1][0] = 4;
    T2[0][0][1] = 1; 
    T2[0][1][1] = 2;
    T2[1][0][1] = 3;
    T2[1][1][1] = 4; 
    T2[0][0][2] = 1; 
    T2[0][1][2] = 2;
    T2[1][0][2] = 3;
    T2[1][1][2] = 4; 

    //output(contract(T2,0,1)); // It works! :D

    // Testing tensor_to_matrix function
    auto M = tensor_to_matrix(T2,std::vector<int>{0,1},std::vector<int>{2});
    copy(M.begin(),M.end(), std::ostream_iterator<int>(std::cout," "));
    std::cout << std::endl;

    return 0;
}
