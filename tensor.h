#include <iostream>
#include <vector>
#include <boost/multi_array.hpp>
#include <array>
#include <functional>
#include <utility>
#include <armadillo>


template<typename T, std::size_t N>
using tensor = boost::multi_array<T,N>;



/*****************************************
 * Output contents of N-component tensor *
 *****************************************/
template<typename T, std::size_t N> 
void output(boost::multi_array<T,N> a)
{
    auto f = [] (int x) { std::cout << x << " "; };
    std::for_each(a.origin(), a.origin()+a.num_elements(), f);
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

template<typename T, std::size_t N, std::size_t M> 
boost::multi_array<T, N+M> 
multiply_tensors(boost::multi_array<T,N> t1, boost::multi_array<T,M> t2)
{
    // declares the components of tensor
    std::vector<int> index; index.reserve(N+M);
    for (auto p = t1.shape(); p != t1.shape()+N; ++p)
        index.push_back(*p);
    for (auto p = t2.shape(); p != t2.shape()+M; ++p)
        index.push_back(*p);

    boost::multi_array<T, N+M> new_tensor(index);

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


/* ****************************************************************** * 
 *  eg. F_abca                                                        *
 *  This method contracts one pair of indices within a single tensor. *
 * ****************************************************************** */
template<typename T, std::size_t N>
boost::multi_array<T,N-2> 
contract(boost::multi_array<T,N> t1, int x, int y){
    /*  eg. F_abca 
     *  This method contracts one pair of indices within a single tensor.
     */

    if ((t1.shape())[x] != (t1.shape())[y]) throw;

    std::vector<int> new_index;
    new_index.reserve(N-2);
    new_index.insert(new_index.end(), t1.shape(), t1.shape()+x);
    new_index.insert(new_index.end(), t1.shape()+x+1,   t1.shape()+x+y);
    new_index.insert(new_index.end(), t1.shape()+x+y+1, t1.shape()+t1.num_dimensions());

    boost::multi_array<T,N-2> ans(new_index);
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

/*****************************
 * Tensor to matrix function *
 *****************************/
template<typename T, std::size_t N>
arma::Mat<T> tensor_to_matrix(boost::multi_array<T,N> t, 
        std::vector<int> first_index, 
        std::vector<int> second_index)
{
    int a = std::accumulate(first_index.begin(), first_index.end(), 1,
            [&t](int x, int y){ return x * (t.shape())[y]; });
    int b = std::accumulate(second_index.begin(), second_index.end(), 1,
            [&t](int x, int y){ return x * (t.shape())[y]; });

    arma::Mat<T> M(a,b);

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


 /**********************************************************************
 * Matrix to Tensor
 *      M_{ (a1a2..) (b1b2..) } -> T_{a1 b1 b2 a2 a3 b3 a4 ...}
 *
 *      the first component of M may not necessarily correspond to the
 *      first few components of the tensor (first_index tells us where
 *      to put the a_i's; first_shape tells us how big the dimension of
 *      each first_index has)
 ***********************************************************************/
template<typename T, std::size_t N>
boost::multi_array<T,N> matrix_to_tensor(arma::Mat<T> A, 
        std::vector<int> first_shape,
        std::vector<int> first_index, 
        std::vector<int> second_shape, 
        std::vector<int> second_index)
{
    // ***_shape == the size of each component
    // ***_index == where in the tensor it is placed
    // Therefore the sizes of these vectors must be equal!
    assert(first_shape.size() == first_index.size());
    assert(second_shape.size() == second_index.size()); 
    assert(N == first_index.size()+second_index.size());

    // Make sure the dimensions of the first and second component of A
    // will indeed bijectively correspond to the proposed tensor indices  
    assert(arma::size(A)[0] == 
        std::accumulate(first_shape.begin(),first_shape.end(), 1, std::multiplies<int>()));
    assert(arma::size(A)[1] == 
        std::accumulate(second_shape.begin(),second_shape.end(), 1, std::multiplies<int>()));

    // Initialize the shape of the tensor
    std::vector<int> total_shape(first_shape.size()+second_shape.size());
    for(int i = 0; i < first_index.size(); ++i)
        total_shape[first_index[i]] = first_shape[i];
    for(int j = 0; j < second_index.size(); ++j)
        total_shape[first_index[j]] = first_shape[j];

    boost::multi_array<T,N> final_tensor(total_shape);

    std::vector<int> tmp1_shape(first_shape.size());
    std::vector<int> tmp2_shape(second_shape.size());
    std::vector<int> tmp_total_shape(first_shape.size()+second_shape.size());

    for(int i = 0; i < arma::size(A)[0]; ++i){
        tmp1_shape = int_to_index(i,first_shape);
        for(int k = 0; k < first_index.size(); ++k)
            tmp_total_shape[first_index[k]] = tmp1_shape[k];
        for (int j = 0; j < arma::size(A)[1]; ++j){
            tmp2_shape = int_to_index(j,second_shape);
            for(int k = 0; k < second_index.size(); ++k)
                tmp_total_shape[second_index[k]] = tmp2_shape[k];
            
            // initialize the final_tensor
            final_tensor(tmp_total_shape) = A(i,j);
        }
    }
    return final_tensor;
}


/*********************
* Left Normalize MPS *
*********************/
template<typename T>
std::vector<tensor<T,3> > left_normalize_MPS (std::vector<tensor<T,3> > MPS)
{
    // Tested? No

    // Follows Section 4.4 
    // (Bringing a matrix product state into canonical form)

    // NOTE: Resizing is not implemented yet.
    // NOTE: MPS[i]     = M^{sg}_{a,b}  ==    M_[a,sg,b]
    //  notation:          paper             tensor implementation
    
    arma::Mat<T> U;
    arma::Col<T> s;
    arma::Mat<T> V;

    for (int i = 0; i < MPS.size()-1; ++i)
    {
        arma::svd(U,s,V,
                tensor_to_matrix(MPS[i],std::vector<int>{0,1},std::vector<int>{2}));
        MPS[i]   = matrix_to_tensor(U,
                      std::vector<int> { MPS[i].shape()[0], MPS[i].shape()[1] },  //first_shape
                      std::vector<int> {         0        ,         1         },  //first_index
                      std::vector<int> {               s.size()               },  //second_shape
                      std::vector<int> {                  2                   }); //second_index
        MPS[i+1] = contract( multiply_tensors( matrix_to_tensor(arma::diagmat(s)*V.t(),
                                                    std::vector<int> { s.size() },  //first_shape
                                                    std::vector<int> {     0    },  //first_index
                                                    std::vector<int> { s.size() },  //second_shape
                                                    std::vector<int> {     1    }), //second_index
                                               MPS[i+1]),
                             1,2); // indices to contract 
        // Explanation:
        // SV_{s,a} * M^{sg}_{a,b} => Tensor_[s,a,a,sg,b] == Tensor_[0,1,2,3,4] 
        // contract -> Tensor_[s,sg,b]
    }
    return MPS;
}

/***********************
* Vidal Representation *
***********************/
template<typename T>
std::vector<tensor<T,3> > MPS_to_Vidal_Form (std::vector<tensor<T,3> > MPS)
{
    // 2L+1 tensors in vidal form
    // odd-indexed (using 0,1,.. convention) are the diagonal elements with
    //      tensor structure T_[a,1,b], effectively a 2-tensor
}
