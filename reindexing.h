
#ifndef MPS_REINDEXING_H
#define MPS_REINDEXING_H

#include <vector>
#include <Eigen/Dense>

int prod(std::vector<int> A, int a, int b)
{
    int ans = 1;
    for(int i = a; i < b; ++i)
        ans *= A[i];
    return ans;
}

std::vector<int> triple(int a, int b, int c)
{
    return std::vector<int>({a,b,c});
}

std::vector<int> vector_subset(std::vector<int> V, int a, int b)
    //returns a vector [a, a+1, ..., b-1, b] *inclusive*
{
    std::vector<int> x(b-a+1);
    for (int i = 0; i < b-a+1; ++i)
        x[i] = V[a+i];
    return x;
}

template<size_t A, size_t B>
std::array<int,B-A+1>
array_range()
{
    std::array<int,B-A+1> x;
    for (int i = 0; i < B-A+1; ++i)
        x[i] = A+i;
    return x;
}

std::vector<int> vector_range(int a, int b)
    //returns a vector [a, a+1, ..., b-1, b] *inclusive*
{
    std::vector<int> x(b-a+1);
    for (int i = 0; i < b-a+1; ++i)
        x[i] = a+i;
    return x;
}

int multi_index_to_number(std::vector<int> index, std::vector<int> shape)
{
    // This method may be opposite to the expected base-10 expansion!!!
    // Example: if the shape is (222)
    //  (ijk)  -->  i + 2j + 4k
    //  
    //  if the shape is (5497)
    //  (pqrs) -->  p + 5q + 20r + 180s
    //
    //  where 20= 5*4
    //       180= 5*4*9

    int prod = 1, sum=0;
    for(int i = 0; i < shape.size(); ++i){
        sum  += index[i]*prod;
        prod *= shape[i];
    }
    return sum;
}

std::vector<int> number_to_multi_index(int alpha, std::vector<int> shape)
{
    std::vector<int> index(shape.size());
    for(int i = 0; i < shape.size(); ++i){
        index[i]  =  (alpha % shape[i]);
        alpha    /=  shape[i];
    }
    return index;
}

template<typename T, typename P>
std::vector<int> reindexing_subset(T container, P index_subset)
{
    std::vector<int> sub_container(index_subset.size());
    //for(auto i : index_subset){
    //    container[i]
    for (int i = 0; i < index_subset.size(); ++i)
        sub_container[i] = container[index_subset[i]];
    return sub_container;
}

Eigen::MatrixXcd
inplace_index_swap_of_underlying_tensor (Eigen::MatrixXcd M, std::vector<int> shape)
{
    /*
     *   Let $M \in Mat(alpha,beta)$, with 
     *     $alpha = shape[0]$
     *     $beta = shape[1]*shape[2]*...*shape[n-1]$. 
     *
     *   Define a new matrix $N \in Mat(gamma,delta)$, with
     *     $gamma = shape[0]*shape[1]$
     *     $delta = shape[2]*...*shape[n-1]$.
     */
    
    int n = shape.size();
    int gamma = prod(shape,0,2); //prod is clopen: [0,2)
    int delta = prod(shape,2,n);

    int alpha = shape[0];
    int beta  = prod(shape,1,n);

    assert(M.rows() == alpha);
    assert(M.cols() ==  beta); // prod is clopen: [1,n)

    Eigen::MatrixXcd N(gamma,delta);

    std::vector<int> I, J;
    int p,q;

    for(int i = 0; i < gamma; ++i){
        for(int j = 0; j < delta; ++j){
            I = number_to_multi_index(i, vector_subset(shape,0,1)); // vector_subset inclusive: [0,1]
            J = number_to_multi_index(j, vector_subset(shape,2,n-1)); // [2,n-1]
            J.insert(J.begin(),I[1]);

            p = I[0];
            q = multi_index_to_number(J, vector_subset(shape,1,n-1)); // [1,n-1]

            N(i,j) = M(p,q);
        }
    }
    return N;
}

#endif //MPS_REINDEXING_H
