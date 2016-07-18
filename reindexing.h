
#ifndef MPS_REINDEXING_H
#define MPS_REINDEXING_H

#include <vector>

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

#endif //MPS_REINDEXING_H
