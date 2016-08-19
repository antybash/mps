
/* 
 * This allows one to run over all possible indices 
 * for a given tensor type.
 */

#ifndef MPS_INCREMENT_INDICES_H
#define MPS_INCREMENT_INDICES_H

#include <vector>
#include <array>
#include <algorithm>
#include <cassert>

/// OUTPUT ///
template<typename T>
void out(T t) {
    copy(t.begin(),t.end(), std::ostream_iterator<int>(std::cout, " "));
    std::cout << std::endl;
}

template<typename T>
void oe(T s)
{
    std::cout << s << std::endl;
}
/// OUTPUT ///


bool check_all_zeros(std::vector<int> index)
{
    for(int i : index)
        if(i > 0) return false;
    return true;
}

template<typename T>
bool check_all_zeros_with_barrier(std::vector<int> index, T k)
{
    std::sort(k.begin(), k.end());
    for(int i = 0; i < index.size(); ++i)
        if( index[i] > 0 && std::find(k.begin(),k.end(),i) == k.end() ) return false;
    return true;
}

template<typename T>
bool check_all_zeros_with_selection(std::vector<int> index, T selection)
{
    std::sort(selection.begin(), selection.end());
    for(int i : selection)
        if( index[i] > 0 ) return false;
    return true;
}

template<typename T>
void increment_index(std::vector<int> &index, T shape)
{
    for (int i = 0; i < index.size(); ++i){
        if(index[i] < shape[i]-1){
            ++index[i];
            break;
        } else {
            index[i] = 0;
        }
    }
}

template<typename T>
void increment_index_with_barrier(std::vector<int> &index, std::vector<int> &shape, T k)
{
    std::sort(k.begin(),k.end());
    for (int i = 0; i < index.size(); ++i){
        if (std::find(k.begin(),k.end(),i) != k.end()) continue;
        if(index[i] < shape[i]-1){
            ++index[i];
            break;
        } else {
            index[i] = 0;
        }
    }
}

template<typename T>
void increment_index_with_selection(std::vector<int> &index, std::vector<int> &shape, T selection)
{
    std::sort(selection.begin(),selection.end());
    for (int i : selection){
        if(index[i] < shape[i]-1){
            ++index[i]; // this part is only encountered once!
            break;
        } else {
            index[i] = 0;
        }
    }
}


template<size_t dA, size_t dB>
void combine(std::vector<int> &res, std::vector<int> indA, std::vector<int> indB, std::array<int,dA> barA, std::array<int,dB> barB)
{
    assert(res.size() == indA.size()+indB.size()-barA.size()-barB.size());
    std::sort(barA.begin(),barA.end());
    std::sort(barB.begin(),barB.end());

    int len = 0;
    for(int i = 0; i < indA.size(); ++i){
        if( std::find(barA.begin(),barA.end(), i) != barA.end() ) 
            continue;
        else {
            res[len] = indA[i];
            ++len;
        }
    }
    for(int i = 0; i < indB.size(); ++i){
        if( std::find(barB.begin(),barB.end(), i) != barB.end() ) 
            continue;
        else {
            res[len] = indB[i];
            ++len;
        }
    }
}

#endif //MPS_INCREMENT_INDICES_H
