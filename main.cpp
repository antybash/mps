#include <iostream>

#include "tensor.h"
#include "DMRG.h"
#include "tests.h"

int check_index = 0;

std::array<int,1> zero      = {0};
std::array<int,1> one       = {1};
std::array<int,1> two       = {2};
std::array<int,1> three     = {3};
std::array<int,1> four      = {4};
std::array<int,1> five      = {5};

std::array<int,2> zeroone   = {0,1};
std::array<int,2> onezero   = {1,0};
std::array<int,2> onetwo    = {1,2};
std::array<int,2> twothree  = {2,3};
std::array<int,2> twoseven  = {2,7};
std::array<int,2> threefour = {3,4};
std::array<int,2> fourfive  = {4,5};
std::array<int,2> threesix  = {3,6};

std::array<int,3> a_123 = {1,2,3};
std::array<int,3> a_135 = {1,3,5};
std::array<int,3> a_345 = {3,4,5};

void check()
{
    std::cout << "main_check_" << check_index << std::endl;
    ++check_index;
}

int main()
{
    mpo::initialize();
    check_all();
    return 0;
}

