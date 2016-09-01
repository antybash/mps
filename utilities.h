
#ifndef UTILITIES_H
#define UTILITIES_H

namespace ar {
    std::array<int,1> zero      = {0};
    std::array<int,1> one       = {1};
    std::array<int,1> two       = {2};
    std::array<int,1> three     = {3};
    std::array<int,1> four      = {4};
    std::array<int,1> five      = {5};

    std::array<int,2> zeroone   = {0,1};
    std::array<int,2> zerotwo   = {0,2};
    std::array<int,2> onezero   = {1,0};
    std::array<int,2> onetwo    = {1,2};
    std::array<int,2> onethree  = {1,3};
    std::array<int,2> twothree  = {2,3};
    std::array<int,2> twoseven  = {2,7};
    std::array<int,2> threefour = {3,4};
    std::array<int,2> fourfive  = {4,5};
    std::array<int,2> threesix  = {3,6};

    std::array<int,3> a_012 = {0,1,2};
    std::array<int,3> a_123 = {1,2,3};
    std::array<int,3> a_135 = {1,3,5};
    std::array<int,3> a_345 = {3,4,5};

    std::array<int,4> a_1234 = {1,2,3,4};

    std::array<int,6> a_123456 = {1,2,3,4,5,6};
}

#endif
