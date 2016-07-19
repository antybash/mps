#include <iostream>

#include "tensor.h"
#include "DMRG.h"
#include "tests.h"

template<typename T>
void out(T t) {
    copy(t.begin(),t.end(), std::ostream_iterator<int>(std::cout, " "));
    std::cout << std::endl;
}

int main()
{
    std::cout << check1() << std::endl;
    std::cout << check2() << std::endl;
    std::cout << check3() << std::endl;
    std::cout << check4() << std::endl;
    std::cout << check5() << std::endl;
    std::cout << check6() << std::endl;
    std::cout << check7() << std::endl;

    typedef std::vector<int> vi;
    tensor<int,1> A(vi({{2}}));
    A[0] = 1; A[1] = 2;
    tensor<int,1> B(vi({{2}}));
    B[0] = 3; B[1] = 4;
    
    output_tensor(contract(A,B,std::array<int,0>(),std::array<int,0>()));

    return 0;
}
