#include <iostream>

#include "tests.h"
#include "tensor.h"

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


    return 0;
}
