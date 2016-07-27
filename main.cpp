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
    // The real test begins.
    // Hamiltonian: ZZIII + IZZII + IIZZI + IIIZZ + IIIII

    typedef std::complex<double> cd;
    typedef std::vector<int> vi;

    Eigen::Matrix2cd Z;
    Z[0][0] = 1; Z[0][1] =  0;
    Z[1][0] = 0; Z[1][1] = -1;

    Eigen::Matrix2cd I;
    I[0][0] = 1; I[0][1] = 0;
    I[1][0] = 0; I[1][1] = 1;

    tensor<cd,4> startH(vi({{1,2,2,3}}));
    for(int sg_up = 0; sg_up < 2; ++sg_up){
        for (int sg_down = 0; sg_down < 2; ++sg_down){
            startH[0][sg_down][sg_up][0] = I[sg_up][sg_down];
            startH[0][sg_down][sg_up][1] = Z[sg_up][sg_down];
            startH[0][sg_down][sg_up][2] = 0;
        }
    }

    tensor<cd,4> middleH(vi({{3,2,2,3}}));

    // lower left 2x2 block == zero
    for(int i = 0; i < 2; ++i)
        for (int j = 1; j < 3; ++j)
            for (int sg_up = 0; sg_up < 2; ++sg_up)
                for (int sg_down = 0; sg_down < 2; ++sg_down)
                    middleH[j][sg_down][sg_up][i] = 0;

    // top-right matrix == zero
    for (int sg_up = 0; sg_up < 2; ++sg_up)
        for (int sg_down = 0; sg_down < 2; ++sg_down)
            middleH[2][sg_down][sg_up][0] = 0;

    // top-left == bottom-right == Identity
    for (int sg_up = 0; sg_up < 2; ++sg_up)
        for (int sg_down = 0; sg_down < 2; ++sg_down){
            middleH[0][sg_down][sg_up][0] = I[sg_up][sg_down];
            middleH[2][sg_down][sg_up][2] = I[sg_up][sg_down];
        }

    // M_(1,2) == M_(2,3) == sigma_Z == Z
    for (int sg_up = 0; sg_up < 2; ++sg_up)
        for (int sg_down = 0; sg_down < 2; ++sg_down){
            middleH[1][sg_down][sg_up][0] = Z[sg_up][sg_down];
            middleH[2][sg_down][sg_up][1] = Z[sg_up][sg_down];
        }

    tensor<cd,4> endH(vi({{3,2,2,1}}));
    for(int sg_up = 0; sg_up < 2; ++sg_up){
        for (int sg_down = 0; sg_down < 2; ++sg_down){
            startH[0][sg_down][sg_up][0] = I[sg_up][sg_down];
            startH[1][sg_down][sg_up][0] = Z[sg_up][sg_down];
            startH[2][sg_down][sg_up][0] = I[sg_up][sg_down];
        }
    }

    // Hamiltonian building blocks are ready!

    tensor<cd,3> start, up, down, end;


    return 0;
}
