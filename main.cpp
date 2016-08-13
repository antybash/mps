#include <iostream>

#include "tensor.h"
#include "DMRG.h"
#include "tests.h"


int main()
{
    // The real test begins.
    // Hamiltonian: ZZIII + IZZII + IIZZI + IIIZZ + IIIII

/*
    Eigen::Matrix2cd Z;
    Z(0,0) = 1; Z(0,1) =  0;
    Z(1,0) = 0; Z(1,1) = -1;

    Eigen::Matrix2cd I;
    I(0,0) = 1; I(0,1) = 0;
    I(1,0) = 0; I(1,1) = 1;

    tensor<cd,4> startH(vi({{1,2,2,3}}));
    for(int sg_up = 0; sg_up < 2; ++sg_up){
        for (int sg_down = 0; sg_down < 2; ++sg_down){
            startH[0][sg_down][sg_up][0] = I(sg_up,sg_down);
            startH[0][sg_down][sg_up][1] = Z(sg_up,sg_down);
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
            middleH[0][sg_down][sg_up][0] = I(sg_up,sg_down);
            middleH[2][sg_down][sg_up][2] = I(sg_up,sg_down);
        }

    // M_(1,2) == M_(2,3) == sigma_Z == Z
    for (int sg_up = 0; sg_up < 2; ++sg_up)
        for (int sg_down = 0; sg_down < 2; ++sg_down){
            middleH[1][sg_down][sg_up][0] = Z(sg_up,sg_down);
            middleH[2][sg_down][sg_up][1] = Z(sg_up,sg_down);
        }

    tensor<cd,4> endH(vi({{3,2,2,1}}));
    for(int sg_up = 0; sg_up < 2; ++sg_up){
        for (int sg_down = 0; sg_down < 2; ++sg_down){
            startH[0][sg_down][sg_up][0] = I(sg_up,sg_down);
            startH[1][sg_down][sg_up][0] = Z(sg_up,sg_down);
            startH[2][sg_down][sg_up][0] = I(sg_up,sg_down);
        }
    }

    // Hamiltonian building blocks are ready!

    tensor<cd,3> start, up, down, end;
*/
    ////////////////////////////////////////////
    // Testing left-normalization MPS creator //
    ////////////////////////////////////////////

    tensor<cd,5> A(vi({ 2,2,2,2,2 }));
    tensor<cd,3> B(vi({ 2,2,2 }));
    tensor<cd,5> C(vi({ 1,2,2,2,1 }));
/*
    A(vi({0, 0, 0, 0, 0})) = 3;
    A(vi({1, 0, 0, 0, 0})) = 2;
    A(vi({0, 1, 0, 0, 0})) = 2;
    A(vi({1, 1, 0, 0, 0})) = 6;
    A(vi({0, 0, 1, 0, 0})) = 9;
    A(vi({1, 0, 1, 0, 0})) = 3;
    A(vi({0, 1, 1, 0, 0})) = 6;
    A(vi({1, 1, 1, 0, 0})) = 1;
    A(vi({0, 0, 0, 1, 0})) = 5;
    A(vi({1, 0, 0, 1, 0})) = 6;
    A(vi({0, 1, 0, 1, 0})) = 7;
    A(vi({1, 1, 0, 1, 0})) = 8;
    A(vi({0, 0, 1, 1, 0})) = 9;
    A(vi({1, 0, 1, 1, 0})) = 3;
    A(vi({0, 1, 1, 1, 0})) = 3;
    A(vi({1, 1, 1, 1, 0})) = 3;
    A(vi({0, 0, 0, 0, 1})) = 3;
    A(vi({1, 0, 0, 0, 1})) = 1;
    A(vi({0, 1, 0, 0, 1})) = 7;
    A(vi({1, 1, 0, 0, 1})) = 4;
    A(vi({0, 0, 1, 0, 1})) = 3;
    A(vi({1, 0, 1, 0, 1})) = 3;
    A(vi({0, 1, 1, 0, 1})) = 10;
    A(vi({1, 1, 1, 0, 1})) = 2;
    A(vi({0, 0, 0, 1, 1})) = 3;
    A(vi({1, 0, 0, 1, 1})) = 1;
    A(vi({0, 1, 0, 1, 1})) = 7;
    A(vi({1, 1, 0, 1, 1})) = 5;
    A(vi({0, 0, 1, 1, 1})) = 9;
    A(vi({1, 0, 1, 1, 1})) = 9;
    A(vi({0, 1, 1, 1, 1})) = 7;
    A(vi({1, 1, 1, 1, 1})) = 9;

    B(vi({0, 0, 0})) = 10;
    B(vi({1, 0, 0})) = 9;
    B(vi({0, 1, 0})) = 7;
    B(vi({1, 1, 0})) = 2;
    B(vi({0, 0, 1})) = 10;
    B(vi({1, 0, 1})) = 5;
    B(vi({0, 1, 1})) = 6;
    B(vi({1, 1, 1})) = 3;
*/

    /*
    C(vi({0, 0, 0, 0, 0})) = 3;
    C(vi({0, 1, 0, 0, 0})) = 6;
    C(vi({0, 0, 1, 0, 0})) = 10;
    C(vi({0, 1, 1, 0, 0})) = 8;
    C(vi({0, 0, 0, 1, 0})) = 8;
    C(vi({0, 1, 0, 1, 0})) = 7;
    C(vi({0, 0, 1, 1, 0})) = 1;
    C(vi({0, 1, 1, 1, 0})) = 4;

    //std::vector<tensor<cd,3> > mpsA = tensor_to_left_normalized_mps(A,10e-4);
    auto mpsC = tensor_to_left_normalized_mps(C,10e-7);

    oe("main.cpp:\n\n\t Final tensors:");
    oe(mpsC.size());
    for(int i = 0; i < mpsC.size(); ++i){
        std::cout << "The " << i << "-th tensor has dimensions: ";
        std::cout << "(" << mpsC[i].shape()[0] << mpsC[i].shape()[1] << mpsC[i].shape()[2] << ")" << std::endl;
        output_tensor(mpsC[i]);
    }

    std::array<int,1> zero  = {0};
    std::array<int,1> two   = {2};
    std::array<int,1> three = {3};

    std::cout << "The 'final' test:" << std::endl;
    output_tensor( contract( contract(mpsC[0],mpsC[1],two,zero), mpsC[2], three, zero) );
    output_tensor( C );
    */

    check_svd_1();
    check_svd_2();
    check_svd_3();

    return 0;
}
