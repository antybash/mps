#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "tensor.h"
#include "DMRG.h"
#include "utilities.h"

int main()
{
    double a = 1.0;
    double b = 1.0; 
    double c = 1.0;

    mpo_heis::initialize(a,b,c);

    std::ofstream outfile;
    std::stringstream in;
    in << "E0_";
    in << std::fixed << std::setprecision(3) << a; in << "_";
    in << std::fixed << std::setprecision(3) << b; in << "_";
    in << std::fixed << std::setprecision(3) << c; in << ".dat";
  
    std::cout << "Writing to: " << in.str() << std::endl;
  
    outfile.open(in.str());

    DMRG<4>  d2(3,(double)1e-3);  d2.sweep(10) ;  d2.output_lowest_energy (outfile);
    DMRG<6>  d4(3,(double)1e-3);  d4.sweep(10) ;  d4.output_lowest_energy (outfile);
    DMRG<8>  d6(3,(double)1e-3);  d6.sweep(10) ;  d6.output_lowest_energy (outfile);
    DMRG<10> d8(3,(double)1e-3);  d8.sweep(10) ;  d8.output_lowest_energy (outfile);
    DMRG<12> d10(3,(double)1e-3); d10.sweep(10) ; d10.output_lowest_energy(outfile);
    DMRG<14> d12(3,(double)1e-3); d12.sweep(10) ; d12.output_lowest_energy(outfile);
    DMRG<16> d14(3,(double)1e-3); d14.sweep(10) ; d14.output_lowest_energy(outfile);
    DMRG<18> d16(3,(double)1e-3); d16.sweep(10) ; d16.output_lowest_energy(outfile);
    DMRG<20> d18(3,(double)1e-3); d18.sweep(10) ; d18.output_lowest_energy(outfile);
    DMRG<22> d20(3,(double)1e-3); d20.sweep(10) ; d20.output_lowest_energy(outfile);
    DMRG<24> d22(3,(double)1e-3); d22.sweep(10) ; d22.output_lowest_energy(outfile);
    DMRG<26> d24(3,(double)1e-3); d24.sweep(10) ; d24.output_lowest_energy(outfile);
    DMRG<28> d26(3,(double)1e-3); d26.sweep(10) ; d26.output_lowest_energy(outfile);
    DMRG<30> d28(3,(double)1e-3); d28.sweep(10) ; d28.output_lowest_energy(outfile);
    DMRG<32> d30(3,(double)1e-3); d30.sweep(10) ; d30.output_lowest_energy(outfile);
    DMRG<34> d32(3,(double)1e-3); d32.sweep(10) ; d32.output_lowest_energy(outfile);
    DMRG<36> d34(3,(double)1e-3); d34.sweep(10) ; d34.output_lowest_energy(outfile);
    DMRG<38> d36(3,(double)1e-3); d36.sweep(10) ; d36.output_lowest_energy(outfile);
    DMRG<40> d38(3,(double)1e-3); d38.sweep(10) ; d38.output_lowest_energy(outfile);
    DMRG<42> d40(3,(double)1e-3); d40.sweep(10) ; d40.output_lowest_energy(outfile);
    DMRG<44> d42(3,(double)1e-3); d42.sweep(10) ; d42.output_lowest_energy(outfile);
    DMRG<46> d44(3,(double)1e-3); d44.sweep(10) ; d44.output_lowest_energy(outfile);
    DMRG<48> d46(3,(double)1e-3); d46.sweep(10) ; d46.output_lowest_energy(outfile);
    DMRG<50> d48(3,(double)1e-3); d48.sweep(10) ; d48.output_lowest_energy(outfile);
    DMRG<52> d50(3,(double)1e-3); d50.sweep(10) ; d50.output_lowest_energy(outfile);
    DMRG<54> d52(3,(double)1e-3); d52.sweep(10) ; d52.output_lowest_energy(outfile);
    DMRG<56> d54(3,(double)1e-3); d54.sweep(10) ; d54.output_lowest_energy(outfile);
    DMRG<58> d56(3,(double)1e-3); d56.sweep(10) ; d56.output_lowest_energy(outfile);
    DMRG<60> d58(3,(double)1e-3); d58.sweep(10) ; d58.output_lowest_energy(outfile);
    DMRG<62> d60(3,(double)1e-3); d60.sweep(10) ; d60.output_lowest_energy(outfile);
    DMRG<64> d62(3,(double)1e-3); d62.sweep(10) ; d62.output_lowest_energy(outfile);
    DMRG<66> d64(3,(double)1e-3); d64.sweep(10) ; d64.output_lowest_energy(outfile);
    DMRG<68> d66(3,(double)1e-3); d66.sweep(10) ; d66.output_lowest_energy(outfile);
    DMRG<70> d68(3,(double)1e-3); d68.sweep(10) ; d68.output_lowest_energy(outfile);
    DMRG<72> d70(3,(double)1e-3); d70.sweep(10) ; d70.output_lowest_energy(outfile);
    DMRG<74> d72(3,(double)1e-3); d72.sweep(10) ; d72.output_lowest_energy(outfile);
    DMRG<76> d74(3,(double)1e-3); d74.sweep(10) ; d74.output_lowest_energy(outfile);
    DMRG<78> d76(3,(double)1e-3); d76.sweep(10) ; d76.output_lowest_energy(outfile);
    DMRG<80> d78(3,(double)1e-3); d78.sweep(10) ; d78.output_lowest_energy(outfile);
    DMRG<82> d80(3,(double)1e-3); d80.sweep(10) ; d80.output_lowest_energy(outfile);
    DMRG<84> d82(3,(double)1e-3); d82.sweep(10) ; d82.output_lowest_energy(outfile);
    DMRG<86> d84(3,(double)1e-3); d84.sweep(10) ; d84.output_lowest_energy(outfile);
    DMRG<88> d86(3,(double)1e-3); d86.sweep(10) ; d86.output_lowest_energy(outfile);
    DMRG<90> d88(3,(double)1e-3); d88.sweep(10) ; d88.output_lowest_energy(outfile);
    DMRG<92> d90(3,(double)1e-3); d90.sweep(10) ; d90.output_lowest_energy(outfile);
    DMRG<94> d92(3,(double)1e-3); d92.sweep(10) ; d92.output_lowest_energy(outfile);
    DMRG<96> d94(3,(double)1e-3); d94.sweep(10) ; d94.output_lowest_energy(outfile);
    DMRG<98> d96(3,(double)1e-3); d96.sweep(10) ; d96.output_lowest_energy(outfile);


    return 0;
}

