# mps

## Dependencies

Linear algebra: [*Eigen*](http://eigen.tuxfamily.org/)

## Running the code

terminal> g++ main.cpp -o dmrg -std=c++11 -g -I /PATH/TO/EIGEN/ 
terminal> ./dmrg

On my system this is how I run the program.

terminal> g++ main.cpp -o dmrg -std=c++11 -g -I /usr/local/include/eigen-eigen-1306d75b4a21 
terminal> ./dmrg

## How it works

The key equation that is being solved iteratively looks like:

    LWR.M = lam.M

The L denotes the *left* side of the triple contraction (psi|H|psi), W
denotes the local Hamiltonian, and R denotes the *right* side of the triple
contraction. M denotes the to-be-determined local MPS state.

A state of the DMRG consists of 4-stacks:

    - The L stack: denoted L3 for *triple* left contraction
    - The left-normalized stack:  denoted L1, consisting of the left-normalized MPS
    - The right-normalized stack: denoted R1, consisting of the right-normalized MPS
    - The R stack: denoted R3 for *triple* right contraction

A left-sweep:

    - Starts with a full L-stack and left-normalized stack and one-by-one
    pops off (or sweeps off) the left-stacks and pushes new updated
    right-normalized states onto the R1 and R3 stacks.

A right-sweep:

    - Sweeps the full R-stacks to and builds the L-stack.

