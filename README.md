# mps

--------------------------------------------------------------------------

## Tue Jun 21, 2016

    - the first code iteration was messy and not organized
    - began second iteration
    - goals: make modular and relatively well-documented 
             (at the very least readable)
    - Interesting Discoveries: multi_array<int,0> is supported and is
    accessed using

            a(vector<int>())
    - Questions: What is a const_iterator?

    - next: testing contract.cpp

--------------------------------------------------------------------------

## Wed Jun 22, 2016

    - contract.cpp seems to be working - no bugs found yet
    - trying to organize code (move some files into header)

--------------------------------------------------------------------------

## Thu Jun 23, 2016

    - contract_check1 and contract_check2 work! (or they were working...)
    - completely mangled the templates and now the errors are unreadable:
    full of instantiation problems.

--------------------------------------------------------------------------

## Wed Jul 13, 2016

    I'm back! Feynman diagram project now done. Time to focus efforts on
    coding for a bit.

    - Finished: explanation of dmrg using mps
        (~/Documents/masters/mps/mps-notes.pdf)
    - Found the reason why the instantiation problems were arising:
        -> added a new "feature": barA, barB, selection could all be
        vector-like containers, however, implemented this with templates
        WITHIN the cpp file (usually templated-functions need to be
                implemented in the header files)
        -> solution: 
            (a) enforce a vector structure for all such 
                index-type containers
            (b) enforce templated-functions ALWAYS go in header file

--------------------------------------------------------------------------

## Thu Jul 14, 2016

    - Question: Just renaming tests.cpp -> test.h fixed a
    multiple-definition problem; why?

    - All tests work. 

    - Question: Why in did I make barA and barB array's!?!
    - Question: Template instantiation requires a lot of overhead (compile
                    time is "supposedly" increased)
                I have not bumped into this, but is this true?
                Is there any way to mitigate this and save on compile time

    - MPS notes say that there are four things to implement:
        1. Full state -> MPS (implement epsilon-cutoff)
        2. Structure for contraction of < psi | psi > and  
            < psi | H | psi >
        3. Reformat tensor to matrix to tensor
        4. Generalized eigenvalue problem
        
    - Decreasing difficulty: 4,2,1,3

    - Began implementing 3; ran into a few problems with Eigen
    	-> solution: take the time to read Eigen's documentation fully
		     it will surely be beneficial later on!!!

--------------------------------------------------------------------------

## Fri Jul 16 2016

    - Finished Reading Eigen documentation

--------------------------------------------------------------------------

## Sun Jul 17 2016

    - Learn to implement operator<< for vector<int>
    - All tests for   tensor_to_matrix,  matrix_to_tensor   
      passed successfully.
    - July 14 - #3 done.

    - Next goal: J14 #4 (generalized eigenvalue problem)
      which is another self-contained problem!
         ----> SOLVED!

         Eigen::GeneralizedSelfAdjointEigenSolver(H,N)
         where $H$ is the Hamiltonian and $N$ is the auxillary matrix
    - July 14 - #4 done.

    - Next goal: J14 #2 -- restructure < psi | psi >

--------------------------------------------------------------------------

## Mon Jul 18 2016
    
    - Implemented J14 #2 -- vector of prefix- and suffix- contractions of 
      < psi | psi > and < psi | H | psi > (all of these are contained in
      DMRG_contractions.h
    - Also implemented an update procedure to the prefix/suffix vector for
      use in the DMRG sweep
    - I have not tested this part of the implementation
    - Next goal: implement the individual site optimization

