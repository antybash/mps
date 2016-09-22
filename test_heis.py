def print_heis_test(n, mat_size=3, num_sweeps=10):
    print("DMRG<" + str(n+2) + "> d" + str(n) + "(" + str(mat_size) + ",(double)1e-3); d" + str(n) + ".sweep(" + str(num_sweeps) + ") ; d" + str(n) + ".output_lowest_energy();")

def test_case(a,b,c):
    for i in range(a,b,c):
        print_heis_test(i)
