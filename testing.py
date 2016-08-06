import scipy
import numpy as np 

def prod ( lst ):
    if (len(lst) == 0):
        return 1;
    return lst[0]* prod(lst[1:])

def print_random_cpp_tensor( index_list, large_limit=10, tensor_name="A" ):
    N = prod ( index_list )
    lst = []
    for i in range(0,N):
        tmp = i
        lst = []
        for k in range(0, len(index_list)):
            lst.append( (tmp / prod( index_list[0:k] ))%index_list[k] );
        print(tensor_name+"(vi({"+(str(lst)[1:])[:-1]+"})) = "+str(scipy.random.random_integers(large_limit))+";");

def random_tensor( shape, large_limit=10 ):
    N = prod ( shape )
    lst = []
    ten = np.ndarray(shape, np.int32)
    for i in range(0,N):
        tmp = i
        lst = []
        for k in range(0, len(shape)):
            lst.append( (tmp / prod( shape[0:k] )) % shape[k] );
        ten[tuple(lst)] = scipy.random.random_integers(large_limit)
    return ten

def number_to_multi_index( num, shape ):
    lst = []
    for k in range(0, len(shape)):
        lst.append( (num / prod( shape[0:k] )) % shape[k] );
    return lst

def tensor_to_matrix( ten , mode=0 ):
    if (mode == 0):
        # (i0i1)(i2...in) = (alpha,beta)
        if (len(ten.shape) < 3):
            return -1;
        N1 = prod(ten.shape[0:2])
        N2 = prod(ten.shape[2: ])
        mat = np.ndarray([N1,N2])
        for i in range(0,N1):
            for j in range(0,N2):
                alpha = number_to_multi_index(i, ten.shape[0:2])
                beta  = number_to_multi_index(j, ten.shape[2: ])
                mat[i,j] = ten[tuple(alpha+beta)]
        return mat

def matrix_to_tensor( mat, row_shape, col_shape ):
    ten = np.ndarray(row_shape+col_shape, mat.dtype);
    for i in range(0, mat.shape[0]):
        for j in range(0,mat.shape[1]):
            alpha = number_to_multi_index(i, row_shape)
            beta  = number_to_multi_index(j, col_shape)
            index = number_to_multi_index(alpha+beta, ten.shape)
            ten[tuple(index)] = mat[i,j]
    return ten
