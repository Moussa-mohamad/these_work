

# import numpy as np
cimport numpy as np
import numpy as np
# import numpy as np
np.import_array()

#include <unordered_set>




def print_hello_pyth(np.ndarray[int, ndim=1] active_faces, np.ndarray[double, ndim=2] blocks, np.ndarray[double, ndim=2] nodes , np.ndarray[int, ndim=1] faces_FEpts, np.ndarray[int, ndim=1] Faces_nodes,  np.ndarray[double, ndim=2] blocks_centroid, np.ndarray[double, ndim=2] local_ref, double lc ):
    
    ### C++ variable construction
    cdef double c_lc = <int> lc

    cdef int* c_Faces_nodes = <int*>Faces_nodes.data

    cdef int* c_active_faces = <int*>active_faces.data
    cdef int active_faces_num = active_faces.shape[0]

    cdef double* c_nodes = <double*>nodes.data
    cdef double* c_blocks = <double*>blocks.data

    cdef double* c_blocks_centroid = <double*>blocks_centroid.data
    cdef double* c_local_ref = <double*>local_ref.data
   
    cdef int nrows = nodes.shape[0]
    cdef int ncols = nodes.shape[1]

    cdef int blocks_num = blocks_centroid.shape[0]
    cdef int centroid_cols = blocks_centroid.shape[1]

    cdef int* c_facesFEpts = <int*>faces_FEpts.data
    
    cdef int faces_num = faces_FEpts.shape[0] -1

    cdef brows = blocks.shape[0]
    cdef bcols = blocks.shape[1]

    ### Call C++ function 

    cdef  Output pointer = print_hello_c(c_blocks, brows, bcols,blocks_num, c_nodes , c_active_faces,  active_faces_num, nrows, ncols, c_facesFEpts, c_Faces_nodes, faces_num, c_blocks_centroid, c_local_ref,  lc   )
  

   
    #res = np.asarray( pointer.ptr1, dtype=np.int32)
    #res = np.asarray(<np.float64_t[:pointer.ptr1]> pointer.ptr2)
    
    ### Convert C++ to python 

    eq_coefs = [pointer.eq_coefs[i] for i in range(pointer.sparse_dim) ]  
    eq_cols = [pointer.eq_cols[i] for i in range(pointer.sparse_dim) ]  
    eq_rows = [pointer.eq_rows[i] for i in range(pointer.sparse_dim) ]  


    mesh_nodes = [[ pointer.mesh_nodes[3*j + i] for i in range(3) ] for j in range(pointer.pts_num) ]

    #TriNodes = [pointer.TriNodes[i] for i in range(pointer.sparse_dim) ]  

    FacesTriNum = [pointer.FacesTriNum[i] for i in range(faces_num + 1) ]
    
    TriNodes= [[[pointer.TriNodes[3*j + i ] for i in range(3) ]  for j in range( FacesTriNum[k],FacesTriNum[k+1] )   ] for k in range(faces_num) ] 

    return  eq_coefs, eq_rows, eq_cols, mesh_nodes, TriNodes, FacesTriNum

