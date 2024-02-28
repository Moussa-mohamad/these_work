

# import numpy as np
cimport numpy as np
import numpy as np
# import numpy as np
np.import_array()

#include <unordered_set>




def print_hello_pyth(np.ndarray[int, ndim=1] active_faces, np.ndarray[double, ndim=2] blocks, np.ndarray[double, ndim=2] nodes , np.ndarray[int, ndim=1] faces_FEpts, np.ndarray[double, ndim=2] blocks_centroid, np.ndarray[double, ndim=2] local_ref ):
    


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

    #cdef int ref_rows = ref.shape[0]
    #cdef int ref_cols = ref.shape[1]

    cdef int* c_facesFEpts = <int*>faces_FEpts.data
    
    cdef int faces_num = faces_FEpts.shape[0] -1

    cdef brows = blocks.shape[0]
    cdef bcols = blocks.shape[1]

    cdef  Pointers pointer = print_hello_c(c_blocks, brows, bcols,blocks_num, c_nodes , c_active_faces,  active_faces_num, nrows, ncols, c_facesFEpts, faces_num, c_blocks_centroid, c_local_ref   )
  
   
    #res = np.asarray( pointer.ptr1, dtype=np.int32)
    #res = np.asarray(<np.float64_t[:nrows]> pointer.ptr1)


   


    #res = np.asarray(<int[:4]> pointer.ptr1).tolist()
    #print(res)

    #np.asarray(<np.float64_t[:cArray.shape[0]]> cArray.data).tolist()
    #res = np.asarray(pointer.ptr2, dtype=np.float64)[:20]

    res = [pointer.ptr2[i] for i in range(pointer.ptr1)]  # Assuming 4 elements in pointer.ptr1

    # af = <np.long[:rows,:cols]> pointer
    
    # anarrray = anarrray.tolist()
    return res
