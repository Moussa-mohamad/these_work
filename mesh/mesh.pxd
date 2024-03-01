
# INCLUDE IMPORTANT MEME SI APPAREMMENT PAS UTILISE !!!
from cython.view cimport array as cvarray
from libc.stdlib cimport malloc, free
# Import the C-level symbols of numpy
# cimport numpy as np
#Import the Python-level symbols of numpy
# import numpy as np



#include <stdio.h>
#include <iostream>
cdef extern from "meshKernel.hxx":
    ctypedef struct Output:
        int sparse_dim ;
        double  *eq_coefs ;
        int  *eq_cols ;
        int  *eq_rows ;
        double  *mesh_nodes ;
        int *TriNodes;
        int *FacesTriNum ;
        int pts_num;


cdef extern from "meshKernel.hxx":
    Output print_hello_c(double* blocks, int brows, int bcols,int blocks_num, double* nodes,int* active_faces,int active_faces_num,  int nrows, int ncols, int* faces_FEpts, int* Faces_nodes, int faces_num, double* blocks_centroid, double* c_local_ref, double lc )

    