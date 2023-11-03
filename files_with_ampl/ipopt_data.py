import os
import pandas as pd  # for pandas.DataFrame objects (https://pandas.pydata.org/)
import numpy as np  # for numpy.matrix objects (https://numpy.org/)
from amplpy import AMPL
from scipy.sparse import csr_matrix

def ipopt_data(Amix,bmix,cmix,length):
    
    length_edge = np.empty((0,2)) 
    b = np.empty((0,2))
    A = np.empty((0,4)) 
    c = np.empty((0,3))

    ind = 1;
    element = Amix

    sparse_matrix = csr_matrix(element)
    # Create a list of vectors for each row
    for row in range(sparse_matrix.shape[0]):

        row_data = sparse_matrix.getrow(row)
        elements = row_data.data
        column_indices = row_data.indices
        for element, col_idx in zip(elements, column_indices):

            A = np.append(A,np.array([[ind, row+1, col_idx, element]]),axis = 0) #row instaed of row +1 for mix
            ind = ind+1;

    element = cmix
    sparse_matrix = csr_matrix(element)

    # Create a list of vectors for each row
    ind = 1;
    for row in range(sparse_matrix.shape[0]):
        row_data = sparse_matrix.getrow(row)
        row_vector = []
        elements = row_data.data
        column_indices = row_data.indices
        for element, col_idx in zip(elements, column_indices):
            c = np.append(c,np.array([[ind,col_idx, element]]),axis = 0)
            ind = ind +1;




    element = bmix
    for i in range(len(element)):
        b= np.append(b,np.array([[i+1] + [element[i]]]),axis=0)
    element = length
    for i in range(len(element)):
        length_edge= np.append(length_edge,np.array([[i+1] + [element[i]]]),axis=0)
    
    
    
    
    A_df = pd.DataFrame(
       A,
        columns=["Aline","ALineInd","AColInd","AVal"],
    ).set_index("Aline")

    # Create a pandas.DataFrame with data for n_min, n_max
    b_df = pd.DataFrame(
        
          b
        ,
        columns=["bline", "BVal"],
    ).set_index("bline")

    c_df = pd.DataFrame(
          c
        ,
        columns=["cline", "CInd", "CVal"],).set_index("cline")
    
    len_df = pd.DataFrame(
          length_edge
        ,
        columns=["lenline", "LenVal" ],).set_index("lenline")
    
    
    return A_df, b_df, c_df, len_df
