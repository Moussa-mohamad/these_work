def mix_data(data ):
    
    import numpy as np

    Astat = np.array(data["Astat"])
    bstat = np.array(data["bstat"])
    Cstat = np.array(data["Gstat"])
    coh = np.array(data["hstat"])
    
    Akin = np.array(data["Akin"])
    bkin = np.array(data["bkin"])
    Gkin = np.array(data["Gkin"])
    zkin = np.array(data["zkin"])
    
    B = data["equilibrium_matrix"]
    N_mat =  data["disp_direction"]
    Live_loadsvect = data["live_loadvect"]
    Dead_loadsvect = data["dead_loadvect"]
    allpolygons = data["allpolygons"]
    
    nb = int(Astat.shape[0]/3)
    A1 = []
    A1 = np.array(Astat) 
    
    A = np.vstack((A1, Cstat))
    
    
    
    firstIndZ = A1.shape[1] # first element's index of the Z vector
    numbZ = Cstat.shape[0]  #number of elements of the vector Z
    A1 = np.hstack((A1, np.zeros((A1.shape[0],Cstat.shape[0]))))

    A1 = np.vstack((A1,np.hstack((np.array(Cstat),np.eye(Cstat.shape[0])))))

    G1 = np.hstack((np.zeros((Cstat.shape[0],Astat.shape[1])),-np.eye(Cstat.shape[0])))

    
    ulb1 = np.append(-bstat,coh)  
    
    firstIndd = A1.shape[1] + 3*len(allpolygons) # first element's index of the d vector
    numbd = Gkin.shape[0]  #number of elements of the vector d
    
    A2 = []
    A2 = np.append(B.T,N_mat.T,axis=1)

    
    A2 = np.array( [np.append(Live_loadsvect,np.zeros(len(zkin)-len(Live_loadsvect))), *A2 ])
    
    G2 = Gkin
    ulb2 = bkin


    def merge_block(matrix1,matrix2):

        result = np.block([[matrix1, np.zeros((matrix1.shape[0], matrix2.shape[1]))],
                           [np.zeros((matrix2.shape[0], matrix1.shape[1])), matrix2]])
        return result

    A = merge_block(A1,A2)
    lb = np.append(-bstat,-np.inf*np.ones(len(coh)))
    ub = np.append(-bstat,coh)
    
    zmix = np.array([*np.zeros(Astat.shape[1] -1),1.0])
    zmix = np.append(zmix,np.zeros(A.shape[1]-len(zmix)))

    G= merge_block(G1,G2)
    ulb = np.append(ulb1,ulb2)
    #ulb = np.append(ulb,-np.dot(Z0.T,d0))

    ub = np.zeros(G.shape[0]) #for inequalities
    lb = -999*np.ones(G.shape[0])  #for inequalities

    data["Amix"] = A
    data["bmix"] = ulb
    data["Gmix"] = G
    data["hmix"] = ub
    data["cmix"] = zmix
    data["numbZ"] = numbZ
    data["nb"] = nb
    print(nb)
    return data
