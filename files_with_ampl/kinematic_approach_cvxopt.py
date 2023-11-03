def kinematic_approach_cvxopt(data,graph_rep ,scale_factor):   
    
    import time
    from cvxopt import solvers, matrix
    import numpy as np
    
    from graphical_rep import graphical_rep
    
    title = "Kinematic Approach Results"
    print("=" * 50)
    print(title.center(50))
    print("=" * 50)
    #Kinematic approach
    
    B = data["equilibrium_matrix"] 
    Dead_loadsvect = data["dead_loadvect"]  
    all_edges_length = data["edges_length"]   
    Live_loadsvect = data["live_loadvect"] 
    coh = data["cohesion_vector"]  
    
    
    Akin = data["Akin"]  
    bkin = data["bkin"]  
    Gkin = data["Gkin"]  
    hkin = data["hkin"]  
    zkin = data["zkin"] 
    

    start_kinematic = time.time()
    kinsolver = solvers.lp(zkin,Gkin,hkin,Akin,bkin)
    end_kinematic = time.time()

    print("Kinematic approach: ", end_kinematic - start_kinematic, "seconds")
    print(kinsolver['primal objective'])
    
    
    Dead_energie = 0
    Dead_energie = sum([Dead_energie  + kinsolver['x'][element]*Dead_loadsvect[element] for element in range(0,B.shape[0])])
    print("Dead loads energie ", + Dead_energie)
    
    Live_energie = 0
    Live_energie = sum([Live_energie  + kinsolver['x'][element]*Live_loadsvect[element] for element in range(0,B.shape[0])])
    
    print("live loads energie ", + Live_energie*kinsolver['y'][0])
    
    print("Internal loads energy",+ np.dot(np.dot(B.T,kinsolver['x'][0:B.shape[0]]).T,kinsolver['y'][1:]))
    
    print("Prm",+ np.dot(kinsolver['x'][-2*B.shape[1]:].T,coh))
    
    
    
    data["solutionvect_kin"] = kinsolver['x']
    data["lagrangevect_kin"] = np.append(np.array(kinsolver['y'][1:]), kinsolver['y'][0])
    
    
    if graph_rep:
        graphical_rep(data,-data["lagrangevect_kin"],data["solutionvect_kin"],data["solutionvect_kin"][B.shape[0]:],[],scale_factor,graph_title = "Kinematic approach" )
    
    return data