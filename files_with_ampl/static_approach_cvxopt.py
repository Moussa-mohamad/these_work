def static_approach_cvxopt(data,graph_rep ,scale_factor):
    import time
    from cvxopt import solvers, matrix
    import numpy as np
    
    
    #Static approach
    title = "Static Approach Results"
    print("=" * 50)
    print(title.center(50))
    print("=" * 50)
    Dead_loadsvect = data["dead_loadvect"]
    Live_loadsvect = data["live_loadvect"]
    C = data["criteria_matrix"]
    B = data["equilibrium_matrix"]
    coh = data["cohesion_vector"]
    testvect = data["testvect"]
    
    Astat =data["Astat"]
    bstat = data["bstat"] 
    Cstat = data["Gstat"] 
    cohstat = data["hstat"] 
    zstat = data["cstat"] 

    start = time.time()
    statsolver = solvers.lp(zstat,Cstat,cohstat,Astat,-bstat)
    end = time.time()
    


    print(statsolver['primal objective'])
    statvect = statsolver['x']

    print("Static approach: ", end - start, "seconds")
    
    #lagrange_stateq = np.array(statsolver['y'])
    lagrange_stateq = []
    for i in statsolver['y']:
        lagrange_stateq.append(i)
        
    
    lagrange_statineq = np.array(statsolver['z'])
   
        
    
    fourth = np.append(lagrange_stateq,lagrange_statineq)
    
    
    if graph_rep:
        from graphical_rep import graphical_rep
        graphical_rep(data,statvect,lagrange_stateq,lagrange_statineq,[],scale_factor,graph_title = "Static approach")
    
    
    Dead_energie = 0
    Dead_energie = sum([Dead_energie  + lagrange_stateq[element]*Dead_loadsvect[element] for element in range(0,len(lagrange_stateq))])
    print("Dead loads energie ", + Dead_energie)
    
    Live_energie = 0
    Live_energie = sum([Live_energie  + lagrange_stateq[element]*Live_loadsvect[element] for element in range(0,len(lagrange_stateq))])
    
    print("live loads energie ", + Live_energie*statsolver['x'][-1:])
    
    print("Internal loads energy",+ np.dot(np.dot(B.T,lagrange_stateq).T,statsolver['x'][0:Astat.size[1]-1]))
    print("Prm",+ np.dot(lagrange_statineq.T,coh))
   
    data["statsolver"] = statsolver
    data["solutionvect_stat"] = statvect
    data["lagrange_stateq"] = lagrange_stateq
    data["lagrange_statineq"] = lagrange_statineq
    data["lagrangevect_stat"] = fourth
    
    
    
    return data
    